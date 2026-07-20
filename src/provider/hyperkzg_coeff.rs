//! This module implements a variant of Nova's [`HyperKZG`](super::hyperkzg) evaluation engine
//! that treats the multilinear polynomial as a **dense monomial-coefficient vector** rather than a
//! Boolean-hypercube evaluation table. It is the "coefficient form" (`HyperKZGCoeffForm`) design.
//!
//! The commitment scheme, powers-of-tau setup, three-point `r, -r, r^2` reduction, polynomial
//! batching, and pairing batching are all identical to evaluation-form HyperKZG (see
//! [`super::hyperkzg`]); the KZG commitment layer is reused verbatim. The only differences are the
//! multilinear-to-univariate reduction rule and the verifier's transition check:
//!
//! * evaluation form binds a variable by interpolating endpoint values: `even + x_i * (odd - even)`;
//! * coefficient form substitutes into the monomial `a_0 + a_1 X_i`: `even + x_i * odd`.
//!
//! Correspondingly the verifier drops the `(1 - x_i)` factor from the transition equation, checking
//! `2 r Y_i == r (v_i^+ + v_i^-) + x_i (v_i^+ - v_i^-)`.
//!
//! Here `poly[k]` is the coefficient of the monomial selected by the binary expansion of `k`
//! (least-significant bit first). To keep the caller-facing interface identical to evaluation-form
//! HyperKZG, this engine reduces variables in the same order relative to `point` as the eval-form
//! engine: round `i` binds coefficient-index bit `i` using `point[ell - 1 - i]`. Equivalently, the
//! evaluation relation proved is
//! `y = sum_k poly[k] * prod_{j : bit_j(k) = 1} point[ell - 1 - j]`.
//! The caller MUST compute `eval` with this coefficient-form relation, not with the evaluation-table
//! MLE helper.
//!
//! The proof shape (`m - 1` intermediate commitments, three KZG witnesses, `3 m` field elements) and
//! all group/pairing costs match evaluation-form HyperKZG. See
//! `jcbase/00-hyperkzg-monomial-coefficient-form.md` for the full design, completeness, and
//! soundness argument (Gemini, eprint 2022/420 §2.4.2/§5; MicroNova §6/Appendix D).
//!
//! Domain separation: [`EvaluationArgument`] here is a distinct type from
//! [`super::hyperkzg::EvaluationArgument`], so a coefficient-form proof cannot be passed to the
//! evaluation-form verifier (or vice versa) at compile time. In addition, `prove`/`verify` inject
//! the Fiat-Shamir domain separator `hyperkzg-monomial-coeff-v1` before deriving any challenge, so
//! the two engines derive independent challenges even for identical inputs.
#![allow(non_snake_case)]
use crate::{
  errors::NovaError,
  provider::{
    hyperkzg::{Commitment, CommitmentEngine, CommitmentKey, ProverKey, VerifierKey},
    traits::{DlogGroup, DlogGroupExt, PairingGroup},
  },
  traits::{
    commitment::CommitmentEngineTrait, evaluation::EvaluationEngineTrait,
    evm_serde::EvmCompatSerde, Engine, TranscriptEngineTrait,
  },
};
use core::{iter, marker::PhantomData, slice};
use ff::Field;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use serde_with::serde_as;

/// Alias to points on G1 that are in preprocessed form
type G1Affine<E> = <<E as Engine>::GE as DlogGroup>::AffineGroupElement;

/// Default number of target chunks used in splitting up polynomial division in the kzg_open closure
const DEFAULT_TARGET_CHUNKS: usize = 1 << 10;

/// Fiat-Shamir domain separator distinguishing coefficient-form HyperKZG from the evaluation form.
const DOMAIN_SEPARATOR: &[u8] = b"hyperkzg-monomial-coeff-v1";

/// Provides an implementation of a coefficient-form polynomial evaluation argument.
///
/// This is a distinct type from [`super::hyperkzg::EvaluationArgument`] so that evaluation-form and
/// coefficient-form proofs cannot be confused at compile time.
#[serde_as]
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(bound = "")]
pub struct EvaluationArgument<E: Engine>
where
  E::GE: PairingGroup,
{
  com: Vec<Commitment<E>>,
  w: [Commitment<E>; 3],
  #[serde_as(as = "Vec<[EvmCompatSerde; 3]>")]
  v: Vec<[E::Scalar; 3]>,
}

impl<E: Engine> EvaluationArgument<E>
where
  E::GE: PairingGroup,
{
  /// Create a new evaluation argument
  pub fn new(com: Vec<Commitment<E>>, w: [Commitment<E>; 3], v: Vec<[E::Scalar; 3]>) -> Self {
    Self { com, w, v }
  }
  /// The KZG commitments to intermediate coefficient polynomials
  pub fn com(&self) -> &[Commitment<E>] {
    &self.com
  }
  /// The KZG witnesses for batch openings
  pub fn w(&self) -> &[Commitment<E>] {
    &self.w
  }
  /// The evaluations of the coefficient polynomials at the challenge points
  pub fn v(&self) -> &[[E::Scalar; 3]] {
    &self.v
  }
}

/// Provides an implementation of a coefficient-form polynomial evaluation engine using KZG.
///
/// It reuses evaluation-form HyperKZG's KZG commitment layer ([`super::hyperkzg::CommitmentEngine`],
/// [`ProverKey`], [`VerifierKey`]) and only changes the multilinear-to-univariate reduction to the
/// monomial-coefficient binding rule. See the [module documentation](self) for the coefficient-index
/// and variable-order conventions.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct EvaluationEngine<E: Engine> {
  _p: PhantomData<E>,
}

/// Convenience alias for the coefficient-form HyperKZG evaluation engine.
pub type HyperKZGCoeffForm<E> = EvaluationEngine<E>;

impl<E: Engine> EvaluationEngine<E>
where
  E::GE: PairingGroup,
{
  // Helper functions that are not part of `EvaluationEngineTrait` but are used to implement it.
  // The transcript labels match evaluation-form HyperKZG; separation from the eval form is provided
  // by the `DOMAIN_SEPARATOR` injected at the start of `prove`/`verify`.
  fn compute_challenge(com: &[Commitment<E>], transcript: &mut <E as Engine>::TE) -> E::Scalar {
    transcript.absorb(b"c", &com.to_vec().as_slice());

    transcript.squeeze(b"c").unwrap()
  }

  // Compute the polynomial-batching challenge q = Hash(vk, C0, ..., C_{k-1}, u0, ..., u_{t-1},
  // (p_i(u_j))_{i,j}).
  fn get_batch_challenge(v: &[[E::Scalar; 3]], transcript: &mut <E as Engine>::TE) -> E::Scalar {
    transcript.absorb(
      b"v",
      &v.iter()
        .flatten()
        .cloned()
        .collect::<Vec<E::Scalar>>()
        .as_slice(),
    );

    transcript.squeeze(b"r").unwrap()
  }

  fn batch_challenge_powers(q: E::Scalar, k: usize) -> Vec<E::Scalar> {
    // Compute powers of q : (1, q, q^2, ..., q^(k-1))
    let mut q_powers = vec![E::Scalar::ONE; k];
    for i in 1..k {
      q_powers[i] = q_powers[i - 1] * q;
    }
    q_powers
  }

  fn verifier_second_challenge(
    W: &[Commitment<E>],
    transcript: &mut <E as Engine>::TE,
  ) -> E::Scalar {
    transcript.absorb(b"W", &W.to_vec().as_slice());

    transcript.squeeze(b"d").unwrap()
  }
}

impl<E> EvaluationEngineTrait<E> for EvaluationEngine<E>
where
  E: Engine<CE = CommitmentEngine<E>>,
  E::GE: PairingGroup,
{
  type EvaluationArgument = EvaluationArgument<E>;
  type ProverKey = ProverKey<E>;
  type VerifierKey = VerifierKey<E>;

  fn setup(
    ck: &<E::CE as CommitmentEngineTrait<E>>::CommitmentKey,
  ) -> Result<(Self::ProverKey, Self::VerifierKey), NovaError> {
    // The KZG keys are basis-independent, so reuse the evaluation-form setup.
    <super::hyperkzg::EvaluationEngine<E> as EvaluationEngineTrait<E>>::setup(ck)
  }

  fn prove(
    ck: &CommitmentKey<E>,
    _pk: &Self::ProverKey,
    transcript: &mut <E as Engine>::TE,
    _C: &Commitment<E>,
    hat_P: &[E::Scalar],
    point: &[E::Scalar],
    _eval: &E::Scalar,
  ) -> Result<Self::EvaluationArgument, NovaError> {
    // Separate the Fiat-Shamir domain from evaluation-form HyperKZG before deriving any challenge.
    transcript.dom_sep(DOMAIN_SEPARATOR);

    let x: Vec<E::Scalar> = point.to_vec();

    //////////////// begin helper closures //////////
    let kzg_open = |f: &[E::Scalar], u: E::Scalar| -> Commitment<E> {
      // Divides polynomial f(x) by (x - u) to obtain the witness polynomial h(x) = f(x)/(x - u)
      // for KZG opening. This is identical to evaluation-form HyperKZG.
      let div_by_monomial =
        |f: &[E::Scalar], u: E::Scalar, target_chunks: usize| -> Vec<E::Scalar> {
          assert!(!f.is_empty());
          let target_chunk_size = f.len() / target_chunks;
          let log2_chunk_size = target_chunk_size.max(1).ilog2();
          let chunk_size = 1 << log2_chunk_size;

          let u_to_the_chunk_size = (0..log2_chunk_size).fold(u, |u_pow, _| u_pow.square());
          let mut result = f.to_vec();
          result
            .par_chunks_mut(chunk_size)
            .zip(f.par_chunks(chunk_size))
            .for_each(|(chunk, f_chunk)| {
              for i in (0..chunk.len() - 1).rev() {
                chunk[i] = f_chunk[i] + u * chunk[i + 1];
              }
            });

          let mut iter = result.chunks_mut(chunk_size).rev();
          if let Some(last_chunk) = iter.next() {
            let mut prev_partial = last_chunk[0];
            for chunk in iter {
              prev_partial = chunk[0] + u_to_the_chunk_size * prev_partial;
              chunk[0] = prev_partial;
            }
          }

          result[1..]
            .par_chunks_exact_mut(chunk_size)
            .rev()
            .for_each(|chunk| {
              let mut prev_partial = chunk[chunk_size - 1];
              for e in chunk.iter_mut().rev().skip(1) {
                prev_partial *= u;
                *e += prev_partial;
              }
            });
          result[1..].to_vec()
        };

      let target_chunks = DEFAULT_TARGET_CHUNKS;
      let h = &div_by_monomial(f, u, target_chunks);

      E::CE::commit(ck, h, &E::Scalar::ZERO)
    };

    let kzg_open_batch = |f: &[Vec<E::Scalar>],
                          u: &[E::Scalar; 3],
                          transcript: &mut <E as Engine>::TE|
     -> (Vec<Commitment<E>>, Vec<[E::Scalar; 3]>) {
      let poly_eval = |f: &[E::Scalar], u: E::Scalar| -> E::Scalar {
        // Horner's method
        let mut acc = E::Scalar::ZERO;
        for &fi in f.iter().rev() {
          acc = acc * u + fi;
        }

        acc
      };

      let scalar_vector_muladd = |a: &mut Vec<E::Scalar>, v: &Vec<E::Scalar>, s: E::Scalar| {
        assert!(a.len() >= v.len());
        a.par_iter_mut().zip(v.par_iter()).for_each(|(a_i, v_i)| {
          *a_i += s * *v_i;
        });
      };

      let kzg_compute_batch_polynomial = |f: &[Vec<E::Scalar>], q: E::Scalar| -> Vec<E::Scalar> {
        let k = f.len(); // Number of polynomials we're batching

        let q_powers = Self::batch_challenge_powers(q, k);

        // Compute B(x) = f[0] + q*f[1] + q^2 * f[2] + ... q^(k-1) * f[k-1]
        let mut B = f[0].clone();
        for i in 1..k {
          scalar_vector_muladd(&mut B, &f[i], q_powers[i]); // B += q_powers[i] * f[i]
        }

        B
      };
      ///////// END kzg_open_batch closure helpers

      let k = f.len();
      // Note: u.len() is always 3.

      // The verifier needs p_i(u_j), so we compute them here (V will compute B(u_j) itself)
      let mut v = vec![[E::Scalar::ZERO; 3]; k];
      v.par_iter_mut().zip_eq(f).for_each(|(v_j, f)| {
        // for each poly f
        v_j.par_iter_mut().enumerate().for_each(|(i, v_ij)| {
          // for each point u
          *v_ij = poly_eval(f, u[i]);
        });
      });

      let q = Self::get_batch_challenge(&v, transcript);
      let B = kzg_compute_batch_polynomial(f, q);

      // Now open B at u0, ..., u_{t-1}
      let w = u
        .into_par_iter()
        .map(|ui| kzg_open(&B, *ui))
        .collect::<Vec<Commitment<E>>>();

      // The prover computes the challenge to keep the transcript in the same
      // state as that of the verifier
      let _d_0 = Self::verifier_second_challenge(&w, transcript);

      (w, v)
    };

    ///// END helper closures //////////

    let ell = x.len();
    let n = hat_P.len();
    assert_eq!(n, 1 << ell); // Below we assume that n is a power of two

    // Phase 1 -- create commitments com_1, ..., com_{ell-1}.
    // We do not compute the final constant polynomial p_ell (and its commitment), as it equals
    // `eval`, which is known to the verifier and can be derived on its side as well.
    let mut polys: Vec<Vec<E::Scalar>> = Vec::new();
    polys.push(hat_P.to_vec());
    for i in 0..ell - 1 {
      let Pi_len = polys[i].len() / 2;
      let mut Pi = vec![E::Scalar::ZERO; Pi_len];

      // Coefficient-form binding: p_{i+1}[j] = even + x_i * odd (no `- even` interpolation term).
      #[allow(clippy::needless_range_loop)]
      Pi.par_iter_mut().enumerate().for_each(|(j, Pi_j)| {
        *Pi_j = polys[i][2 * j] + x[ell - i - 1] * polys[i][2 * j + 1];
      });

      polys.push(Pi);
    }

    // We do not need to commit to the first polynomial as it is already committed.
    // Compute commitments in parallel
    let r = vec![E::Scalar::ZERO; ell - 1];
    let com: Vec<Commitment<E>> = E::CE::batch_commit(ck, &polys[1..], r.as_slice());

    // Phase 2
    // We do not need to add x to the transcript, because in our context x was obtained from the
    // transcript. We also do not need to absorb `C` and `eval` as they are already absorbed by the
    // transcript by the caller.
    let r = Self::compute_challenge(&com, transcript);
    let u = [r, -r, r * r];

    // Phase 3 -- create response
    let (w, v) = kzg_open_batch(&polys, &u, transcript);

    Ok(EvaluationArgument {
      com,
      w: w.try_into().expect("w should have length 3"),
      v,
    })
  }

  /// A method to verify a purported coefficient-form evaluation of a multilinear polynomial.
  fn verify(
    vk: &Self::VerifierKey,
    transcript: &mut <E as Engine>::TE,
    C: &Commitment<E>,
    x: &[E::Scalar],
    y: &E::Scalar,
    pi: &Self::EvaluationArgument,
  ) -> Result<(), NovaError> {
    // Separate the Fiat-Shamir domain from evaluation-form HyperKZG before deriving any challenge.
    transcript.dom_sep(DOMAIN_SEPARATOR);

    let ell = x.len();

    // we do not need to add x to the transcript, because in our context x was
    // obtained from the transcript
    let r = Self::compute_challenge(&pi.com, transcript);

    let u = [r, -r, r * r];

    // Setup vectors (Y, ypos, yneg) from pi.v
    if pi.v.len() != ell || pi.com.len() != ell - 1 {
      return Err(NovaError::ProofVerifyError {
        reason: "Invalid lengths of pi.v".to_string(),
      });
    }

    // Check consistency of (Y, ypos, yneg) via the coefficient-form transition equation
    // 2 r Y == r (ypos + yneg) + x_i (ypos - yneg).
    for i in 0..ell {
      let ypos = pi.v[i][0];
      let yneg = pi.v[i][1];
      let Y = pi.v.get(i + 1).map_or(*y, |v| v[2]);
      if r.double() * Y != r * (ypos + yneg) + x[ell - i - 1] * (ypos - yneg) {
        return Err(NovaError::ProofVerifyError {
          reason: "Inconsistent (Y, ypos, yneg)".to_string(),
        });
      }
      // Note that we don't make any checks about Y[0] here, but our batching
      // check below requires it
    }

    // Check commitments to (Y, ypos, yneg) are valid

    // vk is hashed in transcript already, so we do not add it here

    let q = Self::get_batch_challenge(&pi.v, transcript);

    let d_0 = Self::verifier_second_challenge(&pi.w, transcript);
    let d_1 = d_0.square();

    // We write a special case for t=3, since this what is required for
    // hyperkzg. Following the paper directly, we must compute:
    // let L0 = C_B - vk.G * B_u[0] + W[0] * u[0];
    // let L1 = C_B - vk.G * B_u[1] + W[1] * u[1];
    // let L2 = C_B - vk.G * B_u[2] + W[2] * u[2];
    // let R0 = -W[0];
    // let R1 = -W[1];
    // let R2 = -W[2];
    // let L = L0 + L1*d_0 + L2*d_1;
    // let R = R0 + R1*d_0 + R2*d_1;
    //
    // We group terms to reduce the number of scalar mults (to seven):
    // In Rust, we could use MSMs for these, and speed up verification.
    //
    // Note, that while computing L, the intermediate computation of C_B together with computing
    // L0, L1, L2 can be replaced by single MSM of C with the powers of q multiplied by (1 + d_0 + d_1)
    // with additionally concatenated inputs for scalars/bases.

    let q_power_multiplier = E::Scalar::ONE + d_0 + d_1;

    let q_powers_multiplied: Vec<E::Scalar> =
      iter::successors(Some(q_power_multiplier), |qi| Some(*qi * q))
        .take(ell)
        .collect();

    // Compute the batched openings
    // compute B(u_i) = v[i][0] + q*v[i][1] + ... + q^(t-1) * v[i][t-1]
    let B_u = (0..3)
      .into_par_iter()
      .map(|i| {
        pi.v
          .iter()
          .rev()
          .fold(E::Scalar::ZERO, |acc, v_j| acc * q + v_j[i])
      })
      .collect::<Vec<E::Scalar>>();

    let com_affine: Vec<G1Affine<E>> = pi.com.iter().map(|c| c.into_inner().affine()).collect();
    let w_affine: Vec<G1Affine<E>> = pi.w.iter().map(|c| c.into_inner().affine()).collect();

    let L = E::GE::vartime_multiscalar_mul(
      &[
        &q_powers_multiplied[..],
        &[
          u[0],
          (u[1] * d_0),
          (u[2] * d_1),
          -(B_u[0] + d_0 * B_u[1] + d_1 * B_u[2]),
        ],
      ]
      .concat(),
      &[
        &[C.into_inner().affine()][..],
        &com_affine,
        &w_affine,
        slice::from_ref(&vk.G),
      ]
      .concat(),
    );

    let R0 = pi.w[0].into_inner();
    let R1 = pi.w[1].into_inner();
    let R2 = pi.w[2].into_inner();
    let R = R0 + R1 * d_0 + R2 * d_1;

    // Check that e(L, vk.H) == e(R, vk.tau_H)
    if (E::GE::pairing(&L, &DlogGroup::group(&vk.H)))
      != (E::GE::pairing(&R, &DlogGroup::group(&vk.tau_H)))
    {
      return Err(NovaError::ProofVerifyError {
        reason: "Pairing check failed".to_string(),
      });
    }

    Ok(())
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::provider::{keccak::Keccak256Transcript, Bn256EngineKZG};
  use rand::SeedableRng;

  type E = Bn256EngineKZG;
  type Fr = <E as Engine>::Scalar;

  // Ground-truth coefficient-form evaluation, independent of the prover's fold.
  // Matches the engine's variable order: coefficient-index bit `j` uses `point[ell - 1 - j]`.
  fn coeff_eval_direct(d: &[Fr], point: &[Fr]) -> Fr {
    let ell = point.len();
    assert_eq!(d.len(), 1 << ell);
    let mut acc = Fr::ZERO;
    for (k, &dk) in d.iter().enumerate() {
      let mut term = dk;
      for j in 0..ell {
        if (k >> j) & 1 == 1 {
          term *= point[ell - 1 - j];
        }
      }
      acc += term;
    }
    acc
  }

  // Reference implementation of the coefficient-form fold, matching the prover.
  fn coeff_fold(d: &[Fr], point: &[Fr]) -> Fr {
    let ell = point.len();
    let mut cur = d.to_vec();
    for i in 0..ell {
      let x_i = point[ell - 1 - i];
      let half = cur.len() / 2;
      let mut next = vec![Fr::ZERO; half];
      for j in 0..half {
        next[j] = cur[2 * j] + x_i * cur[2 * j + 1];
      }
      cur = next;
    }
    cur[0]
  }

  #[test]
  fn test_coeff_fold_matches_direct() {
    let mut rng = rand::rngs::StdRng::seed_from_u64(0xC0FFEE);
    for ell in 1..=8 {
      let n = 1usize << ell;
      let d = (0..n).map(|_| Fr::random(&mut rng)).collect::<Vec<_>>();
      let point = (0..ell).map(|_| Fr::random(&mut rng)).collect::<Vec<_>>();
      assert_eq!(coeff_fold(&d, &point), coeff_eval_direct(&d, &point));
    }
  }

  #[test]
  fn test_coeff_basis_vectors() {
    // Property test for bit ordering: f_{e_k}(point) = prod_{j : bit_j(k) = 1} point[ell - 1 - j].
    for ell in 1..=6 {
      let n = 1usize << ell;
      let point = (0..ell).map(|j| Fr::from(3 + j as u64)).collect::<Vec<_>>();
      for k in 0..n {
        let mut d = vec![Fr::ZERO; n];
        d[k] = Fr::ONE;
        let mut expected = Fr::ONE;
        for j in 0..ell {
          if (k >> j) & 1 == 1 {
            expected *= point[ell - 1 - j];
          }
        }
        assert_eq!(coeff_eval_direct(&d, &point), expected);
      }
    }
  }

  #[test]
  fn test_hyperkzg_coeff_eval() {
    // f(X0, X1) = 1 + 2 X0 + 3 X1 + 5 X0 X1 encoded LSB-first as d = (1, 2, 3, 5).
    // With the engine's ordering, X0 <-> point[1], X1 <-> point[0].
    let n = 4;
    let ck: CommitmentKey<E> = CommitmentEngine::setup(b"test", n).unwrap();
    let (pk, vk): (ProverKey<E>, VerifierKey<E>) = EvaluationEngine::setup(&ck).unwrap();

    let poly = vec![Fr::from(1), Fr::from(2), Fr::from(3), Fr::from(5)];
    let C = CommitmentEngine::commit(&ck, &poly, &Fr::ZERO);

    let test_inner = |point: Vec<Fr>, eval: Fr| -> Result<(), NovaError> {
      let mut tr = Keccak256Transcript::new(b"TestEval");
      let proof = EvaluationEngine::prove(&ck, &pk, &mut tr, &C, &poly, &point, &eval).unwrap();
      let mut tr = Keccak256Transcript::new(b"TestEval");
      EvaluationEngine::verify(&vk, &mut tr, &C, &point, &eval, &proof)
    };

    // Correct evaluations must verify. The eval is hard-coded (not computed by any helper) so this
    // test independently pins the coefficient-form relation
    //   y = 1 + 2*point[1] + 3*point[0] + 5*point[0]*point[1].
    for (point, eval) in [
      // point = [p0, p1]                       y = 1 + 2*p1 + 3*p0 + 5*p0*p1
      (vec![Fr::from(0), Fr::from(0)], Fr::from(1)), // 1 + 0 + 0 + 0
      (vec![Fr::from(0), Fr::from(1)], Fr::from(3)), // 1 + 2*1 + 0 + 0
      (vec![Fr::from(1), Fr::from(1)], Fr::from(11)), // 1 + 2 + 3 + 5
      (vec![Fr::from(3), Fr::from(5)], Fr::from(95)), // 1 + 10 + 9 + 75
      (vec![Fr::from(7), Fr::from(2)], Fr::from(96)), // 1 + 4 + 21 + 70
    ] {
      assert!(test_inner(point, eval).is_ok());
    }

    // Incorrect evaluation must fail: at point [3, 5] the correct value is 95, so 96 must be rejected.
    let point = vec![Fr::from(3), Fr::from(5)];
    assert!(test_inner(point, Fr::from(96)).is_err());
  }

  #[test]
  fn test_hyperkzg_coeff_edge_cases() {
    // m = 1 (no intermediate commitments), all-zero, and constant polynomials.
    let cases: Vec<(Vec<Fr>, Vec<Fr>)> = vec![
      (vec![Fr::from(9), Fr::from(4)], vec![Fr::from(6)]), // m = 1
      (
        vec![Fr::ZERO; 8],
        vec![Fr::from(2), Fr::from(3), Fr::from(4)],
      ), // all-zero
      (
        {
          let mut d = vec![Fr::ZERO; 8];
          d[0] = Fr::from(7); // constant polynomial
          d
        },
        vec![Fr::from(2), Fr::from(3), Fr::from(4)],
      ),
    ];

    for (poly, point) in cases {
      let n = poly.len();
      let ck: CommitmentKey<E> = CommitmentEngine::setup(b"test", n).unwrap();
      let (pk, vk) = EvaluationEngine::setup(&ck).unwrap();
      let C = CommitmentEngine::commit(&ck, &poly, &Fr::ZERO);
      let eval = coeff_eval_direct(&poly, &point);

      let mut tr = Keccak256Transcript::new(b"TestEval");
      let proof =
        EvaluationEngine::<E>::prove(&ck, &pk, &mut tr, &C, &poly, &point, &eval).unwrap();
      let mut tr = Keccak256Transcript::new(b"TestEval");
      assert!(EvaluationEngine::verify(&vk, &mut tr, &C, &point, &eval, &proof).is_ok());
    }
  }

  #[test]
  #[cfg(not(feature = "evm"))]
  fn test_hyperkzg_coeff_small() {
    let n = 4;
    let poly = vec![Fr::ONE, Fr::from(2), Fr::from(1), Fr::from(4)];
    let point = vec![Fr::from(4), Fr::from(3)];
    let eval = coeff_eval_direct(&poly, &point);

    let ck: CommitmentKey<E> = CommitmentEngine::setup(b"test", n).unwrap();
    let (pk, vk) = EvaluationEngine::setup(&ck).unwrap();
    let C = CommitmentEngine::commit(&ck, &poly, &Fr::ZERO);

    // prove an evaluation
    let mut prover_transcript = Keccak256Transcript::new(b"TestEval");
    let proof =
      EvaluationEngine::<E>::prove(&ck, &pk, &mut prover_transcript, &C, &poly, &point, &eval)
        .unwrap();
    let post_c_p = prover_transcript.squeeze(b"c").unwrap();

    // verify the evaluation
    let mut verifier_transcript = Keccak256Transcript::new(b"TestEval");
    assert!(
      EvaluationEngine::verify(&vk, &mut verifier_transcript, &C, &point, &eval, &proof).is_ok()
    );
    let post_c_v = verifier_transcript.squeeze(b"c").unwrap();

    // prover and verifier transcripts must be kept in the same state
    assert_eq!(post_c_p, post_c_v);

    // proof size matches evaluation-form HyperKZG (same shape)
    let config = bincode::config::legacy()
      .with_big_endian()
      .with_fixed_int_encoding();
    let proof_bytes =
      bincode::serde::encode_to_vec(&proof, config).expect("Failed to serialize proof");
    assert_eq!(proof_bytes.len(), 336);

    // Change the proof and expect verification to fail
    let mut bad_proof = proof.clone();
    let v1 = bad_proof.v[1];
    bad_proof.v[0].clone_from(&v1);
    let mut verifier_transcript2 = Keccak256Transcript::new(b"TestEval");
    assert!(EvaluationEngine::verify(
      &vk,
      &mut verifier_transcript2,
      &C,
      &point,
      &eval,
      &bad_proof
    )
    .is_err());
  }

  #[test]
  fn test_hyperkzg_coeff_large() {
    for ell in [4, 5, 6] {
      let mut rng = rand::rngs::StdRng::seed_from_u64(ell as u64);
      let n = 1 << ell;

      let poly = (0..n).map(|_| Fr::random(&mut rng)).collect::<Vec<_>>();
      let point = (0..ell).map(|_| Fr::random(&mut rng)).collect::<Vec<_>>();
      let eval = coeff_eval_direct(&poly, &point);

      let ck: CommitmentKey<E> = CommitmentEngine::setup(b"test", n).unwrap();
      let (pk, vk) = EvaluationEngine::setup(&ck).unwrap();
      let C = CommitmentEngine::commit(&ck, &poly, &Fr::ZERO);

      let mut prover_transcript = Keccak256Transcript::new(b"TestEval");
      let proof: EvaluationArgument<E> =
        EvaluationEngine::prove(&ck, &pk, &mut prover_transcript, &C, &poly, &point, &eval)
          .unwrap();

      let mut verifier_tr = Keccak256Transcript::new(b"TestEval");
      assert!(EvaluationEngine::verify(&vk, &mut verifier_tr, &C, &point, &eval, &proof).is_ok());

      // Tampering with an evaluation must fail.
      let mut bad_proof = proof.clone();
      let v1 = bad_proof.v[1];
      bad_proof.v[0].clone_from(&v1);
      let mut verifier_tr2 = Keccak256Transcript::new(b"TestEval");
      assert!(
        EvaluationEngine::verify(&vk, &mut verifier_tr2, &C, &point, &eval, &bad_proof).is_err()
      );
    }
  }

  #[test]
  fn test_hyperkzg_coeff_wrong_point_and_value() {
    let mut rng = rand::rngs::StdRng::seed_from_u64(7);
    let ell = 5;
    let n = 1 << ell;
    let poly = (0..n).map(|_| Fr::random(&mut rng)).collect::<Vec<_>>();
    let point = (0..ell).map(|_| Fr::random(&mut rng)).collect::<Vec<_>>();
    let eval = coeff_eval_direct(&poly, &point);

    let ck: CommitmentKey<E> = CommitmentEngine::setup(b"test", n).unwrap();
    let (pk, vk) = EvaluationEngine::setup(&ck).unwrap();
    let C = CommitmentEngine::commit(&ck, &poly, &Fr::ZERO);

    let mut tr = Keccak256Transcript::new(b"TestEval");
    let proof = EvaluationEngine::<E>::prove(&ck, &pk, &mut tr, &C, &poly, &point, &eval).unwrap();

    // Wrong claimed value.
    let mut tr = Keccak256Transcript::new(b"TestEval");
    assert!(EvaluationEngine::verify(&vk, &mut tr, &C, &point, &(eval + Fr::ONE), &proof).is_err());

    // Perturb a point coordinate that feeds an intermediate (committed) folding round. Even with a
    // consistently recomputed value, the proof (whose intermediate commitments were folded with the
    // original coordinate) must fail an intermediate transition check. Note the last coordinate,
    // `point[0]`, only enters the final linear check, so perturbing it with a matching value would
    // legitimately verify; hence we perturb `point[ell - 1]`, used in round 0.
    let mut bad_point = point.clone();
    bad_point[ell - 1] += Fr::ONE;
    let bad_eval = coeff_eval_direct(&poly, &bad_point);
    let mut tr = Keccak256Transcript::new(b"TestEval");
    assert!(EvaluationEngine::verify(&vk, &mut tr, &C, &bad_point, &bad_eval, &proof).is_err());
  }
}
