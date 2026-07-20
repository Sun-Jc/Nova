//! A representation-changing **adapter** that lets any evaluation-form multilinear PCS
//! (e.g. [`hyperkzg`](super::hyperkzg) or [`mercury`](super::mercury)) prove and verify
//! **monomial-coefficient-form** evaluations *without changing its commitment or its opening
//! algorithm*. The whole change is a point transform applied before the underlying opening.
//!
//! # What it does
//!
//! An evaluation-form multilinear PCS proves, for a committed vector `v` of length `N = 2^m`,
//! the standard multilinear evaluation
//! `y = <v, ⊗_j (1 - r_j, r_j)>` (i.e. `v` is read as a Boolean-hypercube evaluation table).
//!
//! This adapter proves instead the **coefficient-form** evaluation
//! `y = <v, ⊗_j (1, r_j)>` (i.e. `v` is read as the dense monomial-coefficient vector of a
//! multilinear polynomial). See the `hyperkzg_coeff` coefficient-form engine and Gemini
//! (eprint 2022/420 §2.4.2) for the coefficient-form relation.
//!
//! # The transform
//!
//! The two tensor factors are related, per variable, by
//! `(1, r_j) = (1 + r_j) · (1 - r'_j, r'_j)` with `r'_j = r_j / (1 + r_j)`. Taking the tensor
//! product over all `m` variables gives the exact identity
//!
//! ```text
//! <v, ⊗_j (1, r_j)> = ( ∏_j (1 + r_j) ) · <v, ⊗_j (1 - r'_j, r'_j)>
//! ```
//!
//! so a coefficient-form claim `(r, y)` is equivalent to the evaluation-form claim
//! `(r', y / scale)` with `r'_j = r_j / (1 + r_j)` and `scale = ∏_j (1 + r_j)`. Because the
//! transform acts **coordinate-wise** and `scale` is symmetric, it is independent of the engine's
//! internal variable ordering, so the same adapter works for any evaluation-form engine.
//!
//! The correction factor `scale` is public (the verifier computes it from `r` alone) and is folded
//! into the claimed value, so:
//!
//! 1. the commitment algorithm is unchanged — the same `v` is committed with the same KZG commit;
//! 2. the opening algorithm is unchanged — the adapter calls the underlying `prove`/`verify` as-is;
//! 3. opening at `r` reduces to opening the *unmodified* engine at `r'`, and the result equals
//!    reading `v` as monomial coefficients and evaluating at `r`.
//!
//! The transform is undefined when `1 + r_j = 0` for some `j` (i.e. `r_j = -1`); the adapter then
//! returns [`NovaError::ProofVerifyError`]. Under a Fiat-Shamir–sampled point this occurs with
//! negligible probability.
use crate::{
  errors::NovaError,
  traits::{commitment::CommitmentEngineTrait, evaluation::EvaluationEngineTrait, Engine},
};
use core::marker::PhantomData;
use ff::Field;

/// Transform a coefficient-form evaluation point into the evaluation-form point plus the public
/// correction factor.
///
/// Given `point = r`, returns `Some((r', scale))` with `r'_j = r_j / (1 + r_j)` and
/// `scale = ∏_j (1 + r_j)`, satisfying, for every vector `v`,
/// `<v, ⊗_j (1, r_j)> = scale · <v, ⊗_j (1 - r'_j, r'_j)>`.
///
/// Returns `None` if `1 + r_j = 0` for some coordinate (the transform is undefined there).
pub fn coeff_eval_point<F: Field>(point: &[F]) -> Option<(Vec<F>, F)> {
  let mut transformed = Vec::with_capacity(point.len());
  let mut scale = F::ONE;
  for &r_j in point {
    let denom = F::ONE + r_j;
    // `1 + r_j == 0` (i.e. r_j == -1) makes the transform undefined.
    let denom_inv = Option::<F>::from(denom.invert())?;
    transformed.push(r_j * denom_inv);
    scale *= denom;
  }
  Some((transformed, scale))
}

/// An [`EvaluationEngineTrait`] adapter that turns an evaluation-form multilinear PCS `EE` into a
/// coefficient-form one by transforming the opening point (see the [module docs](self)).
///
/// The commitment scheme, prover key, verifier key, and evaluation-argument type are all inherited
/// unchanged from `EE`; only the opening point and claimed value are transformed before delegating
/// to `EE`'s unmodified `prove`/`verify`.
#[derive(Clone, Debug)]
pub struct CoeffEvaluationEngine<E: Engine, EE: EvaluationEngineTrait<E>> {
  // `fn() -> ...` keeps the marker unconditionally `Send + Sync + Clone` without constraining E/EE.
  _p: PhantomData<fn() -> (E, EE)>,
}

/// Coefficient-form HyperKZG obtained by adapting the evaluation-form [`hyperkzg`](super::hyperkzg)
/// engine with the point transform.
pub type HyperKZGCoeffAdapter<E> = CoeffEvaluationEngine<E, super::hyperkzg::EvaluationEngine<E>>;

/// Coefficient-form Mercury obtained by adapting the evaluation-form [`mercury`](super::mercury)
/// engine with the point transform.
pub type MercuryCoeffAdapter<E> = CoeffEvaluationEngine<E, super::mercury::EvaluationEngine<E>>;

impl<E: Engine, EE: EvaluationEngineTrait<E>> EvaluationEngineTrait<E>
  for CoeffEvaluationEngine<E, EE>
{
  type ProverKey = EE::ProverKey;
  type VerifierKey = EE::VerifierKey;
  type EvaluationArgument = EE::EvaluationArgument;

  fn setup(
    ck: &<<E as Engine>::CE as CommitmentEngineTrait<E>>::CommitmentKey,
  ) -> Result<(Self::ProverKey, Self::VerifierKey), NovaError> {
    // The keys are basis-independent; reuse the underlying engine's setup verbatim.
    EE::setup(ck)
  }

  fn prove(
    ck: &<<E as Engine>::CE as CommitmentEngineTrait<E>>::CommitmentKey,
    pk: &Self::ProverKey,
    transcript: &mut E::TE,
    comm: &<<E as Engine>::CE as CommitmentEngineTrait<E>>::Commitment,
    poly: &[E::Scalar],
    point: &[E::Scalar],
    eval: &E::Scalar,
  ) -> Result<Self::EvaluationArgument, NovaError> {
    let (point_eval, scale) = transform_or_err(point)?;
    let eval_eval = *eval * invert_scale(scale);
    // Same commitment, same vector, unchanged opening algorithm — only the point/value change.
    EE::prove(ck, pk, transcript, comm, poly, &point_eval, &eval_eval)
  }

  fn verify(
    vk: &Self::VerifierKey,
    transcript: &mut E::TE,
    comm: &<<E as Engine>::CE as CommitmentEngineTrait<E>>::Commitment,
    point: &[E::Scalar],
    eval: &E::Scalar,
    arg: &Self::EvaluationArgument,
  ) -> Result<(), NovaError> {
    let (point_eval, scale) = transform_or_err(point)?;
    let eval_eval = *eval * invert_scale(scale);
    EE::verify(vk, transcript, comm, &point_eval, &eval_eval, arg)
  }
}

/// Apply [`coeff_eval_point`], mapping the undefined case to a verification error.
fn transform_or_err<F: Field>(point: &[F]) -> Result<(Vec<F>, F), NovaError> {
  coeff_eval_point(point).ok_or_else(|| NovaError::ProofVerifyError {
    reason: "coefficient-form point transform undefined: 1 + r_j == 0".to_string(),
  })
}

/// Invert the correction factor. `scale` is a product of nonzero field elements, hence invertible.
fn invert_scale<F: Field>(scale: F) -> F {
  Option::<F>::from(scale.invert()).expect("scale is a product of nonzero factors, so invertible")
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::{
    provider::{hyperkzg, keccak::Keccak256Transcript, mercury, Bn256EngineKZG},
    spartan::polys::multilinear::MultilinearPolynomial,
    traits::TranscriptEngineTrait,
  };
  use ff::Field;
  use rand::SeedableRng;
  use rand_core::OsRng;

  type E = Bn256EngineKZG;
  type Fr = <E as Engine>::Scalar;
  type CE = <E as Engine>::CE;

  // Ground-truth coefficient-form evaluation using the MSB-first convention that matches
  // `MultilinearPolynomial::evaluate_with` (coordinate `c` binds index bit `m-1-c`).
  fn coeff_eval_msb(v: &[Fr], point: &[Fr]) -> Fr {
    let m = point.len();
    assert_eq!(v.len(), 1 << m);
    let mut acc = Fr::ZERO;
    for (k, &vk) in v.iter().enumerate() {
      let mut term = vk;
      for (c, &pc) in point.iter().enumerate() {
        if (k >> (m - 1 - c)) & 1 == 1 {
          term *= pc;
        }
      }
      acc += term;
    }
    acc
  }

  #[test]
  fn test_point_transform_identity() {
    // coeff_eval_msb(v, r) == scale * evaluate_with(v, r'), for the transform (r', scale).
    let mut rng = rand::rngs::StdRng::seed_from_u64(0xC0FFEE);
    for m in 1..=8 {
      let n = 1usize << m;
      let v = (0..n).map(|_| Fr::random(&mut rng)).collect::<Vec<_>>();
      let r = (0..m).map(|_| Fr::random(&mut rng)).collect::<Vec<_>>();
      let (rp, scale) = coeff_eval_point(&r).unwrap();
      let lhs = coeff_eval_msb(&v, &r);
      let rhs = scale * MultilinearPolynomial::evaluate_with(&v, &rp);
      assert_eq!(lhs, rhs);
    }
  }

  #[test]
  fn test_point_transform_degenerate() {
    // 1 + r_j == 0 (r_j == -1) makes the transform undefined.
    let r = vec![Fr::from(3), -Fr::ONE, Fr::from(5)];
    assert!(coeff_eval_point(&r).is_none());
  }

  // End-to-end: prove/verify a coefficient-form evaluation through the adapter, and independently
  // confirm the adapter's proof is exactly the underlying engine's proof at the transformed point.
  #[test]
  fn test_adapter_hyperkzg() {
    type Adapter = HyperKZGCoeffAdapter<E>;
    for m in [1usize, 2, 4, 5] {
      let n = 1 << m;
      let mut rng = rand::rngs::StdRng::seed_from_u64(m as u64);
      let v = (0..n).map(|_| Fr::random(&mut rng)).collect::<Vec<_>>();
      let r = (0..m).map(|_| Fr::random(&mut rng)).collect::<Vec<_>>();

      let ck = <CE as CommitmentEngineTrait<E>>::CommitmentKey::setup_from_rng(b"test", n, OsRng);
      let (pk, vk) = Adapter::setup(&ck).unwrap();
      let comm = CE::commit(&ck, &v, &Fr::ZERO);
      let eval = coeff_eval_msb(&v, &r); // coefficient-form value at r

      // Adapter prove/verify round-trip.
      let mut tr = Keccak256Transcript::new(b"T");
      let arg = Adapter::prove(&ck, &pk, &mut tr, &comm, &v, &r, &eval).unwrap();
      let mut tr = Keccak256Transcript::new(b"T");
      assert!(Adapter::verify(&vk, &mut tr, &comm, &r, &eval, &arg).is_ok());

      // The adapter's proof is byte-for-byte a plain HyperKZG proof at the transformed point r',
      // with claimed value eval/scale: the raw (unmodified) engine verifies it directly.
      let (rp, scale) = coeff_eval_point(&r).unwrap();
      let eval_eval = eval * scale.invert().unwrap();
      let mut tr = Keccak256Transcript::new(b"T");
      assert!(
        hyperkzg::EvaluationEngine::verify(&vk, &mut tr, &comm, &rp, &eval_eval, &arg).is_ok()
      );

      // Wrong coefficient-form value must be rejected.
      let mut tr = Keccak256Transcript::new(b"T");
      assert!(Adapter::verify(&vk, &mut tr, &comm, &r, &(eval + Fr::ONE), &arg).is_err());
    }
  }

  // Same adapter, applied to Mercury (a different opening algorithm). Mercury requires m > 1.
  #[test]
  fn test_adapter_mercury() {
    type Adapter = MercuryCoeffAdapter<E>;
    for m in [2usize, 4, 5] {
      let n = 1 << m;
      let mut rng = rand::rngs::StdRng::seed_from_u64(100 + m as u64);
      let v = (0..n).map(|_| Fr::random(&mut rng)).collect::<Vec<_>>();
      let r = (0..m).map(|_| Fr::random(&mut rng)).collect::<Vec<_>>();

      let ck = <CE as CommitmentEngineTrait<E>>::CommitmentKey::setup_from_rng(b"test", n, OsRng);
      let (pk, vk) = Adapter::setup(&ck).unwrap();
      let comm = CE::commit(&ck, &v, &Fr::ZERO);
      let eval = coeff_eval_msb(&v, &r);

      let mut tr = Keccak256Transcript::new(b"T");
      let arg = Adapter::prove(&ck, &pk, &mut tr, &comm, &v, &r, &eval).unwrap();
      let mut tr = Keccak256Transcript::new(b"T");
      assert!(Adapter::verify(&vk, &mut tr, &comm, &r, &eval, &arg).is_ok());

      // Cross-check against the raw Mercury verifier at the transformed point.
      let (rp, scale) = coeff_eval_point(&r).unwrap();
      let eval_eval = eval * scale.invert().unwrap();
      let mut tr = Keccak256Transcript::new(b"T");
      assert!(
        mercury::EvaluationEngine::verify(&vk, &mut tr, &comm, &rp, &eval_eval, &arg).is_ok()
      );

      // Wrong value rejected.
      let mut tr = Keccak256Transcript::new(b"T");
      assert!(Adapter::verify(&vk, &mut tr, &comm, &r, &(eval + Fr::ONE), &arg).is_err());
    }
  }

  // Small hard-coded example pinning the coefficient-form semantics end-to-end.
  #[test]
  fn test_adapter_hyperkzg_hardcoded() {
    // f = 1 + 2 X0 + 3 X1 + 5 X0 X1, coefficients LSB-first: v = (1, 2, 3, 5).
    // MSB convention: index bit 1 (high) <-> point[0], index bit 0 (low) <-> point[1].
    // So X0 <-> point[1], X1 <-> point[0]; value = 1 + 2*point[1] + 3*point[0] + 5*point[0]*point[1].
    type Adapter = HyperKZGCoeffAdapter<E>;
    let n = 4;
    let v = vec![Fr::from(1), Fr::from(2), Fr::from(3), Fr::from(5)];
    let ck = <CE as CommitmentEngineTrait<E>>::CommitmentKey::setup_from_rng(b"test", n, OsRng);
    let (pk, vk) = Adapter::setup(&ck).unwrap();
    let comm = CE::commit(&ck, &v, &Fr::ZERO);

    for (point, eval) in [
      // point = [p0, p1]                        value = 1 + 2*p1 + 3*p0 + 5*p0*p1
      (vec![Fr::from(0), Fr::from(0)], Fr::from(1)), // 1
      (vec![Fr::from(1), Fr::from(1)], Fr::from(11)), // 1 + 2 + 3 + 5
      (vec![Fr::from(3), Fr::from(5)], Fr::from(95)), // 1 + 10 + 9 + 75
      (vec![Fr::from(7), Fr::from(2)], Fr::from(96)), // 1 + 4 + 21 + 70
    ] {
      let mut tr = Keccak256Transcript::new(b"T");
      let arg = Adapter::prove(&ck, &pk, &mut tr, &comm, &v, &point, &eval).unwrap();
      let mut tr = Keccak256Transcript::new(b"T");
      assert!(Adapter::verify(&vk, &mut tr, &comm, &point, &eval, &arg).is_ok());
    }

    // Wrong value at [3, 5] (correct is 95) must fail.
    let point = vec![Fr::from(3), Fr::from(5)];
    let mut tr = Keccak256Transcript::new(b"T");
    let arg = Adapter::prove(&ck, &pk, &mut tr, &comm, &v, &point, &Fr::from(96)).unwrap();
    let mut tr = Keccak256Transcript::new(b"T");
    assert!(Adapter::verify(&vk, &mut tr, &comm, &point, &Fr::from(96), &arg).is_err());
  }
}
