//! The Plain Projective SumCheck **prover**: a dense, correctness-first
//! reference over the infinity hypercube `{0,∞}^n`.
//!
//! # Scope — this is a Phase-1 audit oracle, NOT a general prover
//!
//! [`prove_dense_multilinear`] handles exactly *one* narrow shape: the target
//! polynomial `G` handed in fully expanded as a dense table of projective
//! corner values `G[[b]]_D` (one field element per `b ∈ {0,1}^n`), with **every
//! variable multilinear (per-variable degree bound 1)**, so each round message
//! is linear `S_i(T) = a_0 + a_1 T`. It is the design doc's Phase-1 "dense
//! algebraic reference" — the obviously-correct oracle a future general prover
//! is differential-tested against (design §16.4, §18), *not* something that can
//! run real workloads.
//!
//! It deliberately **cannot** express the SumChecks Nova actually needs:
//! - **Higher degree** (Spartan outer `Az·Bz − u·Cz` is degree 3, inner is
//!   degree 2, GKR layers are degree ≥ 2): this prover hard-codes `D = 1` and
//!   emits two-coefficient round messages, with no small-polynomial convolution
//!   round kernel (design §9.3).
//! - **Factorized targets** `Σ_t λ_t ∏_j F_{t,j}` with shared witness factors
//!   and structured factors (`U`, `Eq^∞`): the input is a single flat corner
//!   table, so there is no `VirtualPolynomial`, no per-factor bind state
//!   (design §8.2, §9.1), and no structured-factor machinery (design §5.2,
//!   §7.1).
//! - **GKR**: no layer recursion, no `λ`-power batching, no cross-layer root
//!   binding.
//!
//! A general prover (arbitrary `D`, factorized virtual polynomials, structured
//! factors) is a separate implementation built on those pieces; do not extend
//! this file toward it.
//!
//! The prover is written to satisfy [`super::verifier::verify`] *exactly*: it
//! uses the same transcript labels and the same absorb-then-squeeze order, and
//! emits [`CompressedUniPoly`] round messages compressed with
//! [`UniPoly::compress_projective`]. The verifier is not modified.

use crate::{
  spartan::{
    polys::univariate::{CompressedUniPoly, UniPoly},
    sumcheck::SumcheckProof,
  },
  traits::{Engine, TranscriptEngineTrait},
};
use ff::Field;

/// Output of the plain prover: the proof plus the values a caller (or a
/// differential test) needs to cross-check against the verifier's reduction.
#[derive(Clone, Debug)]
pub struct ProjectiveSumcheckProverOutput<E: Engine> {
  /// The round messages, ready to feed to [`super::verifier::verify`].
  pub proof: SumcheckProof<E>,
  /// The claimed projective sum `C_0 = Σ_b G[[b]]_D` (the verifier's
  /// `initial_claim`).
  pub initial_claim: E::Scalar,
  /// The sampled point `r = (r_0, ..., r_{n-1})`.
  pub point: Vec<E::Scalar>,
  /// The final reduced value `C_n = G(r)`.
  pub final_claim: E::Scalar,
}

/// Runs the Phase-1 dense multilinear Projective SumCheck prover on a corner
/// table.
///
/// **Narrow by design** — this is the `D = 1`, single-flat-table audit oracle,
/// not a general prover; see the module docs for what it cannot do (higher
/// degree, factorized targets, GKR).
///
/// `corners[b]` holds the projective corner value `G[[b]]_D`, indexed by the
/// bit vector `b` with the low bit corresponding to `X_0` (the first variable
/// eliminated). Its length must be `2^num_vars`.
///
/// The prover mirrors the verifier's transcript discipline: for each round it
/// forms the linear round polynomial `S_i`, absorbs it, squeezes the challenge
/// `r_i`, then binds `X_i = r_i` in the monomial basis
/// (`new = a_0 + r_i · a_1`). After `n` rounds the table collapses to the
/// single value `G(r)`.
///
/// # Panics
/// Panics if `corners.len() != 2^num_vars`.
pub fn prove_dense_multilinear<E: Engine>(
  corners: &[E::Scalar],
  num_vars: usize,
  transcript: &mut E::TE,
) -> ProjectiveSumcheckProverOutput<E> {
  assert_eq!(
    corners.len(),
    1 << num_vars,
    "corner table length must be 2^num_vars"
  );

  // C_0 = Σ_b G[[b]]_D, the claimed projective sum.
  let initial_claim: E::Scalar = corners.iter().copied().sum();

  // Working table over surviving corners; the low bit is the current variable
  // (index `suffix` = X_i^0, index `suffix + half` = X_i^∞).
  let mut table = corners.to_vec();
  let mut rounds: Vec<CompressedUniPoly<E::Scalar>> = Vec::with_capacity(num_vars);
  let mut point = Vec::with_capacity(num_vars);

  let mut remaining = num_vars;
  for _ in 0..num_vars {
    let half = 1 << (remaining - 1);

    // Round polynomial: sum the linear slices over all suffixes.
    //   S_i(T) = Σ_suffix (table[0,suffix] + table[1,suffix] · T)
    // At the projective endpoints this is a_0 = S_i(0) and a_1 = [T] S_i.
    let mut a0 = E::Scalar::ZERO;
    let mut a1 = E::Scalar::ZERO;
    for suffix in 0..half {
      a0 += table[suffix];
      a1 += table[suffix + half];
    }

    let poly = UniPoly::<E::Scalar>::from_coeffs_no_trim(vec![a0, a1])
      .expect("two coefficients is a valid univariate polynomial");

    // Absorb the round polynomial (all D+1 coeffs, incl. the linear term),
    // matching the verifier, then squeeze the challenge (never before).
    super::absorb_round::<E>(transcript, &poly);
    let r_i = transcript
      .squeeze(b"projective_sumcheck_challenge")
      .expect("transcript squeeze failed");

    rounds.push(poly.compress_projective());
    point.push(r_i);

    // Bind X_i = r_i in the monomial basis: new = a_0 + r_i · a_1.
    let mut next = vec![E::Scalar::ZERO; half];
    for suffix in 0..half {
      next[suffix] = table[suffix] + r_i * table[suffix + half];
    }
    table = next;
    remaining -= 1;
  }

  // After n rounds the table holds the single reduced value G(r) = C_n.
  let final_claim = table[0];

  ProjectiveSumcheckProverOutput {
    proof: SumcheckProof::new(rounds),
    initial_claim,
    point,
    final_claim,
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::{provider::PallasEngine, spartan::projective_sumcheck::verify};
  use ff::Field;

  type E = PallasEngine;
  type Fr = <E as Engine>::Scalar;

  /// Regression for audit F2: a round whose polynomial is a nonzero constant
  /// (`a_0 != 0`, `a_1 = 0`) previously trimmed to length 1 and panicked in
  /// `compress_projective`. Corner table `[5, 0]` (n=1) makes the single round
  /// `S(T) = 5 + 0·T`.
  #[test]
  fn dense_nonzero_constant_round_does_not_panic() {
    let corners = vec![Fr::from(5), Fr::ZERO];
    let mut ts = <E as Engine>::TE::new(b"projsc_f2");
    let out = prove_dense_multilinear::<E>(&corners, 1, &mut ts);

    // The single round message must carry exactly D = 1 stored coefficient.
    assert_eq!(out.proof.compressed_polys()[0].stored_coeffs().len(), 1);

    let mut ts_ver = <E as Engine>::TE::new(b"projsc_f2");
    let reduction = verify::<E>(out.initial_claim, &[1usize], &out.proof, &mut ts_ver).unwrap();
    assert_eq!(reduction.point, out.point);
    assert_eq!(reduction.final_claim, out.final_claim);
    assert_eq!(out.initial_claim, Fr::from(5));
  }
}
