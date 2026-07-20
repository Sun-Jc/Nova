//! A minimal, correctness-first Projective SumCheck (PoC).
//!
//! Projective SumCheck (ePrint 2026/762) replaces the Boolean summation domain
//! `{0,1}^n` with the *infinity hypercube* `{0,∞}^n`. For a univariate round
//! polynomial `S_i(T) = a_0 + a_1 T + ... + a_{D} T^D` with declared per-round
//! degree bound `D_i`, projective evaluation reads endpoint coefficients:
//! `S_i(0) = a_0` and `S_i(∞) = [T^{D_i}] S_i = a_{D_i}`. Only the per-round
//! consistency identity changes relative to ordinary SumCheck:
//!
//! ```text
//!   C_i  ==  S_i(0) + [T^{D_i}] S_i(T)   (endpoint identity)
//!   C_{i+1} := S_i(r_i)                  (random reduction, unchanged)
//!   C_n  ==  G(r_0, ..., r_{n-1})        (final oracle check, unchanged)
//! ```
//!
//! This module implements only the SumCheck reduction (rounds + final claim
//! `C_n`). Authenticating `G(r)` against a PCS / virtual polynomial is left to
//! the caller — [`verify`] returns the sampled point and the reduced claim
//! `C_n`, and the caller compares `C_n` against its own evaluation of `G(r)`
//! (including any public structured factors such as `U(r)`).
//!
//! # Reuse of existing types
//!
//! Round messages reuse [`CompressedUniPoly`](crate::spartan::polys::univariate::CompressedUniPoly)
//! as their container and the proof reuses
//! [`SumcheckProof`](crate::spartan::sumcheck::SumcheckProof); only the
//! *compression convention* differs. A projective round omits the **constant**
//! term (`a_0 = C_i - a_D`), whereas ordinary Spartan SumCheck omits the
//! **linear** term (`a_1`, under the Boolean identity `S(0)+S(1)`). The two
//! conventions are kept apart purely by method name:
//! `UniPoly::compress_projective` / `CompressedUniPoly::decompress_projective`
//! here, versus `compress` / `decompress` in the Boolean path. Never cross
//! them.
//!
//! Round-polynomial evaluation reuses `UniPoly::evaluate` (basis-agnostic
//! Horner), so there is no bespoke univariate arithmetic in this module.
//!
//! Scope: this is a standalone proof-of-concept, additive to the existing
//! Spartan sumcheck engine.
//!
//! Module layout:
//! - [`prover`] — the Phase-1 dense `D=1` reference oracle
//!   ([`prove_dense_multilinear`]); narrow by design, not a general prover.
//! - [`verifier`] — the round-by-round reduction ([`verify`]).
//! - [`eq_sumcheck`] — Gruen-style projective equality-polynomial build
//!   ([`EqSumCheckInstanceProjective`]); eq-construction parts only.
//! - [`virtual_poly`] — factorized `Σ_t coeff_t ∏_j F_j` prover with an
//!   arbitrary-degree product round kernel ([`VirtualPolynomial`]).
//! - [`batched`] — λ-RLC batched prover folding several instances into one
//!   sumcheck execution ([`prove_batched`]).

pub mod batched;
pub mod eq_sumcheck;
pub mod prover;
pub mod verifier;
pub mod virtual_poly;

pub use batched::{prove_batched, BatchedProverOutput};
pub use eq_sumcheck::EqSumCheckInstanceProjective;
pub use prover::{prove_dense_multilinear, ProjectiveSumcheckProverOutput};
pub use verifier::{verify, ProjectiveSumcheckReduction};
pub use virtual_poly::VirtualPolynomial;

#[cfg(test)]
mod tests {
  use super::*;
  use crate::{
    provider::PallasEngine,
    spartan::{polys::univariate::UniPoly, sumcheck::SumcheckProof},
    traits::{Engine, TranscriptEngineTrait},
  };
  use ff::Field;

  type E = PallasEngine;
  type Fr = <E as Engine>::Scalar;

  #[test]
  fn projective_sumcheck_multilinear_roundtrip() {
    let num_vars = 3usize;
    // Arbitrary corner table.
    let corners: Vec<Fr> = (0..(1 << num_vars))
      .map(|i| Fr::from((7 * i + 3) as u64))
      .collect();

    // Plain prover produces the proof and its own running values.
    let mut prover_ts = <E as Engine>::TE::new(b"projsc_poc_test");
    let out = prove_dense_multilinear::<E>(&corners, num_vars, &mut prover_ts);

    let degree_bounds = vec![1usize; num_vars];
    let mut verifier_ts = <E as Engine>::TE::new(b"projsc_poc_test");
    let reduction = verify::<E>(
      out.initial_claim,
      &degree_bounds,
      &out.proof,
      &mut verifier_ts,
    )
    .unwrap();

    // The prover must satisfy the (unmodified) verifier: identical point and
    // final claim, and the claim must equal the direct sum of corners.
    assert_eq!(reduction.point, out.point);
    assert_eq!(reduction.final_claim, out.final_claim);
    assert_eq!(out.initial_claim, corners.iter().copied().sum());
  }

  #[test]
  fn projective_sumcheck_rejects_wrong_round_count() {
    let num_vars = 2usize;
    let poly = UniPoly::<Fr>::from_coeffs(vec![Fr::ZERO, Fr::ONE]).unwrap();
    let proof = SumcheckProof::<E>::new(vec![poly.compress_projective()]);
    let mut ts = <E as Engine>::TE::new(b"projsc_poc_test");
    let res = verify::<E>(Fr::ZERO, &vec![1usize; num_vars], &proof, &mut ts);
    assert!(res.is_err());
  }

  /// Cross-check that the projective compression round-trips through the
  /// dedicated methods and reconstructs `a_0 = C - a_D`.
  #[test]
  fn projective_compression_roundtrip() {
    // S(T) = 4 + 5T + 6T^2, D = 2. C = S(0) + [T^2]S = 4 + 6 = 10.
    let poly = UniPoly::<Fr>::from_coeffs(vec![Fr::from(4), Fr::from(5), Fr::from(6)]).unwrap();
    let claim = Fr::from(10);
    let restored = poly.compress_projective().decompress_projective(&claim);
    assert_eq!(restored.coeffs(), poly.coeffs());
  }

  /// Exercises `verify_final_claim` on the design §13 coefficient-wise ZeroCheck
  /// virtual polynomial `G = Eq^∞_ρ · (f·g − h·U)`, with `U` built internally
  /// for the homogenizing power. Factor openings are supplied directly (a PCS
  /// stub); the test checks that a correct `C_n` accepts and a tampered one
  /// rejects.
  #[test]
  fn verify_final_claim_zerocheck() {
    // Sampled point r = (r0, r1) and ZeroCheck challenge rho = (rho0, rho1).
    let r0 = Fr::from(5);
    let r1 = Fr::from(7);
    let (rho0, rho1) = (Fr::from(2), Fr::from(3));

    // Multilinear factors evaluated at r (f = X + 2Y, g = 3X + 4Y, h = 3X + 8Y).
    let f_r = r0 + Fr::from(2) * r1;
    let g_r = Fr::from(3) * r0 + Fr::from(4) * r1;
    let h_r = Fr::from(3) * r0 + Fr::from(8) * r1;
    // Projective equality factor Eq^∞_rho(r) = ∏_i ((1 - rho_i) + rho_i r_i).
    let eq_r = ((Fr::ONE - rho0) + rho0 * r0) * ((Fr::ONE - rho1) + rho1 * r1);
    // All-ones factor (computed here only to form the expected value).
    let u_r = (Fr::ONE + r0) * (Fr::ONE + r1);

    // Expected G(r) = Eq · (f·g − h·U).
    let expected = eq_r * (f_r * g_r - h_r * u_r);

    let reduction = ProjectiveSumcheckReduction::<E> {
      point: vec![r0, r1],
      final_claim: expected,
    };

    // evals order: [Eq, f, g, h]; U is NOT passed — it is built internally.
    let evals = vec![eq_r, f_r, g_r, h_r];
    // Terms: +1·(Eq·f·g)  and  −1·(Eq·h).  D = 3, so the second term is
    // homogenized by U^{3-2} = U internally.
    let terms = vec![(Fr::ONE, vec![0usize, 1, 2]), (-Fr::ONE, vec![0usize, 3])];

    assert!(reduction.verify_final_claim(&evals, &terms));

    // Tampered claim must be rejected.
    let bad = ProjectiveSumcheckReduction::<E> {
      point: vec![r0, r1],
      final_claim: expected + Fr::ONE,
    };
    assert!(!bad.verify_final_claim(&evals, &terms));

    // Out-of-range factor index rejects rather than panics.
    let bad_terms = vec![(Fr::ONE, vec![0usize, 99])];
    assert!(!reduction.verify_final_claim(&evals, &bad_terms));
  }
}
