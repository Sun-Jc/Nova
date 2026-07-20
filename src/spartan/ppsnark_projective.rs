//! Projective ppSNARK (Stage 4, in progress).
//!
//! This module is being assembled incrementally. The original eval-basis
//! [`ppsnark`](crate::spartan::ppsnark) is kept intact; this file builds the
//! coefficient-basis (projective) counterpart of its sumchecks so that a full
//! `RelaxedR1CSSNARK` can be composed once every relation is ported.
//!
//! So far it provides the **outer relation** builder: the coefficient-basis
//! ZeroCheck
//! ```text
//!   G_outer(X) = Eq^∞_τ(X) · ( Az(X)·Bz(X) − u·Cz(X)·U(X) − E(X)·U(X) ),
//! ```
//! whose projective sum over `{0,∞}^n` is `Σ_b eq̃(τ,b)·(Az[b]Bz[b] − u·Cz[b] −
//! E[b])` — identical to the eval-basis outer sumcheck's summand (Stage-1
//! relation ledger, §2.1), because the pointwise relation is basis-neutral.
//!
//! The builder returns a [`VirtualPolynomial`] the projective prover/verifier
//! already handle; it does not touch the immutable projective sumcheck verifier.

use crate::{
  spartan::{
    polys::eq_projective::EqPolynomialProjective,
    projective_sumcheck::virtual_poly::VirtualPolynomial,
  },
  traits::Engine,
};
use ff::Field;

/// Builds the projective outer-relation virtual polynomial
/// `Eq^∞_τ · (Az·Bz − u·Cz·U − E·U)` over `num_vars` variables.
///
/// Inputs are the coefficient-basis tables (length `2^num_vars`) of `Az`, `Bz`,
/// `Cz`, `E`; `u` is the relaxed-R1CS scalar; `tau` is the ZeroCheck challenge
/// vector (`|tau| = num_vars`). The all-ones factor `U` used to homogenize the
/// degree-1 `Cz` and `E` terms up to the degree-3 `Az·Bz·Eq^∞` term is added
/// internally by [`VirtualPolynomial::new_homogenized`].
///
/// The resulting virtual polynomial has degree `D = 3` and, when the R1CS
/// instance is satisfied, projective sum `0`.
pub fn build_outer<E: Engine>(
  num_vars: usize,
  az: Vec<E::Scalar>,
  bz: Vec<E::Scalar>,
  cz: Vec<E::Scalar>,
  e: Vec<E::Scalar>,
  u: E::Scalar,
  tau: &[E::Scalar],
) -> VirtualPolynomial<E> {
  assert_eq!(tau.len(), num_vars, "tau must have num_vars entries");
  let eq = EqPolynomialProjective::<E::Scalar>::new(tau.to_vec()).evals();

  // Scale Cz by u once (folded into the factor table). E stays as-is.
  let u_cz: Vec<E::Scalar> = cz.into_iter().map(|c| u * c).collect();

  // Factors: [Eq, Az, Bz, uCz, E].
  //   term +1 · Eq·Az·Bz          (degree 3)
  //   term −1 · Eq·uCz            (degree 2 → homogenized by U)
  //   term −1 · Eq·E              (degree 2 → homogenized by U)
  let factors = vec![eq, az, bz, u_cz, e];
  let terms = vec![
    (E::Scalar::ONE, vec![0usize, 1, 2]),
    (-E::Scalar::ONE, vec![0usize, 3]),
    (-E::Scalar::ONE, vec![0usize, 4]),
  ];

  VirtualPolynomial::new_homogenized(num_vars, factors, terms)
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::{
    provider::PallasEngine,
    spartan::projective_sumcheck::{verify, ProjectiveSumcheckReduction},
    traits::{Engine, TranscriptEngineTrait},
  };
  use ff::Field;

  type E = PallasEngine;
  type Fr = <E as Engine>::Scalar;

  /// End-to-end: a *satisfied* outer relation (Az∘Bz = u·Cz + E pointwise) has
  /// projective sum 0, the projective proof verifies, and verify_final_claim
  /// reconstructs G(r) from the factor openings.
  #[test]
  fn outer_relation_satisfied_reduces_to_zero() {
    let num_vars = 3usize;
    let n = 1usize << num_vars;
    let u = Fr::from(7);

    // Build a satisfying instance: pick Az, Bz, Cz freely, set E = Az∘Bz − u·Cz.
    let az: Vec<Fr> = (0..n).map(|i| Fr::from((2 * i + 1) as u64)).collect();
    let bz: Vec<Fr> = (0..n).map(|i| Fr::from((3 * i + 4) as u64)).collect();
    let cz: Vec<Fr> = (0..n).map(|i| Fr::from((i + 2) as u64)).collect();
    let e: Vec<Fr> = (0..n).map(|i| az[i] * bz[i] - u * cz[i]).collect();

    let tau: Vec<Fr> = (0..num_vars)
      .map(|i| Fr::from((5 * i + 3) as u64))
      .collect();

    let vp = build_outer::<E>(
      num_vars,
      az.clone(),
      bz.clone(),
      cz.clone(),
      e.clone(),
      u,
      &tau,
    );
    assert_eq!(vp.degree(), 3);

    let mut ts = <E as Engine>::TE::new(b"projsc_outer");
    let out = vp.prove(&mut ts);
    // Satisfied ⇒ residual is zero at every corner ⇒ projective sum 0.
    assert_eq!(out.initial_claim, Fr::ZERO);

    // Verify.
    let degree_bounds = vec![3usize; num_vars];
    let mut ts_ver = <E as Engine>::TE::new(b"projsc_outer");
    let reduction = verify::<E>(Fr::ZERO, &degree_bounds, &out.proof, &mut ts_ver).unwrap();
    assert_eq!(reduction.point, out.point);
    assert_eq!(reduction.final_claim, out.final_claim);

    // Final-oracle check from factor openings [Eq, Az, Bz, uCz, E].
    let r = &reduction.point;
    let eval = |t: &[Fr]| -> Fr {
      (0..n)
        .map(|b| {
          let mut acc = t[b];
          for (k, rk) in r.iter().enumerate() {
            let bit = num_vars - 1 - k;
            if (b >> bit) & 1 == 1 {
              acc *= *rk;
            }
          }
          acc
        })
        .sum()
    };
    let eq = EqPolynomialProjective::<Fr>::new(tau.clone()).evals();
    let u_cz: Vec<Fr> = cz.iter().map(|c| u * *c).collect();
    let evals = vec![eval(&eq), eval(&az), eval(&bz), eval(&u_cz), eval(&e)];
    let terms = vec![
      (Fr::ONE, vec![0usize, 1, 2]),
      (-Fr::ONE, vec![0usize, 3]),
      (-Fr::ONE, vec![0usize, 4]),
    ];
    let red = ProjectiveSumcheckReduction::<E> {
      point: reduction.point.clone(),
      final_claim: reduction.final_claim,
    };
    assert!(red.verify_final_claim(&evals, &terms));
  }

  /// An *unsatisfied* instance has nonzero projective sum (the ZeroCheck would
  /// then catch it via a nonzero initial claim).
  #[test]
  fn outer_relation_unsatisfied_is_nonzero() {
    let num_vars = 2usize;
    let n = 1usize << num_vars;
    let u = Fr::from(3);
    let az: Vec<Fr> = (0..n).map(|i| Fr::from((i + 1) as u64)).collect();
    let bz: Vec<Fr> = (0..n).map(|i| Fr::from((2 * i + 1) as u64)).collect();
    let cz: Vec<Fr> = (0..n).map(|i| Fr::from((i + 5) as u64)).collect();
    // E deliberately wrong at one corner.
    let mut e: Vec<Fr> = (0..n).map(|i| az[i] * bz[i] - u * cz[i]).collect();
    e[1] += Fr::ONE;
    let tau: Vec<Fr> = (0..num_vars).map(|i| Fr::from((i + 2) as u64)).collect();

    let vp = build_outer::<E>(num_vars, az, bz, cz, e, u, &tau);
    let mut ts = <E as Engine>::TE::new(b"projsc_outer_bad");
    let out = vp.prove(&mut ts);
    assert_ne!(out.initial_claim, Fr::ZERO);
  }
}
