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

/// Builds the projective **inner ABC** relation `L_row · L_col · val` over
/// `num_vars` variables (Form B — a plain degree-3 product, no eq factor).
///
/// `val` is the caller's folded matrix-value column `val_A + c·val_B + c²·val_C`
/// (public `c`). The projective sum equals `Σ_b L_row[b]·L_col[b]·val[b]`, the
/// eval-basis inner ABC summand.
pub fn build_inner_abc<E: Engine>(
  num_vars: usize,
  l_row: Vec<E::Scalar>,
  l_col: Vec<E::Scalar>,
  val: Vec<E::Scalar>,
) -> VirtualPolynomial<E> {
  let factors = vec![l_row, l_col, val];
  let terms = vec![(E::Scalar::ONE, vec![0usize, 1, 2])];
  VirtualPolynomial::new(num_vars, factors, terms)
}

/// Builds the projective **inner E** relation `Eq^∞_{r_outer} · E` over
/// `num_vars` variables (Form A — rerandomizes E's opening point). The
/// projective sum equals `Σ_b eq̃(r_outer,b)·E[b]`.
pub fn build_inner_e<E: Engine>(
  num_vars: usize,
  e: Vec<E::Scalar>,
  r_outer: &[E::Scalar],
) -> VirtualPolynomial<E> {
  assert_eq!(
    r_outer.len(),
    num_vars,
    "r_outer must have num_vars entries"
  );
  let eq = EqPolynomialProjective::<E::Scalar>::new(r_outer.to_vec()).evals();
  let factors = vec![eq, e];
  let terms = vec![(E::Scalar::ONE, vec![0usize, 1])];
  VirtualPolynomial::new(num_vars, factors, terms)
}

/// Builds the projective **witness-bound** relation `maskedEq^∞_τ · W` over
/// `num_vars` variables (Form D). Certifies the padded tail of `W` is zero:
/// `0 = Σ_{2^m ≤ b} eq̃(τ,b)·W[b]`. `num_masked_vars = m` zeroes the first `2^m`
/// corners of the eq factor.
pub fn build_witness_bound<E: Engine>(
  num_vars: usize,
  w: Vec<E::Scalar>,
  tau: &[E::Scalar],
  num_masked_vars: usize,
) -> VirtualPolynomial<E> {
  assert_eq!(tau.len(), num_vars, "tau must have num_vars entries");
  let masked_eq =
    EqPolynomialProjective::<E::Scalar>::new(tau.to_vec()).masked_evals(num_masked_vars);
  let factors = vec![masked_eq, w];
  let terms = vec![(E::Scalar::ONE, vec![0usize, 1])];
  VirtualPolynomial::new(num_vars, factors, terms)
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

  /// Helper: verify a projective proof at uniform degree D and check the
  /// reduction matches the prover.
  fn check_verifies(
    out: &crate::spartan::projective_sumcheck::ProjectiveSumcheckProverOutput<E>,
    num_vars: usize,
    d: usize,
    label: &'static [u8],
  ) {
    let degree_bounds = vec![d; num_vars];
    let mut ts = <E as Engine>::TE::new(label);
    let reduction = verify::<E>(out.initial_claim, &degree_bounds, &out.proof, &mut ts).unwrap();
    assert_eq!(reduction.point, out.point);
    assert_eq!(reduction.final_claim, out.final_claim);
  }

  /// Inner ABC: projective sum equals Σ_b L_row[b]·L_col[b]·val[b], proof
  /// verifies at D=3.
  #[test]
  fn inner_abc_reduces_and_verifies() {
    let num_vars = 3usize;
    let n = 1usize << num_vars;
    let l_row: Vec<Fr> = (0..n).map(|i| Fr::from((i + 1) as u64)).collect();
    let l_col: Vec<Fr> = (0..n).map(|i| Fr::from((2 * i + 3) as u64)).collect();
    let val: Vec<Fr> = (0..n).map(|i| Fr::from((5 * i + 7) as u64)).collect();
    let expected: Fr = (0..n).map(|i| l_row[i] * l_col[i] * val[i]).sum();

    let vp = build_inner_abc::<E>(num_vars, l_row, l_col, val);
    let mut ts = <E as Engine>::TE::new(b"projsc_abc");
    let out = vp.prove(&mut ts);
    assert_eq!(out.initial_claim, expected);
    check_verifies(&out, num_vars, 3, b"projsc_abc");
  }

  /// Inner E: projective sum equals Σ_b eq̃(r_outer,b)·E[b], proof verifies
  /// at D=2.
  #[test]
  fn inner_e_reduces_and_verifies() {
    let num_vars = 3usize;
    let n = 1usize << num_vars;
    let e: Vec<Fr> = (0..n).map(|i| Fr::from((3 * i + 2) as u64)).collect();
    let r_outer: Vec<Fr> = (0..num_vars).map(|i| Fr::from((i + 4) as u64)).collect();

    // Expected = Σ_b eq̃(r_outer,b)·E[b] = coeff-basis: Eq^∞ corner · E corner.
    let eq = EqPolynomialProjective::<Fr>::new(r_outer.clone()).evals();
    let expected: Fr = (0..n).map(|i| eq[i] * e[i]).sum();

    let vp = build_inner_e::<E>(num_vars, e, &r_outer);
    assert_eq!(vp.degree(), 2);
    let mut ts = <E as Engine>::TE::new(b"projsc_e");
    let out = vp.prove(&mut ts);
    assert_eq!(out.initial_claim, expected);
    check_verifies(&out, num_vars, 2, b"projsc_e");
  }

  /// Witness-bound: with W zero on the unmasked tail, the masked-eq-weighted
  /// projective sum is zero; a nonzero tail entry makes it nonzero.
  #[test]
  fn witness_bound_certifies_zero_tail() {
    let num_vars = 3usize;
    let n = 1usize << num_vars;
    let num_masked_vars = 1usize; // first 2^1 = 2 corners masked out
    let tau: Vec<Fr> = (0..num_vars).map(|i| Fr::from((i + 2) as u64)).collect();

    // W zero everywhere on the unmasked tail (indices >= 2) → sum 0.
    let w_zero_tail: Vec<Fr> = (0..n)
      .map(|i| {
        if i < 2 {
          Fr::from((i + 9) as u64)
        } else {
          Fr::ZERO
        }
      })
      .collect();
    let vp = build_witness_bound::<E>(num_vars, w_zero_tail, &tau, num_masked_vars);
    let mut ts = <E as Engine>::TE::new(b"projsc_wb");
    let out = vp.prove(&mut ts);
    assert_eq!(out.initial_claim, Fr::ZERO);
    check_verifies(&out, num_vars, 2, b"projsc_wb");

    // Nonzero tail entry → nonzero sum.
    let mut w_bad = vec![Fr::ZERO; n];
    w_bad[3] = Fr::ONE;
    let vp2 = build_witness_bound::<E>(num_vars, w_bad, &tau, num_masked_vars);
    let mut ts2 = <E as Engine>::TE::new(b"projsc_wb2");
    let out2 = vp2.prove(&mut ts2);
    assert_ne!(out2.initial_claim, Fr::ZERO);
  }
}
