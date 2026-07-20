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

/// One side (row or col) of the projective **memory / logup** check, as three
/// [`VirtualPolynomial`] sub-instances mirroring ppSNARK's `MemorySumcheckInstance`:
///
/// - **inv** (degree 1, no eq): `t_inv − w_inv`, where `t_inv[i] = TS[i]/(T[i]+r)`
///   and `w_inv[i] = 1/(W[i]+r)`. Projective sum = `Σ (t_inv − w_inv)`.
/// - **T** (degree 3): `Eq^∞_ρ · (t_inv·(T+r) − TS)`. Certifies `t_inv` is the
///   correct fingerprint inverse times read-count.
/// - **W** (degree 3): `Eq^∞_ρ · (w_inv·(W+r) − 1)`. Certifies `w_inv` is the
///   correct fingerprint inverse.
///
/// The auxiliary tables `t_plus_r = T+r`, `w_plus_r = W+r`, `t_inv`, `w_inv`,
/// `ts` are precomputed pointwise by the caller (exactly as
/// `R1CSShapeSparkRepr::compute_oracles`), so each is a single degree-1 factor.
/// `rho` is the memory ZeroCheck challenge (`|rho| = num_vars`). Returns
/// `(inv, t_relation, w_relation)`.
#[allow(clippy::too_many_arguments)]
pub fn build_memory_side<E: Engine>(
  num_vars: usize,
  t_inv: Vec<E::Scalar>,
  w_inv: Vec<E::Scalar>,
  t_plus_r: Vec<E::Scalar>,
  w_plus_r: Vec<E::Scalar>,
  ts: Vec<E::Scalar>,
  rho: &[E::Scalar],
) -> (
  VirtualPolynomial<E>,
  VirtualPolynomial<E>,
  VirtualPolynomial<E>,
) {
  assert_eq!(rho.len(), num_vars, "rho must have num_vars entries");
  let eq = EqPolynomialProjective::<E::Scalar>::new(rho.to_vec()).evals();

  // inv: t_inv − w_inv  (degree 1, no eq).
  let inv = VirtualPolynomial::new(
    num_vars,
    vec![t_inv.clone(), w_inv.clone()],
    vec![(E::Scalar::ONE, vec![0]), (-E::Scalar::ONE, vec![1])],
  );

  // T: Eq·(t_inv·(T+r) − TS).  Factors [Eq, t_inv, t_plus_r, ts].
  //   +Eq·t_inv·t_plus_r  (deg 3)
  //   −Eq·ts              (deg 2 → homogenized by U)
  let t_relation = VirtualPolynomial::new_homogenized(
    num_vars,
    vec![eq.clone(), t_inv, t_plus_r, ts],
    vec![
      (E::Scalar::ONE, vec![0, 1, 2]),
      (-E::Scalar::ONE, vec![0, 3]),
    ],
  );

  // W: Eq·(w_inv·(W+r) − 1).  Factors [Eq, w_inv, w_plus_r].
  //   +Eq·w_inv·w_plus_r  (deg 3)
  //   −Eq·1               (deg 1 → homogenized by U^2, coeff −1 · const 1)
  // The constant term uses the all-ones factor implicitly via homogenization:
  // represent "1" as an explicit all-ones factor table so the term is Eq·ones.
  let n = 1usize << num_vars;
  let ones = vec![E::Scalar::ONE; n];
  let w_relation = VirtualPolynomial::new_homogenized(
    num_vars,
    vec![eq, w_inv, w_plus_r, ones],
    vec![
      (E::Scalar::ONE, vec![0, 1, 2]),
      (-E::Scalar::ONE, vec![0, 3]),
    ],
  );

  (inv, t_relation, w_relation)
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

  /// Memory side with honest inverses: the T and W relations reduce to zero
  /// (t_inv·(T+r) = TS and w_inv·(W+r) = 1 pointwise), and each verifies.
  #[test]
  fn memory_side_honest_inverses_reduce_to_zero() {
    let num_vars = 3usize;
    let n = 1usize << num_vars;

    // Synthetic fingerprints (already include +r): pick nonzero T+r, W+r.
    let t_plus_r: Vec<Fr> = (0..n).map(|i| Fr::from((3 * i + 11) as u64)).collect();
    let w_plus_r: Vec<Fr> = (0..n).map(|i| Fr::from((5 * i + 13) as u64)).collect();
    let ts: Vec<Fr> = (0..n).map(|i| Fr::from((i % 4 + 1) as u64)).collect();

    // Honest inverses: t_inv = TS/(T+r), w_inv = 1/(W+r).
    let t_inv: Vec<Fr> = (0..n)
      .map(|i| ts[i] * t_plus_r[i].invert().unwrap())
      .collect();
    let w_inv: Vec<Fr> = (0..n).map(|i| w_plus_r[i].invert().unwrap()).collect();

    let rho: Vec<Fr> = (0..num_vars).map(|i| Fr::from((i + 6) as u64)).collect();

    let (inv, t_rel, w_rel) =
      build_memory_side::<E>(num_vars, t_inv, w_inv, t_plus_r, w_plus_r, ts, &rho);

    // T and W relations vanish on every corner ⇒ projective sum 0.
    let mut ts_t = <E as Engine>::TE::new(b"projsc_mem_t");
    let out_t = t_rel.prove(&mut ts_t);
    assert_eq!(out_t.initial_claim, Fr::ZERO);
    check_verifies(&out_t, num_vars, 3, b"projsc_mem_t");

    let mut ts_w = <E as Engine>::TE::new(b"projsc_mem_w");
    let out_w = w_rel.prove(&mut ts_w);
    assert_eq!(out_w.initial_claim, Fr::ZERO);
    check_verifies(&out_w, num_vars, 3, b"projsc_mem_w");

    // The inv sub-claim is the logup balance Σ(t_inv − w_inv); it verifies at
    // degree 1 (its value is whatever the multiset balance is, not asserted 0
    // here since these are synthetic fingerprints).
    let mut ts_i = <E as Engine>::TE::new(b"projsc_mem_inv");
    let out_i = inv.prove(&mut ts_i);
    check_verifies(&out_i, num_vars, 1, b"projsc_mem_inv");
  }

  /// A dishonest T-inverse makes the T relation nonzero.
  #[test]
  fn memory_side_dishonest_inverse_is_nonzero() {
    let num_vars = 2usize;
    let n = 1usize << num_vars;
    let t_plus_r: Vec<Fr> = (0..n).map(|i| Fr::from((3 * i + 11) as u64)).collect();
    let w_plus_r: Vec<Fr> = (0..n).map(|i| Fr::from((5 * i + 13) as u64)).collect();
    let ts: Vec<Fr> = (0..n).map(|i| Fr::from((i + 1) as u64)).collect();
    let mut t_inv: Vec<Fr> = (0..n)
      .map(|i| ts[i] * t_plus_r[i].invert().unwrap())
      .collect();
    t_inv[0] += Fr::ONE; // corrupt one entry
    let w_inv: Vec<Fr> = (0..n).map(|i| w_plus_r[i].invert().unwrap()).collect();
    let rho: Vec<Fr> = (0..num_vars).map(|i| Fr::from((i + 6) as u64)).collect();

    let (_, t_rel, _) =
      build_memory_side::<E>(num_vars, t_inv, w_inv, t_plus_r, w_plus_r, ts, &rho);
    let mut ts_t = <E as Engine>::TE::new(b"projsc_mem_bad");
    let out_t = t_rel.prove(&mut ts_t);
    assert_ne!(out_t.initial_claim, Fr::ZERO);
  }
}
