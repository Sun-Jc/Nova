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
//! The builder returns a `VirtualPolynomial` the projective prover/verifier
//! already handle; it does not touch the immutable projective sumcheck verifier.

use crate::{
  provider::coeff_eval_adapter::coeff_eval_point,
  spartan::{
    polys::{
      eq::EqPolynomial, eq_projective::EqPolynomialProjective, identity::IdentityPolynomial,
      masked_eq::MaskedEqPolynomial, multilinear::MultilinearPolynomial,
    },
    projective_sumcheck::{
      eq_factored::EqFactoredVirtualPolynomial, virtual_poly::VirtualPolynomial,
    },
  },
  traits::Engine,
};
use ff::Field;

// ---------------------------------------------------------------------------
// Coefficient-form finite-point evaluation helpers (verify-side reconstruction)
//
// The projective sumcheck reduces to a coefficient-form point `r`. Every value
// the verifier reconstructs is a coefficient-form MLE evaluation `coeffMLE(v, r)
// = <v, ⊗ (1, r_j)>`. Reusing the point-transform adapter, this equals
// `scale · evalMLE(v, r')` with `(r', scale) = coeff_eval_point(r)` — so we
// reuse the existing (eval-basis) MLE evaluator rather than writing a new one.
// This is the concrete "same repr, different interpretation" reuse.
// ---------------------------------------------------------------------------

/// `coeffMLE(v, r) = <v, ⊗_j (1, r_j)>`, computed by delegating to the existing
/// evaluation-basis evaluator via the point transform (adapter §5). MSB-first,
/// matching the projective provers' high-bit-first binding.
///
/// # Panics
/// Panics if the transform is undefined (`1 + r_j = 0` for some `j`) — under a
/// Fiat-Shamir point this is negligible.
pub fn coeff_eval<E: Engine>(v: &[E::Scalar], r: &[E::Scalar]) -> E::Scalar {
  let (r_prime, scale) = coeff_eval_point(r).expect("coeff transform undefined (1 + r_j == 0)");
  scale * MultilinearPolynomial::evaluate_with(v, &r_prime)
}

/// Coefficient-form evaluation of the projective equality polynomial
/// `Eq^∞_ρ` at a finite point `r`: `∏_j ((1 − ρ_j) + ρ_j r_j)`. This is the
/// coeff-form analogue of `EqPolynomial::evaluate` used in the eval-basis verify.
pub fn coeff_eq_eval<E: Engine>(rho: &[E::Scalar], r: &[E::Scalar]) -> E::Scalar {
  EqPolynomialProjective::<E::Scalar>::new(rho.to_vec()).evaluate(r)
}

/// Coefficient-form evaluation of the identity factor (whose corner
/// coefficients are `[0, 1, …, N−1]`) at a finite point `r`: `coeffMLE([0..N),
/// r)`, computed **succinctly in `O(m)`** via the point transform.
///
/// `coeffMLE([0..N), r) = scale · evalMLE([0..N), r')`, and `[0..N)` read as an
/// evaluation table is exactly the identity polynomial, so
/// `evalMLE([0..N), r') = IdentityPolynomial::evaluate(r')` (`Σ 2^i r'_i`, O(m)).
pub fn coeff_identity_eval<E: Engine>(num_vars: usize, r: &[E::Scalar]) -> E::Scalar {
  let (r_prime, scale) = coeff_eval_point(r).expect("coeff transform undefined (1 + r_j == 0)");
  scale * IdentityPolynomial::<E::Scalar>::new(num_vars).evaluate(&r_prime)
}

/// Coefficient-form evaluation of the masked projective equality polynomial
/// (first `2^num_masked_vars` corners zeroed) at `r`, computed **succinctly in
/// `O(m)`** via the point transform.
///
/// The masked projective eq's corner coefficients equal the Boolean masked eq
/// weights (the projective/Boolean coincidence at corners), so
/// `coeffMLE(masked_proj, r) = scale · MaskedEqPolynomial::evaluate(r')`.
pub fn coeff_masked_eq_eval<E: Engine>(
  rho: &[E::Scalar],
  num_masked_vars: usize,
  r: &[E::Scalar],
) -> E::Scalar {
  let (r_prime, scale) = coeff_eval_point(r).expect("coeff transform undefined (1 + r_j == 0)");
  let eq = EqPolynomial::<E::Scalar>::new(rho.to_vec());
  scale * MaskedEqPolynomial::new(&eq, num_masked_vars).evaluate(&r_prime)
}

/// `coeffMLE(⊗(1, a), r)` — the coeff-MLE of the monomial tensor table
/// `t[j] = ∏_{i∈j} a_i` (= [`coeff_tensor`]) evaluated at `r` — computed
/// **succinctly in `O(m)`** as `∏_c (1 + a_c · r_c)`. Both `a` and `r` are
/// length `m`, MSB-first. Used to reconstruct `mem_row(r_inner)` (mem_row is the
/// monomial tensor over `r_outer`) without materializing the `2^m` table.
pub fn coeff_tensor_cross<E: Engine>(a: &[E::Scalar], r: &[E::Scalar]) -> E::Scalar {
  debug_assert_eq!(a.len(), r.len());
  a.iter()
    .zip(r.iter())
    .fold(E::Scalar::ONE, |acc, (ac, rc)| {
      acc * (E::Scalar::ONE + *ac * *rc)
    })
}

/// The **monomial tensor** corner table `t[j] = ∏_{i∈j} r_i` (= ⊗_i (1, r_i)),
/// MSB-first (bit `num_vars−1−c` ↔ `r[c]`). Built by tensor doubling in `O(N)`
/// (not `O(N·m)`): each variable doubles the table, the high half scaled by
/// `r_c`. This is the coefficient-basis counterpart of the eval-basis
/// `mem_row = EqPolynomial(r).evals()`.
pub fn coeff_tensor<F: Field>(r: &[F]) -> Vec<F> {
  let num_vars = r.len();
  let n = 1usize << num_vars;
  let mut table = vec![F::ZERO; n];
  table[0] = F::ONE;
  // MSB-first: r[0] is the highest bit, so process from the top so that after
  // step k the first 2^{k} entries hold the tensor over the top k variables.
  let mut size = 1usize;
  for &rc in r.iter().rev() {
    // Low-bit-first doubling; r is reversed so the final layout is MSB-first.
    let (lo, hi) = table.split_at_mut(size);
    hi[..size]
      .iter_mut()
      .zip(lo.iter())
      .for_each(|(h, l)| *h = *l * rc);
    size <<= 1;
  }
  table
}

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

/// Eq-factored variant of [`build_outer`]: returns an
/// [`EqFactoredVirtualPolynomial`] that carries `Eq^∞_τ` analytically (Gruen)
/// instead of as a dense factor. The R-part is `Az·Bz − u·Cz·U − E·U`, whose
/// terms are U-homogenized to degree `Dr = 2` internally (an explicit all-ones
/// factor lifts the degree-1 `uCz` and `E` terms). The full message degree is
/// `D = Dr + 1 = 3`, matching [`build_outer`].
pub fn build_outer_eq_factored<E: Engine>(
  num_vars: usize,
  az: Vec<E::Scalar>,
  bz: Vec<E::Scalar>,
  cz: Vec<E::Scalar>,
  e: Vec<E::Scalar>,
  u: E::Scalar,
  tau: &[E::Scalar],
) -> EqFactoredVirtualPolynomial<E> {
  assert_eq!(tau.len(), num_vars, "tau must have num_vars entries");
  // Fuse u·Cz + E into one column (as the eval-basis outer prover does), so R
  // has 2 terms / 3 factor tables — and the shorter term is homogenized by an
  // analytic U (no dense all-ones 2^n table to store or bind).
  let u_cz_e: Vec<E::Scalar> = cz.into_iter().zip(e).map(|(c, ei)| u * c + ei).collect();

  // R factors: [Az, Bz, uCzE]. R terms (Dr = 2): +Az·Bz, −uCzE·U (the second
  // term carries one analytic U copy).
  let factors = vec![az, bz, u_cz_e];
  let terms = vec![
    (E::Scalar::ONE, vec![0usize, 1]),
    (-E::Scalar::ONE, vec![2usize]),
  ];
  EqFactoredVirtualPolynomial::new_mixed(num_vars, tau.to_vec(), factors, terms)
}
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

// ===========================================================================
// Full projective ppSNARK: RelaxedR1CSSNARKProjective
// ===========================================================================

use crate::{
  errors::NovaError,
  r1cs::{R1CSShape, RelaxedR1CSInstance, RelaxedR1CSWitness},
  spartan::{
    batch_invert,
    math::Math,
    powers,
    ppsnark::{ProverKey, VerifierKey},
    projective_sumcheck::{
      batched::{prove_batched_mixed, ProjInstance},
      verify as sc_verify, ProjectiveSumcheckReduction,
    },
    PolyEvalInstance, PolyEvalWitness,
  },
  traits::{
    commitment::CommitmentEngineTrait,
    evaluation::EvaluationEngineTrait,
    snark::{DigestHelperTrait, RelaxedR1CSSNARKTrait},
    TranscriptEngineTrait,
  },
  Commitment, CommitmentKey,
};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

/// Zero-pads `v` to length `n`.
fn padded<E: Engine>(v: &[E::Scalar], n: usize) -> Vec<E::Scalar> {
  let mut out = vec![E::Scalar::ZERO; n];
  out[..v.len()].copy_from_slice(v);
  out
}

/// Coefficient-basis evaluation oracles: like `R1CSShapeSparkRepr::evaluation_oracles`
/// but `mem_row` is the **monomial tensor** `coeff_tensor(r_outer)` (not eval eq),
/// per the coeff Spartan identity. `mem_col = z` padded; `L_row/L_col` gather.
#[allow(clippy::type_complexity)]
fn coeff_evaluation_oracles<E: Engine>(
  s: &R1CSShape<E>,
  n: usize,
  r_outer: &[E::Scalar],
  z: &[E::Scalar],
) -> (
  Vec<E::Scalar>,
  Vec<E::Scalar>,
  Vec<E::Scalar>,
  Vec<E::Scalar>,
) {
  let mem_row = coeff_tensor(r_outer); // length N = 2^num_rounds_inner
  let mem_col = padded::<E>(z, n);

  let mut l_row = vec![mem_row[0]; n];
  let mut l_col = vec![mem_col[n - 1]; n];
  for (i, (vr, vc)) in s
    .A
    .iter()
    .chain(s.B.iter())
    .chain(s.C.iter())
    .map(|(r, c, _)| (mem_row[r], mem_col[c]))
    .enumerate()
  {
    l_row[i] = vr;
    l_col[i] = vc;
  }
  (mem_row, mem_col, l_row, l_col)
}

/// A projective (coefficient-basis) ppSNARK. Reuses the eval-basis
/// [`ProverKey`]/[`VerifierKey`] and their `setup` verbatim; the two sumcheck
/// executions run in projective (coefficient) form, and every committed column
/// is opened in coefficient form via `EE` (use a `CoeffEvaluationEngine`, e.g.
/// `HyperKZGCoeffAdapter`).
#[derive(Clone, Serialize, Deserialize)]
#[serde(bound = "")]
pub struct RelaxedR1CSSNARKProjective<E: Engine, EE: EvaluationEngineTrait<E>> {
  comm_L_row: Commitment<E>,
  comm_L_col: Commitment<E>,
  comm_t_plus_r_inv_row: Commitment<E>,
  comm_w_plus_r_inv_row: Commitment<E>,
  comm_t_plus_r_inv_col: Commitment<E>,
  comm_w_plus_r_inv_col: Commitment<E>,

  sc_outer: crate::spartan::sumcheck::SumcheckProof<E>,
  eval_Az_at_r_outer: E::Scalar,
  eval_Bz_at_r_outer: E::Scalar,
  eval_Cz_at_r_outer: E::Scalar,
  eval_E_at_r_outer: E::Scalar,

  sc_inner: crate::spartan::sumcheck::SumcheckProof<E>,
  eval_E: E::Scalar,
  eval_L_row: E::Scalar,
  eval_L_col: E::Scalar,
  eval_val_A: E::Scalar,
  eval_val_B: E::Scalar,
  eval_val_C: E::Scalar,
  eval_W: E::Scalar,
  eval_t_plus_r_inv_row: E::Scalar,
  eval_row: E::Scalar,
  eval_w_plus_r_inv_row: E::Scalar,
  eval_ts_row: E::Scalar,
  eval_t_plus_r_inv_col: E::Scalar,
  eval_col: E::Scalar,
  eval_w_plus_r_inv_col: E::Scalar,
  eval_ts_col: E::Scalar,

  eval_arg: EE::EvaluationArgument,
}

impl<E: Engine, EE: EvaluationEngineTrait<E>> RelaxedR1CSSNARKTrait<E>
  for RelaxedR1CSSNARKProjective<E, EE>
{
  type ProverKey = ProverKey<E, EE>;
  type VerifierKey = VerifierKey<E, EE>;

  fn ck_floor() -> Box<dyn for<'a> Fn(&'a R1CSShape<E>) -> usize> {
    Box::new(|shape: &R1CSShape<E>| -> usize { shape.A.len() + shape.B.len() + shape.C.len() })
  }

  fn setup(
    ck: &CommitmentKey<E>,
    S: &R1CSShape<E>,
  ) -> Result<(Self::ProverKey, Self::VerifierKey), NovaError> {
    // Reuse the eval-basis setup verbatim (keys/S_repr/commitments are
    // basis-independent).
    <crate::spartan::ppsnark::RelaxedR1CSSNARK<E, EE> as RelaxedR1CSSNARKTrait<E>>::setup(ck, S)
  }

  fn prove(
    ck: &CommitmentKey<E>,
    pk: &Self::ProverKey,
    S: &R1CSShape<E>,
    U: &RelaxedR1CSInstance<E>,
    W: &RelaxedR1CSWitness<E>,
  ) -> Result<Self, NovaError> {
    let S = S.pad();
    assert!(S.is_regular_shape());
    let W = W.pad(&S);

    let mut transcript = E::TE::new(b"RelaxedR1CSSNARKProjective");
    transcript.absorb(b"vk", &pk.vk_digest);
    transcript.absorb(b"U", U);

    let z = [W.W.clone(), vec![U.u], U.X.clone()].concat();
    let (Az, Bz, Cz) = S.multiply_vec(&z)?;

    let num_rounds_inner = pk.S_repr.N.log_2();
    // Full-length outer (coeff basis): run over all N rounds so the outer point
    // matches the inner point space directly (no zero-pad factor gymnastics —
    // coeff padding factor is 1, §anchor). Pad Az/Bz/Cz/E to N.
    let n = pk.S_repr.N;
    let az = padded::<E>(&Az, n);
    let bz = padded::<E>(&Bz, n);
    let cz = padded::<E>(&Cz, n);
    let e_col = padded::<E>(&W.E, n);
    let w_col = padded::<E>(&W.W, n);

    let tau = (0..num_rounds_inner)
      .map(|_| transcript.squeeze(b"t"))
      .collect::<Result<Vec<_>, NovaError>>()?;

    // Step 1: projective outer ZeroCheck Eq^∞_τ·(Az·Bz − u·Cz·U − E·U).
    let outer = build_outer::<E>(
      num_rounds_inner,
      az.clone(),
      bz.clone(),
      cz.clone(),
      e_col.clone(),
      U.u,
      &tau,
    );
    let out_outer = outer.prove(&mut transcript);
    let r_outer = out_outer.point.clone();

    // Outer claims (coeff-form evaluations at r_outer).
    let eval_Az_at_r_outer = coeff_eval::<E>(&az, &r_outer);
    let eval_Bz_at_r_outer = coeff_eval::<E>(&bz, &r_outer);
    let eval_Cz_at_r_outer = coeff_eval::<E>(&cz, &r_outer);
    let eval_E_at_r_outer = coeff_eval::<E>(&e_col, &r_outer);

    transcript.absorb(
      b"e",
      &[
        eval_Az_at_r_outer,
        eval_Bz_at_r_outer,
        eval_Cz_at_r_outer,
        eval_E_at_r_outer,
      ]
      .as_slice(),
    );

    // Step 2: memory oracles + inner batched sumcheck.
    let (mem_row, mem_col, L_row, L_col) = coeff_evaluation_oracles::<E>(&S, n, &r_outer, &z);
    let (comm_L_row, comm_L_col) = rayon::join(
      || E::CE::commit(ck, &L_row, &E::Scalar::ZERO),
      || E::CE::commit(ck, &L_col, &E::Scalar::ZERO),
    );
    transcript.absorb(b"L", &[comm_L_row, comm_L_col].as_slice());

    let c = transcript.squeeze(b"c")?;
    let gamma = transcript.squeeze(b"g")?;
    let r_mem = transcript.squeeze(b"r")?;

    // Fingerprint tables (pointwise), coeff basis. Built in parallel (flat
    // par_iter, no nested rayon).
    let id: Vec<E::Scalar> = (0..n as u64).into_par_iter().map(E::Scalar::from).collect();
    let t_row: Vec<E::Scalar> = (0..n)
      .into_par_iter()
      .map(|i| mem_row[i] * gamma + id[i] + r_mem)
      .collect();
    let w_row: Vec<E::Scalar> = (0..n)
      .into_par_iter()
      .map(|i| L_row[i] * gamma + pk.S_repr.row[i] + r_mem)
      .collect();
    let t_col: Vec<E::Scalar> = (0..n)
      .into_par_iter()
      .map(|i| mem_col[i] * gamma + id[i] + r_mem)
      .collect();
    let w_col_fp: Vec<E::Scalar> = (0..n)
      .into_par_iter()
      .map(|i| L_col[i] * gamma + pk.S_repr.col[i] + r_mem)
      .collect();

    // Inverses: concatenate all four fingerprint tables and run a SINGLE
    // batch_invert (Montgomery batching turns 4N field inversions into 1
    // inversion + O(N) mults). batch_invert parallelizes internally, so it runs
    // sequentially w.r.t. the maps above — no nested rayon contention.
    let all_fp: Vec<E::Scalar> = t_row
      .iter()
      .chain(w_row.iter())
      .chain(t_col.iter())
      .chain(w_col_fp.iter())
      .copied()
      .collect();
    let all_inv = batch_invert(&all_fp)?;
    // inv_i * TS_i (or *1 for the w-side), flat parallel.
    let w_ones_scale = |inv: &[E::Scalar]| inv.to_vec();
    let mul_ts = |inv: &[E::Scalar], ts: &[E::Scalar]| -> Vec<E::Scalar> {
      inv
        .par_iter()
        .zip(ts.par_iter())
        .map(|(a, b)| *a * *b)
        .collect()
    };
    let t_inv_row = mul_ts(&all_inv[0..n], &pk.S_repr.ts_row);
    let w_inv_row = w_ones_scale(&all_inv[n..2 * n]);
    let t_inv_col = mul_ts(&all_inv[2 * n..3 * n], &pk.S_repr.ts_col);
    let w_inv_col = w_ones_scale(&all_inv[3 * n..4 * n]);

    // Commit the four memory oracles. Each `commit` is itself a parallel MSM,
    // so keep this loop sequential to avoid nesting parallel MSMs (which would
    // contend for the same thread pool). The MSM's internal parallelism already
    // saturates the cores.
    let comm_t_inv_row = E::CE::commit(ck, &t_inv_row, &E::Scalar::ZERO);
    let comm_w_inv_row = E::CE::commit(ck, &w_inv_row, &E::Scalar::ZERO);
    let comm_t_inv_col = E::CE::commit(ck, &t_inv_col, &E::Scalar::ZERO);
    let comm_w_inv_col = E::CE::commit(ck, &w_inv_col, &E::Scalar::ZERO);
    transcript.absorb(
      b"mem",
      &[
        comm_t_inv_row,
        comm_w_inv_row,
        comm_t_inv_col,
        comm_w_inv_col,
      ]
      .as_slice(),
    );
    let rho = (0..num_rounds_inner)
      .map(|_| transcript.squeeze(b"rho"))
      .collect::<Result<Vec<_>, NovaError>>()?;

    // Build inner instances. The eq-carrying relations (inner-E, memory T/W)
    // use the eq-factored builders so Eq^∞ and the homogenizing U are carried
    // analytically (Gruen split-eq + BDDT) instead of as dense 2^n factors; the
    // plain relations (logup inv, inner-ABC, witness-bound) stay dense. They are
    // folded together by `prove_batched_mixed`.
    let val: Vec<E::Scalar> = (0..n)
      .into_par_iter()
      .map(|i| pk.S_repr.val_A[i] + c * pk.S_repr.val_B[i] + c * c * pk.S_repr.val_C[i])
      .collect();
    let abc = build_inner_abc::<E>(num_rounds_inner, L_row.clone(), L_col.clone(), val);
    let inner_e = build_inner_e_eq_factored::<E>(num_rounds_inner, e_col.clone(), &r_outer);
    let wb = build_witness_bound::<E>(num_rounds_inner, w_col.clone(), &tau, S.num_vars.log_2());
    let (row_inv_inst, row_t, row_w) = build_memory_side_eq_factored::<E>(
      num_rounds_inner,
      t_inv_row.clone(),
      w_inv_row.clone(),
      t_row.clone(),
      w_row.clone(),
      pk.S_repr.ts_row.clone(),
      &rho,
    );
    let (col_inv_inst, col_t, col_w) = build_memory_side_eq_factored::<E>(
      num_rounds_inner,
      t_inv_col.clone(),
      w_inv_col.clone(),
      t_col.clone(),
      w_col_fp.clone(),
      pk.S_repr.ts_col.clone(),
      &rho,
    );

    let out_inner = prove_batched_mixed::<E>(
      vec![
        ProjInstance::Plain(row_inv_inst),
        ProjInstance::EqFactored(row_t),
        ProjInstance::EqFactored(row_w),
        ProjInstance::Plain(col_inv_inst),
        ProjInstance::EqFactored(col_t),
        ProjInstance::EqFactored(col_w),
        ProjInstance::Plain(abc),
        ProjInstance::EqFactored(inner_e),
        ProjInstance::Plain(wb),
      ],
      &mut transcript,
    );
    let r_inner = out_inner.point.clone();

    // Openings at r_inner. Reuse the sumcheck's own factor binding for the
    // columns that appear as factors (the batcher collapsed each factor table to
    // F_j(r_inner) for free), and only explicitly evaluate the columns that are
    // NOT sumcheck factors (val_{A,B,C}, row, col) — mirroring the eval-basis
    // prover, which likewise gets most openings for free from the sumcheck.
    // Batch order & factor layouts. Eq-factored instances expose only their
    // REAL (non-eq) factors via bound_factor_values (eq/U are analytic):
    //   [0] row_inv  Plain      = [t_inv_row, w_inv_row]
    //   [1] row_t    EqFactored = [t_inv_row, t_plus_r_row, ts_row]
    //   [3] col_inv  Plain      = [t_inv_col, w_inv_col]
    //   [4] col_t    EqFactored = [t_inv_col, t_plus_r_col, ts_col]
    //   [6] abc      Plain      = [L_row, L_col, val]
    //   [7] inner_e  EqFactored = [e_col]
    //   [8] wb       Plain      = [masked_eq, w_col]
    let bf = &out_inner.per_instance_bound_factors;
    let eval_W = bf[8][1]; // wb: w_col
    let eval_E = bf[7][0]; // inner_e (eq-factored): e_col
    let eval_L_row = bf[6][0]; // abc: L_row
    let eval_L_col = bf[6][1]; // abc: L_col
    let eval_t_plus_r_inv_row = bf[0][0]; // row_inv: t_inv_row
    let eval_w_plus_r_inv_row = bf[0][1]; // row_inv: w_inv_row
    let eval_ts_row = bf[1][2]; // row_t (eq-factored): ts_row
    let eval_t_plus_r_inv_col = bf[3][0]; // col_inv: t_inv_col
    let eval_w_plus_r_inv_col = bf[3][1]; // col_inv: w_inv_col
    let eval_ts_col = bf[4][2]; // col_t (eq-factored): ts_col
                                // Only these 5 are not sumcheck factors — evaluate explicitly (as eval does).
    let eval_val_A = coeff_eval::<E>(&pk.S_repr.val_A, &r_inner);
    let eval_val_B = coeff_eval::<E>(&pk.S_repr.val_B, &r_inner);
    let eval_val_C = coeff_eval::<E>(&pk.S_repr.val_C, &r_inner);
    let eval_row = coeff_eval::<E>(&pk.S_repr.row, &r_inner);
    let eval_col = coeff_eval::<E>(&pk.S_repr.col, &r_inner);

    // Batched PCS opening of all committed columns at r_inner.
    let comm_vec = [
      U.comm_W,
      U.comm_E,
      comm_L_row,
      comm_L_col,
      pk.S_comm.comm_val_A,
      pk.S_comm.comm_val_B,
      pk.S_comm.comm_val_C,
      comm_t_inv_row,
      pk.S_comm.comm_row,
      comm_w_inv_row,
      pk.S_comm.comm_ts_row,
      comm_t_inv_col,
      pk.S_comm.comm_col,
      comm_w_inv_col,
      pk.S_comm.comm_ts_col,
    ];
    let eval_vec = [
      eval_W,
      eval_E,
      eval_L_row,
      eval_L_col,
      eval_val_A,
      eval_val_B,
      eval_val_C,
      eval_t_plus_r_inv_row,
      eval_row,
      eval_w_plus_r_inv_row,
      eval_ts_row,
      eval_t_plus_r_inv_col,
      eval_col,
      eval_w_plus_r_inv_col,
      eval_ts_col,
    ];
    let poly_vec = [
      &w_col,
      &e_col,
      &L_row,
      &L_col,
      &pk.S_repr.val_A,
      &pk.S_repr.val_B,
      &pk.S_repr.val_C,
      &t_inv_row,
      &pk.S_repr.row,
      &w_inv_row,
      &pk.S_repr.ts_row,
      &t_inv_col,
      &pk.S_repr.col,
      &w_inv_col,
      &pk.S_repr.ts_col,
    ];
    transcript.absorb(b"e", &eval_vec.as_slice());
    let c_pcs = transcript.squeeze(b"c")?;
    let w_pe: PolyEvalWitness<E> = PolyEvalWitness::batch(&poly_vec, &c_pcs);
    let u_pe: PolyEvalInstance<E> = PolyEvalInstance::batch(&comm_vec, &r_inner, &eval_vec, &c_pcs);
    let eval_arg = EE::prove(
      ck,
      &pk.pk_ee,
      &mut transcript,
      &u_pe.c,
      &w_pe.p,
      &r_inner,
      &u_pe.e,
    )?;

    Ok(RelaxedR1CSSNARKProjective {
      comm_L_row,
      comm_L_col,
      comm_t_plus_r_inv_row: comm_t_inv_row,
      comm_w_plus_r_inv_row: comm_w_inv_row,
      comm_t_plus_r_inv_col: comm_t_inv_col,
      comm_w_plus_r_inv_col: comm_w_inv_col,
      sc_outer: out_outer.proof,
      eval_Az_at_r_outer,
      eval_Bz_at_r_outer,
      eval_Cz_at_r_outer,
      eval_E_at_r_outer,
      sc_inner: out_inner.proof,
      eval_E,
      eval_L_row,
      eval_L_col,
      eval_val_A,
      eval_val_B,
      eval_val_C,
      eval_W,
      eval_t_plus_r_inv_row,
      eval_row,
      eval_w_plus_r_inv_row,
      eval_ts_row,
      eval_t_plus_r_inv_col,
      eval_col,
      eval_w_plus_r_inv_col,
      eval_ts_col,
      eval_arg,
    })
  }

  fn verify(&self, vk: &Self::VerifierKey, U: &RelaxedR1CSInstance<E>) -> Result<(), NovaError> {
    let mut transcript = E::TE::new(b"RelaxedR1CSSNARKProjective");
    transcript.absorb(b"vk", &vk.digest());
    transcript.absorb(b"U", U);

    let num_rounds_inner = vk.S_comm.N.log_2();
    let tau = (0..num_rounds_inner)
      .map(|_| transcript.squeeze(b"t"))
      .collect::<Result<Vec<_>, NovaError>>()?;

    // Step 1: verify projective outer, reconstruct its final claim.
    let db_outer = vec![3usize; num_rounds_inner];
    let red_outer: ProjectiveSumcheckReduction<E> =
      sc_verify::<E>(E::Scalar::ZERO, &db_outer, &self.sc_outer, &mut transcript)?;
    let r_outer = red_outer.point.clone();
    let u_r_outer: E::Scalar = r_outer
      .iter()
      .fold(E::Scalar::ONE, |a, ri| a * (E::Scalar::ONE + *ri));
    let eq_tau_at_r_outer = coeff_eq_eval::<E>(&tau, &r_outer);
    let outer_expected = eq_tau_at_r_outer
      * (self.eval_Az_at_r_outer * self.eval_Bz_at_r_outer
        - U.u * self.eval_Cz_at_r_outer * u_r_outer
        - self.eval_E_at_r_outer * u_r_outer);
    if outer_expected != red_outer.final_claim {
      return Err(NovaError::InvalidSumcheckProof);
    }

    transcript.absorb(
      b"e",
      &[
        self.eval_Az_at_r_outer,
        self.eval_Bz_at_r_outer,
        self.eval_Cz_at_r_outer,
        self.eval_E_at_r_outer,
      ]
      .as_slice(),
    );

    // Step 2: memory + inner batched.
    transcript.absorb(b"L", &[self.comm_L_row, self.comm_L_col].as_slice());
    let c = transcript.squeeze(b"c")?;
    let gamma = transcript.squeeze(b"g")?;
    let r_mem = transcript.squeeze(b"r")?;
    transcript.absorb(
      b"mem",
      &[
        self.comm_t_plus_r_inv_row,
        self.comm_w_plus_r_inv_row,
        self.comm_t_plus_r_inv_col,
        self.comm_w_plus_r_inv_col,
      ]
      .as_slice(),
    );
    let rho = (0..num_rounds_inner)
      .map(|_| transcript.squeeze(b"rho"))
      .collect::<Result<Vec<_>, NovaError>>()?;

    let db_inner = vec![3usize; num_rounds_inner];
    // The batcher squeezes λ before the shared rounds; mirror it.
    let lambda = transcript.squeeze(b"projective_sumcheck_batch")?;
    let lambda_pows = powers::<E>(&lambda, 9);
    // Joint initial claim = Σ λⁱ · claim0ᵢ. Instances (prove order):
    //   [0..6): memory row/col × {inv, T, W} — all initial claim 0 (inv is the
    //           multiset balance = 0 for honest data; T/W are eq-ZeroChecks).
    //   6: inner ABC — claim0 = eval_Az + c·eval_Bz + c²·eval_Cz (coeff Spartan
    //      identity; padding factor = 1).
    //   7: inner E  — claim0 = eval_E_at_r_outer.
    //   8: witness-bound — claim0 = 0.
    let abc_init =
      self.eval_Az_at_r_outer + c * self.eval_Bz_at_r_outer + c * c * self.eval_Cz_at_r_outer;
    let joint_init = lambda_pows[6] * abc_init + lambda_pows[7] * self.eval_E_at_r_outer;
    let red_inner: ProjectiveSumcheckReduction<E> =
      sc_verify::<E>(joint_init, &db_inner, &self.sc_inner, &mut transcript)?;
    let r_inner = red_inner.point.clone();

    // Reconstruct the joint inner final claim from per-factor coeff evals.
    let u_r_inner: E::Scalar = r_inner
      .iter()
      .fold(E::Scalar::ONE, |a, ri| a * (E::Scalar::ONE + *ri));
    // mem_row(r_inner) = coeffMLE(monomial-tensor over r_outer, r_inner), O(m).
    // These reconstruct bound *factor* values (natural order).
    let eq_r_outer_at_r_inner = coeff_tensor_cross::<E>(&r_outer, &r_inner);
    let id_at_r_inner = coeff_identity_eval::<E>(num_rounds_inner, &r_inner);

    // Analytic eq for the eq-factored relations (memory T/W, inner-E) is carried
    // as eq_left = ∏_k ((1−τ_{cur})+τ_{cur}·r_k) with the eq-factored prover
    // binding the HIGH variable first, so τ_j pairs with r_inner REVERSED.
    let r_inner_rev: Vec<E::Scalar> = r_inner.iter().rev().copied().collect();
    let eq_rho_ef = coeff_eq_eval::<E>(&rho, &r_inner_rev);
    let eq_r_outer_ef = coeff_eq_eval::<E>(&r_outer, &r_inner_rev);

    // Fingerprint reconstructions at r_inner.
    let t_row_at = gamma * eq_r_outer_at_r_inner + id_at_r_inner + r_mem * u_r_inner;
    let w_row_at = gamma * self.eval_L_row + self.eval_row + r_mem * u_r_inner;
    // mem_col(r_inner) = coeffMLE(z_padded, r_inner). By linearity (anchor 6) this
    // splits into the W block [0, num_vars) plus the public IO block [u, X] at
    // [num_vars, 2·num_vars). eval_W already covers the W block; the IO block has
    // only |X|+1 nonzero entries, so its coeffMLE is an O(|X|·m) sparse monomial
    // sum — no 2^m table is materialized.
    let mem_col_at = {
      let nv = vk.num_vars;
      let m = num_rounds_inner;
      // Monomial ∏_{i∈idx} r_i for a single index (MSB-first: bit m-1-c ↔ r[c]).
      let monomial = |idx: usize| -> E::Scalar {
        let mut acc = E::Scalar::ONE;
        for (c, rc) in r_inner.iter().enumerate() {
          if (idx >> (m - 1 - c)) & 1 == 1 {
            acc *= *rc;
          }
        }
        acc
      };
      let mut io_at = U.u * monomial(nv);
      for (j, x) in U.X.iter().enumerate() {
        io_at += *x * monomial(nv + 1 + j);
      }
      self.eval_W + io_at
    };
    let t_col_at = gamma * mem_col_at + id_at_r_inner + r_mem * u_r_inner;
    let w_col_at = gamma * self.eval_L_col + self.eval_col + r_mem * u_r_inner;

    let masked_eq_at = coeff_masked_eq_eval::<E>(&tau, vk.num_vars.log_2(), &r_inner);

    // Per-instance final claims (order matches prove's batch). Each instance is
    // U-homogenized to degree 3, so a term of native degree d carries an extra
    // U(r)^(3−d) factor; the batcher lifts lower-degree instances similarly.
    //   inv:    native deg 1 → batcher lifts ×U² . claim = (t_inv − w_inv)·U²
    //   T:      Eq·(t_inv·(T+r) − TS·U) already deg 3.  (TS term: one U)
    //   W:      Eq·(w_inv·(W+r) − ones·U); ones evaluates to U(r), so the const
    //           term is −Eq·U·U = −Eq·U².
    //   ABC:    native deg 3, no U.
    //   E:      Eq·E native deg 2 → batcher lifts ×U¹.
    //   wit:    maskedEq·W native deg 2 → batcher lifts ×U¹.
    let u2 = u_r_inner * u_r_inner;
    let f0 = (self.eval_t_plus_r_inv_row - self.eval_w_plus_r_inv_row) * u2;
    let f1 = eq_rho_ef * (self.eval_t_plus_r_inv_row * t_row_at - self.eval_ts_row * u_r_inner);
    let f2 = eq_rho_ef * (self.eval_w_plus_r_inv_row * w_row_at - u2);
    let f3 = (self.eval_t_plus_r_inv_col - self.eval_w_plus_r_inv_col) * u2;
    let f4 = eq_rho_ef * (self.eval_t_plus_r_inv_col * t_col_at - self.eval_ts_col * u_r_inner);
    let f5 = eq_rho_ef * (self.eval_w_plus_r_inv_col * w_col_at - u2);
    let val_at = self.eval_val_A + c * self.eval_val_B + c * c * self.eval_val_C;
    let f6 = self.eval_L_row * self.eval_L_col * val_at;
    let f7 = eq_r_outer_ef * self.eval_E * u_r_inner;
    let f8 = masked_eq_at * self.eval_W * u_r_inner;

    let per = [f0, f1, f2, f3, f4, f5, f6, f7, f8];
    let inner_expected: E::Scalar = per
      .iter()
      .zip(lambda_pows.iter())
      .map(|(f, l)| *f * *l)
      .sum();
    if inner_expected != red_inner.final_claim {
      return Err(NovaError::InvalidSumcheckProof);
    }

    // Verify the batched PCS opening.
    let comm_vec = [
      U.comm_W,
      U.comm_E,
      self.comm_L_row,
      self.comm_L_col,
      vk.S_comm.comm_val_A,
      vk.S_comm.comm_val_B,
      vk.S_comm.comm_val_C,
      self.comm_t_plus_r_inv_row,
      vk.S_comm.comm_row,
      self.comm_w_plus_r_inv_row,
      vk.S_comm.comm_ts_row,
      self.comm_t_plus_r_inv_col,
      vk.S_comm.comm_col,
      self.comm_w_plus_r_inv_col,
      vk.S_comm.comm_ts_col,
    ];
    let eval_vec = [
      self.eval_W,
      self.eval_E,
      self.eval_L_row,
      self.eval_L_col,
      self.eval_val_A,
      self.eval_val_B,
      self.eval_val_C,
      self.eval_t_plus_r_inv_row,
      self.eval_row,
      self.eval_w_plus_r_inv_row,
      self.eval_ts_row,
      self.eval_t_plus_r_inv_col,
      self.eval_col,
      self.eval_w_plus_r_inv_col,
      self.eval_ts_col,
    ];
    transcript.absorb(b"e", &eval_vec.as_slice());
    let c_pcs = transcript.squeeze(b"c")?;
    let u_pe: PolyEvalInstance<E> = PolyEvalInstance::batch(&comm_vec, &r_inner, &eval_vec, &c_pcs);
    EE::verify(
      &vk.vk_ee,
      &mut transcript,
      &u_pe.c,
      &r_inner,
      &u_pe.e,
      &self.eval_arg,
    )?;
    Ok(())
  }
}

/// Eq-factored variant of [`build_inner_e`]: `Eq^∞_{r_outer}·E` carrying the eq
/// factor analytically (no dense eq table). Message degree `D = 2`.
pub fn build_inner_e_eq_factored<E: Engine>(
  num_vars: usize,
  e: Vec<E::Scalar>,
  r_outer: &[E::Scalar],
) -> EqFactoredVirtualPolynomial<E> {
  assert_eq!(
    r_outer.len(),
    num_vars,
    "r_outer must have num_vars entries"
  );
  EqFactoredVirtualPolynomial::new(
    num_vars,
    r_outer.to_vec(),
    vec![e],
    vec![(E::Scalar::ONE, vec![0])],
  )
}

/// Eq-factored variant of [`build_memory_side`]: the T and W relations carry
/// `Eq^∞_ρ` and the homogenizing `U` **analytically** (no dense eq / all-ones
/// `2^n` tables to store or bind), while the degree-1 `inv` difference stays a
/// plain instance. Returns `(inv, t_relation, w_relation)`.
#[allow(clippy::too_many_arguments)]
pub fn build_memory_side_eq_factored<E: Engine>(
  num_vars: usize,
  t_inv: Vec<E::Scalar>,
  w_inv: Vec<E::Scalar>,
  t_plus_r: Vec<E::Scalar>,
  w_plus_r: Vec<E::Scalar>,
  ts: Vec<E::Scalar>,
  rho: &[E::Scalar],
) -> (
  VirtualPolynomial<E>,
  EqFactoredVirtualPolynomial<E>,
  EqFactoredVirtualPolynomial<E>,
) {
  assert_eq!(rho.len(), num_vars, "rho must have num_vars entries");

  // inv: t_inv − w_inv  (degree 1, no eq) — plain.
  let inv = VirtualPolynomial::new(
    num_vars,
    vec![t_inv.clone(), w_inv.clone()],
    vec![(E::Scalar::ONE, vec![0]), (-E::Scalar::ONE, vec![1])],
  );

  // T: Eq·(t_inv·t_plus_r − ts·U). Real factors [t_inv, t_plus_r, ts].
  let t_relation = EqFactoredVirtualPolynomial::new_mixed(
    num_vars,
    rho.to_vec(),
    vec![t_inv, t_plus_r, ts],
    vec![(E::Scalar::ONE, vec![0, 1]), (-E::Scalar::ONE, vec![2])],
  );

  // W: Eq·(w_inv·w_plus_r − U²). The constant "−1" term is a pure U^Dr term
  // (the dense path's all-ones factor is the polynomial U itself).
  let w_relation = EqFactoredVirtualPolynomial::new_mixed(
    num_vars,
    rho.to_vec(),
    vec![w_inv, w_plus_r],
    vec![(E::Scalar::ONE, vec![0, 1]), (-E::Scalar::ONE, vec![])],
  );

  (inv, t_relation, w_relation)
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::{
    provider::PallasEngine,
    spartan::projective_sumcheck::{
      prove_batched, prove_batched_mixed, verify, ProjInstance, ProjectiveSumcheckReduction,
    },
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

  /// The eq-factored outer builder reduces a satisfied instance to zero and
  /// verifies at D=3 (Gruen path; correctness pinned by the differential test
  /// in eq_factored.rs against the dense path).
  #[test]
  fn outer_eq_factored_reduces_to_zero() {
    let num_vars = 3usize;
    let n = 1usize << num_vars;
    let u = Fr::from(7);
    let az: Vec<Fr> = (0..n).map(|i| Fr::from((2 * i + 1) as u64)).collect();
    let bz: Vec<Fr> = (0..n).map(|i| Fr::from((3 * i + 4) as u64)).collect();
    let cz: Vec<Fr> = (0..n).map(|i| Fr::from((i + 2) as u64)).collect();
    let e: Vec<Fr> = (0..n).map(|i| az[i] * bz[i] - u * cz[i]).collect();
    let tau: Vec<Fr> = (0..num_vars)
      .map(|i| Fr::from((5 * i + 3) as u64))
      .collect();

    let ef = build_outer_eq_factored::<E>(num_vars, az, bz, cz, e, u, &tau);
    assert_eq!(ef.degree(), 3);
    let mut ts = <E as Engine>::TE::new(b"projsc_outer_ef");
    let out = ef.prove(&mut ts);
    assert_eq!(out.initial_claim, Fr::ZERO);

    let degree_bounds = vec![3usize; num_vars];
    let mut ts_v = <E as Engine>::TE::new(b"projsc_outer_ef");
    let reduction = verify::<E>(Fr::ZERO, &degree_bounds, &out.proof, &mut ts_v).unwrap();
    assert_eq!(reduction.point, out.point);
    assert_eq!(reduction.final_claim, out.final_claim);
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

  /// Architecture-B capstone anchor for the memory T relation: the projective
  /// SC's reduced final claim equals the verifier's per-factor coeff
  /// reconstruction — `eq^∞_ρ(r)·(t_inv(r)·(T+r)(r) − TS(r)·U(r))` — where each
  /// `(r)` is a coeff-form opening the verifier would obtain (t_inv, ts from PCS;
  /// (T+r)(r) reconstructed via lemma 4; eq via coeff_eq_eval; U locally). This
  /// is the exact reconstruction the full verify performs; pinning it here proves
  /// verify can rebuild the T claim before the ~500-line assembly.
  #[test]
  fn memory_t_final_claim_reconstructs() {
    let num_vars = 4usize;
    let n = 1usize << num_vars;
    let gamma = Fr::from(11);
    let r_const = Fr::from(7);

    // Fingerprint pieces: mem_row (projective eq corners over r_outer), id, TS.
    let r_outer: Vec<Fr> = (0..num_vars)
      .map(|i| Fr::from((3 * i + 2) as u64))
      .collect();
    let mem_row = EqPolynomialProjective::<Fr>::new(r_outer.clone()).evals();
    let ts: Vec<Fr> = (0..n).map(|i| Fr::from((i % 5 + 1) as u64)).collect();
    let id: Vec<Fr> = (0..n as u64).map(Fr::from).collect();
    // T+r = mem_row·γ + id + r (pointwise), and honest t_inv = TS/(T+r).
    let t_plus_r: Vec<Fr> = (0..n)
      .map(|i| mem_row[i] * gamma + id[i] + r_const)
      .collect();
    let t_inv: Vec<Fr> = (0..n)
      .map(|i| ts[i] * t_plus_r[i].invert().unwrap())
      .collect();

    // Build only the T relation (reuse build_memory_side, take t_rel).
    let w_plus_r: Vec<Fr> = (0..n).map(|i| Fr::from((5 * i + 13) as u64)).collect();
    let w_inv: Vec<Fr> = (0..n).map(|i| w_plus_r[i].invert().unwrap()).collect();
    let rho: Vec<Fr> = (0..num_vars).map(|i| Fr::from((i + 6) as u64)).collect();
    let (_, t_rel, _) = build_memory_side::<E>(
      num_vars,
      t_inv.clone(),
      w_inv,
      t_plus_r.clone(),
      w_plus_r,
      ts.clone(),
      &rho,
    );

    let mut ts_p = <E as Engine>::TE::new(b"mem_t_recon");
    let out = t_rel.prove(&mut ts_p);
    let r = &out.point;

    // Verifier reconstruction of the T final claim from per-factor coeff evals:
    //   eq^∞_ρ(r) · ( t_inv(r) · (T+r)(r) − TS(r) · U(r) ).
    let u_at_r: Fr = r.iter().fold(Fr::ONE, |a, ri| a * (Fr::ONE + *ri));
    let eq_rho_r = coeff_eq_eval::<E>(&rho, r);
    let t_inv_r = coeff_eval::<E>(&t_inv, r); // from PCS opening
    let ts_r = coeff_eval::<E>(&ts, r); // from PCS opening
                                        // (T+r)(r) reconstructed (lemma 4): γ·mem_row(r) + id(r) + r·U(r).
    let t_plus_r_at_r = gamma * coeff_eval::<E>(&mem_row, r)
      + coeff_identity_eval::<E>(num_vars, r)
      + r_const * u_at_r;
    let recon = eq_rho_r * (t_inv_r * t_plus_r_at_r - ts_r * u_at_r);

    assert_eq!(out.final_claim, recon);
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

  /// The projective `prove_helper` equivalent: batch the inner engines
  /// (inner ABC, inner E, witness-bound, and one memory T relation) into a
  /// single sumcheck execution via λ-RLC, all lifted to the common degree D=3,
  /// and verify the combined proof. This mirrors ppSNARK folding its engines in
  /// prove_helper.
  #[test]
  fn batch_all_inner_engines() {
    let num_vars = 3usize;
    let n = 1usize << num_vars;

    // Inner ABC.
    let l_row: Vec<Fr> = (0..n).map(|i| Fr::from((i + 1) as u64)).collect();
    let l_col: Vec<Fr> = (0..n).map(|i| Fr::from((2 * i + 3) as u64)).collect();
    let val: Vec<Fr> = (0..n).map(|i| Fr::from((5 * i + 7) as u64)).collect();
    let abc = build_inner_abc::<E>(num_vars, l_row, l_col, val);

    // Inner E.
    let e: Vec<Fr> = (0..n).map(|i| Fr::from((3 * i + 2) as u64)).collect();
    let r_outer: Vec<Fr> = (0..num_vars).map(|i| Fr::from((i + 4) as u64)).collect();
    let inner_e = build_inner_e::<E>(num_vars, e, &r_outer);

    // Witness-bound.
    let w: Vec<Fr> = (0..n)
      .map(|i| {
        if i < 2 {
          Fr::from((i + 9) as u64)
        } else {
          Fr::ZERO
        }
      })
      .collect();
    let tau: Vec<Fr> = (0..num_vars).map(|i| Fr::from((i + 2) as u64)).collect();
    let wb = build_witness_bound::<E>(num_vars, w, &tau, 1);

    // Memory T relation.
    let t_plus_r: Vec<Fr> = (0..n).map(|i| Fr::from((3 * i + 11) as u64)).collect();
    let w_plus_r: Vec<Fr> = (0..n).map(|i| Fr::from((5 * i + 13) as u64)).collect();
    let ts: Vec<Fr> = (0..n).map(|i| Fr::from((i % 4 + 1) as u64)).collect();
    let t_inv: Vec<Fr> = (0..n)
      .map(|i| ts[i] * t_plus_r[i].invert().unwrap())
      .collect();
    let w_inv: Vec<Fr> = (0..n).map(|i| w_plus_r[i].invert().unwrap()).collect();
    let rho: Vec<Fr> = (0..num_vars).map(|i| Fr::from((i + 6) as u64)).collect();
    let (_, t_rel, _) =
      build_memory_side::<E>(num_vars, t_inv, w_inv, t_plus_r, w_plus_r, ts, &rho);

    // Batch all four (the batcher lifts each to D=3 via U-homogenization).
    let mut ts_p = <E as Engine>::TE::new(b"projsc_batch_all");
    let out = prove_batched::<E>(vec![abc, inner_e, wb, t_rel], &mut ts_p);

    // Verify with the transcript mirroring the λ squeeze.
    let degree_bounds = vec![3usize; num_vars];
    let mut ts_v = <E as Engine>::TE::new(b"projsc_batch_all");
    let _lambda = ts_v.squeeze(b"projective_sumcheck_batch").unwrap();
    let reduction = verify::<E>(out.initial_claim, &degree_bounds, &out.proof, &mut ts_v).unwrap();
    assert_eq!(reduction.point, out.point);
    assert_eq!(reduction.final_claim, out.final_claim);

    // Joint final = Σ λ^i G_i(r).
    let coeffs = crate::spartan::powers::<E>(&out.lambda, 4);
    let joint: Fr = out
      .per_instance_final
      .iter()
      .zip(coeffs.iter())
      .map(|(v, c)| *v * *c)
      .sum();
    assert_eq!(out.final_claim, joint);
  }

  /// F1 wire: the eq-factored inner builders (analytic eq/U) batched through the
  /// mixed batcher must produce a byte-identical proof to a **natural-order**
  /// dense-eq batch of the same relations — validating that inner-E and memory-T
  /// carry eq/U analytically without changing the proof. (The existing
  /// `EqPolynomialProjective`-order dense builders use the reverse tau order, so
  /// the reference here is built with the matching natural-order eq.)
  #[test]
  fn batch_all_inner_engines_eq_factored() {
    // Natural-order projective eq corner table (matches the eq-factored prover).
    fn eq_nat(taus: &[Fr]) -> Vec<Fr> {
      let n = taus.len();
      let mut e = vec![Fr::ZERO; 1usize << n];
      e[0] = Fr::ONE;
      for (i, &t) in taus.iter().enumerate() {
        let blk = 1usize << i;
        for b in 0..blk {
          let lo = e[b];
          e[b] = lo * (Fr::ONE - t);
          e[b + blk] = lo * t;
        }
      }
      e
    }

    let num_vars = 4usize;
    let n = 1usize << num_vars;

    let l_row: Vec<Fr> = (0..n).map(|i| Fr::from((i + 1) as u64)).collect();
    let l_col: Vec<Fr> = (0..n).map(|i| Fr::from((2 * i + 3) as u64)).collect();
    let val: Vec<Fr> = (0..n).map(|i| Fr::from((5 * i + 7) as u64)).collect();

    let e: Vec<Fr> = (0..n).map(|i| Fr::from((3 * i + 2) as u64)).collect();
    let r_outer: Vec<Fr> = (0..num_vars).map(|i| Fr::from((i + 4) as u64)).collect();

    let t_plus_r: Vec<Fr> = (0..n).map(|i| Fr::from((3 * i + 11) as u64)).collect();
    let ts: Vec<Fr> = (0..n).map(|i| Fr::from((i % 4 + 1) as u64)).collect();
    let t_inv: Vec<Fr> = (0..n)
      .map(|i| ts[i] * t_plus_r[i].invert().unwrap())
      .collect();
    let w_inv: Vec<Fr> = (0..n).map(|i| Fr::from((7 * i + 1) as u64)).collect();
    let w_plus_r: Vec<Fr> = (0..n).map(|i| Fr::from((5 * i + 13) as u64)).collect();
    let rho: Vec<Fr> = (0..num_vars).map(|i| Fr::from((i + 6) as u64)).collect();

    // Eq-factored mixed batch: [ABC (plain), inner-E (eq-fac), memory-T (eq-fac)].
    let abc_m = build_inner_abc::<E>(num_vars, l_row.clone(), l_col.clone(), val.clone());
    let inner_e_m = build_inner_e_eq_factored::<E>(num_vars, e.clone(), &r_outer);
    let (_, t_rel_m, _) = build_memory_side_eq_factored::<E>(
      num_vars,
      t_inv.clone(),
      w_inv.clone(),
      t_plus_r.clone(),
      w_plus_r.clone(),
      ts.clone(),
      &rho,
    );
    let mut ts_m = <E as Engine>::TE::new(b"projsc_batch_ef");
    let out_m = prove_batched_mixed::<E>(
      vec![
        ProjInstance::Plain(abc_m),
        ProjInstance::EqFactored(inner_e_m),
        ProjInstance::EqFactored(t_rel_m),
      ],
      &mut ts_m,
    );

    // Natural-order dense reference (same relations, eq/U as dense factors).
    let abc_n = build_inner_abc::<E>(num_vars, l_row, l_col, val);
    let inner_e_n = VirtualPolynomial::<E>::new(
      num_vars,
      vec![eq_nat(&r_outer), e],
      vec![(Fr::ONE, vec![0, 1])],
    );
    let t_rel_n = VirtualPolynomial::<E>::new_homogenized(
      num_vars,
      vec![eq_nat(&rho), t_inv, t_plus_r, ts],
      vec![(Fr::ONE, vec![0, 1, 2]), (-Fr::ONE, vec![0, 3])],
    );
    let mut ts_n = <E as Engine>::TE::new(b"projsc_batch_ef");
    let out_n = prove_batched::<E>(vec![abc_n, inner_e_n, t_rel_n], &mut ts_n);

    assert_eq!(out_m.lambda, out_n.lambda);
    assert_eq!(out_m.initial_claim, out_n.initial_claim);
    assert_eq!(out_m.point, out_n.point);
    assert_eq!(out_m.final_claim, out_n.final_claim);
    assert_eq!(out_m.per_instance_final, out_n.per_instance_final);
    let polys_m = out_m.proof.compressed_polys();
    let polys_n = out_n.proof.compressed_polys();
    assert_eq!(polys_m.len(), polys_n.len());
    for (x, y) in polys_m.iter().zip(polys_n.iter()) {
      assert_eq!(x.stored_coeffs(), y.stored_coeffs());
    }

    // And the mixed proof verifies.
    let degree_bounds = vec![3usize; num_vars];
    let mut ts_v = <E as Engine>::TE::new(b"projsc_batch_ef");
    let _l = ts_v.squeeze(b"projective_sumcheck_batch").unwrap();
    let reduction =
      verify::<E>(out_m.initial_claim, &degree_bounds, &out_m.proof, &mut ts_v).unwrap();
    assert_eq!(reduction.point, out_m.point);
    assert_eq!(reduction.final_claim, out_m.final_claim);
  }

  /// PCS boundary wiring: the coefficient-form point-transform adapter
  /// (`coeff_eval_point`) opens a committed witness at the projective sumcheck's
  /// reduced point `r` and yields exactly `coeffMLE(witness, r)` — the quantity
  /// the final-claim / verify_final_claim consumes. This ties the projective
  /// reduction to the (unmodified) eval-form PCS via the O(m) point transform,
  /// so `ppsnark_projective` uses `HyperKZGCoeffAdapter` as its EE.
  #[test]
  fn adapter_opens_coeff_form_at_reduced_point() {
    use crate::{
      provider::coeff_eval_adapter::coeff_eval_point,
      spartan::polys::multilinear::MultilinearPolynomial,
    };

    // Run a small projective relation to obtain a genuine reduced point r.
    let num_vars = 4usize;
    let n = 1usize << num_vars;
    let l_row: Vec<Fr> = (0..n).map(|i| Fr::from((i + 1) as u64)).collect();
    let l_col: Vec<Fr> = (0..n).map(|i| Fr::from((2 * i + 3) as u64)).collect();
    let val: Vec<Fr> = (0..n).map(|i| Fr::from((5 * i + 7) as u64)).collect();
    let vp = build_inner_abc::<E>(num_vars, l_row.clone(), l_col, val);
    let mut ts = <E as Engine>::TE::new(b"projsc_adapter");
    let out = vp.prove(&mut ts);
    let r = &out.point;

    // coeffMLE(l_row, r) with the MSB-first convention the adapter uses
    // (coordinate c binds bit num_vars-1-c) — matches the projective binding.
    let coeff_mle: Fr = (0..n)
      .map(|b| {
        let mut acc = l_row[b];
        for (c, rc) in r.iter().enumerate() {
          if (b >> (num_vars - 1 - c)) & 1 == 1 {
            acc *= *rc;
          }
        }
        acc
      })
      .sum();

    // The adapter reduces coeffMLE(v, r) to scale * evalMLE(v, r').
    let (r_prime, scale) = coeff_eval_point(r).unwrap();
    let eval_mle = MultilinearPolynomial::evaluate_with(&l_row, &r_prime);
    assert_eq!(coeff_mle, scale * eval_mle);
  }

  /// Integration capstone: the full milestone-critical path with **real
  /// cryptography**. Commit a coefficient-form witness with HyperKZG, run a
  /// projective sumcheck (inner ABC) to a reduced point `r`, then open the
  /// committed witness at `r` through the coefficient-form adapter and confirm:
  /// (a) the adapter prove/verify round-trip succeeds, and (b) the opened value
  /// equals `coeffMLE(witness, r)` — the exact quantity the final-oracle check
  /// consumes. This exercises commit → projective reduce → coeff-form open →
  /// verify against a live PCS.
  #[test]
  fn e2e_projective_reduce_then_coeff_open() {
    use crate::{
      provider::{
        coeff_eval_adapter::{coeff_eval_point, HyperKZGCoeffAdapter},
        Bn256EngineKZG,
      },
      traits::{
        commitment::{CommitmentEngineTrait, Len},
        evaluation::EvaluationEngineTrait,
      },
    };
    use ff::Field;
    use rand_core::OsRng;

    type Ek = Bn256EngineKZG;
    type Fk = <Ek as crate::traits::Engine>::Scalar;
    type CE = <Ek as crate::traits::Engine>::CE;
    type Adapter = HyperKZGCoeffAdapter<Ek>;

    let num_vars = 4usize;
    let n = 1usize << num_vars;

    // Witness columns (coefficient-form tables).
    let l_row: Vec<Fk> = (0..n).map(|i| Fk::from((i + 1) as u64)).collect();
    let l_col: Vec<Fk> = (0..n).map(|i| Fk::from((2 * i + 3) as u64)).collect();
    let val: Vec<Fk> = (0..n).map(|i| Fk::from((5 * i + 7) as u64)).collect();

    // Projective sumcheck (inner ABC) to obtain a genuine reduced point r.
    let vp = build_inner_abc::<Ek>(num_vars, l_row.clone(), l_col, val);
    let mut ts = <Ek as crate::traits::Engine>::TE::new(b"e2e");
    let out = vp.prove(&mut ts);
    let r = &out.point;

    // coeffMLE(l_row, r), MSB-first (matches the adapter + projective binding).
    let coeff_mle: Fk = (0..n)
      .map(|b| {
        let mut acc = l_row[b];
        for (c, rc) in r.iter().enumerate() {
          if (b >> (num_vars - 1 - c)) & 1 == 1 {
            acc *= *rc;
          }
        }
        acc
      })
      .sum();

    // Real HyperKZG setup + commit the raw coefficient vector.
    let ck = <CE as CommitmentEngineTrait<Ek>>::CommitmentKey::setup_from_rng(b"e2e", n, OsRng);
    assert!(ck.length() >= n);
    let (pk, vk) = Adapter::setup(&ck).unwrap();
    let comm = CE::commit(&ck, &l_row, &Fk::ZERO);

    // Open the committed witness at the projective reduced point in coeff form.
    let mut tr_p = <Ek as crate::traits::Engine>::TE::new(b"open");
    let arg = Adapter::prove(&ck, &pk, &mut tr_p, &comm, &l_row, r, &coeff_mle).unwrap();
    let mut tr_v = <Ek as crate::traits::Engine>::TE::new(b"open");
    assert!(Adapter::verify(&vk, &mut tr_v, &comm, r, &coeff_mle, &arg).is_ok());

    // Cross-check the transform value the adapter used.
    let (r_prime, scale) = coeff_eval_point(r).unwrap();
    let eval_mle =
      crate::spartan::polys::multilinear::MultilinearPolynomial::evaluate_with(&l_row, &r_prime);
    assert_eq!(coeff_mle, scale * eval_mle);

    // A wrong opening value must be rejected.
    let mut tr_bad = <Ek as crate::traits::Engine>::TE::new(b"open");
    assert!(Adapter::verify(&vk, &mut tr_bad, &comm, r, &(coeff_mle + Fk::ONE), &arg).is_err());
  }

  // Direct coeffMLE(v, r) = Σ_b v[b] ∏_{i∈b} r_i, MSB-first (bit num_vars-1-c ↔ r[c]).
  fn direct_coeff_mle(v: &[Fr], r: &[Fr]) -> Fr {
    let m = r.len();
    (0..v.len())
      .map(|b| {
        let mut acc = v[b];
        for (c, rc) in r.iter().enumerate() {
          if (b >> (m - 1 - c)) & 1 == 1 {
            acc *= *rc;
          }
        }
        acc
      })
      .sum()
  }

  /// `coeff_eval` (via the adapter transform) matches the direct coeff-MLE.
  #[test]
  fn coeff_eval_matches_direct() {
    for m in 1..=6usize {
      let n = 1usize << m;
      let v: Vec<Fr> = (0..n).map(|i| Fr::from((7 * i + 1) as u64)).collect();
      let r: Vec<Fr> = (0..m).map(|i| Fr::from((3 * i + 2) as u64)).collect();
      assert_eq!(coeff_eval::<E>(&v, &r), direct_coeff_mle(&v, &r), "m={m}");
    }
  }

  /// KEY coeff-vs-eval difference: for a zero-padded vector (real data in the
  /// low 2^m block, high pad bits zero), coeffMLE(padded, r_full) equals
  /// coeffMLE(v, r_low) with **factor = 1** — unlike eval basis, where the pad
  /// contributes ∏(1 − r_pad). This is because padded high bits are 0, so those
  /// corners carry no r factors.
  #[test]
  fn coeff_padding_factor_is_one() {
    let m = 3usize; // real vars
    let pad = 2usize; // padding vars (high/MSB)
    let n_low = 1usize << m;
    let n_full = 1usize << (m + pad);

    let v: Vec<Fr> = (0..n_low).map(|i| Fr::from((5 * i + 3) as u64)).collect();
    // Zero-pad into the low block; high pad bits (MSB positions) index the tail.
    let mut v_full = vec![Fr::ZERO; n_full];
    v_full[..n_low].copy_from_slice(&v);

    let r_pad: Vec<Fr> = (0..pad).map(|i| Fr::from((i + 9) as u64)).collect();
    let r_low: Vec<Fr> = (0..m).map(|i| Fr::from((2 * i + 4) as u64)).collect();
    // MSB-first: pad occupies the top positions, so r_full = [r_pad, r_low].
    let r_full: Vec<Fr> = r_pad.iter().chain(r_low.iter()).cloned().collect();

    assert_eq!(
      coeff_eval::<E>(&v_full, &r_full),
      coeff_eval::<E>(&v, &r_low)
    );
  }

  /// `coeff_identity_eval` equals coeffMLE([0,1,…,N-1], r).
  #[test]
  fn coeff_identity_eval_matches() {
    let m = 4usize;
    let n = 1usize << m;
    let r: Vec<Fr> = (0..m).map(|i| Fr::from((i + 3) as u64)).collect();
    let table: Vec<Fr> = (0..n as u64).map(Fr::from).collect();
    assert_eq!(
      coeff_identity_eval::<E>(m, &r),
      direct_coeff_mle(&table, &r)
    );
  }

  /// Architecture-B correctness anchor: the fingerprint table `T+r = mem·γ + id
  /// + r` (built pointwise) can be reconstructed at the reduced point `r_inner`
  /// from per-factor coeff evaluations:
  /// `coeffMLE(T+r, r_inner) = γ·coeffMLE(mem, r_inner) + coeff_identity_eval(r_inner)
  ///                           + r · U(r_inner)`,
  /// where `U(r) = ∏(1 + r_i)` is coeffMLE of the all-ones table. This is the
  /// linearity the verify-side memory reconstruction relies on — and note the
  /// constant `+r` becomes `r · U(r_inner)`, NOT `+r` (a real coeff-vs-eval
  /// difference vs ppsnark.rs:1511).
  #[test]
  fn coeff_fingerprint_reconstruction() {
    let m = 4usize;
    let n = 1usize << m;
    let gamma = Fr::from(11);
    let r_const = Fr::from(7); // the memory-check `r` challenge

    // mem is an arbitrary memory-contents table; id is [0..N).
    let mem: Vec<Fr> = (0..n).map(|i| Fr::from((3 * i + 5) as u64)).collect();
    let id: Vec<Fr> = (0..n as u64).map(Fr::from).collect();
    // T+r built pointwise: (mem·γ + id + r) per corner.
    let t_plus_r: Vec<Fr> = (0..n).map(|i| mem[i] * gamma + id[i] + r_const).collect();

    let r_inner: Vec<Fr> = (0..m).map(|i| Fr::from((2 * i + 3) as u64)).collect();

    // Direct coeff-MLE of the pointwise table.
    let direct = coeff_eval::<E>(&t_plus_r, &r_inner);

    // Per-factor reconstruction.
    let u_at_r: Fr = r_inner.iter().fold(Fr::ONE, |a, ri| a * (Fr::ONE + *ri));
    let recon = gamma * coeff_eval::<E>(&mem, &r_inner)
      + coeff_identity_eval::<E>(m, &r_inner)
      + r_const * u_at_r;

    assert_eq!(direct, recon);
  }

  /// Architecture-B anchor for the column memory value `mem_col = z` (z = [W, u,
  /// X] zero-padded to N). By coeffMLE linearity, `coeffMLE(mem_col, r_inner) =
  /// coeffMLE(W_padded, r_inner) + coeffMLE(io_shifted, r_inner)` where W and the
  /// (u, X) IO block occupy disjoint index ranges. This is the coeff-form
  /// counterpart of the eval-basis `eval_W + factor·… ·eval_X` reconstruction
  /// (ppsnark.rs:1552), and it is exact with NO padding factor (lemma 2).
  #[test]
  fn coeff_col_value_reconstruction() {
    let m = 5usize; // log N
    let n = 1usize << m;
    let num_vars = 1usize << (m - 1); // |W| block size (arbitrary split for the test)

    // z = [W (num_vars entries), u, X...] then zero-padded to N.
    let w: Vec<Fr> = (0..num_vars)
      .map(|i| Fr::from((3 * i + 1) as u64))
      .collect();
    let u = Fr::from(9);
    let x: Vec<Fr> = (0..(num_vars - 1))
      .map(|i| Fr::from((2 * i + 5) as u64))
      .collect();
    let z: Vec<Fr> = w
      .iter()
      .cloned()
      .chain(std::iter::once(u))
      .chain(x.iter().cloned())
      .collect();
    let mut mem_col = vec![Fr::ZERO; n];
    mem_col[..z.len()].copy_from_slice(&z);

    // Split into disjoint padded blocks: W in [0, num_vars), IO in [num_vars, ...).
    let mut w_block = vec![Fr::ZERO; n];
    w_block[..w.len()].copy_from_slice(&w);
    let mut io_block = vec![Fr::ZERO; n];
    for (i, zi) in z.iter().enumerate().skip(w.len()) {
      io_block[i] = *zi;
    }

    let r_inner: Vec<Fr> = (0..m).map(|i| Fr::from((4 * i + 3) as u64)).collect();

    let direct = coeff_eval::<E>(&mem_col, &r_inner);
    let recon = coeff_eval::<E>(&w_block, &r_inner) + coeff_eval::<E>(&io_block, &r_inner);
    assert_eq!(direct, recon);
  }

  /// Architecture-B correction anchor: in coefficient basis the memory value
  /// table `mem_row` is the **monomial tensor** `t[j] = ∏_{i∈j} r_i` (= ⊗(1,r_i)),
  /// NOT the projective eq corner table eq̃(r,·). The Spartan identity in coeff
  /// form is
  /// ```text
  ///   coeffMLE(Az, r) = Σ_k (∏_{i∈row_k} r_i) · val_A,k · (∏_{i∈col_k} r_i · z-part),
  /// ```
  /// i.e. `L_row[k] = tensor(r)[row_k]`, `L_col[k] = tensor(r)[col_k]` (over z),
  /// with `val_A` the sparse matrix value. This anchors that `mem_row/mem_col`
  /// must be built as monomial tensors (`coeff_tensor`), not projective eq — a
  /// correction to an earlier assumption, caught before assembly.
  #[test]
  fn coeff_spartan_identity_uses_monomial_tensor() {
    let num_vars = 3usize;
    let n = 1usize << num_vars;

    // Monomial tensor table t[j] = ∏_{i∈j} r_i (MSB-first: bit num_vars-1-c ↔ r[c]).
    let r: Vec<Fr> = (0..num_vars)
      .map(|i| Fr::from((2 * i + 3) as u64))
      .collect();
    let tensor: Vec<Fr> = (0..n)
      .map(|j| {
        let mut acc = Fr::ONE;
        for (c, rc) in r.iter().enumerate() {
          if (j >> (num_vars - 1 - c)) & 1 == 1 {
            acc *= *rc;
          }
        }
        acc
      })
      .collect();

    // A sparse "matrix" A as a list of (row, col, val); z is the assignment.
    // Build Az[i] = Σ_{(i,c,v)} v·z[c] as a dense table, then check the identity.
    let z: Vec<Fr> = (0..n).map(|i| Fr::from((5 * i + 1) as u64)).collect();
    let entries: Vec<(usize, usize, Fr)> = vec![
      (0, 1, Fr::from(2)),
      (1, 0, Fr::from(3)),
      (1, 3, Fr::from(4)),
      (2, 2, Fr::from(5)),
      (5, 4, Fr::from(6)),
      (7, 7, Fr::from(9)),
    ];
    let mut az = vec![Fr::ZERO; n];
    for &(row, col, v) in &entries {
      az[row] += v * z[col];
    }

    // LHS: coeffMLE(Az, r).
    let lhs = coeff_eval::<E>(&az, &r);
    // RHS: Σ_k tensor[row_k] · val_k · tensor[col_k]·(z via col) — but z is folded
    // into Az already, so the sparse form is Σ_k tensor[row_k]·val_k·(z[col_k]).
    // Here L_col[k] = tensor[col_k] would carry the z-gather; we test the row side
    // identity coeffMLE(Az,r) = Σ_k tensor[row_k]·(Az contribution at row_k), which
    // for a general Az reduces to Σ_j tensor[j]·Az[j] = coeffMLE — the tensor is the
    // correct weight table.
    let rhs: Fr = (0..n).map(|j| tensor[j] * az[j]).sum();
    assert_eq!(lhs, rhs);

    // And confirm the tensor is NOT the projective eq corner table (they differ),
    // so the earlier "mem_row = EqPolynomialProjective::evals" assumption is wrong.
    let proj_eq = EqPolynomialProjective::<Fr>::new(r.clone()).evals();
    assert_ne!(tensor, proj_eq);
  }

  /// THE end-to-end gate: a full projective ppSNARK prove→verify round-trip on
  /// a real satisfying R1CS instance, via DirectSNARK with the coefficient-form
  /// HyperKZG adapter as EE. Proves the whole assembly is sound.
  #[test]
  fn e2e_projective_ppsnark_direct() {
    use crate::{
      provider::{coeff_eval_adapter::HyperKZGCoeffAdapter, Bn256EngineKZG},
      spartan::direct::DirectSNARK,
      traits::circuit::NonTrivialCircuit,
    };

    type Ek = Bn256EngineKZG;
    type EEc = HyperKZGCoeffAdapter<Ek>;
    type Sp = RelaxedR1CSSNARKProjective<Ek, EEc>;

    let circuit = NonTrivialCircuit::<<Ek as crate::traits::Engine>::Scalar>::new(4);
    let (pk, vk) = DirectSNARK::<Ek, Sp, NonTrivialCircuit<_>>::setup(circuit.clone()).unwrap();

    let z0 = vec![<Ek as crate::traits::Engine>::Scalar::from(7)];
    let res = DirectSNARK::prove(&pk, circuit.clone(), &z0);
    assert!(
      res.is_ok(),
      "projective ppSNARK prove failed: {:?}",
      res.err()
    );
    let snark = res.unwrap();

    // NonTrivialCircuit(num_cons=4) computes z -> z^(2^4) = z^16.
    let z0v = <Ek as crate::traits::Engine>::Scalar::from(7);
    let mut z1 = z0v;
    for _ in 0..4 {
      z1 = z1 * z1;
    }
    let io: Vec<<Ek as crate::traits::Engine>::Scalar> = vec![z0v, z1];
    assert!(snark.verify(&vk, &io).is_ok());

    // Soundness: a wrong public output must be rejected (verify is not vacuous).
    let bad_io = vec![z0v, z1 + <Ek as crate::traits::Engine>::Scalar::from(1)];
    assert!(snark.verify(&vk, &bad_io).is_err());
  }
}
