//! Proof objects for the Logup-GKR fractional-sum argument.
//!
//! # One batched tree, not N independent trees
//! The logup instances (`row`, `col`) are padded to a uniform height and folded
//! through a **single** `BatchTree` with **one shared GKR depth**, so the
//! evaluation point is shared *by construction*. There is no per-tree proof and
//! no "assert the two points are equal". Per layer, the instances' column pairs
//! are batched into one sumcheck via a fresh `λ`. Sumcheck round polynomials use
//! Nova's `CompressedUniPoly`.

use crate::spartan::logup_gkr::fraction::Fraction;
use crate::spartan::polys::univariate::CompressedUniPoly;
use crate::traits::Engine;
use serde::{Deserialize, Serialize};

/// `rlc(a, b, r) = a + r·(b - a)` — the two-to-one fold of split claims.
#[inline(always)]
fn rlc<F: ff::Field>(a: F, b: F, r: F) -> F {
  a + r * (b - a)
}

/// A claim about one instance's column pair after a layer: the fraction
/// `num/den`.
pub type LayerClaim<E> = Fraction<<E as Engine>::Scalar>;

/// The split (even/odd) final claim of one instance at one layer: the `left`
/// and `right` child fractions `(nL,dL)` and `(nR,dR)`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(bound = "")]
pub struct LayerFinalClaim<E: Engine> {
  /// Left child `(nL, dL)`.
  pub left: Fraction<E::Scalar>,
  /// Right child `(nR, dR)`.
  pub right: Fraction<E::Scalar>,
}

impl<E: Engine> LayerFinalClaim<E> {
  /// Builds from the four folded sumcheck evaluations `nL, nR, dL, dR`.
  ///
  /// WARNING: the layer sumcheck's `bound_evals` come out **interleaved** as
  /// `[nL, dR, nR, dL]` (matching the flattened virtual-polynomial order
  /// `numerator_left, denominator_right, numerator_right, denominator_left`).
  /// Map them explicitly — `new(evals[0], evals[2], evals[3], evals[1])` —
  /// never slice `evals[0..4]` into the parameters positionally.
  pub fn new(nL: E::Scalar, nR: E::Scalar, dL: E::Scalar, dR: E::Scalar) -> Self {
    Self {
      left: Fraction::new(nL, dL),
      right: Fraction::new(nR, dR),
    }
  }

  /// Folds the split claim into the next layer's claim via `rlc` at `r`.
  pub fn fold_into_next_claim(&self, r: E::Scalar) -> LayerClaim<E> {
    Fraction::new(
      rlc(self.left.num, self.right.num, r),
      rlc(self.left.den, self.right.den, r),
    )
  }

  /// The gate output `left + right` (projective fraction add).
  pub fn compute_gate(&self) -> Fraction<E::Scalar> {
    self.left + self.right
  }
}

/// Sumcheck transcript for one batched tree layer: one compressed round
/// polynomial per variable (the instances are batched into this single
/// sumcheck via `λ`).
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(bound = "")]
pub struct LayerSumcheck<E: Engine> {
  /// Compressed round polynomials, one per sumcheck round of this layer.
  pub round_polys: Vec<CompressedUniPoly<E::Scalar>>,
}

/// A single **batched** Logup-GKR proof over all instances (row, col).
///
/// Instance-indexed vectors are ordered `[row, col]`. `initial_claims` are the
/// per-instance output-layer fractions (observed first). `final_claims[layer]`
/// holds one split claim per instance; `sumchecks[layer]` is the one batched
/// sumcheck for that layer. Ordering is output→input, and the top transition
/// (0-variable layer) carries no sumcheck, so
/// `sumchecks.len() + 1 == final_claims.len()`.
///
/// The per-layer batching challenge `λ` is **not** stored here: it is a
/// Fiat-Shamir challenge the verifier re-samples fresh at each layer (reusing
/// one `λ` across layers would let the prover adaptively forge each layer's
/// claims). Likewise `initial_claims` are not bound to committed data by this
/// proof alone; soundness closes at the host's reconcile step, where the
/// returned `openings` must match the fractions the host recomputes from its
/// real `L_row`/`L_col` openings.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(bound = "")]
pub struct LogupGkrProof<E: Engine> {
  /// Per-instance output-layer claims (the root fractions before fold-down).
  pub initial_claims: Vec<LayerClaim<E>>,
  /// Per-layer, per-instance split final claims, output→input.
  pub final_claims: Vec<Vec<LayerFinalClaim<E>>>,
  /// Per-layer batched sumcheck (one fewer than `final_claims`).
  pub sumchecks: Vec<LayerSumcheck<E>>,
}

/// Verifier output — a **continuation token** carrying the shared opening
/// claim, with its point field named Nova's `eval_point`.
///
/// The host reads only [`Self::eval_point`] to build its merged PCS opening set;
/// `openings` (the per-instance input-layer fractions the GKR reduced to, order
/// `[row, col]`) is handed to the host's reconcile step, which is where the
/// `0/den` zero-sum check lives — **not** inside the GKR verifier. See the host
/// contract in [`crate::spartan::logup_gkr`].
#[derive(Clone, Debug)]
pub struct LogupGkrOpeningClaim<E: Engine> {
  eval_point: Vec<E::Scalar>,
  openings: Vec<Fraction<E::Scalar>>,
}

impl<E: Engine> LogupGkrOpeningClaim<E> {
  /// Constructs the token (only the GKR verifier should call this).
  pub fn new(eval_point: Vec<E::Scalar>, openings: Vec<Fraction<E::Scalar>>) -> Self {
    Self {
      eval_point,
      openings,
    }
  }

  /// The single shared evaluation point the host opens its columns at.
  pub fn eval_point(&self) -> &[E::Scalar] {
    &self.eval_point
  }

  /// The reduced per-instance input-layer fractions (`[row, col]`), which the
  /// host's reconcile step recomputes from its own `L_row`/`L_col` openings and
  /// compares against.
  pub fn openings(&self) -> &[Fraction<E::Scalar>] {
    &self.openings
  }
}
