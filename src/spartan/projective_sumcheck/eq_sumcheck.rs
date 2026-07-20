//! `EqSumCheckInstanceProjective`: the eq-build portion of a Gruen-style
//! ([eprint 2024/108, §3]) equality-polynomial sumcheck instance, adapted to
//! the **projective** (coefficient-basis) equality polynomial
//! `Eq^∞_ρ(X) = ∏_i ((1 - ρ_i) + ρ_i X_i)`.
//!
//! # Scope of this file
//!
//! This implements only the parts **closely related to the equality
//! polynomial's construction** — the pieces that mirror
//! [`EqPolynomial`](crate::spartan::polys::eq::EqPolynomial) /
//! [`EqPolynomialProjective`](crate::spartan::polys::eq_projective::EqPolynomialProjective):
//!
//! - the **split-half prefix eq tables** (`poly_eq_left` / `poly_eq_right`),
//!   built in parallel with `rayon::join`, so the eq factor is stored as two
//!   `2^{n/2}` halves instead of one dense `2^n` table (Gruen's core space win);
//! - the per-variable **endpoint cache** of the projective evaluation points,
//!   and the running eq scalar updated by [`bound`](Self::bound).
//!
//! The round-polynomial kernels (the `evaluation_points_*` methods that combine
//! these halves against instance-specific witness columns) are **out of scope
//! here** and are left to the eventual projective sumcheck round logic.
//!
//! # Projective endpoints: two, not three
//!
//! Gruen's Boolean instance caches three points `{0, ∞, -1}` per `τ`. The
//! projective per-variable factor `(1 - τ) + τ·X` is read at only the two
//! projective endpoints:
//! ```text
//!   Eq^∞ at X=0  = constant coeff   = 1 - τ,
//!   Eq^∞ at X=∞  = leading coeff    = τ.
//! ```
//! So the endpoint cache stores just `(1 - τ, τ)`.

use crate::traits::Engine;
use ff::Field;
use rayon::prelude::*;

/// The eq-construction state of a projective equality-polynomial sumcheck
/// instance: split-half prefix tables plus cached projective endpoints.
///
/// Round-polynomial computation is intentionally not part of this struct yet;
/// see the module docs.
pub struct EqSumCheckInstanceProjective<E: Engine> {
  /// Number of variables `n` at construction.
  init_num_vars: usize,
  /// Size of the first (left) half, `⌊n/2⌋`.
  first_half: usize,
  /// Size of the second (right) half, `n - ⌊n/2⌋`.
  second_half: usize,
  /// 1-based round counter (starts at 1, like the Boolean instance).
  round: usize,
  /// The challenge vector `ρ` (the "taus").
  taus: Vec<E::Scalar>,
  /// Running eq scalar for the already-bound prefix (starts at `1`).
  eval_eq_left: E::Scalar,
  /// Prefix eq tables for the left half (excluding `ρ_0`, folded via
  /// `eval_eq_left`): `poly_eq_left[k]` has `2^k` entries.
  poly_eq_left: Vec<Vec<E::Scalar>>,
  /// Prefix eq tables for the right half: `poly_eq_right[k]` has `2^k` entries.
  poly_eq_right: Vec<Vec<E::Scalar>>,
  /// Per-variable projective endpoints `(Eq^∞ at 0, Eq^∞ at ∞) = (1 - τ, τ)`.
  eq_tau_0_inf: Vec<(E::Scalar, E::Scalar)>,
}

impl<E: Engine> EqSumCheckInstanceProjective<E> {
  /// Builds the projective eq instance from the challenge vector `ρ` (`taus`).
  ///
  /// Mirrors the Boolean Gruen constructor: split `ρ` into halves, build each
  /// half's prefix eq tables in parallel, and cache the projective endpoints.
  /// The corner-coefficient recurrence is identical to
  /// [`EqPolynomialProjective::evals_from_points`](crate::spartan::polys::eq_projective::EqPolynomialProjective::evals_from_points)
  /// (split each accumulated coefficient into its `×τ` and `×(1-τ)` branches);
  /// only the endpoint cache differs (two projective points instead of three).
  pub fn new(taus: Vec<E::Scalar>) -> Self {
    let l = taus.len();
    let first_half = l / 2;

    // Prefix eq tables: result[i] is the eq corner-coefficient table over the
    // first `i` taus. Each step appends the `×τ` branch and turns the original
    // half into the `×(1-τ)` branch — division-free, parallel.
    let compute_eq_polynomials = |taus: Vec<&E::Scalar>| -> Vec<Vec<E::Scalar>> {
      let len = taus.len();
      let mut result = Vec::with_capacity(len + 1);
      result.push(vec![E::Scalar::ONE]);

      for (i, &tau) in taus.iter().enumerate() {
        let prev = &result[i];
        let mut v_next = prev.to_vec();
        v_next.par_extend(prev.par_iter().map(|v| *v * tau));
        let (first, last) = v_next.split_at_mut(prev.len());
        first.par_iter_mut().zip(last).for_each(|(a, b)| *a -= *b);
        result.push(v_next);
      }

      result
    };

    // Left half drops ρ_0 (folded into `eval_eq_left` by `bound` at round 1);
    // both halves are consumed most-significant-first, matching the Boolean
    // instance's split.
    let (left_taus, right_taus) = taus.split_at(first_half);
    let left_taus = left_taus.iter().skip(1).rev().collect::<Vec<_>>();
    let right_taus = right_taus.iter().rev().collect::<Vec<_>>();

    let (poly_eq_left, poly_eq_right) = rayon::join(
      || compute_eq_polynomials(left_taus),
      || compute_eq_polynomials(right_taus),
    );

    // Projective endpoints: (Eq^∞ at 0, Eq^∞ at ∞) = (1 - τ, τ).
    let eq_tau_0_inf = taus
      .par_iter()
      .map(|tau| (E::Scalar::ONE - tau, *tau))
      .collect::<Vec<_>>();

    Self {
      init_num_vars: l,
      first_half,
      second_half: l - first_half,
      round: 1,
      taus,
      eval_eq_left: E::Scalar::ONE,
      poly_eq_left,
      poly_eq_right,
      eq_tau_0_inf,
    }
  }

  /// Folds the just-bound variable's challenge `r` into the running eq scalar
  /// `eval_eq_left` and advances the round.
  ///
  /// Unlike the Boolean Gruen instance (whose factor is
  /// `(1-τ)(1-r) + τ·r`), the **projective** per-variable factor is
  /// `(1 - τ) + τ·r` — the `(1 - X)` term is gone, exactly as in
  /// [`EqPolynomialProjective::evaluate`](crate::spartan::polys::eq_projective::EqPolynomialProjective::evaluate).
  /// So the running eq prefix is `∏_i ((1 - τ_i) + τ_i · r_i)`.
  pub fn bound(&mut self, r: &E::Scalar) {
    let tau = self.taus[self.round - 1];
    self.eval_eq_left *= (E::Scalar::ONE - tau) + tau * r;
    self.round += 1;
  }

  /// The number of variables `n`.
  pub fn num_vars(&self) -> usize {
    self.init_num_vars
  }

  /// The size of the first (left) half, `⌊n/2⌋`.
  pub fn first_half(&self) -> usize {
    self.first_half
  }

  /// The size of the second (right) half, `n - ⌊n/2⌋`.
  pub fn second_half(&self) -> usize {
    self.second_half
  }

  /// The running eq scalar for the bound prefix (`1` before any `bound`).
  pub fn eval_eq_left(&self) -> E::Scalar {
    self.eval_eq_left
  }

  /// The cached projective endpoints `(1 - τ_i, τ_i)` per variable.
  pub fn endpoints(&self) -> &[(E::Scalar, E::Scalar)] {
    &self.eq_tau_0_inf
  }

  /// The full right-half eq corner table (`2^{second_half}` entries).
  pub fn poly_eq_right_full(&self) -> &[E::Scalar] {
    &self.poly_eq_right[self.second_half]
  }

  /// The full left-half eq corner table (`2^{first_half - 1}` entries; empty
  /// prefix when `first_half == 0`).
  pub fn poly_eq_left_full(&self) -> &[E::Scalar] {
    self.poly_eq_left.last().unwrap()
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::{
    provider::PallasEngine, spartan::polys::eq_projective::EqPolynomialProjective, traits::Engine,
  };

  type E = PallasEngine;
  type Fr = <E as Engine>::Scalar;

  #[test]
  fn build_matches_eq_projective() {
    // n = 5 → first_half = 2, second_half = 3.
    let taus: Vec<Fr> = (0..5).map(|i| Fr::from((3 * i + 1) as u64)).collect();
    let inst = EqSumCheckInstanceProjective::<E>::new(taus.clone());

    // Right half table = Eq^∞ corner coefficients over taus[first_half..].
    let right_expected =
      EqPolynomialProjective::<Fr>::new(taus[inst.first_half..].to_vec()).evals();
    assert_eq!(inst.poly_eq_right_full(), right_expected.as_slice());

    // Left half table = Eq^∞ corner coefficients over taus[1..first_half]
    // (ρ_0 is folded into eval_eq_left, not the table).
    let left_expected =
      EqPolynomialProjective::<Fr>::new(taus[1..inst.first_half].to_vec()).evals();
    assert_eq!(inst.poly_eq_left_full(), left_expected.as_slice());

    // Endpoint cache holds the two projective points (1 - τ, τ).
    for (i, &tau) in taus.iter().enumerate() {
      assert_eq!(inst.endpoints()[i], (Fr::ONE - tau, tau));
    }
  }

  #[test]
  fn bound_tracks_finite_eq_prefix() {
    let taus: Vec<Fr> = vec![Fr::from(2), Fr::from(3), Fr::from(5)];
    let mut inst = EqSumCheckInstanceProjective::<E>::new(taus.clone());
    assert_eq!(inst.eval_eq_left(), Fr::ONE);

    // Bind two challenges; eval_eq_left must equal ∏ ((1-τ_i) + τ_i r_i).
    let rs = [Fr::from(7), Fr::from(11)];
    inst.bound(&rs[0]);
    inst.bound(&rs[1]);

    let expected =
      ((Fr::ONE - taus[0]) + taus[0] * rs[0]) * ((Fr::ONE - taus[1]) + taus[1] * rs[1]);
    assert_eq!(inst.eval_eq_left(), expected);
  }
}
