//! `EqPolynomialProjective`: the **coefficient-form (projective)** equality
//! polynomial, the drop-in analogue of [`EqPolynomial`](super::eq::EqPolynomial)
//! for Projective SumCheck over the infinity hypercube `{0,∞}^n`.
//!
//! Where the ordinary equality MLE uses per-variable factors
//! `(1 - ρ_i)(1 - X_i) + ρ_i X_i` (so that its value at a Boolean *point* `b`
//! is `eq̃(ρ, b)`), the projective equality polynomial drops the `(1 - X_i)`
//! term:
//! $$
//! \operatorname{Eq}^{\infty}_{\rho}(X) = \prod_{i} \big((1 - \rho_i) + \rho_i X_i\big).
//! $$
//! Its per-variable factor has constant coefficient `1 - ρ_i` and linear
//! coefficient `ρ_i`, so reading the projective corner coefficient at
//! `b ∈ {0,1}^n` gives
//! $$
//! \operatorname{Eq}^{\infty}_{\rho}(b^\infty)
//!   = \prod_i \begin{cases} 1 - \rho_i, & b_i = 0 \\ \rho_i, & b_i = 1 \end{cases}
//!   = \widetilde{\operatorname{eq}}(\rho, b),
//! $$
//! i.e. the *same* Boolean equality weight the ordinary eq MLE produces at its
//! corners. That coincidence is exactly what lets this replace `EqPolynomial`
//! as the ZeroCheck weight in the coefficient basis.
//!
//! Finite-point evaluation, in contrast, differs from the ordinary eq (no
//! `1 - X_i` term):
//! $$
//! \operatorname{Eq}^{\infty}_{\rho}(r) = \prod_i \big((1 - \rho_i) + \rho_i r_i\big).
//! $$
//!
//! This is a structured factor: it is described by `O(n)` data (`ρ`) and never
//! needs to be materialized as a dense `2^n` table for its finite-point value.
//! [`evals`](EqPolynomialProjective::evals) is provided for reference/testing;
//! a factorized prover should consume the structured form instead.

use ff::PrimeField;
use rayon::prelude::*;

/// The coefficient-form (projective) equality polynomial
/// `Eq^∞_ρ(X) = ∏_i ((1 - ρ_i) + ρ_i X_i)`.
///
/// Stores the challenge vector `ρ`. See the module docs for how its projective
/// corner coefficients recover the Boolean equality weight `eq̃(ρ, b)`.
#[derive(Debug)]
pub struct EqPolynomialProjective<Scalar: PrimeField> {
  /// The challenge vector `ρ` (one component per variable).
  pub r: Vec<Scalar>,
}

impl<Scalar: PrimeField> EqPolynomialProjective<Scalar> {
  /// Creates a new `EqPolynomialProjective` from the challenge vector `ρ`.
  pub const fn new(r: Vec<Scalar>) -> Self {
    EqPolynomialProjective { r }
  }

  /// The number of variables `n = |ρ|`.
  pub fn num_vars(&self) -> usize {
    self.r.len()
  }

  /// Evaluates `Eq^∞_ρ` at a finite point `rx`:
  /// `∏_i ((1 - ρ_i) + ρ_i · rx_i)`.
  ///
  /// Panics if `rx` and `ρ` have different lengths.
  pub fn evaluate(&self, rx: &[Scalar]) -> Scalar {
    assert_eq!(self.r.len(), rx.len());
    self
      .r
      .iter()
      .zip(rx.iter())
      .map(|(&rho_i, &rx_i)| (Scalar::ONE - rho_i) + rho_i * rx_i)
      .fold(Scalar::ONE, |acc, item| acc * item)
  }

  /// Returns the `2^n` projective corner coefficients
  /// `Eq^∞_ρ(b^∞) = eq̃(ρ, b)` for every `b ∈ {0,1}^n`.
  ///
  /// Provided for reference and differential testing; a structured factorized
  /// prover should not materialize this dense table.
  pub fn evals(&self) -> Vec<Scalar> {
    Self::evals_from_points(&self.r)
  }

  /// Computes the `2^|ρ|` projective corner coefficients from `ρ` without an
  /// intermediate polynomial.
  ///
  /// Each variable splits an accumulated coefficient into its `X_i^0` branch
  /// (multiply by `1 - ρ_i`) and its `X_i^1` branch (multiply by `ρ_i`); the
  /// low bit corresponds to `X_0`. The resulting table equals the Boolean
  /// equality weights `eq̃(ρ, b)` — division-free, valid even when `ρ_i` or
  /// `1 - ρ_i` is zero.
  pub fn evals_from_points(r: &[Scalar]) -> Vec<Scalar> {
    let ell = r.len();
    let mut evals: Vec<Scalar> = vec![Scalar::ZERO; (2_usize).pow(ell as u32)];
    let mut size = 1;
    evals[0] = Scalar::ONE;

    for &rho in r.iter().rev() {
      let (evals_left, evals_right) = evals.split_at_mut(size);
      let (evals_right, _) = evals_right.split_at_mut(size);

      // X_i^1 branch = ρ_i · prev; X_i^0 branch = (1 - ρ_i) · prev.
      zip_with_for_each!(par_iter_mut, (evals_left, evals_right), |x, y| {
        *y = *x * rho;
        *x -= &*y;
      });

      size *= 2;
    }

    evals
  }

  /// Projective corner coefficients of the **masked** projective equality
  /// polynomial: identical to [`evals`](Self::evals) but with the first
  /// `2^num_masked_vars` corner coefficients set to zero.
  ///
  /// This is the coefficient-basis analogue of `MaskedEqPolynomial`, used by
  /// ppSNARK's witness-bound sumcheck to certify that the padded tail of a
  /// witness is zero: `0 = Σ_{2^m ≤ b < 2^n} eq̃(ρ,b)·W[b]`. As a structured
  /// factor it may be materialized (Gruen's split does not apply to the masked
  /// variant).
  pub fn masked_evals(&self, num_masked_vars: usize) -> Vec<Scalar> {
    let mut evals = self.evals();
    let masked = 1usize << num_masked_vars;
    evals[..masked].iter_mut().for_each(|e| *e = Scalar::ZERO);
    evals
  }
}

impl<Scalar: PrimeField> FromIterator<Scalar> for EqPolynomialProjective<Scalar> {
  fn from_iter<I: IntoIterator<Item = Scalar>>(iter: I) -> Self {
    let r: Vec<_> = iter.into_iter().collect();
    EqPolynomialProjective { r }
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::provider::{bn256_grumpkin::bn256, pasta::pallas, secp_secq::secp256k1};
  use crate::spartan::polys::eq::EqPolynomial;

  /// The projective corner coefficients must equal the Boolean equality
  /// weights `eq̃(ρ, b)` — i.e. match the ordinary `EqPolynomial::evals`.
  fn corners_match_boolean_eq_with<F: PrimeField>() {
    let rho = vec![F::from(2), F::from(3), F::from(5)];
    let proj = EqPolynomialProjective::<F>::new(rho.clone()).evals();
    let ord = EqPolynomial::<F>::new(rho).evals();
    assert_eq!(proj, ord);
  }

  /// Finite-point evaluation follows `∏_i ((1 - ρ_i) + ρ_i r_i)`.
  ///
  /// Note this does **not** equal the ordinary eq MLE at a Boolean point: the
  /// projective/Boolean coincidence is at the *corner coefficient* level
  /// (tested above), not at finite-point evaluation (the `1 - X_i` term is
  /// gone). E.g. for ρ = (2, 3) at b = (1, 0): projective gives
  /// `1 · (1 - 3) = -2`, whereas ordinary eq gives `2 · (-2) = -4`.
  fn evaluate_with<F: PrimeField>() {
    let rho = vec![F::from(2), F::from(3)];
    let poly = EqPolynomialProjective::<F>::new(rho.clone());

    // Direct product formula at an arbitrary point.
    let r = vec![F::from(5), F::from(7)];
    let expected = ((F::ONE - rho[0]) + rho[0] * r[0]) * ((F::ONE - rho[1]) + rho[1] * r[1]);
    assert_eq!(poly.evaluate(&r), expected);
  }

  /// masked_evals zeroes exactly the first 2^m corners and keeps the rest.
  fn masked_evals_with<F: PrimeField>() {
    let rho = vec![F::from(2), F::from(3), F::from(5)];
    let full = EqPolynomialProjective::<F>::new(rho.clone()).evals();
    let masked = EqPolynomialProjective::<F>::new(rho).masked_evals(1);
    assert_eq!(masked[0], F::ZERO);
    assert_eq!(masked[1], F::ZERO);
    assert_eq!(masked[2..], full[2..]);
  }

  #[test]
  fn test_eq_projective() {
    corners_match_boolean_eq_with::<pallas::Scalar>();
    corners_match_boolean_eq_with::<bn256::Scalar>();
    corners_match_boolean_eq_with::<secp256k1::Scalar>();
    evaluate_with::<pallas::Scalar>();
    evaluate_with::<bn256::Scalar>();
    evaluate_with::<secp256k1::Scalar>();
    masked_evals_with::<pallas::Scalar>();
    masked_evals_with::<bn256::Scalar>();
    masked_evals_with::<secp256k1::Scalar>();
  }
}
