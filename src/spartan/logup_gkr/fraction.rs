//! Projective fraction arithmetic for the Logup-GKR fractional-sum tree.
//!
//! A fraction `num/den` is kept in projective form so the GKR circuit never
//! performs a field inversion: the fraction-add gate combines two children
//! `(n_l, d_l)` and `(n_r, d_r)` into `(n_l·d_r + n_r·d_l, d_l·d_r)`. This is the
//! root of why Logup-GKR avoids the inverse-polynomial commitments that the
//! current ppSNARK memory-check pays for.

use ff::Field;
use serde::{Deserialize, Serialize};

/// A projective fraction `numerator / denominator` over the engine scalar field.
///
/// Equality of the represented rationals is cross-multiplicative
/// (`a/b == c/d` iff `a·d == c·b`); this type does not normalize.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(bound = "F: Serialize + for<'a> Deserialize<'a>")]
pub struct Fraction<F: Field> {
  /// Numerator.
  pub num: F,
  /// Denominator.
  pub den: F,
}

impl<F: Field> Fraction<F> {
  /// Creates a new projective fraction `num/den`.
  pub fn new(num: F, den: F) -> Self {
    Self { num, den }
  }

  /// The additive identity `0/1`, used as the padding cell for the tree.
  pub fn zero() -> Self {
    Self {
      num: F::ZERO,
      den: F::ONE,
    }
  }

  /// Projective fraction addition: `a/b + c/d = (a·d + c·b)/(b·d)`.
  ///
  /// This is the 2-to-1 gate applied at every internal node of the tree.
  pub fn add(&self, other: &Self) -> Self {
    Self {
      num: self.num * other.den + other.num * self.den,
      den: self.den * other.den,
    }
  }
}
