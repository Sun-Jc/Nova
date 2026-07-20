//! Identity polynomial implementation.

use core::marker::PhantomData;
use ff::PrimeField;

/// A polynomial that evaluates to the identity of its input
pub struct IdentityPolynomial<Scalar: PrimeField> {
  ell: usize,
  _p: PhantomData<Scalar>,
}

impl<Scalar: PrimeField> IdentityPolynomial<Scalar> {
  /// Creates a new identity polynomial with the given number of variables
  pub fn new(ell: usize) -> Self {
    IdentityPolynomial {
      ell,
      _p: PhantomData,
    }
  }

  /// Evaluates the polynomial at the given point
  pub fn evaluate(&self, r: &[Scalar]) -> Scalar {
    assert_eq!(self.ell, r.len());
    let mut power_of_two = 1_u64;
    (0..self.ell)
      .rev()
      .map(|i| {
        let result = Scalar::from(power_of_two) * r[i];
        power_of_two *= 2;
        result
      })
      .fold(Scalar::ZERO, |acc, item| acc + item)
  }

  /// Returns the **coefficient-basis (projective)** dense corner table of the
  /// identity: the length-`2^ell` vector `[0, 1, 2, ..., 2^ell − 1]`.
  ///
  /// In the evaluation basis the identity is `O(ell)` structured (its *value*
  /// at Boolean point `b` is the index `b`, `id(r) = Σ 2^i r_i`). In the
  /// coefficient basis a factor whose *corner coefficient* equals the index
  /// requires the dense table `[0, 1, …, N−1]`. This is a public setup-time
  /// constant (used by ppSNARK's memory fingerprint `T = mem·γ + i`), not a
  /// per-proof witness, but it does degrade from `O(ell)` to `O(2^ell)`; a
  /// structured coefficient-basis representation, if found, is a future
  /// optimization.
  pub fn projective_corner_table(&self) -> Vec<Scalar> {
    (0..(1u64 << self.ell)).map(Scalar::from).collect()
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::provider::pasta::pallas;

  #[test]
  fn projective_corner_table_is_index_sequence() {
    let id = IdentityPolynomial::<pallas::Scalar>::new(3);
    let table = id.projective_corner_table();
    assert_eq!(table.len(), 8);
    for (i, &v) in table.iter().enumerate() {
      assert_eq!(v, pallas::Scalar::from(i as u64));
    }
  }
}
