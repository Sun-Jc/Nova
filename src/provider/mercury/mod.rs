use ff::PrimeField;

use crate::spartan::polys::univariate::UniPoly;

mod ipa;
mod kzg;
mod split_polynomial;

#[cfg(test)]
mod tests;

impl<Scalar: PrimeField> UniPoly<Scalar> {
  /// Remove (trailing) zero coefficients of high degree monomials
  fn trim(&mut self) {
    while self.coeffs.last().unwrap() == &Scalar::ZERO {
      self.coeffs.pop();
    }
  }

  fn raise(&mut self, n: usize) {
    self.coeffs.resize(n, Scalar::ZERO);
  }

  fn log_n(&self) -> u32 {
    self.coeffs.len().next_power_of_two().ilog2()
  }
}
