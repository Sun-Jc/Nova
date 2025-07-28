#![allow(dead_code)]

use std::cmp::min;

use ff::PrimeField;
use rayon::iter::{
  IndexedParallelIterator, IntoParallelRefIterator, IntoParallelRefMutIterator, ParallelIterator,
};

use crate::spartan::polys::univariate::{gaussian_elimination, UniPoly};

mod degree_check;
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

  /// Compute f(x) / (x - alpha)
  /// Using Horner's Method
  /// Returns (quotient_polynomial, remainder)
  fn into_div_by_deg_one_polynomial(self, alpha: &Scalar) -> (UniPoly<Scalar>, Scalar) {
    let n = self.coeffs.len();
    let mut res = self.coeffs;
    for i in (0..n - 1).rev() {
      let tmp = res[i + 1] * alpha;
      res[i] += tmp;
    }
    let (remainder, coeffs) = res.split_first().unwrap();

    let coeffs = coeffs.to_owned();

    let mut res = UniPoly { coeffs };
    res.trim();

    (res, *remainder)
  }

  // Multiply by (x + zeta)
  fn mul_by_deg_one_polynomial(&self, zeta: &Scalar) -> Self {
    let n = self.coeffs.len();
    let rhs = self.coeffs.par_iter().skip(1).map(|c| *c * *zeta);
    let res = self
      .coeffs
      .par_iter()
      .take(n - 1)
      .zip(rhs)
      .map(|(c, r)| *c + r)
      .collect::<Vec<_>>();

    let last = [self.coeffs[n - 1]];
    let first = [self.coeffs[0] * *zeta];

    let coeffs = first.into_iter().chain(res).chain(last).collect();

    UniPoly { coeffs }
  }

  fn into_sub_by_polynomial(self, rhs: &UniPoly<Scalar>) -> Self {
    let mut res = self;
    if rhs.coeffs.len() > res.coeffs.len() {
      res.raise(rhs.coeffs.len());
    }
    let n = min(rhs.coeffs.len(), res.coeffs.len());
    res
      .coeffs
      .par_iter_mut()
      .take(n)
      .zip(rhs.coeffs.par_iter().take(n))
      .for_each(|(l, r)| {
        *l -= *r;
      });

    res.trim();

    res
  }

  // Only linear or quadratic polynomials are supported
  // Adapted from `UniPoly::from_evals`
  fn from_evals_with_xs(xs: &[Scalar], evals: &[Scalar]) -> Self {
    let n = evals.len();

    let mut matrix: Vec<Vec<Scalar>> = Vec::with_capacity(n);
    for i in 0..n {
      let mut row = Vec::with_capacity(n);
      let x = xs[i];
      row.push(Scalar::ONE);
      row.push(x);
      for j in 2..n {
        row.push(row[j - 1] * x);
      }
      row.push(evals[i]);
      matrix.push(row);
    }

    let coeffs = gaussian_elimination(&mut matrix);
    Self { coeffs }
  }
}
