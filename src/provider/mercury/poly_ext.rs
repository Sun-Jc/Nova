#![allow(dead_code)]

use std::cmp::min;

use ff::PrimeField;
use rand_core::OsRng;
use rayon::iter::{
  IndexedParallelIterator, IntoParallelRefIterator, IntoParallelRefMutIterator, ParallelIterator,
};

use crate::spartan::polys::univariate::{gaussian_elimination, UniPoly};

impl<Scalar: PrimeField> UniPoly<Scalar> {
  /// Remove (trailing) zero coefficients of high degree monomials
  pub fn trim(&mut self) {
    while !self.coeffs.is_empty() && self.coeffs.last().unwrap() == &Scalar::ZERO {
      self.coeffs.pop();
    }
  }

  pub fn raise(&mut self, n: usize) {
    assert!(n >= self.coeffs.len());
    self.coeffs.resize(n, Scalar::ZERO);
  }

  pub fn log_n(&self) -> u32 {
    self.coeffs.len().next_power_of_two().ilog2()
  }

  /// Compute f(x) / (x - alpha)
  /// Using Horner's Method
  /// Returns (quotient_polynomial, remainder)
  pub fn into_div_by_deg_one_polynomial(self, alpha: &Scalar) -> (UniPoly<Scalar>, Scalar) {
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
  pub fn mul_by_deg_one_polynomial(&self, zeta: &Scalar) -> Self {
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

  pub fn into_add_by_polynomial(self, rhs: &UniPoly<Scalar>) -> Self {
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
        *l += *r;
      });

    res.trim();

    res
  }

  pub fn into_sub_by_polynomial(self, rhs: &UniPoly<Scalar>) -> Self {
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
  pub fn from_evals_with_xs(xs: &[Scalar], evals: &[Scalar]) -> Self {
    if evals.len() == 1 {
      return Self {
        coeffs: vec![evals[0]],
      };
    }

    assert_eq!(xs.len(), evals.len());
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

  // Only linear or quadratic polynomials are supported
  // Adapted from `UniPoly::from_evals`
  pub fn from_evals_with_xs2(xs: &[&Scalar], evals: &[Scalar]) -> Self {
    if evals.len() == 1 {
      return Self {
        coeffs: vec![evals[0]],
      };
    }

    assert_eq!(xs.len(), evals.len());
    let n = evals.len();

    let mut matrix: Vec<Vec<Scalar>> = Vec::with_capacity(n);
    for i in 0..n {
      let mut row = Vec::with_capacity(n);
      let x = xs[i];
      row.push(Scalar::ONE);
      row.push(*x);
      for j in 2..n {
        row.push(row[j - 1] * x);
      }
      row.push(evals[i]);
      matrix.push(row);
    }

    let coeffs = gaussian_elimination(&mut matrix);
    Self { coeffs }
  }

  // prod{ (x - alpha_i) }
  pub fn make_vanishing_poly(alphas: &[Scalar]) -> Self {
    let mut res = Self {
      coeffs: vec![-alphas[0], Scalar::ONE],
    };

    for alpha in alphas.iter().skip(1) {
      res = res.mul_by_deg_one_polynomial(&-*alpha);
    }

    res
  }

  pub fn scale(&mut self, s: &Scalar) {
    self.coeffs.par_iter_mut().for_each(|c| *c *= *s);
  }

  pub fn add_with(&mut self, polynomials: &[UniPoly<Scalar>], scalars: &[Scalar]) {
    self.coeffs.par_iter_mut().enumerate().for_each(|(i, c)| {
      *c += polynomials
        .par_iter()
        .zip(scalars.par_iter())
        .map(|(p, s)| *p.coeffs.get(i).unwrap_or(&Scalar::ZERO) * s)
        .sum::<Scalar>();
    });
  }
}

/// Parallelized transpose of a vector treated as a matrix
/// Input: vec, num_row, num_col
/// Output: transposed vector
pub fn transpose_parallel<T: Copy + Send + Sync>(
  vec: &[T],
  num_row: usize,
  num_col: usize,
) -> Vec<T> {
  assert_eq!(vec.len(), num_row * num_col);

  let mut result = vec![vec[0].clone(); vec.len()];

  result.par_iter_mut().enumerate().for_each(|(i, elem)| {
    let row = i / num_row;
    let col = i % num_row;
    let original_index = col * num_col + row;
    *elem = vec[original_index];
  });

  result
}

/// Parallel polynomial addition
/// Adds a list of polynomials together in parallel
pub fn parallel_poly_add<Scalar: PrimeField>(polys: &[UniPoly<Scalar>]) -> UniPoly<Scalar> {
  if polys.is_empty() {
    return UniPoly { coeffs: vec![] };
  }

  if polys.len() == 1 {
    return polys[0].clone();
  }

  // Find the maximum degree
  let max_len = polys.iter().map(|p| p.coeffs.len()).max().unwrap_or(0);

  // Initialize result with zeros
  let mut result_coeffs = vec![Scalar::ZERO; max_len];

  // Parallel addition of all polynomials
  result_coeffs
    .par_iter_mut()
    .enumerate()
    .for_each(|(i, coeff)| {
      *coeff = polys
        .iter()
        .map(|poly| poly.coeffs.get(i).copied().unwrap_or(Scalar::ZERO))
        .sum();
    });

  let mut result = UniPoly {
    coeffs: result_coeffs,
  };
  result.trim();
  result
}

/// Test method for polynomial linear combination
/// Input: (polynomials, scalars), generates random r, evaluates and linear combines, asserts zero
pub fn test_poly_linear_combination_zero<Scalar: PrimeField>(
  polys: &[UniPoly<Scalar>],
  scalars: &[Scalar],
) -> bool {
  assert_eq!(polys.len(), scalars.len());

  if polys.is_empty() {
    return true;
  }

  // Generate random evaluation point
  let r = Scalar::random(OsRng);

  // Evaluate each polynomial at r and compute linear combination
  let linear_combination: Scalar = polys
    .iter()
    .zip(scalars.iter())
    .map(|(poly, scalar)| poly.evaluate(&r) * scalar)
    .sum();

  // Check if the linear combination evaluates to zero
  linear_combination == Scalar::ZERO
}

#[cfg(test)]
mod tests {
  use super::*;
  use ff::Field;

  // Mock scalar type for testing
  type TestScalar = halo2curves::bn256::Fr;

  #[test]
  fn test_transpose_parallel() {
    let vec = vec![1, 2, 3, 4, 5, 6];
    let transposed = transpose_parallel(&vec, 2, 3);
    // Original: [[1, 2, 3], [4, 5, 6]]
    // Transposed: [[1, 4], [2, 5], [3, 6]]
    assert_eq!(transposed, vec![1, 4, 2, 5, 3, 6]);
  }

  #[test]
  fn test_parallel_poly_add() {
    let poly1 = UniPoly {
      coeffs: vec![TestScalar::ONE, TestScalar::from(2u64)],
    };
    let poly2 = UniPoly {
      coeffs: vec![TestScalar::from(3u64), TestScalar::from(4u64)],
    };

    let result = parallel_poly_add(&[poly1, poly2]);
    assert_eq!(result.coeffs.len(), 2);
    assert_eq!(result.coeffs[0], TestScalar::from(4u64)); // 1 + 3
    assert_eq!(result.coeffs[1], TestScalar::from(6u64)); // 2 + 4
  }

  #[test]
  fn test_poly_linear_combination_zero_functionality() {
    // Create two polynomials that should cancel out
    let poly1 = UniPoly {
      coeffs: vec![TestScalar::ONE, TestScalar::from(2u64)],
    };
    let poly2 = UniPoly {
      coeffs: vec![-TestScalar::ONE, -TestScalar::from(2u64)],
    };

    let scalars = vec![TestScalar::ONE, TestScalar::ONE];
    let polys = vec![poly1, poly2];

    // This should be true since poly1 + poly2 = 0
    assert!(test_poly_linear_combination_zero(&polys, &scalars));
  }
}
