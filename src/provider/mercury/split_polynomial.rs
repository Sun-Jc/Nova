#![allow(dead_code)]

use ff::PrimeField;
use rayon::prelude::*;

use crate::{
  provider::mercury::kzg::div_polynomial_by_deg_one, spartan::polys::univariate::UniPoly,
};

/// Split a univariate polynomial into column polynomials by viewing coefficients as a 2D matrix
pub fn split_polynomial<Scalar: PrimeField>(
  polynomial: &UniPoly<Scalar>,
  log_n: u32,
) -> Vec<UniPoly<Scalar>> {
  assert!(log_n % 2 == 0, "log_n must be even (log_n = 2*b)");

  let b = 1 << (log_n / 2);

  assert_eq!(
    polynomial.coeffs.len(),
    1 << log_n,
    "Polynomial must have exactly 2^log_n coefficients"
  );

  (0..b)
    .into_par_iter()
    .map(|col| {
      let column_coeffs = (0..b)
        .into_par_iter()
        .map(|row| polynomial.coeffs[row * b + col])
        .collect();

      UniPoly {
        coeffs: column_coeffs,
      }
    })
    .collect()
}

/// Divide a polynomial by x^b - alpha
/// Returns (quotient, remainder)
pub fn divide_polynomial_by_x_b_alpha<Scalar: PrimeField>(
  polynomial: &UniPoly<Scalar>,
  log_n: u32,
  alpha: &Scalar,
) -> (UniPoly<Scalar>, UniPoly<Scalar>) {
  assert!(log_n % 2 == 0, "log_n must be even (log_n = 2*b)");

  let b = 1 << (log_n / 2);

  let split_res = split_polynomial(polynomial, log_n);

  let (quotients, remainder_coefficients): (Vec<Vec<Scalar>>, Vec<Scalar>) = split_res
    .par_iter()
    .map(|p| {
      let (mut quotient, remainder) = div_polynomial_by_deg_one(p, alpha);
      quotient.raise(b);
      (quotient.coeffs, remainder)
    })
    .unzip();

  let untransposed_quotient_coefficients = quotients.into_iter().flatten().collect::<Vec<_>>();
  let quotient_coefficients: Vec<Scalar> = (0..b * b)
    .into_par_iter()
    .map(|i| {
      let col = i / b;
      let row = i % b;
      untransposed_quotient_coefficients[row * b + col]
    })
    .collect();

  (
    UniPoly {
      coeffs: quotient_coefficients,
    },
    UniPoly {
      coeffs: remainder_coefficients,
    },
  )
}
