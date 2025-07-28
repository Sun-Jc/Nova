#![allow(dead_code)]

use ff::PrimeField;

use crate::spartan::polys::univariate::UniPoly;

/// Compute f(x) / (x - alpha)
/// Using Horner's Method
/// Returns (quotient_polynomial, remainder)
pub fn div_polynomial_by_deg_one<Scalar: PrimeField>(
  dividend: &UniPoly<Scalar>,
  alpha: &Scalar,
) -> (UniPoly<Scalar>, Scalar) {
  let n = dividend.coeffs.len();
  let mut res = dividend.coeffs.clone();
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
