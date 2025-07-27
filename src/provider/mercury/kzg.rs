#![allow(dead_code)]

use ff::PrimeField;

use crate::spartan::polys::univariate::UniPoly;

/// Compute f(x) / (x - alpha)
/// Using Horner's Method
pub fn compute_quotient<Scalar: PrimeField>(
  dividend: &UniPoly<Scalar>,
  alpha: &Scalar,
) -> UniPoly<Scalar> {
  let n = dividend.coeffs.len();
  let mut q = vec![Scalar::ZERO; n];
  for i in (1..n).rev() {
    q[i - 1] = dividend.coeffs[i] + q[i] * alpha;
  }
  let mut res = UniPoly { coeffs: q };
  res.trim();
  res
}
