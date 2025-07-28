#![allow(dead_code)]

use ff::PrimeField;

use crate::spartan::polys::univariate::UniPoly;

pub fn d_polynomial<Scalar: PrimeField>(
  polynomial: &UniPoly<Scalar>,
  log_n: u32,
) -> UniPoly<Scalar> {
  let mut d_poly = polynomial.clone();
  d_poly.raise(1 << log_n);
  d_poly.coeffs.reverse();
  d_poly.trim();
  d_poly
}
