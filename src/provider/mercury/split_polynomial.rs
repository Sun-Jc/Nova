#![allow(dead_code)]

use ff::PrimeField;
use rayon::prelude::*;

use crate::{provider::mercury::poly_ext::transpose_parallel, spartan::polys::univariate::UniPoly};

/// Split a univariate polynomial into column polynomials by viewing coefficients as a 2D matrix
pub fn split_polynomial<Scalar: PrimeField>(
  polynomial: &UniPoly<Scalar>,
  b: usize,
) -> Vec<UniPoly<Scalar>> {
  let coeffs = transpose_parallel(&polynomial.coeffs, b, b);

  coeffs
    .par_chunks(b)
    .map(|c| UniPoly {
      coeffs: c.to_owned(),
    })
    .collect()
}

/// Divide a polynomial by x^b - alpha
/// Returns (quotient, remainder)
pub fn divide_by_binomial<Scalar: PrimeField>(
  split_res: &[UniPoly<Scalar>],
  alpha: &Scalar,
) -> (UniPoly<Scalar>, UniPoly<Scalar>) {
  let num_col = split_res.len();
  let num_row = split_res[0].coeffs.len();

  let (quotients, remainder_coefficients): (Vec<Vec<Scalar>>, Vec<Scalar>) = split_res
    .into_par_iter()
    .map(|p| {
      let (mut quotient, remainder) = p.clone().into_div_by_deg_one_polynomial(alpha);
      quotient.raise(num_row);
      (quotient.coeffs, remainder)
    })
    .unzip();

  assert_eq!(quotients.len(), num_col);
  assert_eq!(quotients[0].len(), num_row);

  let untransposed_quotient_coefficients = quotients.into_iter().flatten().collect::<Vec<_>>();
  let quotient_coefficients =
    transpose_parallel(&untransposed_quotient_coefficients, num_row, num_col);

  assert_eq!(quotient_coefficients.len(), num_col * num_row);

  (
    UniPoly {
      coeffs: quotient_coefficients,
    },
    UniPoly {
      coeffs: remainder_coefficients,
    },
  )
}
