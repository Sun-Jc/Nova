#![allow(dead_code)]

use ff::PrimeField;
use halo2curves::fft::best_fft;
use rayon::iter::{
  IndexedParallelIterator, IntoParallelRefIterator, IntoParallelRefMutIterator, ParallelIterator,
};

use crate::spartan::polys::univariate::UniPoly;

pub(crate) fn omega<Scalar: PrimeField>(log_n: u32) -> Scalar {
  Scalar::ROOT_OF_UNITY.pow([1_u64 << (Scalar::S - log_n)])
}

#[derive(Clone)]
pub struct InputPolynomials<Scalar: PrimeField> {
  pub lhs: Vec<UniPoly<Scalar>>,
  pub rhs: Vec<UniPoly<Scalar>>,
}

pub struct IPAWitness<Scalar: PrimeField> {
  pub s_polynomial: UniPoly<Scalar>,
  pub products: Vec<Scalar>,
}

fn fft<const IFFT: bool, Scalar: PrimeField>(v: &mut Vec<Scalar>, log_n: u32) {
  let omega = if IFFT {
    omega::<Scalar>(log_n).invert().unwrap()
  } else {
    omega::<Scalar>(log_n)
  };

  v.resize(1 << log_n, Scalar::ZERO);
  best_fft(v, omega, log_n);
}

impl<Scalar: PrimeField> UniPoly<Scalar> {
  pub fn into_evaluations(self) -> Vec<Scalar> {
    let log_n = self.log_n();
    let mut res = self.coeffs;
    fft::<false, _>(&mut res, log_n);
    res
  }

  pub fn from_evaluations(evals: Vec<Scalar>) -> Self {
    let mut res = UniPoly { coeffs: evals };
    let log_n = res.log_n();

    fft::<true, _>(&mut res.coeffs, log_n);
    res.trim();

    let n_inv = Scalar::from(1 << log_n as u64).invert().unwrap();
    res.coeffs.par_iter_mut().for_each(|x| *x *= n_inv);

    res
  }
}

pub fn make_s_polynomial<Scalar: PrimeField>(
  a_polys: Vec<UniPoly<Scalar>>,
  b_polys: Vec<UniPoly<Scalar>>,
  log_n: u32,
  gamma: &Scalar,
) -> UniPoly<Scalar> {
  assert!(log_n >= 1);

  assert!(log_n >= a_polys[0].log_n());
  assert!(log_n >= b_polys[0].log_n());

  let mut a_polys = a_polys;
  let mut b_polys = b_polys;

  a_polys.par_iter_mut().for_each(|p| p.trim());
  b_polys.par_iter_mut().for_each(|p| p.trim());

  // corner case: all constants
  if a_polys.iter().all(|p| p.coeffs.len() == 1) && b_polys.iter().all(|p| p.coeffs.len() == 1) {
    return UniPoly {
      coeffs: vec![Scalar::ZERO],
    };
  }

  let n = 1 << log_n;
  let n2 = n * 2;

  debug_assert!(a_polys
    .iter()
    .zip(b_polys.iter())
    .all(|(a, b)| a.coeffs.len() <= n && b.coeffs.len() <= n));

  // x^{n-1} * [ sum { (a(x) * b(1/x) + a(1/x) * b(x)) * gamma^j } ]
  let mut lhs_evals = vec![Scalar::ZERO; n2];
  let mut gamma_pow = Scalar::ONE;

  a_polys.into_iter().zip(b_polys).for_each(|(a, b)| {
    let mut a = a;
    let mut b = b;
    a.raise(n2);
    b.raise(n2);

    let (a_evals, b_evals) = rayon::join(|| a.into_evaluations(), || b.into_evaluations());

    let lhs0 = a_evals[0] * b_evals[0];
    lhs_evals[0] += (lhs0 + lhs0) * gamma_pow;

    lhs_evals[1..]
      .par_iter_mut()
      .enumerate()
      .for_each(|(i, v)| {
        let i = i + 1;
        let tmp = a_evals[i] * b_evals[n2 - i] + a_evals[n2 - i] * b_evals[i];
        *v += tmp * gamma_pow;
      });

    gamma_pow *= gamma;
  });

  // dbg!(&lhs_evals);

  // x^{n-1}
  let omega_n_1 = omega::<Scalar>(log_n + 1).pow([n as u64 - 1]);
  let mut x_pow = omega_n_1;
  for v in lhs_evals.iter_mut().take(n2).skip(1) {
    *v *= x_pow;
    x_pow *= omega_n_1;
  }

  dbg!(&lhs_evals);

  let mut coeffs = UniPoly::from_evaluations(lhs_evals).coeffs;

  assert!(coeffs.len() < n2);

  coeffs.drain(..n);

  UniPoly { coeffs }
}

impl<Scalar: PrimeField> IPAWitness<Scalar> {
  pub fn generate(log_n: u32, input_polynomials: InputPolynomials<Scalar>, gamma: &Scalar) -> Self {
    assert_eq!(input_polynomials.lhs.len(), input_polynomials.rhs.len());
    assert!(!input_polynomials.lhs.is_empty());

    assert!(log_n < Scalar::S);

    let compute_inner_product = |(a, b): (&UniPoly<Scalar>, &UniPoly<Scalar>)| {
      a.coeffs
        .par_iter()
        .zip(b.coeffs.par_iter())
        .map(|(ca, cb)| *ca * *cb)
        .sum::<Scalar>()
    };

    let products = input_polynomials
      .lhs
      .iter()
      .zip(input_polynomials.rhs.iter())
      .map(compute_inner_product)
      .collect();

    let s_polynomial =
      make_s_polynomial(input_polynomials.lhs, input_polynomials.rhs, log_n, gamma);

    Self {
      s_polynomial,
      products,
    }
  }
}
