use ff::{Field, PrimeField};
use rand_core::OsRng;
use rayon::iter::{IntoParallelIterator, IntoParallelRefIterator, ParallelIterator};

use crate::{
  provider::mercury::{
    ipa::{omega, IPAWitness, InputPolynomials},
    kzg::compute_quotient,
  },
  spartan::polys::univariate::UniPoly,
};

type F = halo2curves::bn256::Fr;

fn make_fft_domain<Scalar: PrimeField>(log_n: u32) -> Vec<Scalar> {
  let length = 1 << log_n;
  let mut domain = vec![Scalar::ONE; length];

  let omega = omega::<Scalar>(log_n);

  for i in 1..length {
    domain[i] = omega * domain[i - 1];
  }

  {
    for e in domain.iter().take(length).skip(1) {
      let mut res = Scalar::ONE;
      for _ in 0..length {
        res *= e;
      }
      assert_eq!(res, Scalar::ONE);
    }
  }

  domain
}

fn make_random_poly<Scalar: PrimeField>(log_n: u32) -> UniPoly<Scalar> {
  UniPoly {
    coeffs: (0..1 << log_n)
      .into_par_iter()
      .map(|_| Scalar::random(OsRng))
      .collect::<Vec<_>>(),
  }
}

fn make_poly<Scalar: PrimeField>(coeffs: &[u64]) -> UniPoly<Scalar> {
  UniPoly {
    coeffs: coeffs.iter().map(|c| Scalar::from(*c)).collect::<Vec<_>>(),
  }
}

fn verify_ipa_witness<Scalar: PrimeField>(
  input_polynomials: &InputPolynomials<Scalar>,
  witness: &IPAWitness<Scalar>,
  gamma: &Scalar,
  r: &Scalar,
) {
  let s_poly = witness.s_polynomial.clone();

  let mut gamma_pow = Scalar::ONE;

  let r_inv = r.invert().unwrap();
  let s_r = s_poly.evaluate(r);
  let s_r_inv = s_poly.evaluate(&r_inv);

  let mut lhs = Scalar::ZERO;
  let mut rhs = *r * s_r + r_inv * s_r_inv;

  let num_products = input_polynomials.lhs.len();

  for i in 0..num_products {
    let v = witness.products[i];

    let tmp = v * gamma_pow;
    rhs += tmp + tmp;

    let a_poly = &input_polynomials.lhs[i];
    let b_poly = &input_polynomials.rhs[i];

    let a_r = a_poly.evaluate(r);
    let a_r_inv = a_poly.evaluate(&r_inv);

    let b_r = b_poly.evaluate(r);
    let b_r_inv = b_poly.evaluate(&r_inv);

    let tmp = a_r * b_r_inv + a_r_inv * b_r;

    lhs += tmp * gamma_pow;

    gamma_pow *= gamma;
  }

  assert_eq!(lhs, rhs);
}

#[test]
fn test_fft() {
  let log_n = 10;
  let poly = make_random_poly::<F>(log_n);

  let domain = make_fft_domain::<F>(log_n);

  let real_evals = domain
    .par_iter()
    .map(|x| poly.evaluate(x))
    .collect::<Vec<_>>();

  let fft_evals = poly.clone().into_evaluations();

  assert_eq!(real_evals, fft_evals);

  let poly2 = UniPoly::from_evaluations(fft_evals.clone());

  assert_eq!(poly, poly2);
}

#[test]
fn test_mercury_ipa_non_constant() {
  let log_n = 11;

  let poly_a = make_random_poly(log_n);
  let poly_b = make_random_poly(log_n);
  let poly_c = make_random_poly(log_n);

  let gamma = F::random(OsRng);
  let r = F::random(OsRng);

  let input = InputPolynomials {
    lhs: vec![poly_a.clone(), poly_b.clone(), poly_c.clone()],
    rhs: vec![poly_b.clone(), poly_c.clone(), poly_c.clone()],
  };

  let witness = IPAWitness::generate(log_n, input.clone(), &gamma);

  verify_ipa_witness(&input, &witness, &gamma, &r);
}

#[test]
fn test_mercury_ipa_constant() {
  let log_n = 1;

  let poly_a = make_poly(&[rand::random::<u64>()]);
  let poly_b = make_poly(&[rand::random::<u64>()]);
  let poly_c = make_poly(&[rand::random::<u64>()]);

  let gamma = F::random(OsRng);
  let r = F::random(OsRng);

  let input = InputPolynomials {
    lhs: vec![poly_a.clone(), poly_b.clone(), poly_c.clone()],
    rhs: vec![poly_b.clone(), poly_c.clone(), poly_c.clone()],
  };

  let witness = IPAWitness::generate(log_n, input.clone(), &gamma);

  verify_ipa_witness(&input, &witness, &gamma, &r);
}

#[test]
fn test_div() {
  // (x - 5) * (x - 2)
  // = x^2 - 7x + 10
  let poly = UniPoly {
    coeffs: vec![F::from(10), -F::from(7), F::from(1)],
  };

  let alpha = F::from(5);

  let quotient = compute_quotient(&poly, &alpha);

  assert_eq!(quotient.coeffs, vec![-F::from(2), F::from(1)]);
}
