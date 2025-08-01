use std::time::Instant;

use ff::{Field, PrimeField};
use halo2curves::bn256;
use rand_core::OsRng;
use rayon::iter::{IntoParallelIterator, IntoParallelRefIterator, ParallelIterator};
use serde::Serialize;

use crate::provider::hyperkzg;
use crate::provider::keccak::Keccak256Transcript;
use crate::spartan::polys::eq::EqPolynomial;
use crate::spartan::polys::multilinear::MultilinearPolynomial;
use crate::traits::commitment::CommitmentEngineTrait;
use crate::traits::evaluation::EvaluationEngineTrait;
use crate::traits::{Engine, TranscriptEngineTrait};
use crate::{
  provider::{
    hyperkzg::{CommitmentKey, VerifierKey},
    mercury::{
      self,
      degree_check::d_polynomial,
      ipa::{omega, IPAWitness, InputPolynomials},
      kzg::EvaluationProcess,
      split_polynomial::split_polynomial,
    },
    traits::DlogGroup,
    Bn256EngineKZG,
  },
  spartan::polys::univariate::UniPoly,
};

type F = halo2curves::bn256::Fr;
type E = Bn256EngineKZG;
type EE = mercury::engine::EvaluationEngine<E>;

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
fn test_div_deg_one() {
  // (x - 5) * (x - 2) + 3
  // = x^2 - 7x + 13
  let poly = UniPoly {
    coeffs: vec![F::from(13), -F::from(7), F::from(1)],
  };

  let alpha = F::from(5);

  let (quotient, remainder) = poly.into_div_by_deg_one_polynomial(&alpha);

  assert_eq!(quotient.coeffs, vec![-F::from(2), F::from(1)]);
  assert_eq!(remainder, F::from(3));
}

#[test]
fn test_div_deg_one_zero_remainder() {
  // (x - 5) * (x - 2)
  // = x^2 - 7x + 10
  let poly = UniPoly {
    coeffs: vec![F::from(10), -F::from(7), F::from(1)],
  };

  let alpha = F::from(5);

  let (quotient, remainder) = poly.into_div_by_deg_one_polynomial(&alpha);

  assert_eq!(quotient.coeffs, vec![-F::from(2), F::from(1)]);
  assert_eq!(remainder, F::from(0));
}

// #[test]
// fn test_split_polynomial() {
//   let log_n = 22;
//   let poly = make_random_poly::<F>(log_n);

//   let split_res = split_polynomial(&poly, log_n);

//   let b = 1 << (log_n / 2);

//   assert_eq!(split_res.len(), b);

//   // f(x) = sum_i { x^i f_i(x^b) }
//   let r = F::random(OsRng);
//   let lhs = poly.evaluate(&r);
//   let r_b = r.pow([b as u64]);
//   let rhs = split_res
//     .iter()
//     .enumerate()
//     .map(|(i, p)| r.pow([i as u64]) * p.evaluate(&r_b))
//     .sum::<F>();

//   assert_eq!(lhs, rhs);
// }

// #[test]
// fn test_div_x_b_alpha() {
//   // f(X) = (X^b − α) * q(X) + g(X).
//   let log_n = 10;
//   let poly = make_random_poly::<F>(log_n);

//   let alpha = F::random(OsRng);

//   let (quotient, remainder) = divide_polynomial_by_x_b_alpha(&poly, log_n, &alpha);

//   let r = F::random(OsRng);

//   let lhs = poly.evaluate(&r);

//   let b = 1 << (log_n / 2);
//   let r_b = r.pow([b as u64]);
//   let rhs = (r_b - alpha) * quotient.evaluate(&r) + remainder.evaluate(&r);

//   assert_eq!(lhs, rhs);
// }

#[test]
fn test_degree_check() {
  let log_n = 12;
  let poly = make_random_poly::<F>(log_n);

  // d(x) = X^{1<<log_n - 1} * f(1/X)
  let d_poly = d_polynomial(&poly, log_n);

  let r = F::random(OsRng);

  let r_pow = r.pow([((1 << log_n) - 1) as u64]);
  let r_inv = r.invert().unwrap();

  let lhs = d_poly.evaluate(&r);
  let rhs = r_pow * poly.evaluate(&r_inv);

  assert_eq!(lhs, rhs);
}

#[test]
fn test_mul_by_deg_one_polynomial() {
  // (x + 1)(x + 2) * (x - 2)(x - 3)
  // = (x^2 + 3x + 2) * ..
  let poly = make_poly(&[2, 3, 1]);

  let zeta = &-F::from(2);

  let res = poly.mul_by_deg_one_polynomial(zeta);

  let zeta = &-F::from(3);
  let res = res.mul_by_deg_one_polynomial(zeta);

  assert_eq!(
    res.coeffs,
    vec![
      F::from(12),
      F::from(8),
      -F::from(7),
      -F::from(2),
      F::from(1)
    ]
  );
}

#[test]
fn test_sub_by_deg_one_polynomial() {
  // (x + 1)(x + 2) - (x + 2)
  // = (x^2 + 3x + 2) - (x + 2)
  // = x^2 + 2x
  let poly = make_poly::<F>(&[2, 3, 1]);

  let lhs = make_poly(&[2, 1]);

  let res = poly.into_sub_by_polynomial(&lhs);

  assert_eq!(res.coeffs, vec![F::ZERO, F::from(2), F::from(1)]);
}

#[test]
fn test_from_evals_with_xs() {
  let xs = vec![F::from(0), F::from(1), F::from(5)];
  let evals = vec![F::from(1), F::from(1), F::from(3)];

  let poly = UniPoly::from_evals_with_xs(&xs, &evals);

  let tenth = F::from(10).invert().unwrap();

  assert_eq!(poly.coeffs, vec![F::ONE, -tenth, tenth]);
}

#[test]
fn test_batch_kzg() {
  let log_n = 12;

  let g_poly = make_random_poly::<F>(log_n);
  let h_poly = make_random_poly::<F>(log_n);
  let s_poly = make_random_poly::<F>(log_n);
  let d_poly = make_random_poly::<F>(log_n);

  let alpha = F::random(OsRng);
  let zeta = F::random(OsRng);

  let ck = CommitmentKey::<E>::setup_from_rng(b"test", 1 << (log_n * 2), OsRng);
  let vk = VerifierKey {
    G: bn256::G1::gen().affine(),
    H: bn256::G2::gen().affine(),
    tau_H: ck.tau_H.clone(),
  };

  let start = Instant::now();
  let (witness, eval) =
    EvaluationProcess::init_eval_0(&alpha, &zeta, &g_poly, &h_poly, &s_poly, &d_poly);

  let dur0 = start.elapsed();
  let start = Instant::now();

  let mut witness = witness;

  witness.commit_1(&ck);

  let dur1 = start.elapsed();
  let start = Instant::now();

  // let beta = witness.sample_beta_2();
  let beta = F::random(OsRng);

  witness.open_phase_3(&beta);

  let dur3 = start.elapsed();
  let start = Instant::now();

  witness.commit_quot_m_4(&ck);

  let dur4 = start.elapsed();
  let start = Instant::now();

  // let z = witness.sample_z_5();
  let z = F::random(OsRng);

  witness.open_phase_6(&z);

  let dur6 = start.elapsed();
  let start = Instant::now();

  witness.commit_quot_l_7(&ck);

  let dur7 = start.elapsed();

  let proof = witness.into_proof_8();

  let start = Instant::now();
  proof.verify(&vk, &eval);
  let dur8 = start.elapsed();

  dbg!(&dur0, &dur1, &dur3, &dur4, &dur6, &dur7, &dur8);
}

#[test]
fn test_mercury_ee() {
  let log_n = 11;
  let poly = make_random_poly::<F>(log_n);
  let point = (0..log_n).map(|_| F::random(OsRng)).collect::<Vec<_>>();

  let ck = <<E as Engine>::CE as CommitmentEngineTrait<E>>::CommitmentKey::setup_from_rng(
    b"test",
    1 << log_n,
    OsRng,
  );

  let (pk, vk) = EE::setup(&ck);

  let eval = MultilinearPolynomial::new(poly.coeffs.clone()).evaluate(&point);

  let mut transcript = <E as Engine>::TE::new(b"test");

  let comm = <E as Engine>::CE::commit(&ck, &poly.coeffs, &F::ZERO);

  let start = Instant::now();
  let arg = EE::prove(
    &ck,
    &pk,
    &mut transcript,
    &comm,
    &poly.coeffs,
    &point,
    &eval,
  )
  .unwrap();
  let dur = start.elapsed();

  println!("Mercury: {:?}", &dur);

  {
    let (pk, vk) = hyperkzg::EvaluationEngine::setup(&ck);

    let mut transcript = <E as Engine>::TE::new(b"test");

    let start = Instant::now();
    let arg = hyperkzg::EvaluationEngine::<E>::prove(
      &ck,
      &pk,
      &mut transcript,
      &comm,
      &poly.coeffs,
      &point,
      &eval,
    )
    .unwrap();
    let dur = start.elapsed();

    println!("HyperKZG: {:?}", &dur);

    let mut transcript = <E as Engine>::TE::new(b"test");

    hyperkzg::EvaluationEngine::<E>::verify(&vk, &mut transcript, &comm, &point, &eval, &arg)
      .unwrap();
  }

  let mut transcript = <E as Engine>::TE::new(b"test");

  EE::verify(&vk, &mut transcript, &comm, &point, &eval, &arg).unwrap();
}
