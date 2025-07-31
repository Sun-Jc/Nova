#![allow(dead_code)]
//! # KZG Batch Proof
//! [BCGH20]
use std::fmt::Debug;

use ff::{Field, PrimeField};
use rayon::iter::{IndexedParallelIterator, IntoParallelIterator, ParallelIterator};

use crate::{
  provider::{
    hyperkzg::{Commitment, CommitmentEngine, CommitmentKey, VerifierKey},
    traits::{DlogGroup, PairingGroup},
  },
  spartan::polys::univariate::UniPoly,
  traits::{commitment::CommitmentEngineTrait, Engine},
};

fn empty_poly<Scalar: PrimeField>() -> UniPoly<Scalar> {
  UniPoly {
    coeffs: Default::default(),
  }
}

pub struct BatchEvaluation<Scalar: PrimeField> {
  pub zeta: Scalar,
  pub alpha: Scalar,

  pub g_zeta: Scalar,
  pub g_zeta_inv: Scalar,

  pub h_zeta: Scalar,
  pub h_zeta_inv: Scalar,
  pub h_alpha: Scalar,

  pub s_zeta: Scalar,
  pub s_zeta_inv: Scalar,

  pub d_zeta: Scalar,
}

#[derive(Debug)]
pub struct BatchKZGProof<E: Engine>
where
  E::GE: PairingGroup,
{
  pub comm_g: Commitment<E>,
  pub comm_h: Commitment<E>,
  pub comm_s: Commitment<E>,
  pub comm_d: Commitment<E>,

  pub beta: E::Scalar,
  pub z: E::Scalar,

  pub comm_quot_m: Commitment<E>,
  pub comm_quot_l: Commitment<E>,
}

#[derive(Debug)]
pub struct EvaluationProcess<E: Engine>
where
  E::GE: PairingGroup,
{
  pub zeta: E::Scalar,
  pub zeta_inv: E::Scalar,
  pub alpha: E::Scalar,

  pub comm_g: Commitment<E>,
  pub comm_h: Commitment<E>,
  pub comm_s: Commitment<E>,
  pub comm_d: Commitment<E>,

  pub beta: E::Scalar,
  pub z: E::Scalar,

  pub g_poly: UniPoly<E::Scalar>,
  pub h_poly: UniPoly<E::Scalar>,
  pub s_poly: UniPoly<E::Scalar>,
  pub d_poly: UniPoly<E::Scalar>,

  pub g_star: UniPoly<E::Scalar>,
  pub h_star: UniPoly<E::Scalar>,
  pub s_star: UniPoly<E::Scalar>,
  pub d_star: UniPoly<E::Scalar>,

  pub g_minus_g_star: UniPoly<E::Scalar>,
  pub h_minus_h_star: UniPoly<E::Scalar>,
  pub s_minus_s_star: UniPoly<E::Scalar>,
  pub d_minus_d_star: UniPoly<E::Scalar>,

  pub quot_m_poly: UniPoly<E::Scalar>,
  pub comm_qm: Commitment<E>,

  pub quot_l_poly: UniPoly<E::Scalar>,
  pub comm_ql: Commitment<E>,
}

impl<E: Engine> EvaluationProcess<E>
where
  E::GE: PairingGroup,
  E: Engine<CE = CommitmentEngine<E>>,
{
  pub fn init_eval_0(
    alpha: &E::Scalar,
    zeta: &E::Scalar,
    g_poly: &UniPoly<E::Scalar>,
    h_poly: &UniPoly<E::Scalar>,
    s_poly: &UniPoly<E::Scalar>,
    d_poly: &UniPoly<E::Scalar>,
  ) -> (Self, BatchEvaluation<E::Scalar>) {
    let zeta_inv = &zeta.invert().unwrap();

    let xss = vec![
      vec![*zeta, *zeta_inv],
      vec![*zeta, *zeta_inv, *alpha],
      vec![*zeta, *zeta_inv],
      vec![*zeta],
    ];

    let [g_zeta, g_zeta_inv, h_zeta, h_zeta_inv, h_alpha, s_zeta, s_zeta_inv, d_zeta] = {
      // Compute Evaluations
      // g | zeta
      //   |        zeta_inv
      // h | zeta
      //   |        zeta_inv
      //   |                    alpha
      // s | zeta
      //   |        zeta_inv
      // d | zeta

      let polys = vec![
        g_poly, g_poly, h_poly, h_poly, h_poly, s_poly, s_poly, d_poly,
      ];
      let xs = xss.iter().flatten().collect::<Vec<_>>();

      polys
        .into_par_iter()
        .zip(xs)
        .map(|(poly, x)| poly.evaluate(x))
        .collect::<Vec<_>>()
        .try_into()
        .unwrap()
    };

    let [g_star, h_star, s_star, d_star] = {
      // Interpolate
      // g_star | zeta, zeta_inv
      // h_star | zeta, zeta_inv, alpha
      // s_star | zeta, zeta_inv
      // d_star | zeta

      let evals = vec![
        vec![g_zeta, g_zeta_inv],
        vec![h_zeta, h_zeta_inv, h_alpha],
        vec![s_zeta, s_zeta_inv],
        vec![d_zeta],
      ];

      xss
        .clone()
        .into_par_iter()
        .zip(evals)
        .map(|(xs, evals)| UniPoly::from_evals_with_xs(&xs, &evals))
        .collect::<Vec<_>>()
        .try_into()
        .unwrap()
    };

    let [g_minus_g_star, h_minus_h_star, s_minus_s_star, d_minus_d_star] = {
      // Compute Diff Polynomial
      // g(x) - g_star(x)
      // h(x) - h_star(x)
      // s(x) - s_star(x)
      // d(x) - d_star(x)

      let lhs = [
        g_poly.clone(),
        h_poly.clone(),
        s_poly.clone(),
        d_poly.clone(),
      ];
      let rhs = [&g_star, &h_star, &s_star, &d_star];

      lhs
        .into_par_iter()
        .zip(rhs)
        .map(|(poly, poly_star)| poly.into_sub_by_polynomial(poly_star))
        .collect::<Vec<_>>()
        .try_into()
        .unwrap()
    };

    let witness = Self {
      zeta: *zeta,
      zeta_inv: *zeta_inv,
      alpha: *alpha,

      g_poly: g_poly.clone(),
      h_poly: h_poly.clone(),
      s_poly: s_poly.clone(),
      d_poly: d_poly.clone(),

      g_star,
      h_star,
      s_star,
      d_star,

      g_minus_g_star,
      h_minus_h_star,
      s_minus_s_star,
      d_minus_d_star,

      comm_g: Default::default(),
      comm_h: Default::default(),
      comm_s: Default::default(),
      comm_d: Default::default(),
      comm_qm: Default::default(),
      comm_ql: Default::default(),

      z: Default::default(),
      beta: Default::default(),
      quot_m_poly: empty_poly(),
      quot_l_poly: empty_poly(),
    };

    let eval = BatchEvaluation {
      zeta: *zeta,
      alpha: *alpha,
      g_zeta,
      g_zeta_inv,
      s_zeta,
      s_zeta_inv,
      h_zeta,
      h_zeta_inv,
      h_alpha,
      d_zeta,
    };

    (witness, eval)
  }

  pub fn commit_1(&mut self, ck: &CommitmentKey<E>) {
    let vs = vec![
      &self.g_poly.coeffs,
      &self.h_poly.coeffs,
      &self.s_poly.coeffs,
      &self.d_poly.coeffs,
    ];

    let [comm_g, comm_h, comm_s, comm_d] = vs
      .into_par_iter()
      .map(|v| E::CE::commit(ck, v, &E::Scalar::ZERO))
      .collect::<Vec<_>>()
      .try_into()
      .unwrap();

    self.comm_g = comm_g;
    self.comm_h = comm_h;
    self.comm_s = comm_s;
    self.comm_d = comm_d;
  }

  pub fn sample_beta_2(&self) -> E::Scalar {
    todo!()
  }

  pub fn open_phase_3(&mut self, beta: &E::Scalar) {
    self.beta = *beta;

    let beta_2 = *beta * *beta;
    let beta_3 = beta_2 * beta;

    let m_poly = {
      // compute m(X)
      //             z_poly_t_s1 * g_diff
      // + beta   * (z_poly_t_s2 * h_diff)
      // + beta^2 * (z_poly_t_s3 * s_diff)
      // + beta^3 * (z_poly_t_s4 * d_diff)

      let lhs = vec![
        &self.g_minus_g_star,
        &self.h_minus_h_star,
        &self.s_minus_s_star,
        &self.d_minus_d_star,
      ];

      let rhs = vec![
        vec![-self.alpha],
        vec![],
        vec![-self.alpha],
        vec![-self.alpha, -self.zeta_inv],
      ];

      let [g_prod, h_prod, s_prod, d_prod] = lhs
        .into_par_iter()
        .zip(rhs)
        .map(|(poly, rhs)| {
          let mut poly = poly.clone();
          rhs
            .into_iter()
            .for_each(|c| poly = poly.mul_by_deg_one_polynomial(&c));

          poly
        })
        .collect::<Vec<_>>()
        .try_into()
        .unwrap();

      let mut h_prod = h_prod;
      let mut s_prod = s_prod;
      let mut d_prod = d_prod;

      h_prod.scale(&beta);
      s_prod.scale(&beta_2);
      d_prod.scale(&beta_3);

      let m_poly = g_prod
        .clone()
        .into_add_by_polynomial(&h_prod)
        .into_add_by_polynomial(&s_prod)
        .into_add_by_polynomial(&d_prod);

      m_poly
    };

    let q_m_poly = {
      // Compute Quotient Polynomial
      // m(X) / (X - alpha) / (X - zeta) / (X - zeta_inv)
      let (q_m_poly, rem) = m_poly.into_div_by_deg_one_polynomial(&self.alpha);
      assert_eq!(rem, E::Scalar::ZERO);
      let (q_m_poly, rem) = q_m_poly.into_div_by_deg_one_polynomial(&self.zeta);
      assert_eq!(rem, E::Scalar::ZERO);
      let (q_m_poly, rem) = q_m_poly.into_div_by_deg_one_polynomial(&self.zeta_inv);
      assert_eq!(rem, E::Scalar::ZERO);

      q_m_poly
    };

    self.quot_m_poly = q_m_poly;
  }

  pub fn commit_quot_m_4(&mut self, ck: &CommitmentKey<E>) {
    let comm_quot_m = E::CE::commit(ck, &self.quot_m_poly.coeffs, &E::Scalar::ZERO);
    self.comm_qm = comm_quot_m;
  }

  pub fn sample_z_5(&self) -> E::Scalar {
    todo!()
  }

  pub fn open_phase_6(&mut self, z: &E::Scalar) {
    self.z = *z;
    // L(X) = m_z(X) - Z_T(z) \cdot q_m(X)

    let beta_2 = self.beta * self.beta;
    let beta_3 = beta_2 * self.beta;

    let t_s1_eval_at_z = *z - self.alpha;
    let t_s2_eval_at_z = E::Scalar::ONE;
    let t_s3_eval_at_z = t_s1_eval_at_z;
    let t_s4_eval_at_z = (*z - self.zeta_inv) * t_s1_eval_at_z;
    let t_eval_at_z = t_s4_eval_at_z * (*z - self.zeta);

    let [g_poly_z, h_poly_z, s_poly_z, d_poly_z] = {
      // Compute mz(X)
      // =          z_poly_t_s1(z) * (g_poly(X) - g_star(z))
      // + beta   * z_poly_t_s2(z) * (h_poly(X) - h_star(z))
      // + beta^2 * z_poly_t_s3(z) * (s_poly(X) - s_star(z))
      // + beta^3 * z_poly_t_s4(z) * (d_poly(X) - d_star(z))
      //   lh s     mhs               rhs

      let lhs = vec![E::Scalar::ONE, self.beta, beta_2, beta_3];

      let mhs = vec![
        t_s1_eval_at_z,
        t_s2_eval_at_z,
        t_s3_eval_at_z,
        t_s4_eval_at_z,
      ];

      let star_eval = [&self.g_star, &self.h_star, &self.s_star, &self.d_star]
        .into_par_iter()
        .map(|poly| poly.evaluate(z))
        .collect::<Vec<_>>();

      let rhs = [&self.g_poly, &self.h_poly, &self.s_poly, &self.d_poly]
        .into_par_iter()
        .zip(star_eval)
        .map(|(poly, star_eval)| {
          let mut res = poly.clone();
          res.coeffs[0] -= star_eval;
          res
        })
        .collect::<Vec<_>>();

      lhs
        .into_par_iter()
        .zip(mhs)
        .zip(rhs)
        .map(|((lhs, mhs), rhs)| {
          let mut res = rhs.clone();
          res.scale(&(mhs * lhs));
          res
        })
        .collect::<Vec<_>>()
        .try_into()
        .unwrap()
    };

    // TODO: can be parallelized
    let mz_poly = g_poly_z
      .into_add_by_polynomial(&h_poly_z)
      .into_add_by_polynomial(&s_poly_z)
      .into_add_by_polynomial(&d_poly_z);

    let l_poly = {
      let mut q_m = self.quot_m_poly.clone();
      q_m.scale(&t_eval_at_z);

      mz_poly.clone().into_sub_by_polynomial(&q_m)
    };

    let (quot_l_poly, rem) = l_poly.into_div_by_deg_one_polynomial(z);
    assert_eq!(rem, E::Scalar::ZERO);

    self.quot_l_poly = quot_l_poly;
  }

  pub fn commit_quot_l_7(&mut self, ck: &CommitmentKey<E>) {
    let comm_quot_l = E::CE::commit(ck, &self.quot_l_poly.coeffs, &E::Scalar::ZERO);
    self.comm_ql = comm_quot_l;
  }

  pub fn into_proof_8(self) -> BatchKZGProof<E> {
    BatchKZGProof {
      comm_g: self.comm_g,
      comm_h: self.comm_h,
      comm_s: self.comm_s,
      comm_d: self.comm_d,

      beta: self.beta,
      z: self.z,
      comm_quot_m: self.comm_qm,
      comm_quot_l: self.comm_ql,
    }
  }
}

impl<E: Engine> BatchKZGProof<E>
where
  E: Engine<CE = CommitmentEngine<E>>,
  E::GE: PairingGroup,
{
  pub fn verify(&self, vk: &VerifierKey<E>, eval: &BatchEvaluation<E::Scalar>) {
    let zeta = eval.zeta;
    let zeta_inv = zeta.invert().unwrap();
    let alpha = eval.alpha;

    // TODO: To Sample Locally
    let z = self.z;
    let beta = self.beta;

    let g_zeta = eval.g_zeta;
    let g_zeta_inv = eval.g_zeta_inv;
    let h_zeta = eval.h_zeta;
    let h_zeta_inv = eval.h_zeta_inv;
    let h_alpha = eval.h_alpha;
    let s_zeta = eval.s_zeta;
    let s_zeta_inv = eval.s_zeta_inv;
    let d_zeta = eval.d_zeta;

    let g_star = UniPoly::from_evals_with_xs(&[zeta, zeta_inv], &[g_zeta, g_zeta_inv]);
    let h_star =
      UniPoly::from_evals_with_xs(&[zeta, zeta_inv, alpha], &[h_zeta, h_zeta_inv, h_alpha]);
    let s_star = UniPoly::from_evals_with_xs(&[zeta, zeta_inv], &[s_zeta, s_zeta_inv]);
    let d_star = UniPoly::from_evals_with_xs(&[zeta], &[d_zeta]);

    let g_star_eval = g_star.evaluate(&z);
    let h_star_eval = h_star.evaluate(&z);
    let s_star_eval = s_star.evaluate(&z);
    let d_star_eval = d_star.evaluate(&z);

    let van_zeta = z - zeta;
    let van_zeta_inv = z - zeta_inv;
    let van_alpha = z - alpha;

    let z_eval_t_s1 = van_alpha;
    let z_eval_t_s2 = E::Scalar::ONE;
    let z_eval_t_s3 = van_alpha;
    let z_eval_t_s4 = van_zeta_inv * van_alpha;
    let z_eval_t = z_eval_t_s4 * van_zeta;

    let comm_g = self.comm_g;
    let comm_h = self.comm_h;
    let comm_s = self.comm_s;
    let comm_d = self.comm_d;

    let comm_quot_m = self.comm_quot_m;
    let comm_quot_l = self.comm_quot_l;

    let beta_2 = beta * beta;
    let beta_3 = beta_2 * beta;

    let mut f = comm_g * z_eval_t_s1;

    f = f + comm_h * (beta * z_eval_t_s2);

    f = f + comm_s * (beta_2 * z_eval_t_s3);

    f = f + comm_d * (beta_3 * z_eval_t_s4);

    f = f + comm_quot_m * (-z_eval_t);

    let scalar = z_eval_t_s1 * g_star_eval
      + beta * z_eval_t_s2 * h_star_eval
      + beta_2 * z_eval_t_s3 * s_star_eval
      + beta_3 * z_eval_t_s4 * d_star_eval;

    f = f + Commitment::new(DlogGroup::group(&vk.G)) * (-scalar);

    let ll = f + comm_quot_l * z;
    let ll = ll.into_inner();
    let lr = DlogGroup::group(&vk.H);

    let rl = comm_quot_l.into_inner();
    let rr = DlogGroup::group(&vk.tau_H);

    let pairing_l = E::GE::pairing(&ll, &lr);
    let pairing_r = E::GE::pairing(&rl, &rr);

    assert!(pairing_l == pairing_r);
  }
}
