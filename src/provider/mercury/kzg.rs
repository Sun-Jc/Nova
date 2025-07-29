#![allow(dead_code)]
//! # KZG Batch Proof
//! [BCGH20]
use std::fmt::Debug;

use ff::{Field, PrimeField};
use rand_core::OsRng;
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
pub struct EvaluationEngine<E: Engine>
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

  //   pub z_poly_t_s1: UniPoly<E::Scalar>,
  //   pub z_poly_t_s2: UniPoly<E::Scalar>,
  //   pub z_poly_t_s3: UniPoly<E::Scalar>,
  //   pub z_poly_t_s4: UniPoly<E::Scalar>,
  pub m_poly: UniPoly<E::Scalar>,
  pub quot_m_poly: UniPoly<E::Scalar>,
  pub comm_qm: Commitment<E>,

  pub mz_poly: UniPoly<E::Scalar>,

  pub l_poly: UniPoly<E::Scalar>,
  pub quot_l_poly: UniPoly<E::Scalar>,
  pub comm_ql: Commitment<E>,
}

impl<E: Engine> EvaluationEngine<E>
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

    // let z_poly_t_s1 = UniPoly::make_vanishing_poly(&[*alpha]);
    // let z_poly_t_s2 = UniPoly::make_vanishing_poly(&[E::Scalar::ONE]);
    // let z_poly_t_s3 = UniPoly::make_vanishing_poly(&[*alpha]);
    // let z_poly_t_s4 = UniPoly::make_vanishing_poly(&[*zeta_inv, *alpha]);

    {
      // Check Interpolation
      let polys = vec![
        g_poly, g_poly, h_poly, h_poly, h_poly, s_poly, s_poly, d_poly,
      ];

      let xs = xss.iter().flatten().collect::<Vec<_>>();

      let evals1 = polys
        .into_par_iter()
        .zip(xs.clone())
        .map(|(poly, x)| poly.evaluate(x))
        .collect::<Vec<_>>();

      let polys = vec![
        &g_star, &g_star, &h_star, &h_star, &h_star, &s_star, &s_star, &d_star,
      ];

      let evals2 = polys
        .into_par_iter()
        .zip(xs)
        .map(|(poly, x)| poly.evaluate(x))
        .collect::<Vec<_>>();

      assert_eq!(evals1, evals2);

      let polys = vec![
        &g_minus_g_star,
        &g_minus_g_star,
        &h_minus_h_star,
        &h_minus_h_star,
        &s_minus_s_star,
        &s_minus_s_star,
        &d_minus_d_star,
      ];

      let evals3 = polys
        .into_par_iter()
        .map(|poly| poly.evaluate(zeta))
        .collect::<Vec<_>>();

      for e in evals3 {
        assert_eq!(e, E::Scalar::ZERO);
      }
    }

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

      //   z_poly_t_s1,
      //   z_poly_t_s2,
      //   z_poly_t_s3,
      //   z_poly_t_s4,
      comm_g: Default::default(),
      comm_h: Default::default(),
      comm_s: Default::default(),
      comm_d: Default::default(),
      comm_qm: Default::default(),
      comm_ql: Default::default(),

      z: Default::default(),
      beta: Default::default(),
      m_poly: empty_poly(),
      quot_m_poly: empty_poly(),
      mz_poly: empty_poly(),
      l_poly: empty_poly(),
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
      self.g_poly.coeffs.clone(),
      self.h_poly.coeffs.clone(),
      self.s_poly.coeffs.clone(),
      self.d_poly.coeffs.clone(),
    ];

    let comms = E::CE::batch_commit(ck, &vs, &[E::Scalar::ZERO; 4]);
    let [comm_g, comm_h, comm_s, comm_d] = comms.try_into().unwrap();

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
          rhs.into_iter().for_each(|c| {
            let x = poly.mul_by_deg_one_polynomial(&c);

            {
              let r = E::Scalar::random(OsRng);
              let eval1 = poly.evaluate(&r) * (r + c);
              let eval2 = x.evaluate(&r);

              assert_eq!(eval1, eval2);
            }

            poly = x;
          });

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

      {
        let r = E::Scalar::random(OsRng);
        let v1 = g_prod.evaluate(&r);
        let v2 = h_prod.evaluate(&r);
        let v3 = s_prod.evaluate(&r);
        let v4 = d_prod.evaluate(&r);

        let v = m_poly.evaluate(&r);

        assert_eq!(v, v1 + v2 + v3 + v4);
      }

      m_poly
    };

    self.m_poly = m_poly.clone();

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

    {
      let mut res = E::Scalar::ZERO;

      let eval1 = self.g_minus_g_star.evaluate(z) * (*z - self.alpha);
      let eval2 = g_poly_z.evaluate(z);

      res += eval1;

      assert_eq!(eval1, eval2);

      let eval1 = self.h_minus_h_star.evaluate(z) * self.beta;
      let eval2 = h_poly_z.evaluate(z);

      res += eval1;

      assert_eq!(eval1, eval2);

      let eval1 = self.s_minus_s_star.evaluate(z) * beta_2 * (*z - self.alpha);
      let eval2 = s_poly_z.evaluate(z);

      res += eval1;

      assert_eq!(eval1, eval2);

      let eval1 =
        self.d_minus_d_star.evaluate(z) * beta_3 * (*z - self.alpha) * (*z - self.zeta_inv);
      let eval2 = d_poly_z.evaluate(z);

      res += eval1;

      assert_eq!(eval1, eval2);

      let mz_poly = g_poly_z
        .clone()
        .into_add_by_polynomial(&h_poly_z)
        .into_add_by_polynomial(&s_poly_z)
        .into_add_by_polynomial(&d_poly_z);

      let eval1 = mz_poly.evaluate(&z);
      let eval2 = self.m_poly.evaluate(&z);

      assert_eq!(eval1, res);

      assert_eq!(res, eval2);
    }

    // TODO: can be parallelized
    let mz_poly = g_poly_z
      .into_add_by_polynomial(&h_poly_z)
      .into_add_by_polynomial(&s_poly_z)
      .into_add_by_polynomial(&d_poly_z);

    self.mz_poly = mz_poly.clone();

    {
      let eval1 = mz_poly.evaluate(&z);
      let eval2 = self.m_poly.evaluate(&z);

      assert_eq!(eval1, eval2);
    }

    {
      let r = E::Scalar::random(OsRng);

      let lhs = mz_poly.evaluate(&r);

      let rhs1 = (*z - self.alpha) * (self.g_poly.evaluate(&r) - self.g_star.evaluate(&z));
      let rhs2 = self.h_poly.evaluate(&r) - self.h_star.evaluate(&z);
      let rhs3 = (*z - self.alpha) * (self.s_poly.evaluate(&r) - self.s_star.evaluate(&z));
      let rhs4 = (*z - self.alpha)
        * (*z - self.zeta_inv)
        * (self.d_poly.evaluate(&r) - self.d_star.evaluate(&z));

      let rhs = rhs1 + rhs2 * self.beta + rhs3 * beta_2 + rhs4 * beta_3;

      assert_eq!(lhs, rhs);
    }

    let l_poly = {
      let mut q_m = self.quot_m_poly.clone();
      q_m.scale(&t_eval_at_z);

      mz_poly.clone().into_sub_by_polynomial(&q_m)
    };

    self.l_poly = l_poly.clone();

    {
      let r = E::Scalar::random(OsRng);
      assert_eq!(
        l_poly.evaluate(&r),
        mz_poly.evaluate(&r) - self.quot_m_poly.evaluate(&r) * t_eval_at_z
      );

      assert_eq!(l_poly.evaluate(&z), E::Scalar::ZERO);
      assert_eq!(
        mz_poly.evaluate(&z) - self.quot_m_poly.evaluate(&z) * t_eval_at_z,
        E::Scalar::ZERO
      );
    }

    let (quot_l_poly, rem) = l_poly.into_div_by_deg_one_polynomial(z);
    assert_eq!(rem, E::Scalar::ZERO);

    self.quot_l_poly = quot_l_poly;
  }

  pub fn commit_quot_l_7(&mut self, ck: &CommitmentKey<E>) {
    let comm_quot_l = E::CE::commit(ck, &self.quot_l_poly.coeffs, &E::Scalar::ZERO);
    self.comm_ql = comm_quot_l;
  }

  pub fn check_7(&mut self, ck: &CommitmentKey<E>, vk: &VerifierKey<E>) {
    {
      let r = E::Scalar::random(OsRng);
      let ql_eval = self.quot_l_poly.evaluate(&r);
      let l_eval = if self.l_poly.coeffs.len() == 1 {
        self.l_poly.coeffs[0]
      } else {
        self.l_poly.evaluate(&r)
      };

      assert_eq!(ql_eval * (r - self.z), l_eval);
    }

    {
      let comm_l = E::CE::commit(ck, &self.l_poly.coeffs, &E::Scalar::ZERO);
      let z = self.z;
      let comm_ql = self.comm_ql;
      let ll = comm_l + comm_ql * z;
      let lr = DlogGroup::group(&vk.H);
      let rl = comm_ql;
      let rr = DlogGroup::group(&vk.tau_H);

      let pairing_l = E::GE::pairing(&ll.into_inner(), &lr);
      let pairing_r = E::GE::pairing(&rl.into_inner(), &rr);

      assert!(pairing_l == pairing_r);

      // q_m(tau) * (tau-z) = m(tau)
      // c_qm * (tau - z_0) (tau - z_1)...  == c_m
      // c_l + z * c_ql == c_ql * tau
    }
  }

  pub fn into_proof_8(self) -> BatchKZGProof<E> {
    {
      assert_ne!(self.comm_g, Default::default());
      assert_ne!(self.comm_h, Default::default());
      assert_ne!(self.comm_s, Default::default());
      assert_ne!(self.comm_d, Default::default());
      assert_ne!(self.comm_qm, Default::default());
      assert_ne!(self.comm_ql, Default::default());
    }
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
  <E::GE as PairingGroup>::GT: Debug,
{
  pub fn verify(
    &self,
    vk: &VerifierKey<E>,
    eval: &BatchEvaluation<E::Scalar>,
  ) -> (<E::GE as PairingGroup>::GT, <E::GE as PairingGroup>::GT) {
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

    assert_eq!(ll, rl);

    // assert_eq!(pairing_l, pairing_r);
    (pairing_l, pairing_r)
  }
}
