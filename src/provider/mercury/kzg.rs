#![allow(dead_code)]

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
pub struct PolyCommitment<E: Engine>
where
  E::GE: PairingGroup,
{
  pub comm_g: Commitment<E>,
  pub comm_h: Commitment<E>,
  pub comm_s: Commitment<E>,
  pub comm_d: Commitment<E>,
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
  pub beta: E::Scalar,
  pub z: E::Scalar,

  pub comm_quot_m: Commitment<E>,
  pub comm_quot_l: Commitment<E>,
}

impl<E: Engine> PolyCommitment<E>
where
  E::GE: PairingGroup,
  E: Engine<CE = CommitmentEngine<E>>,
{
  pub fn commit_to(
    ck: &CommitmentKey<E>,
    g_poly: &UniPoly<E::Scalar>,
    h_poly: &UniPoly<E::Scalar>,
    s_poly: &UniPoly<E::Scalar>,
    d_poly: &UniPoly<E::Scalar>,
  ) -> Self {
    let vs = vec![
      g_poly.coeffs.clone(),
      h_poly.coeffs.clone(),
      s_poly.coeffs.clone(),
      d_poly.coeffs.clone(),
    ];

    let comms = E::CE::batch_commit(ck, &vs, &[E::Scalar::ONE; 4]);
    let [comm_g, comm_h, comm_s, comm_d] = comms.try_into().unwrap();

    Self {
      comm_g,
      comm_h,
      comm_s,
      comm_d,
    }
  }
}

#[derive(Debug)]
pub struct BatchKZGWitness<Scalar: PrimeField> {
  pub zeta: Scalar,
  pub zeta_inv: Scalar,
  pub alpha: Scalar,
  pub beta: Scalar,

  pub z: Scalar,

  pub g_poly: UniPoly<Scalar>,
  pub h_poly: UniPoly<Scalar>,
  pub s_poly: UniPoly<Scalar>,
  pub d_poly: UniPoly<Scalar>,

  pub g_star: UniPoly<Scalar>,
  pub h_star: UniPoly<Scalar>,
  pub s_star: UniPoly<Scalar>,
  pub d_star: UniPoly<Scalar>,

  pub g_minus_g_star: UniPoly<Scalar>,
  pub h_minus_h_star: UniPoly<Scalar>,
  pub s_minus_s_star: UniPoly<Scalar>,
  pub d_minus_d_star: UniPoly<Scalar>,

  pub z_poly_t_s1: UniPoly<Scalar>,
  pub z_poly_t_s2: UniPoly<Scalar>,
  pub z_poly_t_s3: UniPoly<Scalar>,
  pub z_poly_t_s4: UniPoly<Scalar>,

  pub m_poly: UniPoly<Scalar>,
  pub quot_m_poly: UniPoly<Scalar>,

  pub mz_poly: UniPoly<Scalar>,
  pub l_poly: UniPoly<Scalar>,
  pub quot_l_poly: UniPoly<Scalar>,
}

impl<Scalar: PrimeField> BatchKZGWitness<Scalar> {
  pub fn init_eval(
    alpha: &Scalar,
    zeta: &Scalar,
    g_poly: &UniPoly<Scalar>,
    h_poly: &UniPoly<Scalar>,
    s_poly: &UniPoly<Scalar>,
    d_poly: &UniPoly<Scalar>,
  ) -> (Self, BatchEvaluation<Scalar>) {
    let zeta_inv = &zeta.invert().unwrap();

    let polys = vec![
      g_poly, g_poly, h_poly, h_poly, h_poly, s_poly, s_poly, d_poly,
    ];
    let xs = [zeta, zeta_inv, zeta, zeta_inv, alpha, zeta, zeta_inv, zeta];

    let evals = polys
      .into_par_iter()
      .zip(xs)
      .map(|(poly, x)| poly.evaluate(x))
      .collect::<Vec<_>>();

    let [g_zeta, g_zeta_inv, h_zeta, h_zeta_inv, h_alpha, s_zeta, s_zeta_inv, d_zeta] =
      evals.try_into().unwrap();

    let xss = vec![
      vec![*zeta, *zeta_inv],
      vec![*zeta, *zeta_inv, *alpha],
      vec![*zeta, *zeta_inv],
      vec![*zeta],
    ];
    let evals = vec![
      vec![g_zeta, g_zeta_inv],
      vec![h_zeta, h_zeta_inv, h_alpha],
      vec![s_zeta, s_zeta_inv],
      vec![d_zeta],
    ];

    let polys = xss
      .into_par_iter()
      .zip(evals)
      .map(|(xs, evals)| UniPoly::from_evals_with_xs(&xs, &evals))
      .collect::<Vec<_>>();

    let [g_star, h_star, s_star, d_star]: [UniPoly<Scalar>; 4] = polys.try_into().unwrap();

    let diff_polys = [
      g_poly.clone(),
      h_poly.clone(),
      s_poly.clone(),
      d_poly.clone(),
    ]
    .into_par_iter()
    .zip([&g_star, &h_star, &s_star, &d_star])
    .map(|(poly, poly_star)| poly.into_sub_by_polynomial(poly_star))
    .collect::<Vec<_>>();

    let [g_minus_g_star, h_minus_h_star, s_minus_s_star, d_minus_d_star]: [UniPoly<Scalar>; 4] =
      diff_polys.try_into().unwrap();

    let z_poly_t_s1 = UniPoly::make_vanishing_poly(&[*zeta, *zeta_inv]);
    let z_poly_t_s2 = UniPoly::make_vanishing_poly(&[*zeta, *zeta_inv, *alpha]);
    let z_poly_t_s3 = UniPoly::make_vanishing_poly(&[*zeta, *zeta_inv]);
    let z_poly_t_s4 = UniPoly::make_vanishing_poly(&[*zeta]);

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

      z_poly_t_s1,
      z_poly_t_s2,
      z_poly_t_s3,
      z_poly_t_s4,

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

  pub fn open_phase_1(&mut self, beta: &Scalar) {
    // g_diff * (t/s1) (alpha)
    // h_diff * (t/s2) (1)
    // s_diff * (t/s3) (alpha)
    // d_diff * (t/s4) (1/zeta, alpha)

    let rhs = vec![
      vec![self.alpha],
      vec![self.alpha],
      vec![self.zeta_inv, self.alpha],
    ];

    let lhs = vec![
      &self.g_minus_g_star,
      &self.s_minus_s_star,
      &self.d_minus_d_star,
    ];

    let [g_prod, s_prod, d_prod] = lhs
      .into_par_iter()
      .zip(rhs)
      .map(|(poly, rhs)| {
        let mut poly = poly.clone();
        rhs
          .into_iter()
          .for_each(|alpha| poly = poly.mul_by_deg_one_polynomial(&-alpha));

        poly
      })
      .collect::<Vec<_>>()
      .try_into()
      .unwrap();

    let mut g_prod = g_prod;
    let mut h_prod = self.h_minus_h_star.clone();
    let mut s_prod = s_prod;
    let mut d_prod = d_prod;

    let mut beta_acc = *beta;

    g_prod.scale(&beta_acc);
    beta_acc *= *beta;

    h_prod.scale(&beta_acc);
    beta_acc *= *beta;

    s_prod.scale(&beta_acc);
    beta_acc *= *beta;

    d_prod.scale(&beta_acc);

    let mut m_poly = g_prod;
    m_poly = m_poly.into_add_by_polynomial(&h_prod);
    m_poly = m_poly.into_add_by_polynomial(&s_prod);
    m_poly = m_poly.into_add_by_polynomial(&d_prod);

    self.m_poly = m_poly.clone();

    let (q_m_poly, rem) = m_poly.into_div_by_deg_one_polynomial(&self.alpha);
    assert_eq!(rem, Scalar::ZERO);
    let (q_m_poly, rem) = q_m_poly.into_div_by_deg_one_polynomial(&self.zeta);
    assert_eq!(rem, Scalar::ZERO);
    let (q_m_poly, rem) = q_m_poly.into_div_by_deg_one_polynomial(&self.zeta_inv);
    assert_eq!(rem, Scalar::ZERO);

    self.quot_m_poly = q_m_poly;
  }

  pub fn open_phase_2(&mut self, z: &Scalar) {
    // L(X) = m_z(X) - Z_T(z) \cdot q_m(X)

    let vanishing_terms = vec![
      vec![self.alpha],
      vec![],
      vec![self.alpha],
      vec![self.zeta_inv, self.alpha],
    ];

    let [g_poly_res, h_poly_res, s_poly_res, d_poly_res] =
      vec![&self.g_poly, &self.h_poly, &self.s_poly, &self.d_poly]
        .into_par_iter()
        .zip(vanishing_terms)
        .map(|(poly, terms)| {
          let eval = poly.evaluate(z);
          let mut z_eval = Scalar::ONE;
          terms.into_iter().for_each(|term| {
            z_eval *= term;
          });

          let mut res = self
            .g_poly
            .clone()
            .into_sub_by_polynomial(&UniPoly { coeffs: vec![eval] });
          res.scale(&z_eval);

          res
        })
        .collect::<Vec<_>>()
        .try_into()
        .unwrap();

    let mut h_poly_res = h_poly_res;
    let mut s_poly_res = s_poly_res;
    let mut d_poly_res = d_poly_res;

    let mut beta_acc = self.beta;

    h_poly_res.scale(&beta_acc);
    beta_acc *= self.beta;

    h_poly_res.scale(&beta_acc);
    beta_acc *= self.beta;

    s_poly_res.scale(&beta_acc);
    beta_acc *= self.beta;

    d_poly_res.scale(&beta_acc);

    let mut mz_poly = g_poly_res;
    mz_poly = mz_poly.into_add_by_polynomial(&h_poly_res);
    mz_poly = mz_poly.into_add_by_polynomial(&s_poly_res);
    mz_poly = mz_poly.into_add_by_polynomial(&d_poly_res);

    self.mz_poly = mz_poly.clone();

    let (quot_l_poly, rem) = mz_poly.into_div_by_deg_one_polynomial(z);
    assert_eq!(rem, Scalar::ZERO);

    self.quot_l_poly = quot_l_poly;
  }

  pub fn finalize<E>(
    self,
    ck: &<E::CE as CommitmentEngineTrait<E>>::CommitmentKey,
  ) -> BatchKZGProof<E>
  where
    E::GE: PairingGroup,
    E: Engine<CE = CommitmentEngine<E>, Scalar = Scalar>,
  {
    let vs = vec![self.quot_m_poly.coeffs, self.quot_l_poly.coeffs];

    let comms = E::CE::batch_commit(ck, &vs, &[Scalar::ONE; 2]);
    let [comm_quot_m, comm_quot_l] = comms.try_into().unwrap();

    BatchKZGProof {
      beta: self.beta,
      z: self.z,
      comm_quot_m,
      comm_quot_l,
    }
  }
}

impl<E: Engine> BatchKZGProof<E>
where
  E: Engine<CE = CommitmentEngine<E>>,
  E::GE: PairingGroup,
{
  pub fn verify(
    &self,
    vk: &VerifierKey<E>,
    comms: &PolyCommitment<E>,
    eval: &BatchEvaluation<E::Scalar>,
  ) {
    let zeta = eval.zeta;
    let zeta_inv = zeta.invert().unwrap();
    let alpha = eval.alpha;
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

    let comm_g = comms.comm_g;
    let comm_h = comms.comm_h;
    let comm_s = comms.comm_s;
    let comm_d = comms.comm_d;

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

    f = f + Commitment::new(<E::GE as DlogGroup>::gen()) * scalar;

    let ll = f + comm_quot_l * z;
    let ll = ll.into_inner();
    let lr = <E::GE as PairingGroup>::G2::gen();

    let rl = comm_quot_l.into_inner();
    let rr = DlogGroup::group(&vk.tau_H);

    let pairing_l = E::GE::pairing(&ll, &lr);
    let pairing_r = E::GE::pairing(&rl, &rr);

    assert!(pairing_l == pairing_r);
  }
}
