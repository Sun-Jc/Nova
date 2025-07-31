#![allow(dead_code)]

use std::marker::PhantomData;

use crate::{
  errors::NovaError,
  provider::{
    hyperkzg::{Commitment, CommitmentEngine},
    mercury::{
      ipa::{make_s_polynomial, IPAWitness, InputPolynomials},
      kzg::{BatchEvaluation, BatchKZGProof, EvaluationProcess},
      split_polynomial::{
        divide_polynomial_by_x_b_alpha, divide_polynomial_by_x_b_alpha2, split_polynomial,
      },
    },
    traits::{DlogGroup, PairingGroup},
  },
  spartan::polys::{eq::EqPolynomial, univariate::UniPoly},
  traits::{
    commitment::CommitmentEngineTrait, evaluation::EvaluationEngineTrait, Engine,
    TranscriptEngineTrait,
  },
};

use ff::{Field, PrimeField};
use rand_core::OsRng;
use rayon::iter::{
  IndexedParallelIterator, IntoParallelIterator, IntoParallelRefMutIterator, ParallelIterator,
};
use serde::{Deserialize, Serialize};

/// Alias to points on G1 that are in preprocessed form
type G1Affine<E> = <<E as Engine>::GE as DlogGroup>::AffineGroupElement;

/// Alias to points on G1 that are in preprocessed form
type G2Affine<E> = <<<E as Engine>::GE as PairingGroup>::G2 as DlogGroup>::AffineGroupElement;

/// Provides an implementation of generators for proving evaluations
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(bound = "")]
pub struct ProverKey<E: Engine> {
  _p: PhantomData<E>,
}

/// Provides an implementation of a polynomial evaluation engine using KZG
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct EvaluationEngine<E: Engine> {
  _p: PhantomData<E>,
}

/// A verifier key
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(bound = "")]
pub struct VerifierKey<E: Engine>
where
  E::GE: PairingGroup,
{
  pub(crate) g: G1Affine<E>,
  pub(crate) h: G2Affine<E>,
  pub(crate) tau_h: G2Affine<E>,
}

/// Provides an implementation of a polynomial evaluation argument
/// 8 G + 8 F
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(bound = "")]
pub struct EvaluationArgument<E: Engine>
where
  E::GE: PairingGroup,
{
  comm_h: Commitment<E>,
  comm_g: Commitment<E>,
  comm_q: Commitment<E>,
  comm_s: Commitment<E>,
  comm_d: Commitment<E>,

  comm_quot_m: Commitment<E>,
  comm_quot_l: Commitment<E>,

  comm_quot_f_x_zeta: Commitment<E>,

  g_zeta: E::Scalar,
  g_zeta_inv: E::Scalar,

  h_zeta: E::Scalar,
  h_zeta_inv: E::Scalar,
  h_alpha: E::Scalar,

  s_zeta: E::Scalar,
  s_zeta_inv: E::Scalar,

  d_zeta: E::Scalar,

  // ! TODO: remove these
  zeta: E::Scalar,
  alpha: E::Scalar,
  beta: E::Scalar,
  z: E::Scalar,
  gamma: E::Scalar,
}

// P_u(X) = prod(uiX^{2^i} + 1 - u_i)
pub fn make_pu_poly<Scalar: PrimeField>(u: &[Scalar]) -> (UniPoly<Scalar>, Vec<Scalar>) {
  // let log_n = u.len();

  // let mv_eq_poly = EqPolynomial::new(u.to_vec());
  // let xs = 0..1 << log_n as usize;

  // let to_bits = |i: usize| {
  //   (0..log_n as usize)
  //     .map(|j| {
  //       if (i >> j) & 1 == 1 {
  //         Scalar::ONE
  //       } else {
  //         Scalar::ZERO
  //       }
  //     })
  //     .rev()
  //     .collect::<Vec<_>>()
  // };

  // let evals = xs
  //   .into_par_iter()
  //   .map(|i| mv_eq_poly.evaluate(&to_bits(i)))
  //   .collect::<Vec<_>>();

  // (
  //   UniPoly {
  //     coeffs: evals.clone(),
  //   },
  //   evals,
  // )

  let mut evals = EqPolynomial::new(u.to_vec()).evals();
  // transpose evals
  let n_evals = evals.clone();
  let n = n_evals.len().isqrt();
  for i in 0..n {
    for j in 0..n {
      evals[i * n + j] = n_evals[j * n + i];
    }
  }

  (
    UniPoly {
      coeffs: evals.clone(),
    },
    evals,
  )
}

impl<E: Engine> EvaluationEngineTrait<E> for EvaluationEngine<E>
where
  E: Engine<CE = CommitmentEngine<E>>,
  E::GE: PairingGroup,
{
  type ProverKey = ProverKey<E>;
  type VerifierKey = VerifierKey<E>;
  type EvaluationArgument = EvaluationArgument<E>;

  /// A method to perform any additional setup needed to produce proofs of evaluations
  fn setup(
    ck: &<<E as Engine>::CE as CommitmentEngineTrait<E>>::CommitmentKey,
  ) -> (Self::ProverKey, Self::VerifierKey) {
    (
      ProverKey { _p: PhantomData },
      VerifierKey {
        g: E::GE::gen().affine(),
        h: <<E::GE as PairingGroup>::G2 as DlogGroup>::gen().affine(),
        tau_h: ck.tau_H,
      },
    )
  }

  /// A method to prove the evaluation of a multilinear polynomial
  fn prove(
    ck: &<<E as Engine>::CE as CommitmentEngineTrait<E>>::CommitmentKey,
    _pk: &Self::ProverKey,
    transcript: &mut E::TE,
    comm: &<<E as Engine>::CE as CommitmentEngineTrait<E>>::Commitment,
    poly: &[E::Scalar],
    point: &[E::Scalar],
    _eval: &E::Scalar,
  ) -> Result<Self::EvaluationArgument, NovaError> {
    let comm_f = comm;

    let mut point = point.to_vec();
    // point.reverse();

    let (f_poly, log_n, point_l, point_r) = {
      // Round 0: Prepare
      let mut f = UniPoly {
        coeffs: poly.to_vec(),
      };

      // transpose f_poly
      // let n = f.coeffs.len().isqrt();
      // for i in 0..n {
      //   for j in 0..n {
      //     f.coeffs[i * n + j] = poly[j * n + i];
      //   }
      // }

      let mut log_n = f.log_n();
      if log_n % 2 == 1 {
        log_n += 1;
      }

      f.raise(1 << log_n as usize);

      let (point_l, point_r) = point.split_at(log_n as usize / 2);

      (f, log_n, point_l.to_vec(), point_r.to_vec())
    };

    // let mut point_l = point_l.to_vec();
    // point_l.reverse();
    // let mut point_r = point_r.to_vec();
    // point_r.reverse();

    // let point_l = point_r.to_vec();
    // let point_r = point_l.to_vec();

    let b = 1 << (log_n / 2);

    let (p_poly1, eq_evals1) = make_pu_poly(&point_l);
    let (p_poly2, eq_evals2) = make_pu_poly(&point_r);

    {
      dbg!(&point_l);
      dbg!(&point_r);
      dbg!(&p_poly1.coeffs);
      dbg!(&p_poly2.coeffs);
      let r = E::Scalar::from(12);
      dbg!(&p_poly1.evaluate(&r));
      dbg!(&p_poly2.evaluate(&r));
    }

    {
      let eval_pu = |u: &[E::Scalar], x: E::Scalar| {
        let mut res = E::Scalar::ONE;
        for (i, &ui) in u.iter().rev().enumerate() {
          res *= x.pow_vartime(&[1 << i as u64]) * ui + E::Scalar::ONE - ui;
        }
        res
      };

      let r = E::Scalar::random(OsRng);
      // assert_eq!(eval_pu(&point_l, r), p_poly1.evaluate(&r));

      // assert_eq!(eval_pu(&point_r, r), p_poly2.evaluate(&r));
    }

    // Get h_poly
    // let h_poly = {
    //   let mut sub_fs = f_poly
    //     .coeffs
    //     .chunks(b)
    //     .map(|c| UniPoly { coeffs: c.to_vec() })
    //     .collect::<Vec<_>>();

    //   // let mut sub_fs = split_polynomial(&f_poly, log_n);

    //   // sub_fs.par_iter_mut().enumerate().for_each(|(i, sub_f)| {
    //   //   let lhs = eq_evals1[i];
    //   //   sub_f.scale(&lhs);
    //   // });

    //   for i in 0..sub_fs.len() {
    //     let lhs = eq_evals1[i];
    //     sub_fs[i].scale(&lhs);
    //   }

    //   // TODO: parallelize this
    //   let mut h_poly = sub_fs[0].clone();
    //   for sub_f in sub_fs.iter().skip(1) {
    //     h_poly = h_poly.into_add_by_polynomial(sub_f);
    //   }
    //   h_poly
    // };

    let h_poly = {
      let mut coeffs = vec![E::Scalar::ZERO; b];
      for j in 0..b {
        for i in 0..b {
          coeffs[j] += f_poly.coeffs[j * b + i] * eq_evals1[i];
        }
      }
      UniPoly { coeffs }
    };

    {
      dbg!(&h_poly.coeffs);
      // todo!();
    }

    {
      // check h_poly
      let mut e2 = E::Scalar::ZERO;
      for (i, h) in h_poly.coeffs.iter().enumerate() {
        e2 += *h * eq_evals2[i];
      }

      dbg!(&eq_evals1);
      dbg!(&eq_evals2);

      dbg!(&h_poly.coeffs);

      // assert_eq!(e2, *_eval);
      dbg!(&e2);
    }

    let comm_h = E::CE::commit(ck, &h_poly.coeffs, &E::Scalar::ZERO);

    transcript.absorb(b"f", &[*comm_f, comm_h].to_vec().as_slice());
    // let alpha = transcript.squeeze(b"alpha")?;
    let alpha = E::Scalar::from(1);

    // Get q(X) and g(X)
    let (q_poly, g_poly) = divide_polynomial_by_x_b_alpha(&f_poly, log_n, &alpha);

    {
      // Check q,g
      // f(X) = (X^b - alpha) * q(X) + g(X)
      let r = E::Scalar::random(OsRng);
      let lhs = f_poly.evaluate(&r);
      let rhs =
        (r.pow([1 << (log_n / 2) as u64]) - alpha) * q_poly.evaluate(&r) + g_poly.evaluate(&r);
      assert_eq!(lhs, rhs);

      dbg!(&f_poly.coeffs);
      dbg!(&q_poly.coeffs);
      dbg!(&g_poly.coeffs);
    }

    {
      // gx is col
      let split_res = split_polynomial(&f_poly, log_n);
      for i in 0..b {
        assert_eq!(split_res[i].evaluate(&alpha), g_poly.coeffs[i]);
      }
    }

    {
      // let mut g_poly_coeffs = g_poly.coeffs.clone();
      // let split_res = f_poly
      //   .coeffs
      //   .chunks(b)
      //   .map(|c| UniPoly { coeffs: c.to_vec() })
      //   .collect::<Vec<_>>();
      // for i in 0..b {
      //   g_poly_coeffs[i] = split_res[i].evaluate(&alpha);
      // }

      // let g_poly = UniPoly {
      //   coeffs: g_poly_coeffs,
      // };

      // check g and h_alpha
      assert_eq!(g_poly.coeffs.len(), b);

      let res = eq_evals1
        .iter()
        .zip(g_poly.coeffs.iter())
        .map(|a| *a.0 * *a.1)
        .sum::<E::Scalar>();

      dbg!(&g_poly.coeffs.len());
      dbg!(&h_poly.evaluate(&alpha));
      assert_eq!(res, h_poly.evaluate(&alpha));
    }

    // let mut q_poly = q_poly;
    // q_poly.trim();

    dbg!(&q_poly.coeffs.len());
    dbg!(&g_poly.coeffs.len());

    // assert_eq!(q_poly.log_n(), log_n / 2);
    // assert_eq!(g_poly.log_n(), log_n / 2);

    let comm_q = E::CE::commit(ck, &q_poly.coeffs, &E::Scalar::ZERO);
    let comm_g = E::CE::commit(ck, &g_poly.coeffs, &E::Scalar::ZERO);

    transcript.absorb(b"g", &[comm_g, comm_q].to_vec().as_slice());

    let h_alpha = h_poly.evaluate(&alpha);

    transcript.absorb(b"h_alpha", &h_alpha);

    let gamma = transcript.squeeze(b"gamma")?;

    // Get s(X) for ipa
    let s_poly = {
      dbg!(&g_poly.coeffs.len());
      dbg!(&h_poly.coeffs.len());
      dbg!(&p_poly1.coeffs.len());
      dbg!(&p_poly2.coeffs.len());

      // let p_poly1 = UniPoly { coeffs: eq_evals1 };

      let left_polys = vec![p_poly1.clone(), p_poly2.clone()];
      let right_polys = vec![g_poly.clone(), h_poly.clone()];

      assert_eq!(p_poly1.coeffs.len(), g_poly.coeffs.len());

      {
        let wit = IPAWitness::generate(
          log_n / 2,
          InputPolynomials {
            lhs: left_polys.clone(),
            rhs: right_polys.clone(),
          },
          &gamma,
        );

        dbg!(log_n);
        dbg!(&wit.s_polynomial.coeffs.len());

        assert_eq!(wit.products[0], h_alpha);
        assert_eq!(wit.products[1], *_eval);
      }

      make_s_polynomial(left_polys, right_polys, log_n / 2, &gamma)
    };

    {
      // Check s
      let r = E::Scalar::from(12);
      let r_inv = r.invert().unwrap();
      let lhs_l = g_poly.evaluate(&r) * p_poly1.evaluate(&r_inv)
        + g_poly.evaluate(&r_inv) * p_poly1.evaluate(&r);
      let lhs_r = h_poly.evaluate(&r) * p_poly2.evaluate(&r_inv)
        + h_poly.evaluate(&r_inv) * p_poly2.evaluate(&r);

      let lhs = lhs_l + lhs_r * gamma;

      let rhs = h_poly.evaluate(&alpha) + *_eval * gamma;
      let rhs = rhs + rhs;
      let rhs = rhs + r * s_poly.evaluate(&r) + r_inv * s_poly.evaluate(&r_inv);

      assert_eq!(lhs, rhs);
    }

    dbg!(&s_poly.coeffs.len());

    // Get d(X) for degree check
    let d_poly = {
      let mut d_poly = g_poly.clone();
      assert_eq!(d_poly.coeffs.len(), 1 << (log_n / 2));
      d_poly.coeffs.reverse();
      d_poly
    };

    {
      // check d_poly
      let r = E::Scalar::from(12);
      let r_inv = r.invert().unwrap();
      let rhs = d_poly.evaluate(&r);
      // assert_eq!(1_u64 << (log_n / 2) - 1, g_poly.coeffs.len() as u64 - 1);
      let lhs = g_poly.evaluate(&r_inv) * r.pow_vartime(&[(1_u64 << (log_n / 2)) - 1]);
      dbg!(&g_poly.coeffs);
      dbg!(&d_poly.coeffs);
      assert_eq!(lhs, rhs);
    }

    let comm_s = E::CE::commit(ck, &s_poly.coeffs, &E::Scalar::ZERO);
    let comm_d = E::CE::commit(ck, &d_poly.coeffs, &E::Scalar::ZERO);

    transcript.absorb(b"d", &[comm_s, comm_d].to_vec().as_slice());

    let zeta = transcript.squeeze(b"zeta")?;

    let (process, batch_eval) =
      EvaluationProcess::<E>::init_eval_0(&alpha, &zeta, &g_poly, &h_poly, &s_poly, &d_poly);

    let (quot_poly, rem) = {
      let g_zeta = batch_eval.g_zeta;

      // let zeta_b = zeta.pow_vartime(&[1_u64 << (log_n / 2)]);

      let zeta_b = zeta.pow(&[b as u64]);

      dbg!(&[1_u64 << (log_n / 2)]);
      let zeta_b_alpha = zeta_b - alpha;

      let mut new_q = q_poly.clone();
      new_q.scale(&zeta_b_alpha);

      // f_poly
      //   .clone()
      //   .into_sub_by_polynomial(&new_q)
      //   .into_sub_by_polynomial(&UniPoly {
      //     coeffs: vec![g_zeta],
      //   })
      //   .into_div_by_deg_one_polynomial(&zeta)

      let mut tmp = f_poly.clone().into_sub_by_polynomial(&new_q);

      tmp.coeffs[0] -= g_zeta;

      tmp.into_div_by_deg_one_polynomial(&zeta)
    };

    let mut process = process;

    assert_eq!(rem, E::Scalar::ZERO);

    let comm_quot_f_x_zeta = E::CE::commit(ck, &quot_poly.coeffs, &E::Scalar::ZERO);

    transcript.absorb(b"quot_f_x_zeta", &[comm_quot_f_x_zeta].to_vec().as_slice());

    let beta = transcript.squeeze(b"beta")?;

    process.open_phase_3(&beta);

    let comm_quot_m = E::CE::commit(ck, &process.quot_m_poly.coeffs, &E::Scalar::ZERO);

    transcript.absorb(b"quot_m", &[comm_quot_m].to_vec().as_slice());

    let z = transcript.squeeze(b"z")?;

    process.open_phase_6(&z);

    let comm_quot_l = E::CE::commit(ck, &process.quot_l_poly.coeffs, &E::Scalar::ZERO);

    Ok(EvaluationArgument {
      g_zeta: batch_eval.g_zeta,
      g_zeta_inv: batch_eval.g_zeta_inv,
      h_zeta: batch_eval.h_zeta,
      h_zeta_inv: batch_eval.h_zeta_inv,
      h_alpha,
      s_zeta: batch_eval.s_zeta,
      s_zeta_inv: batch_eval.s_zeta_inv,
      d_zeta: batch_eval.d_zeta,

      comm_h,
      comm_g,
      comm_q,
      comm_s,
      comm_d,
      comm_quot_m,
      comm_quot_l,
      comm_quot_f_x_zeta,

      zeta,
      alpha,
      beta,
      z,
      gamma,
    })
  }

  /// A method to verify the purported evaluation of a multilinear polynomials
  fn verify(
    vk: &Self::VerifierKey,
    _transcript: &mut E::TE,
    comm: &<<E as Engine>::CE as CommitmentEngineTrait<E>>::Commitment,
    point: &[E::Scalar],
    eval: &E::Scalar,
    arg: &Self::EvaluationArgument,
  ) -> Result<(), NovaError> {
    let log_n = point.len();
    assert_eq!(log_n % 2, 0);
    let (point_l, point_r) = point.split_at(log_n / 2);
    let (point_l, point_r) = (point_l.to_vec(), point_r.to_vec());

    let (eq_poly_1, _) = make_pu_poly(&point_l);
    let (eq_poly_2, _) = make_pu_poly(&point_r);

    let zeta = arg.zeta;
    let zeta_inv = zeta.invert().unwrap();

    let pu1_zeta = eq_poly_1.evaluate(&zeta);
    let pu2_zeta = eq_poly_2.evaluate(&zeta);
    let pu1_zeta_inv = eq_poly_1.evaluate(&zeta_inv);
    let pu2_zeta_inv = eq_poly_2.evaluate(&zeta_inv);

    // Check IPA
    let g_zeta = arg.g_zeta;
    let h_zeta = arg.h_zeta;
    let g_zeta_inv = arg.g_zeta_inv;
    let h_zeta_inv = arg.h_zeta_inv;
    let s_zeta = arg.s_zeta;
    let s_zeta_inv = arg.s_zeta_inv;
    let gamma = arg.gamma;
    let h_alpha = arg.h_alpha;

    // dbg!(eval_p(&point_l, E::Scalar::from(2)));
    dbg!(&eq_poly_1.evaluate(&E::Scalar::from(2)));

    let mut lhs = g_zeta * pu1_zeta_inv + g_zeta_inv * pu1_zeta;
    lhs += gamma * (h_zeta * pu2_zeta_inv + h_zeta_inv * pu2_zeta);

    let mut rhs = h_alpha + gamma * *eval;
    rhs += rhs;
    rhs += zeta * s_zeta + zeta_inv * s_zeta_inv;

    assert_eq!(lhs, rhs);

    if lhs != rhs {
      return Err(NovaError::ProofVerifyError {
        reason: "IPA check failed".to_string(),
      });
    }

    // Check degree
    let d_zeta = arg.d_zeta;
    let g_zeta_inv = arg.g_zeta_inv;
    let zeta_b_one = zeta.pow_vartime(&[(1_u64 << (log_n / 2)) - 1]);
    dbg!(&zeta_b_one);
    dbg!(&g_zeta_inv);
    dbg!(&d_zeta);
    if d_zeta != zeta_b_one * g_zeta_inv {
      return Err(NovaError::ProofVerifyError {
        reason: "Degree check failed".to_string(),
      });
    }

    // Check g(X)
    let alpha = arg.alpha;

    let comm_f = comm;
    let comm_g = arg.comm_g;
    let comm_h = arg.comm_h;
    let comm_s = arg.comm_s;
    let comm_d = arg.comm_d;
    let comm_q = arg.comm_q;
    let comm_quot_f_x_zeta = arg.comm_quot_f_x_zeta;

    let zeta_b = zeta_b_one * zeta;
    let zeta_b_alpha = zeta_b - alpha;
    let g1 = Commitment::new(DlogGroup::group(&vk.g));

    let g2 = <<E as Engine>::GE as PairingGroup>::G2::gen();
    let tau2 = <E::GE as PairingGroup>::G2::group(&vk.tau_h);

    let ll = *comm_f + comm_q * (-zeta_b_alpha) + g1 * (-g_zeta);
    let lr = g2.clone();
    let rl = comm_quot_f_x_zeta;
    let rr = g2 * (-zeta) + tau2;

    let pairing_l = E::GE::pairing(&ll.into_inner(), &lr);
    let pairing_r = E::GE::pairing(&rl.into_inner(), &rr);
    if pairing_l != pairing_r {
      return Err(NovaError::ProofVerifyError {
        reason: "g Check Paring check failed".to_string(),
      });
    }

    // Check KZG
    {
      let z = arg.z;
      let beta = arg.beta;

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

      let comm_quot_m = arg.comm_quot_m;
      let comm_quot_l = arg.comm_quot_l;

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

      f = f + Commitment::new(DlogGroup::group(&vk.g)) * (-scalar);

      let ll = f + comm_quot_l * z;
      let ll = ll.into_inner();
      let lr = DlogGroup::group(&vk.h);

      let rl = comm_quot_l.into_inner();
      let rr = DlogGroup::group(&vk.tau_h);

      let pairing_l = E::GE::pairing(&ll, &lr);
      let pairing_r = E::GE::pairing(&rl, &rr);

      if pairing_l != pairing_r {
        return Err(NovaError::ProofVerifyError {
          reason: "KZG Pairing check failed".to_string(),
        });
      }
    }

    Ok(())
  }
}
