#![allow(dead_code)]

use std::marker::PhantomData;

use crate::{
  errors::NovaError,
  provider::{
    hyperkzg::{Commitment, CommitmentEngine},
    mercury::{
      ipa::make_s_polynomial,
      kzg::EvaluationProcess,
      split_polynomial::{divide_by_binomial, split_polynomial},
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
fn make_pu_poly<Scalar: PrimeField>(u: &[Scalar]) -> UniPoly<Scalar> {
  let mut u = u.to_vec();
  u.reverse();

  let evals = EqPolynomial::new(u).evals();

  UniPoly { coeffs: evals }
}

fn eval_pu_poly<Scalar: PrimeField>(u: &[Scalar], r: &Scalar) -> Scalar {
  let mut res = Scalar::ONE;
  for (i, u_i) in u.iter().enumerate() {
    res *= *u_i * r.pow([1 << i]) + Scalar::ONE - u_i;
  }
  res
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

    let point = point.to_vec();

    let mut point = point.to_vec();
    point.reverse();

    let (f_poly, log_n, b, point_l, point_r) = {
      // Round 0: Prepare
      let mut f = UniPoly {
        coeffs: poly.to_vec(),
      };

      let mut log_n = point.len();
      if log_n % 2 != 0 {
        log_n += 1;

        point.push(E::Scalar::ZERO);
        f.raise(1 << log_n);
      }

      let b = 1 << (log_n / 2);

      f.raise(1 << log_n);

      let (point_l, point_r) = point.split_at(log_n / 2);

      (f, log_n, b, point_l.to_vec(), point_r.to_vec())
    };

    let p_poly1 = make_pu_poly(&point_l);
    let p_poly2 = make_pu_poly(&point_r);

    assert_eq!(p_poly1.coeffs.len(), b);
    assert_eq!(p_poly2.coeffs.len(), b);

    #[cfg(debug_assertions)]
    {
      // Check pu1(X), pu2(X)
      use rand_core::OsRng;

      let r = E::Scalar::random(OsRng);
      let pu1_eval_expected = eval_pu_poly(&point_l, &r);
      let pu2_eval_expected = eval_pu_poly(&point_r, &r);

      let pu1_eval_actual = p_poly1.evaluate(&r);
      let pu2_eval_actual = p_poly2.evaluate(&r);

      assert_eq!(pu1_eval_expected, pu1_eval_actual);
      assert_eq!(pu2_eval_expected, pu2_eval_actual);
    }

    let h_poly = {
      let mut coeffs = vec![E::Scalar::ZERO; b];

      for (j, coeff) in coeffs.iter_mut().enumerate() {
        for (i, eq_eval) in p_poly1.coeffs.iter().enumerate().take(b) {
          *coeff += f_poly.coeffs[j * b + i] * eq_eval;
        }
      }
      UniPoly { coeffs }
    };

    #[cfg(debug_assertions)]
    {
      // Check h, eq_evals2 ipa vs eval

      let eq_evals2 = p_poly2.coeffs.clone();

      assert_eq!(eq_evals2.len(), h_poly.coeffs.len());

      let ip = eq_evals2
        .iter()
        .zip(h_poly.coeffs.iter())
        .map(|(a, b)| *a * *b)
        .sum::<E::Scalar>();

      assert_eq!(ip, *_eval);
    }

    let comm_h = E::CE::commit(ck, &h_poly.coeffs, &E::Scalar::ZERO);

    transcript.absorb(b"f", &[*comm_f, comm_h].to_vec().as_slice());
    let alpha = transcript.squeeze(b"alpha")?;

    // Get q(X) and g(X)
    let (q_poly, g_poly) = {
      let f_col_sub_poly = split_polynomial(&f_poly, b);

      assert_eq!(f_col_sub_poly.len(), b);
      assert_eq!(f_col_sub_poly[0].coeffs.len(), b);

      divide_by_binomial(&f_col_sub_poly, &alpha)
    };

    assert_eq!(g_poly.coeffs.len(), b);

    #[cfg(debug_assertions)]
    {
      // Check g
      let f_col_sub_poly = split_polynomial(&f_poly, b);
      let f_alphas = f_col_sub_poly
        .iter()
        .map(|p| p.evaluate(&alpha))
        .collect::<Vec<_>>();

      f_alphas
        .iter()
        .zip(g_poly.coeffs.iter())
        .for_each(|(a, b)| {
          assert_eq!(*b, *a);
        });
    }

    #[cfg(debug_assertions)]
    {
      // Check q, g
      // f(r) = (X^b - alpha) * q(r) + g(r)

      use rand_core::OsRng;
      let r = E::Scalar::random(OsRng);

      let f_r = f_poly.evaluate(&r);
      let q_r = q_poly.evaluate(&r);
      let g_r = g_poly.evaluate(&r);

      assert_eq!(f_r, (r.pow([b as u64]) - alpha) * q_r + g_r);
    }

    let mut q_poly = q_poly;
    q_poly.trim();

    let comm_q = E::CE::commit(ck, &q_poly.coeffs, &E::Scalar::ZERO);
    let comm_g = E::CE::commit(ck, &g_poly.coeffs, &E::Scalar::ZERO);

    transcript.absorb(b"g", &[comm_g, comm_q].to_vec().as_slice());

    let h_alpha = h_poly.evaluate(&alpha);

    #[cfg(debug_assertions)]
    {
      // Check g eq_evals1 ipa vs h_alpha
      let ip = p_poly1
        .coeffs
        .iter()
        .zip(g_poly.coeffs.iter())
        .map(|(a, b)| *a * *b)
        .sum::<E::Scalar>();

      assert_eq!(ip, h_alpha);
    }

    transcript.absorb(b"h_alpha", &h_alpha);

    let gamma = transcript.squeeze(b"gamma")?;

    // Get s(X) for ipa
    let s_poly = {
      let left_polys = vec![p_poly1.clone(), p_poly2.clone()];
      let right_polys = vec![g_poly.clone(), h_poly.clone()];

      assert_eq!(p_poly1.coeffs.len(), g_poly.coeffs.len());
      assert_eq!(p_poly2.coeffs.len(), h_poly.coeffs.len());

      make_s_polynomial(left_polys, right_polys, log_n as u32 / 2, &gamma)
    };

    #[cfg(debug_assertions)]
    {
      // Check s_poly IPA proof

      use rand_core::OsRng;
      let r = E::Scalar::random(OsRng);
      let r_inv = r.invert().unwrap();

      let g_r = g_poly.evaluate(&r);
      let g_r_inv = g_poly.evaluate(&r_inv);
      let h_r = h_poly.evaluate(&r);
      let h_r_inv = h_poly.evaluate(&r_inv);
      let pu1_r = eval_pu_poly(&point_l, &r);
      let pu1_r_inv = eval_pu_poly(&point_l, &r_inv);
      let pu2_r = eval_pu_poly(&point_r, &r);
      let pu2_r_inv = eval_pu_poly(&point_r, &r_inv);

      let s_r = s_poly.evaluate(&r);
      let s_r_inv = s_poly.evaluate(&r_inv);

      let mut lhs = g_r * pu1_r_inv + g_r_inv * pu1_r;
      lhs += gamma * (h_r * pu2_r_inv + h_r_inv * pu2_r);

      let mut rhs = h_alpha + gamma * _eval;
      rhs += rhs;
      rhs += r * s_r + r_inv * s_r_inv;

      assert_eq!(lhs, rhs);
    }

    // Get d(X) for degree check
    let d_poly = {
      let mut d_poly = g_poly.clone();

      assert_eq!(d_poly.coeffs.len(), b);

      d_poly.coeffs.reverse();
      d_poly
    };

    let comm_s = E::CE::commit(ck, &s_poly.coeffs, &E::Scalar::ZERO);
    let comm_d = E::CE::commit(ck, &d_poly.coeffs, &E::Scalar::ZERO);

    transcript.absorb(b"d", &[comm_s, comm_d].to_vec().as_slice());

    let zeta = transcript.squeeze(b"zeta")?;

    let (process, batch_eval) =
      EvaluationProcess::<E>::init_eval_0(&alpha, &zeta, &g_poly, &h_poly, &s_poly, &d_poly);

    let (quot_poly, rem) = {
      let g_zeta = batch_eval.g_zeta;

      let zeta_b = zeta.pow([b as u64]);

      let zeta_b_alpha = zeta_b - alpha;

      let mut new_q = q_poly.clone();
      new_q.scale(&zeta_b_alpha);

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
    transcript: &mut E::TE,
    comm: &<<E as Engine>::CE as CommitmentEngineTrait<E>>::Commitment,
    point: &[E::Scalar],
    eval: &E::Scalar,
    arg: &Self::EvaluationArgument,
  ) -> Result<(), NovaError> {
    let mut point = point.to_vec();
    point.reverse();

    let mut log_n = point.len();
    if log_n % 2 != 0 {
      log_n += 1;

      point.push(E::Scalar::ZERO);
    }

    let (point_l, point_r) = point.split_at(log_n / 2);
    let (point_l, point_r) = (point_l.to_vec(), point_r.to_vec());

    let zeta = arg.zeta;
    let zeta_inv = zeta.invert().unwrap();

    let pu1_zeta = eval_pu_poly(&point_l, &zeta);
    let pu1_zeta_inv = eval_pu_poly(&point_l, &zeta_inv);
    let pu2_zeta = eval_pu_poly(&point_r, &zeta);
    let pu2_zeta_inv = eval_pu_poly(&point_r, &zeta_inv);

    let g_zeta = arg.g_zeta;
    let h_zeta = arg.h_zeta;
    let g_zeta_inv = arg.g_zeta_inv;
    let h_zeta_inv = arg.h_zeta_inv;
    let s_zeta = arg.s_zeta;
    let s_zeta_inv = arg.s_zeta_inv;
    let gamma = arg.gamma;
    let h_alpha = arg.h_alpha;
    let d_zeta = arg.d_zeta;

    let alpha = arg.alpha;

    let comm_f = comm;
    let comm_g = arg.comm_g;
    let comm_h = arg.comm_h;
    let comm_s = arg.comm_s;
    let comm_d = arg.comm_d;
    let comm_q = arg.comm_q;
    let comm_quot_f_x_zeta = arg.comm_quot_f_x_zeta;

    // 1. Check IPA
    {
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
    }

    // 2. Check degree
    let zeta_b_one = zeta.pow_vartime([(1_u64 << (log_n / 2)) - 1]);
    {
      if d_zeta != zeta_b_one * g_zeta_inv {
        return Err(NovaError::ProofVerifyError {
          reason: "Degree check failed".to_string(),
        });
      }
    }

    // 3. Check f(X) / (X^b - alpha) = (q(X), g(x))
    {
      let zeta_b = zeta_b_one * zeta;
      let zeta_b_alpha = zeta_b - alpha;
      let g1 = Commitment::new(DlogGroup::group(&vk.g));

      let g2 = <<E as Engine>::GE as PairingGroup>::G2::gen();
      let tau2 = <E::GE as PairingGroup>::G2::group(&vk.tau_h);

      let ll = *comm_f + comm_q * (-zeta_b_alpha) + g1 * (-g_zeta);
      let lr = g2;
      let rl = comm_quot_f_x_zeta;
      let rr = g2 * (-zeta) + tau2;

      let pairing_l = E::GE::pairing(&ll.into_inner(), &lr);
      let pairing_r = E::GE::pairing(&rl.into_inner(), &rr);

      if pairing_l != pairing_r {
        return Err(NovaError::ProofVerifyError {
          reason: "g Check Paring check failed".to_string(),
        });
      }

      {
        // alternative: ll + \zeta CH == CH, tau
        let ll1 = ll + comm_quot_f_x_zeta * zeta;
        let rl1 = comm_quot_f_x_zeta;
        let rr = tau2;

        let z = arg.z;
        let beta = arg.beta;
        let comm_quot_m = arg.comm_quot_m;
        let comm_quot_l = arg.comm_quot_l;

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

        let ll2 = f + comm_quot_l * z;
        let lr = DlogGroup::group(&vk.h);

        let rl2 = comm_quot_l;
        let rr = DlogGroup::group(&vk.tau_h);

        let d = E::Scalar::random(OsRng);

        let ll = ll1 + ll2 * d;
        let rl = rl1 + rl2 * d;

        let pairing_l = E::GE::pairing(&ll.into_inner(), &lr);
        let pairing_r = E::GE::pairing(&rl.into_inner(), &rr);

        if pairing_l != pairing_r {
          return Err(NovaError::ProofVerifyError {
            reason: "alternative: g Check Paring check failed".to_string(),
          });
        }
      }
    }

    // 4. Check KZG
    let z = arg.z;
    let beta = arg.beta;
    let comm_quot_m = arg.comm_quot_m;
    let comm_quot_l = arg.comm_quot_l;
    {
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

    {
      // Check FS
      transcript.absorb(b"f", &[*comm_f, comm_h].to_vec().as_slice());
      if alpha != transcript.squeeze(b"alpha")? {
        return Err(NovaError::ProofVerifyError {
          reason: "alpha sample failed".to_string(),
        });
      }
      transcript.absorb(b"g", &[comm_g, comm_q].to_vec().as_slice());
      transcript.absorb(b"h_alpha", &h_alpha);
      if gamma != transcript.squeeze(b"gamma")? {
        return Err(NovaError::ProofVerifyError {
          reason: "gamma sample failed".to_string(),
        });
      }

      transcript.absorb(b"d", &[comm_s, comm_d].to_vec().as_slice());

      let zeta = transcript.squeeze(b"zeta")?;
      if zeta != arg.zeta {
        return Err(NovaError::ProofVerifyError {
          reason: "zeta sample failed".to_string(),
        });
      }

      transcript.absorb(b"quot_f_x_zeta", &[comm_quot_f_x_zeta].to_vec().as_slice());

      if beta != transcript.squeeze(b"beta")? {
        return Err(NovaError::ProofVerifyError {
          reason: "beta sample failed".to_string(),
        });
      }

      transcript.absorb(b"quot_m", &[comm_quot_m].to_vec().as_slice());

      if z != transcript.squeeze(b"z")? {
        return Err(NovaError::ProofVerifyError {
          reason: "z sample failed".to_string(),
        });
      }
    }

    Ok(())
  }
}
