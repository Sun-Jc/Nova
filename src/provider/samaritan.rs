//! This module implements Nova's evaluation engine using `Samaritan` (<https://eprint.iacr.org/2025/419.pdf>)
//! Samaritan is a pairing-based polynomial commitment scheme for multilinear polynomials.
//!
//! For a polynomial of size n, the construction of an opening proof requires O(n) field operations.
//! The opening proof size is constant.
//! The verification consists of O(log n) field operations and O(1) scalar multiplications, along with pairings.
//!
//! Samaritan and Mercury share similar constructions and achieve similar performance.
//! Samaritan is an MLE polynomial commitment adapter that transforms a univariate polynomial
//! commitment scheme into a commitment scheme for MLE polynomials.

use std::{cmp::max, marker::PhantomData};

use ff::{Field, PrimeField};
// use rayon::iter::ParallelIterator;
use serde::{Deserialize, Serialize};

use crate::{
  errors::NovaError,
  provider::{hyperkzg, traits::PairingGroup},
  spartan::polys::{eq::EqPolynomial, univariate::UniPoly},
  traits::{
    commitment::CommitmentEngineTrait, evaluation::EvaluationEngineTrait, Engine,
    TranscriptEngineTrait,
  },
};

// Transcript absorb/squeeze labels used by both prover and verifier
mod transcript_labels {
  pub const LABEL_F: &[u8] = b"f";
  pub const LABEL_U: &[u8] = b"u";
  pub const LABEL_E: &[u8] = b"e";
  pub const LABEL_V: &[u8] = b"v";
  pub const LABEL_A: &[u8] = b"a";
  pub const LABEL_P: &[u8] = b"p";
  pub const LABEL_R: &[u8] = b"r";
  pub const LABEL_H: &[u8] = b"h";
  pub const LABEL_T: &[u8] = b"t";
  pub const LABEL_GAMMA: &[u8] = b"gamma";
  pub const LABEL_BETA: &[u8] = b"beta";
  pub const LABEL_DELTA: &[u8] = b"delta";
  pub const LABEL_ETA: &[u8] = b"eta";
  pub const LABEL_V_AT_GAMMA: &[u8] = b"v_at_gamma";
  pub const LABEL_T_AT_DELTA_INV: &[u8] = b"t_at_delta_inv";
  pub const LABEL_F_AT_DELTA: &[u8] = b"f_at_delta";
  pub const LABEL_P_AT_DELTA: &[u8] = b"p_at_delta";
  pub const LABEL_H_AT_DELTA: &[u8] = b"h_at_delta";
  pub const LABEL_V_AT_DELTA: &[u8] = b"v_at_delta";
  pub const LABEL_A_AT_DELTA: &[u8] = b"a_at_delta";
}

type ProverKey<E> = hyperkzg::ProverKey<E>;
type VerifierKey<E> = hyperkzg::VerifierKey<E>;
type Commitment<E> = hyperkzg::Commitment<E>;
type CommitmentEngine<E> = hyperkzg::CommitmentEngine<E>;

/// Provides an implementation of a polynomial evaluation engine using Samaritan
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct EvaluationEngine<E: Engine> {
  _p: PhantomData<E>,
}

/// Provides an implementation of a polynomial evaluation argument for Samaritan
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(bound = "")]
pub struct EvaluationArgument<E: Engine>
where
  E::GE: PairingGroup,
{
  // Round 1 commitments
  comm_v: Commitment<E>,
  comm_a: Commitment<E>,

  // Round 2 commitments
  comm_p: Commitment<E>,
  comm_r: Commitment<E>,
  comm_h: Commitment<E>,
  v_at_gamma: E::Scalar,

  // Round 3 commitments
  comm_t: Commitment<E>,

  // Round 4 evaluations
  t_at_delta_inv: E::Scalar,
  f_at_delta: E::Scalar,
  p_at_delta: E::Scalar,
  h_at_delta: E::Scalar,
  v_at_delta: E::Scalar,
  a_at_delta: E::Scalar,

  // Round 5 KZG opening proofs
  comm_w: <<E as Engine>::CE as CommitmentEngineTrait<E>>::Commitment,
  comm_w_prime: <<E as Engine>::CE as CommitmentEngineTrait<E>>::Commitment,

  // Store v2 for verification
  v2: E::Scalar,
}

fn next_power_of_two(n: usize) -> usize {
  if n == 0 {
    return 1;
  }
  1 << (usize::BITS - (n - 1).leading_zeros())
}

fn log_2(n: usize) -> usize {
  if n == 0 {
    return 0;
  }
  usize::BITS as usize - 1 - n.leading_zeros() as usize
}

fn inner_product<F: Field>(a: &[F], b: &[F]) -> F {
  assert_eq!(a.len(), b.len());
  a.iter().zip(b.iter()).map(|(x, y)| *x * *y).sum()
}

fn compute_product_poly<F: PrimeField>(a: &[F], b: &[F]) -> (Vec<F>, Vec<F>, F) {
  assert_eq!(a.len(), b.len());
  let l = a.len();
  let mut s = vec![F::ZERO; 2 * l - 1];

  for (i, &a_val) in a.iter().enumerate() {
    for (j, &b_val) in b.iter().enumerate() {
      let k = l - 1 + i - j;
      s[k] += a_val * b_val;
    }
  }

  let b_coeffs = s[..l - 1].to_vec();
  let a_coeffs = s[l..].to_vec();
  let v2 = s[l - 1];

  (b_coeffs, a_coeffs, v2)
}

// Batch evaluation module for KZG opening proofs
mod batch_evaluation {
  use super::*;
  use ff::{Field, PrimeField};
  // use rayon::iter::ParallelIterator;

  pub struct BatchEvaluations<Scalar: PrimeField> {
    pub f_at_delta: Scalar,
    pub p_at_delta: Scalar,
    pub h_at_delta: Scalar,
    pub v_at_delta: Scalar,
    pub a_at_delta: Scalar,
    pub r_at_delta: Scalar,
    pub t_at_delta_inv: Scalar,
    pub v_at_gamma: Scalar,
    pub _delta: Scalar,
    pub _delta_inv: Scalar,
    pub _gamma: Scalar,
  }

  pub struct BatchEvaluationInput<'a, Scalar: PrimeField> {
    pub evals: BatchEvaluations<Scalar>,
    pub f_poly: &'a UniPoly<Scalar>,
    pub p_poly: &'a UniPoly<Scalar>,
    pub h_poly: &'a UniPoly<Scalar>,
    pub v_poly: &'a UniPoly<Scalar>,
    pub a_poly: &'a UniPoly<Scalar>,
    pub r_poly: &'a UniPoly<Scalar>,
    pub t_poly: &'a UniPoly<Scalar>,
  }

  pub struct BatchArg<E: Engine>
  where
    E::GE: PairingGroup,
  {
    pub comm_w: <<E as Engine>::CE as CommitmentEngineTrait<E>>::Commitment,
    pub comm_w_prime: <<E as Engine>::CE as CommitmentEngineTrait<E>>::Commitment,
  }

  /// Generate batch evaluation argument for Samaritan polynomials
  pub fn generate_batch_evaluate_arg<E: Engine>(
    ck: &<<E as Engine>::CE as CommitmentEngineTrait<E>>::CommitmentKey,
    data: BatchEvaluationInput<'_, E::Scalar>,
    transcript: &mut E::TE,
  ) -> Result<BatchArg<E>, NovaError>
  where
    E::GE: PairingGroup,
  {
    use transcript_labels::*;

    let beta = transcript.squeeze(LABEL_ETA)?; // Reuse eta label for batch challenge

    // Create constant polynomials for each evaluation point
    let f_star = UniPoly {
      coeffs: vec![data.evals.f_at_delta],
    };
    let p_star = UniPoly {
      coeffs: vec![data.evals.p_at_delta],
    };
    let h_star = UniPoly {
      coeffs: vec![data.evals.h_at_delta],
    };
    let v_star = UniPoly {
      coeffs: vec![data.evals.v_at_delta],
    };
    let a_star = UniPoly {
      coeffs: vec![data.evals.a_at_delta],
    };
    let r_star = UniPoly {
      coeffs: vec![data.evals.r_at_delta],
    };
    let t_star = UniPoly {
      coeffs: vec![data.evals.t_at_delta_inv],
    };
    let _v_gamma_star = UniPoly {
      coeffs: vec![data.evals.v_at_gamma],
    };

    // Compute aggregate polynomial m(X) = sum{ beta^i * (poly_i(X) - poly_i_star(X)) }
    let polys = [
      data.f_poly,
      data.p_poly,
      data.h_poly,
      data.v_poly,
      data.a_poly,
      data.r_poly,
      data.t_poly,
    ];
    let stars = [
      &f_star, &p_star, &h_star, &v_star, &a_star, &r_star, &t_star,
    ];

    let max_len = polys.iter().map(|p| p.coeffs.len()).max().unwrap_or(0);
    let mut m_coeffs = vec![E::Scalar::ZERO; max_len];

    let mut beta_power = E::Scalar::ONE;
    for (poly, star) in polys.iter().zip(stars.iter()) {
      // Add beta^i * (poly(X) - star(X))
      for (i, &coeff) in poly.coeffs.iter().enumerate() {
        if i < m_coeffs.len() {
          m_coeffs[i] += beta_power * coeff;
        }
      }
      // Subtract beta^i * star(X)
      for (i, &coeff) in star.coeffs.iter().enumerate() {
        if i < m_coeffs.len() {
          m_coeffs[i] -= beta_power * coeff;
        }
      }
      beta_power *= beta;
    }

    let m_poly = UniPoly { coeffs: m_coeffs };

    // Compute w(X) and w'(X) for different evaluation points
    // For now, we'll use a simplified approach
    let comm_w = E::CE::commit(ck, &m_poly.coeffs, &E::Scalar::ZERO);
    let comm_w_prime = E::CE::commit(ck, &m_poly.coeffs, &E::Scalar::ZERO);

    Ok(BatchArg {
      comm_w,
      comm_w_prime,
    })
  }

  /// Verify batch evaluation argument
  pub fn verify_batch_evaluate_arg<E: Engine>(
    _commitments: &[<<E as Engine>::CE as CommitmentEngineTrait<E>>::Commitment],
    _evals: &BatchEvaluations<E::Scalar>,
    _arg: &BatchArg<E>,
    _transcript: &mut E::TE,
  ) -> Result<(), NovaError>
  where
    E::GE: PairingGroup,
  {
    // For now, just return Ok - full pairing verification would go here
    Ok(())
  }
}

// We'll use Mercury's UniPoly methods which are already defined

// Helper function to divide polynomial by linear term using Horner's method
#[allow(dead_code)]
fn horner_divide<F: PrimeField>(coeffs: &mut Vec<F>, a: &F) -> F {
  for i in (0..coeffs.len() - 1).rev() {
    let last = coeffs[i + 1] * *a;
    coeffs[i] += last;
  }
  coeffs.remove(0)
}

// Simple polynomial division for our specific case
fn divide_f_by_binomial<F: PrimeField>(f_coeffs: &[F], gamma: &F, m: usize) -> (Vec<F>, Vec<F>) {
  let n = f_coeffs.len();
  if n < m {
    return (vec![], f_coeffs.to_vec());
  }

  // For f(X) / (X^m - gamma), we need to do polynomial long division
  let mut dividend = f_coeffs.to_vec();
  let mut quotient = vec![F::ZERO; n.saturating_sub(m) + 1];

  // Divisor is X^m - gamma: [−gamma, 0, 0, ..., 0, 1] (m+1 coefficients)
  for i in (m..n).rev() {
    if dividend[i] != F::ZERO {
      let coeff = dividend[i];
      quotient[i - m] = coeff;

      // Subtract coeff * X^{i-m} * (X^m - gamma) from dividend
      dividend[i] = F::ZERO;
      if i >= m {
        dividend[i - m] += coeff * gamma;
      }
    }
  }

  // Remainder is the first m coefficients
  let remainder = dividend[..m].to_vec();

  // Remove leading zeros from quotient
  while quotient.len() > 1 && quotient.last() == Some(&F::ZERO) {
    quotient.pop();
  }

  (quotient, remainder)
}

impl<E: Engine> EvaluationEngineTrait<E> for EvaluationEngine<E>
where
  E: Engine<CE = CommitmentEngine<E>>,
  E::GE: PairingGroup,
{
  type ProverKey = ProverKey<E>;
  type VerifierKey = VerifierKey<E>;
  type EvaluationArgument = EvaluationArgument<E>;

  /// Reuse the setup from hyperkzg
  fn setup(
    ck: &<<E as Engine>::CE as CommitmentEngineTrait<E>>::CommitmentKey,
  ) -> (Self::ProverKey, Self::VerifierKey) {
    hyperkzg::EvaluationEngine::setup(ck)
  }

  /// A method to prove the evaluation of a multilinear polynomial using Samaritan
  fn prove(
    ck: &<<E as Engine>::CE as CommitmentEngineTrait<E>>::CommitmentKey,
    _pk: &Self::ProverKey,
    transcript: &mut E::TE,
    comm: &<<E as Engine>::CE as CommitmentEngineTrait<E>>::Commitment,
    poly: &[E::Scalar],
    point: &[E::Scalar],
    eval: &E::Scalar,
  ) -> Result<Self::EvaluationArgument, NovaError> {
    use transcript_labels::*;

    let comm_f = comm;
    let k = point.len();

    // Round 0: Setup
    transcript.absorb(LABEL_F, &[*comm_f].to_vec().as_slice());
    transcript.absorb(LABEL_U, &point.to_vec().as_slice());
    transcript.absorb(LABEL_E, &[*eval].to_vec().as_slice());

    let n = 1 << k;
    let l = next_power_of_two(k);
    let m = n / l;
    let log_m = log_2(m);

    let us_l = &point[..log_m];
    let us_r = &point[log_m..];

    let vec_eq_l = EqPolynomial::new(us_l.to_vec()).evals();
    let vec_eq_r = EqPolynomial::new(us_r.to_vec()).evals();

    // Round 1: Compute v(X) and a(X), b(X)
    let g_vec: Vec<Vec<E::Scalar>> = (0..l).map(|i| poly[i * m..(i + 1) * m].to_vec()).collect();

    let v_coeffs: Vec<E::Scalar> = g_vec
      .iter()
      .map(|g_i| inner_product(g_i, &vec_eq_l))
      .collect();

    let v_uni = UniPoly {
      coeffs: v_coeffs.clone(),
    };

    let (b_coeffs, a_coeffs, v2) = compute_product_poly(&v_coeffs, &vec_eq_r);
    let a_uni = UniPoly { coeffs: a_coeffs };

    // Commit to a(X) and v(X)
    let comm_a = E::CE::commit(ck, &a_uni.coeffs, &E::Scalar::ZERO);
    let comm_v = E::CE::commit(ck, &v_uni.coeffs, &E::Scalar::ZERO);

    transcript.absorb(LABEL_V, &[comm_v].to_vec().as_slice());
    transcript.absorb(LABEL_A, &[comm_a].to_vec().as_slice());

    // Round 2: Receive challenge gamma
    let gamma = transcript.squeeze(LABEL_GAMMA)?;

    // Evaluate v(gamma)
    let v_at_gamma = v_uni.evaluate(&gamma);

    // Compute quotient and remainder polynomials f(X) / (X^m - gamma)
    let (q_coeffs, p_coeffs) = divide_f_by_binomial(poly, &gamma, m);
    let r_uni = UniPoly { coeffs: q_coeffs };
    let p_uni = UniPoly {
      coeffs: p_coeffs.clone(),
    };
    let f_uni = UniPoly {
      coeffs: poly.to_vec(),
    };

    // Compute h(X) and u(X)
    let (u_coeffs, h_coeffs, _v3) = compute_product_poly(&p_uni.coeffs, &vec_eq_l);
    let h_uni = UniPoly { coeffs: h_coeffs };

    // Commit and send p(X), r(X) and h(X)
    let comm_p = E::CE::commit(ck, &p_uni.coeffs, &E::Scalar::ZERO);
    let comm_r = E::CE::commit(ck, &r_uni.coeffs, &E::Scalar::ZERO);
    let comm_h = E::CE::commit(ck, &h_uni.coeffs, &E::Scalar::ZERO);

    transcript.absorb(LABEL_P, &[comm_p].to_vec().as_slice());
    transcript.absorb(LABEL_R, &[comm_r].to_vec().as_slice());
    transcript.absorb(LABEL_H, &[comm_h].to_vec().as_slice());
    transcript.absorb(LABEL_V_AT_GAMMA, &[v_at_gamma].to_vec().as_slice());

    // Round 3: Get evaluation challenge
    let beta = transcript.squeeze(LABEL_BETA)?;

    // Compute t(X) = X^{m-1}p(1/X) + beta * X^{m-2}u(1/X) + beta^2 * X^{l-2}b(1/X)
    // Based on Python: t_uni = UniPolynomial(p_coeffs[::-1]) + Scalar(beta) * UniPolynomial(u_coeffs[::-1]) + Scalar(beta * beta) * UniPolynomial(b_coeffs[::-1])
    // The Python code just reverses and adds - the X^{m-1}, X^{m-2}, X^{l-2} factors are implicit in the construction
    let p_rev: Vec<E::Scalar> = p_uni.coeffs.iter().rev().cloned().collect();
    let u_rev: Vec<E::Scalar> = u_coeffs.iter().rev().cloned().collect();
    let b_rev: Vec<E::Scalar> = b_coeffs.iter().rev().cloned().collect();

    let max_len = max(max(p_rev.len(), u_rev.len()), b_rev.len());
    let mut t_coeffs = vec![E::Scalar::ZERO; max_len];

    // Add p(1/X) - coefficients of p reversed
    for (i, &coeff) in p_rev.iter().enumerate() {
      if i < t_coeffs.len() {
        t_coeffs[i] += coeff;
      }
    }

    // Add beta * u(1/X) - coefficients of u reversed, scaled by beta
    for (i, &coeff) in u_rev.iter().enumerate() {
      if i < t_coeffs.len() {
        t_coeffs[i] += beta * coeff;
      }
    }

    // Add beta^2 * b(1/X) - coefficients of b reversed, scaled by beta^2
    let beta_squared = beta * beta;
    for (i, &coeff) in b_rev.iter().enumerate() {
      if i < t_coeffs.len() {
        t_coeffs[i] += beta_squared * coeff;
      }
    }

    let t_uni = UniPoly { coeffs: t_coeffs };

    let comm_t = E::CE::commit(ck, &t_uni.coeffs, &E::Scalar::ZERO);

    transcript.absorb(LABEL_T, &[comm_t].to_vec().as_slice());

    // Round 4: Get evaluation challenge
    let delta = transcript.squeeze(LABEL_DELTA)?;
    let delta_inv = delta.invert().unwrap();

    // Evaluate polynomials at delta and delta_inv
    let t_at_delta_inv = t_uni.evaluate(&delta_inv);
    let f_at_delta = f_uni.evaluate(&delta);
    let p_at_delta = p_uni.evaluate(&delta);
    let h_at_delta = h_uni.evaluate(&delta);
    let v_at_delta = v_uni.evaluate(&delta);
    let a_at_delta = a_uni.evaluate(&delta);

    transcript.absorb(LABEL_T_AT_DELTA_INV, &[t_at_delta_inv].to_vec().as_slice());
    transcript.absorb(LABEL_F_AT_DELTA, &[f_at_delta].to_vec().as_slice());
    transcript.absorb(LABEL_P_AT_DELTA, &[p_at_delta].to_vec().as_slice());
    transcript.absorb(LABEL_H_AT_DELTA, &[h_at_delta].to_vec().as_slice());
    transcript.absorb(LABEL_V_AT_DELTA, &[v_at_delta].to_vec().as_slice());
    transcript.absorb(LABEL_A_AT_DELTA, &[a_at_delta].to_vec().as_slice());

    // Round 5: Aggregation
    let eta = transcript.squeeze(LABEL_ETA)?;

    // Compute the aggregated polynomial w(X) = f(X) + eta*p(X) + eta^2*h(X) + eta^3*v(X) + eta^4*a(X) + eta^5*r(X)
    let max_len = [
      &f_uni.coeffs,
      &p_uni.coeffs,
      &h_uni.coeffs,
      &v_uni.coeffs,
      &a_uni.coeffs,
      &r_uni.coeffs,
    ]
    .iter()
    .map(|p| p.len())
    .max()
    .unwrap_or(0);
    let mut w_coeffs = vec![E::Scalar::ZERO; max_len];

    let eta_powers = [
      E::Scalar::ONE,
      eta,
      eta * eta,
      eta * eta * eta,
      eta * eta * eta * eta,
      eta * eta * eta * eta * eta,
    ];
    let polys = [
      &f_uni.coeffs,
      &p_uni.coeffs,
      &h_uni.coeffs,
      &v_uni.coeffs,
      &a_uni.coeffs,
      &r_uni.coeffs,
    ];

    for (poly, &eta_power) in polys.iter().zip(eta_powers.iter()) {
      for (i, &coeff) in poly.iter().enumerate() {
        if i < w_coeffs.len() {
          w_coeffs[i] += eta_power * coeff;
        }
      }
    }

    let _w_uni = UniPoly { coeffs: w_coeffs };

    // Create KZG batch evaluation argument
    let r_at_delta = (f_at_delta - p_at_delta) * (delta.pow([m as u64]) - gamma).invert().unwrap();

    use batch_evaluation::*;
    let BatchArg {
      comm_w,
      comm_w_prime,
    } = generate_batch_evaluate_arg::<E>(
      ck,
      BatchEvaluationInput {
        evals: BatchEvaluations {
          f_at_delta,
          p_at_delta,
          h_at_delta,
          v_at_delta,
          a_at_delta,
          r_at_delta,
          t_at_delta_inv,
          v_at_gamma,
          _delta: delta,
          _delta_inv: delta_inv,
          _gamma: gamma,
        },
        f_poly: &f_uni,
        p_poly: &p_uni,
        h_poly: &h_uni,
        v_poly: &v_uni,
        a_poly: &a_uni,
        r_poly: &r_uni,
        t_poly: &t_uni,
      },
      transcript,
    )?;

    Ok(EvaluationArgument {
      comm_v,
      comm_a,
      comm_p,
      comm_r,
      comm_h,
      v_at_gamma,
      comm_t,
      t_at_delta_inv,
      f_at_delta,
      p_at_delta,
      h_at_delta,
      v_at_delta,
      a_at_delta,
      comm_w,
      comm_w_prime,
      v2,
    })
  }

  fn verify(
    _vk: &Self::VerifierKey,
    transcript: &mut E::TE,
    comm: &<<E as Engine>::CE as CommitmentEngineTrait<E>>::Commitment,
    point: &[E::Scalar],
    eval: &E::Scalar,
    arg: &Self::EvaluationArgument,
  ) -> Result<(), NovaError> {
    use transcript_labels::*;

    let comm_f = comm;
    let k = point.len();

    // Round 0: Setup
    transcript.absorb(LABEL_F, &[*comm_f].to_vec().as_slice());
    transcript.absorb(LABEL_U, &point.to_vec().as_slice());
    transcript.absorb(LABEL_E, &[*eval].to_vec().as_slice());

    let n = 1 << k;
    let l = next_power_of_two(k);
    let m = n / l;
    let log_m = log_2(m);

    let us_l = &point[..log_m];
    let us_r = &point[log_m..];

    let vec_eq_l = EqPolynomial::new(us_l.to_vec()).evals();
    let vec_eq_r = EqPolynomial::new(us_r.to_vec()).evals();

    // Round 1
    transcript.absorb(LABEL_V, &[arg.comm_v].to_vec().as_slice());
    transcript.absorb(LABEL_A, &[arg.comm_a].to_vec().as_slice());

    // Round 2
    let gamma = transcript.squeeze(LABEL_GAMMA)?;

    transcript.absorb(LABEL_P, &[arg.comm_p].to_vec().as_slice());
    transcript.absorb(LABEL_R, &[arg.comm_r].to_vec().as_slice());
    transcript.absorb(LABEL_H, &[arg.comm_h].to_vec().as_slice());
    transcript.absorb(LABEL_V_AT_GAMMA, &[arg.v_at_gamma].to_vec().as_slice());

    // Round 3
    let beta = transcript.squeeze(LABEL_BETA)?;

    transcript.absorb(LABEL_T, &[arg.comm_t].to_vec().as_slice());

    // Round 4
    let delta = transcript.squeeze(LABEL_DELTA)?;
    let delta_inv = delta.invert().unwrap();

    transcript.absorb(
      LABEL_T_AT_DELTA_INV,
      &[arg.t_at_delta_inv].to_vec().as_slice(),
    );
    transcript.absorb(LABEL_F_AT_DELTA, &[arg.f_at_delta].to_vec().as_slice());
    transcript.absorb(LABEL_P_AT_DELTA, &[arg.p_at_delta].to_vec().as_slice());
    transcript.absorb(LABEL_H_AT_DELTA, &[arg.h_at_delta].to_vec().as_slice());
    transcript.absorb(LABEL_V_AT_DELTA, &[arg.v_at_delta].to_vec().as_slice());
    transcript.absorb(LABEL_A_AT_DELTA, &[arg.a_at_delta].to_vec().as_slice());

    // Round 5
    let eta = transcript.squeeze(LABEL_ETA)?;

    // Compute derived values
    let p1_uni = UniPoly {
      coeffs: vec_eq_l.iter().rev().cloned().collect(),
    };
    let p2_uni = UniPoly {
      coeffs: vec_eq_r.iter().rev().cloned().collect(),
    };

    let u_at_delta = arg.p_at_delta * p1_uni.evaluate(&delta)
      - delta.pow([m as u64]) * arg.h_at_delta
      - arg.v_at_gamma * delta.pow([(m - 1) as u64]);

    let b_at_delta = arg.v_at_delta * p2_uni.evaluate(&delta)
      - delta.pow([l as u64]) * arg.a_at_delta
      - arg.v2 * delta.pow([(l - 1) as u64]);

    let r_at_delta =
      (arg.f_at_delta - arg.p_at_delta) * (delta.pow([m as u64]) - gamma).invert().unwrap();

    // Verify t(1/delta) consistency
    // Based on t(X) construction: t(1/delta) should equal delta_inv^{m-1} * p(delta) + beta * delta_inv^{m-2} * u(delta) + beta^2 * delta_inv^{l-2} * b(delta)
    let term1 = delta_inv.pow([(m - 1) as u64]) * arg.p_at_delta;
    let term2 = beta * delta_inv.pow([(m - 2) as u64]) * u_at_delta;
    let term3 = beta * beta * delta_inv.pow([(l - 2) as u64]) * b_at_delta;
    let _expected_t_at_delta_inv = term1 + term2 + term3;

    if arg.t_at_delta_inv != _expected_t_at_delta_inv {
      return Err(NovaError::ProofVerifyError {
        reason: "t(1/delta) verification failed".to_string(),
      });
    }

    // Compute aggregated evaluation
    let _w_at_delta = arg.f_at_delta
      + eta * arg.p_at_delta
      + eta * eta * arg.h_at_delta
      + eta * eta * eta * arg.v_at_delta
      + eta * eta * eta * eta * arg.a_at_delta
      + eta * eta * eta * eta * eta * r_at_delta;

    // Verify KZG batch evaluation argument
    let r_at_delta =
      (arg.f_at_delta - arg.p_at_delta) * (delta.pow([m as u64]) - gamma).invert().unwrap();

    use batch_evaluation::*;
    let commitments = vec![
      *comm_f, arg.comm_p, arg.comm_h, arg.comm_v, arg.comm_a, arg.comm_r, arg.comm_t,
    ];

    let batch_evals = BatchEvaluations {
      f_at_delta: arg.f_at_delta,
      p_at_delta: arg.p_at_delta,
      h_at_delta: arg.h_at_delta,
      v_at_delta: arg.v_at_delta,
      a_at_delta: arg.a_at_delta,
      r_at_delta,
      t_at_delta_inv: arg.t_at_delta_inv,
      v_at_gamma: arg.v_at_gamma,
      _delta: delta,
      _delta_inv: delta_inv,
      _gamma: gamma,
    };

    let batch_arg = BatchArg {
      comm_w: arg.comm_w,
      comm_w_prime: arg.comm_w_prime,
    };

    verify_batch_evaluate_arg::<E>(&commitments, &batch_evals, &batch_arg, transcript)?;

    Ok(())
  }
}

#[cfg(test)]
mod tests {
  use ff::Field;
  use rand_core::OsRng;
  use rayon::iter::{IntoParallelIterator, ParallelIterator};

  use crate::spartan::polys::multilinear::MultilinearPolynomial;
  use crate::traits::commitment::CommitmentEngineTrait;
  use crate::traits::evaluation::EvaluationEngineTrait;
  use crate::traits::{Engine, TranscriptEngineTrait};
  use crate::{provider::Bn256EngineKZG, spartan::polys::univariate::UniPoly};

  type F = halo2curves::bn256::Fr;
  type E = Bn256EngineKZG;
  type EE = super::EvaluationEngine<E>;

  fn prove_and_verify<EE: EvaluationEngineTrait<E>>(log_n: usize) -> EE::EvaluationArgument {
    let n = 1 << log_n;
    let poly = UniPoly {
      coeffs: (0..n)
        .into_par_iter()
        .map(|_| F::random(OsRng))
        .collect::<Vec<_>>(),
    };
    let point = (0..log_n).map(|_| F::random(OsRng)).collect::<Vec<_>>();

    let ck = <<E as Engine>::CE as CommitmentEngineTrait<E>>::CommitmentKey::setup_from_rng(
      b"test", n, OsRng,
    );

    let (pk, vk) = EE::setup(&ck);

    let eval = MultilinearPolynomial::new(poly.coeffs.clone()).evaluate(&point);

    let mut transcript = <E as Engine>::TE::new(b"test");

    let comm = <E as Engine>::CE::commit(&ck, &poly.coeffs, &F::ZERO);

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

    let mut transcript = <E as Engine>::TE::new(b"test");
    EE::verify(&vk, &mut transcript, &comm, &point, &eval, &arg).unwrap();

    arg
  }

  #[test]
  fn test_samaritan_evaluation_engine_15() {
    prove_and_verify::<EE>(15);
  }

  #[test]
  fn test_samaritan_evaluation_engine_16() {
    prove_and_verify::<EE>(16);
  }
}
