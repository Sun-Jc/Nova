#![allow(unused)]
use ff::Field;
use halo2curves::ff_ext::quadratic;
use itertools::Itertools;
use rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator, ParallelIterator};

use crate::{
  errors::NovaError,
  spartan::{
    polys::{eq::EqPolynomial, multilinear::MultilinearPolynomial, univariate::UniPoly},
    sumcheck::{eq_sumcheck::EqSumCheckInstance, SumcheckProof},
  },
  traits::{Engine, TranscriptEngineTrait},
};

fn make_gkr_trace<F: Field>(mut p: Vec<F>, mut q: Vec<F>) -> (Vec<Vec<F>>, Vec<Vec<F>>) {
  let len = p.len().ilog2() as usize;
  let mut ps = Vec::with_capacity(len);
  let mut qs = Vec::with_capacity(len);

  ps.push(p.clone());
  qs.push(q.clone());

  while ps.last().unwrap().len() > 1 {
    let p = ps.last().unwrap();
    let q = qs.last().unwrap();

    let n = p.len() >> 1;
    // let mut p_new = Vec::with_capacity(n);
    // let mut q_new = Vec::with_capacity(n);

    let (p_first, p_second) = p.split_at(n);
    let (q_first, q_second) = q.split_at(n);

    let p_first = p_first.par_iter();
    let p_second = p_second.par_iter();
    let q_first = q_first.par_iter();
    let q_second = q_second.par_iter();

    let p_new = p_first
      .zip_eq(p_second)
      .zip_eq(q_first.clone())
      .zip_eq(q_second.clone())
      .map(|((((p_first, p_second), q_first), q_second))| {
        *p_first * *q_second + *p_second * *q_first
      })
      .collect::<Vec<_>>();

    let q_new = q_first
      .zip_eq(q_second)
      .map(|((q_first, q_second))| *q_first * *q_second)
      .collect::<Vec<_>>();

    ps.push(p_new);
    qs.push(q_new);
  }

  ps.reverse();
  qs.reverse();

  (ps, qs)
}

// P:
// p0, q0
// p0', p1', q0', q1'
// Claim: p0 + q0 * b
// Final Eval: (p0', q0', p1', q1')
// ---> Check eq * ( p0', q0') == final eval

// p[X] = sum eq(i, X)
// (p_prev[0, i] * q_prev[1, i] + p_prev[1, i] * q_prev[0, i])
// q[X] = sum eq(i, X)
// (q_prev[0, i] * q_prev[1, i])

#[derive(Debug, Default)]
struct LogupGKREval<Scalar: Field> {
  pub r: Vec<Scalar>,
  pub p0_val: Scalar,
  pub p1_val: Scalar,
  pub q0_val: Scalar,
  pub q1_val: Scalar,
}

#[derive(Debug, Default, Clone)]
struct LogupGKRLayerClaim<Scalar: Field> {
  pub p0: Scalar,
  pub p1: Scalar,
  pub q0: Scalar,
  pub q1: Scalar,
}

#[derive(Debug, Default, Clone)]
struct LogupGKRProof<E: Engine> {
  pub sumcheck_proofs: Vec<SumcheckProof<E>>,
  pub layer_claims: Vec<LogupGKRLayerClaim<E::Scalar>>,
}

impl<E: Engine> LogupGKRProof<E> {
  pub fn final_eval(&self) -> LogupGKRLayerClaim<E::Scalar> {
    self.layer_claims.last().unwrap().clone()
  }

  pub fn verify(&self, transcript: &mut E::TE) -> Result<Vec<E::Scalar>, NovaError> {
    let claims = &self.layer_claims;
    let proofs = &self.sumcheck_proofs;

    let mut tau: Vec<E::Scalar> = vec![];
    let mut p_claim = claims[0].p0;
    let mut q_claim = claims[0].q0;

    for layer in 0..proofs.len() {
      // rlc of p and q
      let rlc_coeff = transcript.squeeze(b"r")?;

      let claim = p_claim + q_claim * rlc_coeff;

      let proof = &proofs[layer];

      let (actual_final_eval, rs) = proof.verify(claim, layer, 3, transcript)?;

      let poly_eq = EqPolynomial::new(tau.clone());
      let eval_eq = poly_eq.evaluate(&rs);

      let LogupGKRLayerClaim { p0, p1, q0, q1 } = claims[layer + 1];

      let expected_final_eval = eval_eq * (p0 * q1 + p1 * q0 + rlc_coeff * q0 * q1);

      if actual_final_eval != expected_final_eval {
        return Err(NovaError::InvalidSumcheckProof);
      }

      if layer == proofs.len() - 1 {
        return Ok(rs);
      }

      let r0 = transcript.squeeze(b"r")?;

      tau = [vec![r0], rs].concat();

      let one = E::Scalar::ONE;
      p_claim = p0 * (one - r0) + p1 * r0;
      q_claim = q0 * (one - r0) + q1 * r0;
    }

    unreachable!()
  }

  /// Returns 1) proof, 2) tau, 3) (root_p, root_q)
  pub fn prove(
    ps: Vec<E::Scalar>,
    qs: Vec<E::Scalar>,
    transcript: &mut E::TE,
  ) -> Result<(Self, Vec<E::Scalar>, (E::Scalar, E::Scalar)), NovaError> {
    let (p_trace, q_trace) = make_gkr_trace(ps, qs);

    let root_p = *p_trace[0].last().unwrap();
    let root_q = *q_trace[0].last().unwrap();

    let mut p_claim = root_p;
    let mut q_claim = root_q;

    let mut tau = vec![];

    let mut proofs = vec![];
    let mut claims = vec![LogupGKRLayerClaim {
      p0: p_claim,
      p1: p_claim,
      q0: q_claim,
      q1: q_claim,
    }];

    let mut n = 1;

    for layer in 1..p_trace.len() {
      let rlc_coeff = transcript.squeeze(b"r")?;

      let mut claim = p_claim + q_claim * rlc_coeff;

      {
        let p_poly = p_trace[layer - 1].clone();
        let q_poly = q_trace[layer - 1].clone();
        let p_poly = MultilinearPolynomial::new(p_poly);
        let q_poly = MultilinearPolynomial::new(q_poly);
        let p_eval = p_poly.evaluate(&tau);
        let q_eval = q_poly.evaluate(&tau);
        assert_eq!(p_eval, p_claim);
        assert_eq!(q_eval, q_claim);
      }

      let (p0, p1) = p_trace[layer].split_at(n);
      let (q0, q1) = q_trace[layer].split_at(n);
      let mut p0_polynomial = MultilinearPolynomial::new(p0.to_vec());
      let mut p1_polynomial = MultilinearPolynomial::new(p1.to_vec());
      let mut q0_polynomial = MultilinearPolynomial::new(q0.to_vec());
      let mut q1_polynomial = MultilinearPolynomial::new(q1.to_vec());

      let mut sc_inst = EqSumCheckInstance::<E>::new(tau.clone());

      assert_eq!(1 << tau.len(), p0.len());

      let mut rs = vec![];
      let mut polys = vec![];
      for round in 0..n.ilog2() as usize {
        let (eval_point_h_0, eval_point_h_2, eval_point_h_3) = sc_inst
          .evaluation_points_cubic_with_four_inputs(
            &p0_polynomial,
            &p1_polynomial,
            &q0_polynomial,
            &q1_polynomial,
            &rlc_coeff,
          );

        let eval_point_h_1 = claim - eval_point_h_0;

        let poly = UniPoly::from_evals(&[
          eval_point_h_0,
          eval_point_h_1,
          eval_point_h_2,
          eval_point_h_3,
        ]);

        transcript.absorb(b"p", &poly);

        let r = transcript.squeeze(b"c")?;

        rs.push(r);
        polys.push(poly.compress());

        rayon::join(
          || rayon::join(|| sc_inst.bound(&r), || p0_polynomial.bind_poly_var_top(&r)),
          || {
            rayon::join(
              || p1_polynomial.bind_poly_var_top(&r),
              || {
                rayon::join(
                  || q0_polynomial.bind_poly_var_top(&r),
                  || q1_polynomial.bind_poly_var_top(&r),
                )
              },
            )
          },
        );

        claim = poly.evaluate(&r);
      }

      let (p0, p1, q0, q1) = (
        p0_polynomial[0],
        p1_polynomial[0],
        q0_polynomial[0],
        q1_polynomial[0],
      );

      let r0 = transcript.squeeze(b"r")?;
      tau = [vec![r0], rs].concat();

      proofs.push(SumcheckProof::new(polys));
      claims.push(LogupGKRLayerClaim { p0, p1, q0, q1 });

      let one = E::Scalar::ONE;
      p_claim = p0 * (one - r0) + p1 * r0;
      q_claim = q0 * (one - r0) + q1 * r0;

      n *= 2;
    }

    #[cfg(debug_assertions)]
    {
      let p_poly = p_trace.last().unwrap().clone();
      let p_poly = MultilinearPolynomial::new(p_poly.to_vec());
      let q_poly = MultilinearPolynomial::new(q_trace.last().unwrap().to_vec());

      let p_eval = p_poly.evaluate(&tau);
      let q_eval = q_poly.evaluate(&tau);

      assert_eq!(p_claim, p_eval);
      assert_eq!(q_claim, q_eval);
    }

    Ok((
      Self {
        sumcheck_proofs: proofs,
        layer_claims: claims,
      },
      tau[1..].to_vec(),
      (root_p, root_q),
    ))
  }
}

#[cfg(test)]
mod tests {
  use crate::provider::Bn256EngineIPA;

  use super::*;

  type E = Bn256EngineIPA;
  type Fr = <E as Engine>::Scalar;

  // * P0 = 1, -P1 = m
  // * Q0 = V + alpha * Addr, -Q1 = T + alpha * Addr
  // * Claim = 0

  #[test]
  fn test_logup_gkr() {
    let log_n = 10;
    let n = 1 << log_n;
    let p0 = (1..=n).map(|i| Fr::from(i as u64)).collect::<Vec<_>>();
    let q0 = (1..=n).map(|i| Fr::from(i as u64)).collect::<Vec<_>>();
    let p1 = (1..=n).map(|i| Fr::from(i as u64)).collect::<Vec<_>>();
    let q1 = (1..=n).map(|i| Fr::from(i as u64)).collect::<Vec<_>>();

    let ps: Vec<Fr> = [p0.clone(), p1.clone()].concat();
    let qs: Vec<Fr> = [q0.clone(), q1.clone()].concat();

    let root_expected = ps
      .iter()
      .zip(qs.iter())
      .map(|(p, q)| p * q.invert().unwrap())
      .sum::<Fr>();

    let mut transcript = <E as Engine>::TE::new(b"Test");

    let start_time = std::time::Instant::now();
    let (proof, tau_prover, (root_p, root_q)) =
      LogupGKRProof::<E>::prove(ps, qs, &mut transcript).unwrap();
    let duration = start_time.elapsed();
    println!("Time taken: {:?}", duration);

    {
      let root_actual = root_p * root_q.invert().unwrap();
      assert_eq!(root_expected, root_actual);
    }

    let mut transcript = <E as Engine>::TE::new(b"Test");
    let rs = proof.verify(&mut transcript).unwrap();

    assert_eq!(rs, tau_prover);

    let p0_poly = MultilinearPolynomial::new(p0.clone());
    let p1_poly = MultilinearPolynomial::new(p1.clone());
    let q0_poly = MultilinearPolynomial::new(q0.clone());
    let q1_poly = MultilinearPolynomial::new(q1.clone());

    let p0_eval = p0_poly.evaluate(&rs);
    let p1_eval = p1_poly.evaluate(&rs);
    let q0_eval = q0_poly.evaluate(&rs);
    let q1_eval = q1_poly.evaluate(&rs);

    let LogupGKRLayerClaim { p0, p1, q0, q1 } = proof.final_eval();

    assert_eq!(p0_eval, p0);
    assert_eq!(p1_eval, p1);
    assert_eq!(q0_eval, q0);
    assert_eq!(q1_eval, q1);
  }
}

struct LogupGKRMemoryCheckProof<E: Engine> {
  pub p1: Vec<E::Scalar>,
  pub q0: Vec<E::Scalar>,
  pub q1: Vec<E::Scalar>,
  pub proof: LogupGKRProof<E>,
  pub tau: Vec<E::Scalar>,
  pub root_q: E::Scalar,
}

impl<E: Engine> LogupGKRMemoryCheckProof<E> {
  pub fn prove_memory_check(
    addr: Vec<usize>,
    table: Vec<E::Scalar>,
    transcript: &mut E::TE,
  ) -> Result<Self, NovaError> {
    assert!(addr.len() >= table.len());
    let q0 = addr.par_iter().map(|addr| table[*addr]).collect::<Vec<_>>();
    let mut q1 = table;
    q1.resize(q0.len(), E::Scalar::ZERO);
    let p0 = vec![E::Scalar::ONE; q0.len()];
    // count addr by group by par
    let mut p1 = vec![0usize; q0.len()];
    for a in addr.iter() {
      p1[*a] += 1;
    }
    let p1 = p1
      .par_iter()
      .map(|p| E::Scalar::from(*p as u64))
      .collect::<Vec<_>>();

    // let ps = [p0.clone(), p1.clone()].concat();
    // let qs = [q0.clone(), q1.clone()].concat();

    let (proof, tau, (root_p, root_q)) = LogupGKRProof::<E>::prove(ps, qs, transcript)?;

    Ok(LogupGKRMemoryCheckProof {
      p1,
      q0,
      q1,
      proof,
      tau,
      root_q,
    })
  }

  pub fn prove_memory_check(
    &self,
    addr: Vec<usize>,
    table: Vec<E::Scalar>,
    transcript: &mut E::TE,
  ) -> Result<(), NovaError> {
  }
}
