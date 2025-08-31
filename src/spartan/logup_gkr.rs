#![allow(unused)]
use ff::Field;
use halo2curves::ff_ext::quadratic;

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

  while p.len() > 1 {
    let n = p.len() >> 1;
    let mut p_new = Vec::with_capacity(n);
    let mut q_new = Vec::with_capacity(n);

    for i in 0..n {
      let first = i;
      let second = n + i;

      p_new.push(p[first] * q[second] + p[second] * q[first]);
      q_new.push(q[first] * q[second]);
    }

    ps.push(p_new.clone());
    qs.push(q_new.clone());

    (p, q) = (p_new, q_new);
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
struct GKRLeafEval<Scalar: Field> {
  pub r: Vec<Scalar>,
  pub rlc_coeffs: Vec<Scalar>,
  pub eval: Scalar,
}

#[derive(Debug, Default)]
struct LogupGKREval<Scalar: Field> {
  pub r: Vec<Scalar>,
  pub p0_val: Scalar,
  pub p1_val: Scalar,
  pub q0_val: Scalar,
  pub q1_val: Scalar,
}

#[derive(Debug, Default)]
struct LogupGKRLayerClaim<Scalar: Field> {
  pub p0: Scalar,
  pub p1: Scalar,
  pub q0: Scalar,
  pub q1: Scalar,
}

fn verify_gkr<E: Engine>(
  claims: Vec<LogupGKRLayerClaim<E::Scalar>>,
  proofs: Vec<SumcheckProof<E>>,
  transcript: &mut E::TE,
) -> Result<(), NovaError> {
  let mut tau: Vec<E::Scalar> = vec![];
  let mut p = claims[0].p0;
  let mut q = claims[0].q0;

  dbg!(claims.len());
  dbg!(proofs.len());

  for layer in 0..proofs.len() {
    let proof = &proofs[layer];

    // rlc of p and q
    let rlc_coeff = transcript.squeeze(b"r")?;

    let claim = p + q * rlc_coeff;

    let (actual_final_eval, r) = proof.verify(claim, layer, 3, transcript)?;

    let eq_poly = EqPolynomial::new(tau.clone());

    let eval_eq = eq_poly.evaluate(&r);

    let LogupGKRLayerClaim { p0, p1, q0, q1 } = claims[layer + 1];

    let expected_final_eval = eval_eq * (p0 * q1 + p1 * q0 + rlc_coeff * (q0 * q1));

    println!(
      "verify layer: {}, rlc_coeff: {:?},  eq: {:?},actual_final_eval: {:?}, expected_final_eval: {:?}",
      layer + 1,
      rlc_coeff,
      eval_eq,
      actual_final_eval,
      expected_final_eval
    );

    assert_eq!(actual_final_eval, expected_final_eval);

    // if actual_final_eval != expected_final_eval {
    //   return Err(NovaError::InvalidSumcheckProof);
    // }

    let r0 = transcript.squeeze(b"r")?;

    tau = [vec![r0], r].concat();
    let one = E::Scalar::ONE;
    p = p0 * (one - r0) + p1 * r0;
    q = q0 * (one - r0) + q1 * r0;
  }

  Ok(())
}

fn prove_gkr<E: Engine>(
  p_trace: Vec<Vec<E::Scalar>>,
  q_trace: Vec<Vec<E::Scalar>>,
  transcript: &mut E::TE,
) -> Result<
  (
    Vec<SumcheckProof<E>>,
    Vec<E::Scalar>,
    Vec<LogupGKRLayerClaim<E::Scalar>>,
  ),
  NovaError,
> {
  let mut p_claim = *p_trace[0].last().unwrap();
  let mut q_claim = *q_trace[0].last().unwrap();

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

    println!("layer: {}, rlc_coeff: {:?}", layer, rlc_coeff);

    let claim = p_claim + q_claim * rlc_coeff;

    // if layer > 1 {
    //   let LogupGKRLayerClaim { p0, p1, q0, q1 } = claims.last().unwrap();
    //   let e_claim = *p0 * *q1 + *p1 * *q0 + rlc_coeff * (*q0 * *q1);
    //   assert_eq!(claim, e_claim);
    // }

    let (p0, p1) = p_trace[layer].split_at(n);
    let (q0, q1) = q_trace[layer].split_at(n);
    let mut p0_polynomial = MultilinearPolynomial::new(p0.to_vec());
    let mut p1_polynomial = MultilinearPolynomial::new(p1.to_vec());
    let mut q0_polynomial = MultilinearPolynomial::new(q0.to_vec());
    let mut q1_polynomial = MultilinearPolynomial::new(q1.to_vec());

    let mut sc_inst = EqSumCheckInstance::<E>::new(tau.clone());

    let mut rs = vec![];
    let mut polys = vec![];

    dbg!(n.ilog2() as usize);

    for round in 0..n.ilog2() as usize {
      dbg!(round);
      let (eval_point_h_0_0, eval_point_h_2_0, eval_point_h_3_0) =
        sc_inst.evaluation_points_cubic_cross_term(&p0_polynomial, &q1_polynomial);

      let (eval_point_h_0_1, eval_point_h_2_1, eval_point_h_3_1) =
        sc_inst.evaluation_points_cubic_cross_term(&p1_polynomial, &q0_polynomial);

      let (eval_point_h_0_2, eval_point_h_2_2, eval_point_h_3_2) =
        sc_inst.evaluation_points_cubic_cross_term(&q0_polynomial, &q1_polynomial);

      let eval_point_h_0 = eval_point_h_0_0 + eval_point_h_0_1 + eval_point_h_0_2 * rlc_coeff;
      let eval_point_h_2 = eval_point_h_2_0 + eval_point_h_2_1 + eval_point_h_2_2 * rlc_coeff;
      let eval_point_h_3 = eval_point_h_3_0 + eval_point_h_3_1 + eval_point_h_3_2 * rlc_coeff;
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

      sc_inst.bound(&r);
      p0_polynomial.bind_poly_var_top(&r);
      p1_polynomial.bind_poly_var_top(&r);
      q0_polynomial.bind_poly_var_top(&r);
      q1_polynomial.bind_poly_var_top(&r);
    }

    let r0 = transcript.squeeze(b"r")?;
    tau = [vec![r0], rs].concat();

    {
      assert_eq!(p0_polynomial.len(), 1);
      assert_eq!(p1_polynomial.len(), 1);
      assert_eq!(q0_polynomial.len(), 1);
      assert_eq!(q1_polynomial.len(), 1);
    }

    let (p0, p1, q0, q1) = (
      p0_polynomial[0],
      p1_polynomial[0],
      q0_polynomial[0],
      q1_polynomial[0],
    );

    proofs.push(SumcheckProof::new(polys));
    claims.push(LogupGKRLayerClaim { p0, p1, q0, q1 });

    let one = E::Scalar::ONE;
    p_claim = p0_polynomial[0] * (one - r0) + p1_polynomial[0] * r0;
    q_claim = q0_polynomial[0] * (one - r0) + q1_polynomial[0] * r0;

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

  Ok((proofs, tau, claims))
}

#[cfg(test)]
mod tests {
  use crate::provider::Bn256EngineIPA;

  use super::*;

  type E = Bn256EngineIPA;
  type Fr = <E as Engine>::Scalar;

  #[test]
  fn test_logup_gkr() {
    let log_n = 5;
    let n = 1 << log_n;
    let p0 = (1..=n).map(|i| Fr::from(i as u64)).collect::<Vec<_>>();
    let q0 = (1..=n).map(|i| Fr::from(i as u64)).collect::<Vec<_>>();
    let p1 = (1..=n).map(|i| Fr::from(i as u64)).collect::<Vec<_>>();
    let q1 = (1..=n).map(|i| Fr::from(i as u64)).collect::<Vec<_>>();

    let p: Vec<Fr> = [p0.clone(), p1.clone()].concat();
    let q: Vec<Fr> = [q0.clone(), q1.clone()].concat();

    let (p_trace, q_trace) = make_gkr_trace(p, q);

    {
      let mut res = Fr::ZERO;
      for i in 0..n {
        res += p0[i] * q0[i].invert().unwrap();
        res += p1[i] * q1[i].invert().unwrap();
      }
      assert_eq!(res, p_trace[0][0] * q_trace[0][0].invert().unwrap());
    }

    let mut transcript = <E as Engine>::TE::new(b"Test");

    let (proofs, tau, claims) = prove_gkr::<E>(p_trace, q_trace, &mut transcript).unwrap();

    let mut transcript = <E as Engine>::TE::new(b"Test");
    verify_gkr(claims, proofs, &mut transcript).unwrap();
  }
}
