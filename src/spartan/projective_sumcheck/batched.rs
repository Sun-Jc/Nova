//! Batched Projective SumCheck prover (Stage 3): folds several
//! [`VirtualPolynomial`] instances into a **single** sumcheck execution via a
//! random-linear-combination over powers of a challenge `λ`, mirroring
//! ppSNARK's `prove_helper` (`ppsnark.rs`).
//!
//! Given instances `G_0, ..., G_{k-1}` over the same number of variables `n`
//! **and the same declared degree `D`** (U-homogenize lower-degree instances
//! first, since zero-padding would break their projective endpoint), the
//! batcher squeezes `λ`, forms `coeffs = [1, λ, λ², ...]`, and proves the joint
//! claim `Σ_i λ^i · C_0^{(i)}` for the combined polynomial `Σ_i λ^i · G_i`.
//! Each round it combines the per-instance round coefficient vectors by `λ^i`,
//! absorbs the single combined message, squeezes the shared challenge `r_i`,
//! and binds every instance at `r_i`. All instances share one point `r`.
//!
//! The combined message is a standard projective round polynomial, so the
//! **unmodified** [`super::verifier::verify`] verifies the batched proof at
//! degree `D = max_i degree_i`.

use crate::{
  spartan::{polys::univariate::UniPoly, sumcheck::SumcheckProof},
  traits::{Engine, TranscriptEngineTrait},
};
use ff::Field;

use super::{prover::ProjectiveSumcheckProverOutput, virtual_poly::VirtualPolynomial};

/// Output of the batched prover: the combined proof, the batching challenge
/// `λ`, the sampled point, and the per-instance reduced claims `G_i(r)` (the
/// caller checks each against its own factor openings).
#[derive(Clone, Debug)]
pub struct BatchedProverOutput<E: Engine> {
  /// The combined proof, verifiable at degree `D = max_i degree_i`.
  pub proof: SumcheckProof<E>,
  /// The batching challenge `λ`.
  pub lambda: E::Scalar,
  /// The joint initial claim `Σ_i λ^i C_0^{(i)}`.
  pub initial_claim: E::Scalar,
  /// The sampled point `r`.
  pub point: Vec<E::Scalar>,
  /// The joint reduced claim `Σ_i λ^i G_i(r)`.
  pub final_claim: E::Scalar,
  /// Per-instance reduced claims `G_i(r)`, in input order.
  pub per_instance_final: Vec<E::Scalar>,
  /// Per-instance bound factor values `F_j(r)` (factor order), one inner vector
  /// per instance. Lets the caller reuse the sumcheck's own binding as each
  /// factor's opening at `r`, avoiding an O(N) MLE re-evaluation per column.
  pub per_instance_bound_factors: Vec<Vec<E::Scalar>>,
}

/// Proves a batch of [`VirtualPolynomial`] instances in one sumcheck execution.
///
/// Preconditions (asserted): the batch is non-empty and all instances have the
/// same `num_vars`. The instances' per-instance initial claims are assumed to
/// be already absorbed by the caller (matching ppSNARK's `prove_helper`
/// contract); the batcher squeezes `λ` immediately.
///
/// **Transcript contract for verification.** The batcher squeezes the batching
/// challenge `λ` under label `b"projective_sumcheck_batch"` *before* the shared
/// rounds. The unmodified [`super::verifier::verify`] does not know about `λ`,
/// so the caller/verifier MUST squeeze the same challenge from its transcript
/// (same label) immediately before calling `verify` — exactly as ppSNARK's
/// verifier mirrors `prove_helper`'s challenge squeeze. After that squeeze, the
/// combined proof is a standard projective proof at degree `D = max_i degree_i`.
pub fn prove_batched<E: Engine>(
  mut instances: Vec<VirtualPolynomial<E>>,
  transcript: &mut E::TE,
) -> BatchedProverOutput<E> {
  assert!(!instances.is_empty(), "batched prover needs >= 1 instance");
  let num_vars = instances[0].num_vars();
  for inst in &instances {
    assert_eq!(inst.num_vars(), num_vars, "all instances share num_vars");
  }
  // Bring every instance to the common degree D = max_i D_i by U-homogenization
  // (NOT zero-padding, which would break a lower-degree instance's projective
  // endpoint: under declared D, [T^D]S = 0 gives a_0 + 0, not the true claim
  // a_0 + a_{D_i}; design §5, §17.3). This matches ppSNARK's prove_helper, which
  // requires equal degree across the folded engines.
  let degree = instances.iter().map(|i| i.degree()).max().unwrap();
  for inst in &mut instances {
    inst.lift_to_degree(degree);
  }

  // λ and its powers, one per instance.
  let lambda = transcript
    .squeeze(b"projective_sumcheck_batch")
    .expect("transcript squeeze failed");
  let coeffs = crate::spartan::powers::<E>(&lambda, instances.len());

  // Joint initial claim Σ_i λ^i C_0^{(i)}.
  let initial_claim: E::Scalar = instances
    .iter()
    .zip(coeffs.iter())
    .map(|(inst, c)| inst.claim0() * *c)
    .sum();

  let mut rounds = Vec::with_capacity(num_vars);
  let mut point = Vec::with_capacity(num_vars);

  let mut remaining = num_vars;
  for _ in 0..num_vars {
    // Combine each instance's round coefficients (padded to D) by λ^i.
    let mut combined = vec![E::Scalar::ZERO; degree + 1];
    for (inst, c) in instances.iter().zip(coeffs.iter()) {
      let round = inst.round_full_coeffs(remaining);
      for (dst, src) in combined.iter_mut().zip(round.iter()) {
        *dst += *src * *c;
      }
    }

    let poly = UniPoly::<E::Scalar>::from_coeffs_no_trim(combined)
      .expect("combined round polynomial has D+1 >= 2 coefficients");

    transcript.absorb(b"projective_sumcheck_round", &poly);
    let r_i = transcript
      .squeeze(b"projective_sumcheck_challenge")
      .expect("transcript squeeze failed");

    rounds.push(poly.compress_projective());
    point.push(r_i);

    for inst in &mut instances {
      inst.bind(remaining, r_i);
    }
    remaining -= 1;
  }

  let per_instance_final: Vec<E::Scalar> = instances.iter().map(|i| i.reduced_claim()).collect();
  let per_instance_bound_factors: Vec<Vec<E::Scalar>> =
    instances.iter().map(|i| i.bound_factor_values()).collect();
  let final_claim: E::Scalar = per_instance_final
    .iter()
    .zip(coeffs.iter())
    .map(|(f, c)| *f * *c)
    .sum();

  BatchedProverOutput {
    proof: SumcheckProof::new(rounds),
    lambda,
    initial_claim,
    point,
    final_claim,
    per_instance_final,
    per_instance_bound_factors,
  }
}

impl<E: Engine> BatchedProverOutput<E> {
  /// Convenience view of the joint reduction as a single-instance output.
  pub fn as_single(&self) -> ProjectiveSumcheckProverOutput<E> {
    ProjectiveSumcheckProverOutput {
      proof: self.proof.clone(),
      initial_claim: self.initial_claim,
      point: self.point.clone(),
      final_claim: self.final_claim,
    }
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::{provider::PallasEngine, spartan::projective_sumcheck::verify, traits::Engine};
  use ff::Field;

  type E = PallasEngine;
  type Fr = <E as Engine>::Scalar;

  /// Batch two instances over the same variables, both at degree D=3: a plain
  /// f·g·h product, and a p·q product U-homogenized up to degree 3. The combined
  /// proof must verify at D=3, and the verifier's reduced claim must equal the
  /// batcher's joint final claim (= Σ_i λ^i G_i(r)).
  #[test]
  fn batch_mixed_degree_instances() {
    let num_vars = 3usize;
    let n = 1usize << num_vars;

    let f: Vec<Fr> = (0..n).map(|i| Fr::from((2 * i + 1) as u64)).collect();
    let g: Vec<Fr> = (0..n).map(|i| Fr::from((3 * i + 5) as u64)).collect();
    let h: Vec<Fr> = (0..n).map(|i| Fr::from((7 * i + 2) as u64)).collect();
    let p: Vec<Fr> = (0..n).map(|i| Fr::from((i + 9) as u64)).collect();
    let q: Vec<Fr> = (0..n).map(|i| Fr::from((5 * i + 4) as u64)).collect();

    let inst0 =
      VirtualPolynomial::<E>::new(num_vars, vec![f, g, h], vec![(Fr::ONE, vec![0, 1, 2])]); // D=3
    let inst1 = VirtualPolynomial::<E>::new(num_vars, vec![p, q], vec![(Fr::ONE, vec![0, 1])]); // D=2; batcher lifts to 3

    let mut ts = <E as Engine>::TE::new(b"projsc_batch");
    let out = prove_batched::<E>(vec![inst0, inst1], &mut ts);

    // Joint final = Σ λ^i G_i(r).
    let coeffs = crate::spartan::powers::<E>(&out.lambda, 2);
    let joint: Fr = out
      .per_instance_final
      .iter()
      .zip(coeffs.iter())
      .map(|(v, c)| *v * *c)
      .sum();
    assert_eq!(out.final_claim, joint);

    // The batched proof verifies at D = 3. The verifier transcript must mirror
    // the prover's λ squeeze (the caller owns this batching challenge, exactly
    // as ppSNARK's prove_helper squeezes it before the shared rounds).
    let degree_bounds = vec![3usize; num_vars];
    let mut ts_ver = <E as Engine>::TE::new(b"projsc_batch");
    let lambda_ver = ts_ver.squeeze(b"projective_sumcheck_batch").unwrap();
    assert_eq!(lambda_ver, out.lambda);
    let reduction =
      verify::<E>(out.initial_claim, &degree_bounds, &out.proof, &mut ts_ver).unwrap();
    assert_eq!(reduction.point, out.point);
    assert_eq!(reduction.final_claim, out.final_claim);
  }
}
