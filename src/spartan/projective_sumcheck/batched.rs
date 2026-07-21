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

use super::{
  eq_factored::EqFactoredVirtualPolynomial, prover::ProjectiveSumcheckProverOutput,
  virtual_poly::VirtualPolynomial,
};

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

    super::absorb_round::<E>(transcript, &poly);
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

/// A heterogeneous projective sumcheck instance for the eq-factored batched
/// prover: either a plain product ([`VirtualPolynomial`], no eq factor) or an
/// eq-factored relation ([`EqFactoredVirtualPolynomial`], analytic eq + U +
/// BDDT). Both expose the same stepwise interface (`degree` / `claim0` /
/// `round_coeffs` / `bind` / `reduced_claim`), letting the batcher fold them
/// into one sumcheck without materializing eq or U as dense factors.
pub enum ProjInstance<E: Engine> {
  /// A plain sum-of-products with no eq factor (e.g. inner-ABC, logup inv).
  Plain(VirtualPolynomial<E>),
  /// An eq-factored relation carrying eq/U analytically (e.g. inner-E, memory).
  EqFactored(EqFactoredVirtualPolynomial<E>),
}

impl<E: Engine> ProjInstance<E> {
  fn num_vars(&self) -> usize {
    match self {
      ProjInstance::Plain(v) => v.num_vars(),
      ProjInstance::EqFactored(v) => v.num_vars(),
    }
  }

  fn degree(&self) -> usize {
    match self {
      ProjInstance::Plain(v) => v.degree(),
      ProjInstance::EqFactored(v) => v.degree(),
    }
  }

  fn claim0(&self) -> E::Scalar {
    match self {
      ProjInstance::Plain(v) => v.claim0(),
      ProjInstance::EqFactored(v) => v.claim0(),
    }
  }

  fn lift_to_degree(&mut self, target: usize) {
    match self {
      ProjInstance::Plain(v) => v.lift_to_degree(target),
      ProjInstance::EqFactored(v) => v.lift_to_degree(target),
    }
  }

  /// This instance's own degree-`D` round message `[a_0, …, a_D]`.
  fn round_coeffs(&self, remaining: usize) -> Vec<E::Scalar> {
    match self {
      ProjInstance::Plain(v) => v.round_full_coeffs(remaining),
      ProjInstance::EqFactored(v) => v.round_message(),
    }
  }

  /// Binds the current variable to `r`; `msg` is this instance's own
  /// `round_coeffs` output (used by the eq-factored path to update its claim).
  fn bind(&mut self, remaining: usize, r: E::Scalar, msg: &[E::Scalar]) {
    match self {
      ProjInstance::Plain(v) => v.bind(remaining, r),
      ProjInstance::EqFactored(v) => v.bind_round(r, msg),
    }
  }

  fn reduced_claim(&self) -> E::Scalar {
    match self {
      ProjInstance::Plain(v) => v.reduced_claim(),
      ProjInstance::EqFactored(v) => v.reduced_claim(),
    }
  }

  fn bound_factor_values(&self) -> Vec<E::Scalar> {
    match self {
      ProjInstance::Plain(v) => v.bound_factor_values(),
      ProjInstance::EqFactored(v) => v.bound_factor_values(),
    }
  }
}

/// Batched projective sumcheck over a **heterogeneous** set of instances,
/// folding them via a λ-power RLC into one execution — the eq-factored analogue
/// of [`prove_batched`]. Eq/U-carrying instances contribute their analytic
/// eq-factored round messages (with BDDT), so no dense eq/U table is bound.
///
/// Preconditions (asserted): non-empty, all instances share `num_vars`. Instances
/// are U-homogenized to the common degree `D = max_i degree_i` (raising an
/// eq-factored instance's degree adds analytic `U` copies; a plain instance is
/// lifted with its existing dense-`U` `lift_to_degree`). Same transcript contract
/// as [`prove_batched`]: the caller/verifier must squeeze `λ` under
/// `b"projective_sumcheck_batch"` immediately before [`super::verifier::verify`].
pub fn prove_batched_mixed<E: Engine>(
  mut instances: Vec<ProjInstance<E>>,
  transcript: &mut E::TE,
) -> BatchedProverOutput<E> {
  assert!(!instances.is_empty(), "batched prover needs >= 1 instance");
  let num_vars = instances[0].num_vars();
  for inst in &instances {
    assert_eq!(inst.num_vars(), num_vars, "all instances share num_vars");
  }

  let degree = instances.iter().map(|i| i.degree()).max().unwrap();
  for inst in &mut instances {
    inst.lift_to_degree(degree);
  }

  let lambda = transcript
    .squeeze(b"projective_sumcheck_batch")
    .expect("transcript squeeze failed");
  let coeffs = crate::spartan::powers::<E>(&lambda, instances.len());

  let initial_claim: E::Scalar = instances
    .iter()
    .zip(coeffs.iter())
    .map(|(inst, c)| inst.claim0() * *c)
    .sum();

  let mut rounds = Vec::with_capacity(num_vars);
  let mut point = Vec::with_capacity(num_vars);

  let mut remaining = num_vars;
  for _ in 0..num_vars {
    // Each instance's own round message, combined by λ^i.
    let msgs: Vec<Vec<E::Scalar>> = instances
      .iter()
      .map(|inst| inst.round_coeffs(remaining))
      .collect();

    let mut combined = vec![E::Scalar::ZERO; degree + 1];
    for (m, c) in msgs.iter().zip(coeffs.iter()) {
      for (dst, src) in combined.iter_mut().zip(m.iter()) {
        *dst += *src * *c;
      }
    }

    let poly = UniPoly::<E::Scalar>::from_coeffs_no_trim(combined)
      .expect("combined round polynomial has D+1 >= 2 coefficients");

    super::absorb_round::<E>(transcript, &poly);
    let r_i = transcript
      .squeeze(b"projective_sumcheck_challenge")
      .expect("transcript squeeze failed");

    rounds.push(poly.compress_projective());
    point.push(r_i);

    for (inst, m) in instances.iter_mut().zip(msgs.iter()) {
      inst.bind(remaining, r_i, m);
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

  /// Natural-order projective eq corner table (matching the eq-factored prover's
  /// convention) for building the naive reference.
  fn eq_nat(taus: &[Fr]) -> Vec<Fr> {
    let n = taus.len();
    let mut e = vec![Fr::ZERO; 1usize << n];
    e[0] = Fr::ONE;
    for (i, &t) in taus.iter().enumerate() {
      let blk = 1usize << i;
      for b in 0..blk {
        let lo = e[b];
        e[b] = lo * (Fr::ONE - t);
        e[b + blk] = lo * t;
      }
    }
    e
  }

  /// Differential test for the heterogeneous eq-factored batcher: a mixed batch
  /// of one eq-factored instance `Eq^∞_τ·(f·g − h·U)` (analytic eq + U + BDDT)
  /// and one plain product `p·q·s` must produce a **byte-identical** proof to the
  /// naive batch that carries eq and U as dense factors.
  #[test]
  fn prove_batched_mixed_matches_naive() {
    use crate::spartan::projective_sumcheck::eq_factored::EqFactoredVirtualPolynomial;

    let num_vars = 5usize;
    let n = 1usize << num_vars;
    let f: Vec<Fr> = (0..n).map(|i| Fr::from((2 * i + 1) as u64)).collect();
    let g: Vec<Fr> = (0..n).map(|i| Fr::from((3 * i + 5) as u64)).collect();
    let h: Vec<Fr> = (0..n).map(|i| Fr::from((7 * i + 2) as u64)).collect();
    let p: Vec<Fr> = (0..n).map(|i| Fr::from((i + 9) as u64)).collect();
    let q: Vec<Fr> = (0..n).map(|i| Fr::from((5 * i + 4) as u64)).collect();
    let s: Vec<Fr> = (0..n).map(|i| Fr::from((4 * i + 3) as u64)).collect();
    let taus: Vec<Fr> = (0..num_vars)
      .map(|i| Fr::from((3 * i + 2) as u64))
      .collect();

    // Mixed: eq-factored A = Eq·(f·g − h·U), plain B = p·q·s.
    let a = EqFactoredVirtualPolynomial::<E>::new_mixed(
      num_vars,
      taus.clone(),
      vec![f.clone(), g.clone(), h.clone()],
      vec![(Fr::ONE, vec![0, 1]), (-Fr::ONE, vec![2])],
    );
    let b = VirtualPolynomial::<E>::new(
      num_vars,
      vec![p.clone(), q.clone(), s.clone()],
      vec![(Fr::ONE, vec![0, 1, 2])],
    );
    let mut ts_m = <E as Engine>::TE::new(b"projsc_mixed");
    let out_m = prove_batched_mixed::<E>(
      vec![ProjInstance::EqFactored(a), ProjInstance::Plain(b)],
      &mut ts_m,
    );

    // Naive: both instances as dense VirtualPolynomials (eq + U dense).
    let eq = eq_nat(&taus);
    let naive_a = VirtualPolynomial::<E>::new_homogenized(
      num_vars,
      vec![eq, f, g, h],
      vec![(Fr::ONE, vec![0, 1, 2]), (-Fr::ONE, vec![0, 3])],
    );
    let naive_b =
      VirtualPolynomial::<E>::new(num_vars, vec![p, q, s], vec![(Fr::ONE, vec![0, 1, 2])]);
    let mut ts_n = <E as Engine>::TE::new(b"projsc_mixed");
    let out_n = prove_batched::<E>(vec![naive_a, naive_b], &mut ts_n);

    assert_eq!(out_m.lambda, out_n.lambda);
    assert_eq!(out_m.initial_claim, out_n.initial_claim);
    assert_eq!(out_m.point, out_n.point);
    assert_eq!(out_m.final_claim, out_n.final_claim);
    assert_eq!(out_m.per_instance_final, out_n.per_instance_final);
    let polys_m = out_m.proof.compressed_polys();
    let polys_n = out_n.proof.compressed_polys();
    assert_eq!(polys_m.len(), polys_n.len());
    for (x, y) in polys_m.iter().zip(polys_n.iter()) {
      assert_eq!(x.stored_coeffs(), y.stored_coeffs());
    }

    // And the mixed proof verifies at D = 3.
    let degree_bounds = vec![3usize; num_vars];
    let mut ts_v = <E as Engine>::TE::new(b"projsc_mixed");
    let _lambda = ts_v.squeeze(b"projective_sumcheck_batch").unwrap();
    let reduction =
      verify::<E>(out_m.initial_claim, &degree_bounds, &out_m.proof, &mut ts_v).unwrap();
    assert_eq!(reduction.point, out_m.point);
    assert_eq!(reduction.final_claim, out_m.final_claim);
  }
}
