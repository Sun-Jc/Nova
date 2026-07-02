//! Prover for the Logup-GKR fractional-sum argument.
//!
//! **This prover satisfies the verifier.** The protocol — transcript order,
//! challenge labels, and the per-layer sumcheck contract — is defined in
//! `verifier.rs`; this file produces a proof the verifier accepts, importing the
//! verifier's `spec` labels and `absorb_fraction` rather than restating them.
//!
//! It builds one batched tree over all logup instances, folds it leaf→root, and
//! per internal layer runs one transparent cubic sumcheck
//! ([`prove_layer_sumcheck`]) reducing the merged fraction-sum claim to the next
//! layer.
//!
//! ## Per-layer gate (the verifier's contract)
//! A layer proves, for all instances batched by a per-layer `λ`,
//! `Σ_i (v_p_i + λ v_q_i) = Σ_x eq(τ,x) · Σ_i [nL·dR + nR·dL + λ·dL·dR]_i`.
//! The `eq(τ, ·)` factor is carried as an explicit MLE so the sumcheck's final
//! value reconciles as the verifier demands: `eq(τ,r)·G(r)`. Round polynomials
//! are degree 3 (`eq` deg 1 × gate deg 2 = `spec::LAYER_SC_DEGREE`).
//!
//! ## Endianness
//! The tree folds MSB-first (`Layer::fold_up`, halves `i`/`i+n`), and the
//! sumcheck binds the top variable first, so the sumcheck point is directly the
//! fraction evaluation point of the next layer.

use crate::errors::NovaError;
use crate::spartan::logup_gkr::layer::Layer;
use crate::spartan::logup_gkr::proof::{
  LayerClaim, LayerFinalClaim, LayerSumcheck, LogupGkrOpeningClaim, LogupGkrProof,
};
use crate::spartan::polys::multilinear::MultilinearPolynomial;
use crate::traits::{Engine, TranscriptEngineTrait};
use ff::Field;

// The prover satisfies the protocol defined by the verifier: it uses the
// verifier's transcript labels (`spec`) and fraction-absorb order, and produces
// a sumcheck whose final value reconciles as `eq(τ,r)·G(r)` (the verifier's
// contract).
use crate::spartan::logup_gkr::verifier::{absorb_fraction, spec};

/// One instance's four child half-MLEs at a layer: `nL, nR, dL, dR`, each of
/// length `n` (the halves of the child layer of length `2n`).
struct Halves<E: Engine> {
  nl: MultilinearPolynomial<E::Scalar>,
  nr: MultilinearPolynomial<E::Scalar>,
  dl: MultilinearPolynomial<E::Scalar>,
  dr: MultilinearPolynomial<E::Scalar>,
}

/// Transparent cubic sumcheck for one GKR layer, proving
/// `claim = Σ_x eq(τ, x) · Σ_i [ nLᵢ·dRᵢ + nRᵢ·dLᵢ + λ·dLᵢ·dRᵢ ](x)`
/// over `num_rounds = |τ|` variables. Unlike `prove_batched_cubic`, the `eq`
/// factor here is an explicit MLE (`EqPolynomial(τ).evals()`), so the verifier
/// reconciles the final evaluation as the transparent `eq(τ,r) · G(r)` — the
/// same shape `prove_cubic_with_three_inputs` uses (ppsnark verify).
///
/// Returns the compressed round polynomials, the sumcheck point `r`, and each
/// instance's `(nL, nR, dL, dR)` evaluated at `r`.
#[allow(clippy::type_complexity)]
fn prove_layer_sumcheck<E: Engine>(
  claim: E::Scalar,
  taus: &[E::Scalar],
  halves: &mut [Halves<E>],
  lambda: E::Scalar,
  transcript: &mut E::TE,
) -> Result<
  (
    Vec<crate::spartan::polys::univariate::CompressedUniPoly<E::Scalar>>,
    Vec<E::Scalar>,
    Vec<[E::Scalar; 4]>,
  ),
  NovaError,
> {
  use crate::spartan::polys::eq::EqPolynomial;
  use crate::spartan::polys::univariate::{CompressedUniPoly, UniPoly};

  let num_rounds = taus.len();
  let mut eq = MultilinearPolynomial::new(EqPolynomial::new(taus.to_vec()).evals());

  let mut r: Vec<E::Scalar> = Vec::with_capacity(num_rounds);
  let mut polys: Vec<CompressedUniPoly<E::Scalar>> = Vec::with_capacity(num_rounds);
  let mut claim_per_round = claim;

  for _ in 0..num_rounds {
    let len = eq.len();
    let half = len / 2;

    // Evaluate the round polynomial P(t) = Σ_x eq_t(x)·G_t(x) at t = 0,1,2,3,
    // where at parameter t each MLE m contributes m0 + t·(m1 − m0) (m0 = low
    // half, m1 = high half) — the MSB-first bind direction.
    let mut p = [E::Scalar::ZERO; 4];
    for x in 0..half {
      let e0 = eq.Z[x];
      let e1 = eq.Z[x + half];
      let e_step = e1 - e0;
      // per-instance gate contributions
      let mut g = [E::Scalar::ZERO; 4];
      for h in halves.iter() {
        let nl0 = h.nl.Z[x];
        let nl1 = h.nl.Z[x + half];
        let nr0 = h.nr.Z[x];
        let nr1 = h.nr.Z[x + half];
        let dl0 = h.dl.Z[x];
        let dl1 = h.dl.Z[x + half];
        let dr0 = h.dr.Z[x];
        let dr1 = h.dr.Z[x + half];
        for (k, &t) in [
          E::Scalar::ZERO,
          E::Scalar::ONE,
          E::Scalar::from(2u64),
          E::Scalar::from(3u64),
        ]
        .iter()
        .enumerate()
        {
          let nl = nl0 + t * (nl1 - nl0);
          let nr = nr0 + t * (nr1 - nr0);
          let dl = dl0 + t * (dl1 - dl0);
          let dr = dr0 + t * (dr1 - dr0);
          g[k] += nl * dr + nr * dl + lambda * (dl * dr);
        }
      }
      for k in 0..4 {
        let t = E::Scalar::from(k as u64);
        let e = e0 + t * e_step;
        p[k] += e * g[k];
      }
    }

    let poly = UniPoly::from_evals(&p);
    transcript.absorb(spec::ROUND_POLY, &poly);
    let r_i = transcript.squeeze(spec::ROUND_CHALLENGE)?;
    r.push(r_i);
    polys.push(poly.compress());
    claim_per_round = poly.evaluate(&r_i);

    // Bind the top variable of eq and every half-MLE.
    eq.bind_poly_var_top(&r_i);
    for h in halves.iter_mut() {
      h.nl.bind_poly_var_top(&r_i);
      h.nr.bind_poly_var_top(&r_i);
      h.dl.bind_poly_var_top(&r_i);
      h.dr.bind_poly_var_top(&r_i);
    }
  }

  let _ = claim_per_round;
  let finals: Vec<[E::Scalar; 4]> = halves
    .iter()
    .map(|h| [h.nl.Z[0], h.nr.Z[0], h.dl.Z[0], h.dr.Z[0]])
    .collect();
  Ok((polys, r, finals))
}

/// Proves the fractional-sum identity `Σ p/q = root` for all instances in one
/// batched tree (`inputs` holds one input `Layer` per instance, `[row, col]`),
/// returning the proof and the shared opening claim.
///
/// All instances must have the same height (power of two, ≥ 2). The soundness
/// red line — absorbing every leaf/root/claim before sampling a challenge — is
/// enforced here: `initial_claims` and each layer's `final_claims` are absorbed
/// before the next challenge is drawn.
pub fn prove<E: Engine>(
  inputs: Vec<Layer<E>>,
  transcript: &mut E::TE,
) -> Result<(LogupGkrProof<E>, LogupGkrOpeningClaim<E>), NovaError> {
  let m = inputs.len();
  if m == 0 {
    return Err(NovaError::InvalidNumInstances);
  }
  let num_vars = inputs[0].num_vars();
  if num_vars == 0 {
    // A single-cell input has nothing to fold; the argument needs height ≥ 2.
    return Err(NovaError::InvalidSumcheckProof);
  }
  for inp in &inputs {
    if inp.num_vars() != num_vars {
      return Err(NovaError::InvalidSumcheckProof);
    }
  }

  // Build every instance's tree, leaf→root. trees[instance][layer], where
  // trees[i][j] has (num_vars - j) variables; [0] = input, [num_vars] = root.
  let trees: Vec<Vec<Layer<E>>> = inputs.into_iter().map(|l| l.build_tree()).collect();

  // initial_claims = each instance's output (root) fraction, absorbed first.
  let initial_claims: Vec<LayerClaim<E>> = trees
    .iter()
    .map(|t| {
      let (n, d) = t[num_vars].output_fraction();
      LayerClaim::<E>::new(n, d)
    })
    .collect();
  for c in &initial_claims {
    absorb_fraction::<E>(transcript, c.num, c.den);
  }

  // Running per-instance claims (v_p_i, v_q_i) about `layers[j]` at `eval_point`.
  // Start at the root (j = num_vars, empty point).
  let mut running: Vec<(E::Scalar, E::Scalar)> =
    initial_claims.iter().map(|c| (c.num, c.den)).collect();
  let mut eval_point: Vec<E::Scalar> = Vec::new();

  // Ordered output→input as we go: step j reduces a claim about layers[j] to
  // layers[j-1], for j = num_vars, num_vars-1, ..., 1.
  let mut sumchecks: Vec<LayerSumcheck<E>> = Vec::with_capacity(num_vars.saturating_sub(1));
  let mut final_claims_by_layer: Vec<Vec<LayerFinalClaim<E>>> = Vec::with_capacity(num_vars);

  for j in (1..=num_vars).rev() {
    // Fresh λ per layer, after the previous claims were absorbed.
    let lambda = transcript.squeeze(spec::LAMBDA)?;

    // Children are the two halves of layers[j-1] (which has num_vars-j+1 vars,
    // so each half has num_vars-j vars). nL/nR = num halves, dL/dR = den halves.
    let child = |i: usize| &trees[i][j - 1];
    let child_len = 1usize << (num_vars - j + 1);
    let n = child_len / 2;

    // final claims (nL, nR, dL, dR) per instance for this layer.
    let mut layer_finals: Vec<LayerFinalClaim<E>> = Vec::with_capacity(m);

    if num_vars - j == 0 {
      // Base case (j = num_vars, root reduction): the child layer has exactly
      // two cells; read the split directly, no sumcheck.
      for i in 0..m {
        let c = child(i);
        layer_finals.push(LayerFinalClaim::<E>::new(
          c.num.Z[0], // nL
          c.num.Z[1], // nR
          c.den.Z[0], // dL
          c.den.Z[1], // dR
        ));
      }
      let _ = (lambda, n); // λ unused at the base (single term, no batching)
    } else {
      // Sumcheck over (num_vars - j) variables reduces
      //   Σ_i (v_p_i + λ v_q_i) = Σ_x eq(point, x)·Σ_i [nL·dR + nR·dL + λ·dL·dR]_i
      let claim: E::Scalar = running.iter().map(|(p, q)| *p + lambda * *q).sum();

      let mut halves: Vec<Halves<E>> = Vec::with_capacity(m);
      for i in 0..m {
        let c = child(i);
        halves.push(Halves {
          nl: MultilinearPolynomial::new(c.num.Z[..n].to_vec()),
          nr: MultilinearPolynomial::new(c.num.Z[n..child_len].to_vec()),
          dl: MultilinearPolynomial::new(c.den.Z[..n].to_vec()),
          dr: MultilinearPolynomial::new(c.den.Z[n..child_len].to_vec()),
        });
      }

      let (round_polys, r, finals) =
        prove_layer_sumcheck::<E>(claim, &eval_point, &mut halves, lambda, transcript)?;

      for f in &finals {
        layer_finals.push(LayerFinalClaim::<E>::new(f[0], f[1], f[2], f[3]));
      }
      sumchecks.push(LayerSumcheck { round_polys });
      eval_point = r; // sumcheck point (length num_vars - j)
    }

    // Absorb final claims, sample fold challenge, update running claims + point.
    for fc in &layer_finals {
      absorb_fraction::<E>(transcript, fc.left.num, fc.left.den);
      absorb_fraction::<E>(transcript, fc.right.num, fc.right.den);
    }
    let fold_r = transcript.squeeze(spec::FOLD)?;

    running = layer_finals
      .iter()
      .map(|fc| {
        let c = fc.fold_into_next_claim(fold_r);
        (c.num, c.den)
      })
      .collect();
    // New point for layers[j-1] is [fold_r, ...sumcheck_point] (fold_r as MSB).
    let mut next_point = Vec::with_capacity(eval_point.len() + 1);
    next_point.push(fold_r);
    next_point.extend_from_slice(&eval_point);
    eval_point = next_point;

    final_claims_by_layer.push(layer_finals);
  }

  // Proof is ordered output→input (the order we produced), matching hp.
  // openings = each instance's input-layer fraction at the final eval_point
  // (length num_vars); equals the last running claim by construction.
  let openings: Vec<LayerClaim<E>> = running
    .iter()
    .map(|(n, d)| LayerClaim::<E>::new(*n, *d))
    .collect();

  let proof = LogupGkrProof {
    initial_claims,
    final_claims: final_claims_by_layer,
    sumchecks,
  };
  let claim = LogupGkrOpeningClaim::new(eval_point, openings);
  Ok((proof, claim))
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::spartan::logup_gkr::verifier;
  use crate::spartan::polys::multilinear::MultilinearPolynomial;
  use crate::traits::TranscriptEngineTrait;

  type E = crate::provider::Bn256EngineKZG;
  type Fr = <E as Engine>::Scalar;

  fn mle(v: Vec<u64>) -> MultilinearPolynomial<Fr> {
    MultilinearPolynomial::new(v.into_iter().map(Fr::from).collect())
  }

  fn frac_eq(n0: Fr, d0: Fr, n1: Fr, d1: Fr) -> bool {
    n0 * d1 == n1 * d0
  }

  // Full prove→verify round trip: the verifier must accept, and its returned
  // openings must equal the prover's actual input-layer fractions at the shared
  // evaluation point.
  fn round_trip(inputs: Vec<(Vec<u64>, Vec<u64>)>) {
    let layers: Vec<Layer<E>> = inputs
      .iter()
      .map(|(n, d)| Layer::<E> {
        num: mle(n.clone()),
        den: mle(d.clone()),
      })
      .collect();
    // Keep copies to independently evaluate at the eval_point.
    let raw: Vec<(MultilinearPolynomial<Fr>, MultilinearPolynomial<Fr>)> = inputs
      .iter()
      .map(|(n, d)| (mle(n.clone()), mle(d.clone())))
      .collect();

    let mut tr_p = <E as Engine>::TE::new(b"gkr-test");
    let (proof, claim) = prove::<E>(layers, &mut tr_p).expect("prove");

    let mut tr_v = <E as Engine>::TE::new(b"gkr-test");
    let vclaim = verifier::verify::<E>(&proof, &mut tr_v).expect("verify");

    // Verifier's eval_point and openings must match the prover's.
    assert_eq!(
      vclaim.eval_point(),
      claim.eval_point(),
      "eval_point mismatch"
    );
    let pt = vclaim.eval_point();
    for (i, (n, d)) in raw.iter().enumerate() {
      let en = n.evaluate(pt);
      let ed = d.evaluate(pt);
      let op = vclaim.openings()[i];
      assert!(
        frac_eq(en, ed, op.num, op.den),
        "opening[{i}] must equal input layer fraction at eval_point"
      );
    }
  }

  // NOTE: the three round_trip tests below currently FAIL and are #[ignore]d.
  // Blocker: `SumcheckProof::prove_batched_cubic`'s final evaluation does not
  // reconcile against the naive `eq(taus,r)·Σα·A·B·C` (verified by probe: even
  // 1-variable does not match), whereas `prove_cubic_with_three_inputs` does.
  // The prover/verifier framework is otherwise complete and hp-aligned; the fix
  // is to run the layer sumcheck on a transparent-eq backend. See
  // jcbase/stage23-progress.md. Do NOT delete these — they are the acceptance
  // test for the fix.
  #[test]
  fn round_trip_single_instance_n4() {
    round_trip(vec![(vec![1, 2, 3, 4], vec![5, 6, 7, 8])]);
  }

  #[test]
  fn round_trip_two_instances_n4() {
    round_trip(vec![
      (vec![1, 2, 3, 4], vec![5, 6, 7, 8]),
      (vec![9, 8, 7, 6], vec![2, 3, 4, 5]),
    ]);
  }

  #[test]
  fn round_trip_two_instances_n8() {
    round_trip(vec![
      (vec![3, 1, 4, 1, 5, 9, 2, 6], vec![2, 7, 1, 8, 2, 8, 1, 8]),
      (vec![1, 1, 1, 1, 1, 1, 1, 1], vec![3, 1, 4, 1, 5, 9, 2, 6]),
    ]);
  }

  #[test]
  fn verify_rejects_tampered_final_claim() {
    let layers = vec![Layer::<E> {
      num: mle(vec![1, 2, 3, 4]),
      den: mle(vec![5, 6, 7, 8]),
    }];
    let mut tr_p = <E as Engine>::TE::new(b"gkr-test");
    let (mut proof, _) = prove::<E>(layers, &mut tr_p).unwrap();
    // Corrupt one final claim.
    proof.final_claims[0][0].left.num += Fr::from(1);
    let mut tr_v = <E as Engine>::TE::new(b"gkr-test");
    assert!(verifier::verify::<E>(&proof, &mut tr_v).is_err());
  }
}
