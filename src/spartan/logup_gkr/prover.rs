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
use crate::spartan::logup_gkr::fraction::Fraction;
use crate::spartan::logup_gkr::layer::Layer;
use crate::spartan::logup_gkr::proof::{
  LayerClaim, LayerFinalClaim, LayerSumcheck, LogupGkrOpeningClaim, LogupGkrProof,
};
use crate::spartan::polys::multilinear::MultilinearPolynomial;
use crate::traits::{Engine, TranscriptEngineTrait};
use ff::Field;
use rayon::prelude::*;

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

// Returns: (n0, d0, n_lead, d_lead)
fn compute_eval_gate<Scalar: Field>(
  nl0: Scalar,
  nl1: Scalar,
  nr0: Scalar,
  nr1: Scalar,
  dl0: Scalar,
  dl1: Scalar,
  dr0: Scalar,
  dr1: Scalar,
) -> ((Scalar, Scalar), (Scalar, Scalar)) {
  let n0 = nl0 * dr0 + nr0 * dl0;
  let d0 = dl0 * dr0;

  let nl_lead = nl1 - nl0;
  let dr_lead = dr1 - dr0;
  let nr_lead = nr1 - nr0;
  let dl_lead = dl1 - dl0;
  let n_lead = nl_lead * dr_lead + nr_lead * dl_lead;
  let d_lead = dl_lead * dr_lead;

  ((n0, d0), (n_lead, d_lead))
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
  let mut running_claim = claim;

  // Per-instance λ-powers (num gets λ^{2i}, den λ^{2i+1}) are hoisted out of
  // the x-loop, and the sum over x is parallelized (each x is independent; the
  // per-x [4] contributions reduce by elementwise add).
  // λ^{2i} accumulates by multiplying λ² each step, replacing a per-instance
  // pow_vartime (O(n·log) field muls) with a single running product (O(n)).
  let mut w_num = E::Scalar::ONE;
  let weights: Vec<(E::Scalar, E::Scalar)> = (0..halves.len())
    .map(|_| {
      let w_den = w_num * lambda;
      let pair = (w_num, w_den);
      w_num = w_den * lambda;
      pair
    })
    .collect();

  for _ in 0..num_rounds {
    let len = eq.len();
    let half = len / 2;

    let mut acc_0 = E::Scalar::ZERO;
    let mut acc_lead = E::Scalar::ZERO;

    let (eq_first, eq_last) = eq.Z.split_at(half);

    let mut e = E::Scalar::ZERO;
    let mut f = E::Scalar::ZERO;
    for (eq0, eq1) in eq_first.iter().zip(eq_last.iter()) {
      e += *eq0;
      f += *eq1 - *eq0;
    }

    for (h, &(w_num, w_den)) in halves.iter().zip(weights.iter()) {
      let (nl_first, nl_last) = h.nl.Z.split_at(half);
      let (dl_first, dl_last) = h.dl.Z.split_at(half);
      let (nr_first, nr_last) = h.nr.Z.split_at(half);
      let (dr_first, dr_last) = h.dr.Z.split_at(half);

      for (((((((nl0, dl0), nr0), dr0), nl1), dl1), nr1), dr1) in nl_first
        .iter()
        .zip(dl_first.iter())
        .zip(nr_first.iter())
        .zip(dr_first.iter())
        .zip(nl_last.iter())
        .zip(dl_last.iter())
        .zip(nr_last.iter())
        .zip(dr_last.iter())
      {
        let ((n0, d0), (n_lead, d_lead)) =
          compute_eval_gate(*nl0, *nl1, *nr0, *nr1, *dl0, *dl1, *dr0, *dr1);

        acc_0 += n0 * w_num + d0 * w_den;
        acc_lead += n_lead * w_num + d_lead * w_den;
      }
    }
    // (E + F x) * (A x^2 + B x + C)
    let c = acc_0;
    let a = acc_lead;
    let full_eval_1 = running_claim - c * e;
    let b = full_eval_1 * (f + e).invert().unwrap() - a - c;

    let coeff_3 = a * f;
    let coeff_2 = b * f + a * e;
    let coeff_1 = f * c + e * b;
    let coeff_0 = c * e;
    let uni_poly = UniPoly::from_coeffs(vec![coeff_0, coeff_1, coeff_2, coeff_3]).unwrap();

    // x3: A * F
    // x2: B F + A E
    // x: F C + E B
    // 1: C E

    // A, C, E, F
    // B?
    // (E + F)  ( A + B + C ) = K
    // B = K / (E + F) - A - C

    // Evaluate the round polynomial P(t) = Σ_x eq_t(x)·G_t(x) at t = 0,1,2,3,
    // where at parameter t each MLE m contributes m0 + t·(m1 - m0) (m0 = low
    // half, m1 = high half) — the MSB-first bind direction.
    //

    let ts = [
      E::Scalar::ZERO,
      E::Scalar::ONE,
      E::Scalar::from(2u64),
      E::Scalar::from(3u64),
    ];
    let p = (0..half)
      .into_par_iter()
      .map(|x| {
        let e0 = eq.Z[x];
        let e_step = eq.Z[x + half] - e0;
        // per-instance gate contributions, batched by distinct powers of λ.
        let mut g = [E::Scalar::ZERO; 4];
        for (h, &(w_num, w_den)) in halves.iter().zip(weights.iter()) {
          let nl0 = h.nl.Z[x];
          let nl1 = h.nl.Z[x + half];
          let nr0 = h.nr.Z[x];
          let nr1 = h.nr.Z[x + half];
          let dl0 = h.dl.Z[x];
          let dl1 = h.dl.Z[x + half];
          let dr0 = h.dr.Z[x];
          let dr1 = h.dr.Z[x + half];

          let mut g_tmp = [E::Scalar::ZERO; 4];

          for (k, &t) in ts.iter().enumerate() {
            let nl = nl0 + t * (nl1 - nl0);
            let nr = nr0 + t * (nr1 - nr0);
            let dl = dl0 + t * (dl1 - dl0);
            let dr = dr0 + t * (dr1 - dr0);
            // gate = fraction-add of the two children; batched value =
            // λ^{2i}·gate.num + λ^{2i+1}·gate.den.
            let gate = Fraction::new(nl, dl) + Fraction::new(nr, dr);
            g[k] += w_num * gate.num + w_den * gate.den;

            g_tmp[k] = w_num * gate.num + w_den * gate.den;
          }

          {
            let nl_delta = nl1 - nl0;
            let nr_delta = nr1 - nr0;
            let dl_delta = dl1 - dl0;
            let dr_delta = dr1 - dr0;
            let nl2 = nl1 + nl_delta;
            let nl3 = nl2 + nl_delta;
            let nr2 = nr1 + nr_delta;
            let nr3 = nr2 + nr_delta;
            let dl2 = dl1 + dl_delta;
            let dl3 = dl2 + dl_delta;
            let dr2 = dr1 + dr_delta;
            let dr3 = dr2 + dr_delta;
            let mut f_rec: Vec<Fraction<E::Scalar>> = Vec::with_capacity(4);
            f_rec.push(Fraction::new(nl0, dl0) + Fraction::new(nr0, dr0));
            f_rec.push(Fraction::new(nl1, dl1) + Fraction::new(nr1, dr1));
            f_rec.push(Fraction::new(nl2, dl2) + Fraction::new(nr2, dr2));
            f_rec.push(Fraction::new(nl3, dl3) + Fraction::new(nr3, dr3));
            let mut g_rec = [E::Scalar::ZERO; 4];
            g_rec[0] = f_rec[0].num * w_num + f_rec[0].den * w_den;
            g_rec[1] = f_rec[1].num * w_num + f_rec[1].den * w_den;
            g_rec[2] = f_rec[2].num * w_num + f_rec[2].den * w_den;
            g_rec[3] = f_rec[3].num * w_num + f_rec[3].den * w_den;
            assert_eq!(g_rec, g_tmp);
          }
        }
        let mut px = [E::Scalar::ZERO; 4];
        for (k, &t) in ts.iter().enumerate() {
          px[k] = (e0 + t * e_step) * g[k];
        }
        px
      })
      .reduce(
        || [E::Scalar::ZERO; 4],
        |mut acc, px| {
          for k in 0..4 {
            acc[k] += px[k];
          }
          acc
        },
      );

    let poly = UniPoly::from_evals(&p);

    {
      assert_eq!(poly, uni_poly);
    }

    transcript.absorb(spec::ROUND_POLY, &poly);
    let r_i = transcript.squeeze(spec::ROUND_CHALLENGE)?;
    r.push(r_i);
    polys.push(poly.compress());
    running_claim = poly.evaluate(&r_i);

    // Bind the top variable of eq and every half-MLE.
    eq.bind_poly_var_top(&r_i);
    for h in halves.iter_mut() {
      h.nl.bind_poly_var_top(&r_i);
      h.nr.bind_poly_var_top(&r_i);
      h.dl.bind_poly_var_top(&r_i);
      h.dr.bind_poly_var_top(&r_i);
    }
  }

  let _ = running_claim;
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
    absorb_fraction::<E>(transcript, *c);
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
      // Sumcheck over (num_vars - j) variables. The 2m sub-claims are batched by
      // distinct powers of λ (Horner over [p_0,q_0,p_1,q_1,...]) — see verifier.
      let claim: E::Scalar = {
        let mut acc = E::Scalar::ZERO;
        let mut pw = E::Scalar::ONE;
        for (p, q) in &running {
          acc += pw * *p;
          pw *= lambda;
          acc += pw * *q;
          pw *= lambda;
        }
        acc
      };

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
      absorb_fraction::<E>(transcript, fc.left);
      absorb_fraction::<E>(transcript, fc.right);
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

  // Full prove -> verify round trips over the transparent-eq layer sumcheck.
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

  // Tamper-rejection with m >= 2: mutating per-instance final claims makes the
  // verifier reject. NOTE: this is a general tamper test, NOT a discriminating
  // C1 test — the offset below is also rejected under the old linear-lambda
  // batching (caught by transcript binding + the gate's bilinearity). The
  // lambda-batching was still changed to Horner (distinct powers per instance)
  // to align with hp and remove the non-injective RLC; see
  // jcbase/audit-verifier-vs-hp.md.
  #[test]
  fn verify_rejects_complementary_instance_offset() {
    let inputs = [
      (vec![1u64, 2, 3, 4], vec![5u64, 6, 7, 8]),
      (vec![9u64, 8, 7, 6], vec![2u64, 3, 4, 5]),
    ];
    let layers: Vec<Layer<E>> = inputs
      .iter()
      .map(|(n, d)| Layer::<E> {
        num: mle(n.clone()),
        den: mle(d.clone()),
      })
      .collect();
    let mut tr_p = <E as Engine>::TE::new(b"gkr-test");
    let (proof_ok, _) = prove::<E>(layers, &mut tr_p).unwrap();

    // Sanity: the honest proof verifies.
    let mut tr_v0 = <E as Engine>::TE::new(b"gkr-test");
    assert!(verifier::verify::<E>(&proof_ok, &mut tr_v0).is_ok());

    // Target a SUMCHECK layer (final_claims index >= 1; index 0 is the base
    // case whose per-instance cross-mult check is not batched). For num_vars=2,
    // index 1 is the input-layer reduction. Apply a complementary offset to the
    // two instances' left numerators.
    let mut proof = proof_ok.clone();
    assert!(proof.final_claims.len() >= 2 && proof.final_claims[1].len() == 2);
    let delta = Fr::from(7);
    proof.final_claims[1][0].left.num += delta;
    proof.final_claims[1][1].left.num -= delta;

    let mut tr_v = <E as Engine>::TE::new(b"gkr-test");
    assert!(
      verifier::verify::<E>(&proof, &mut tr_v).is_err(),
      "tampered per-instance final claims must be rejected"
    );
  }
}
