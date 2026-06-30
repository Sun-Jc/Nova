//! Verifier for the Logup-GKR fractional-sum argument.
//!
//! Replays the fold-down of the **single batched tree** (all instances share
//! one GKR depth and one OOD point — audit D1) and reduces it to an eval
//! request on the input layer ([`LogupGkrOpeningClaim`]). It opens no
//! commitments and does **not** perform the `0/den` zero-sum check: per the
//! primary reference (hp `fractional_gkr/verifier.rs`), that check belongs to
//! the host's reconcile step. The GKR verifier only checks internal consistency
//! (output claims vs `final_claims[0]` gates, and each layer's batched sumcheck
//! against the merged gate claim) and returns the shared OOD point plus the
//! per-instance reduced fractions.

use crate::errors::NovaError;
use crate::spartan::logup_gkr::fraction::Fraction;
use crate::spartan::logup_gkr::proof::{LogupGkrOpeningClaim, LogupGkrProof};
use crate::traits::Engine;

/// Degree bound of each layer's batched sumcheck round polynomial under Nova's
/// `prove_batched_cubic` (`eq · (p0·q1 + p1·q0 + λ·q0·q1)`, eq factored via
/// Gruen). To be confirmed against Nova's emission when the prover lands; hp's
/// own `max_degree = 2` does **not** transfer (different eq handling).
const LAYER_SC_DEGREE: usize = 3;

/// Verifies the batched memory-check proof (all logup instances in one tree) and
/// returns the shared opening claim. The host then rerandomizes the openings
/// into the inner sumcheck and runs the zero-sum reconcile.
///
/// # Stage 1
/// Fold-down skeleton with the correct single-tree structure and check
/// placement; the exact sumcheck/transcript wiring is finalized with the prover.
///
/// Structure (hp `verifier.rs:21-99`): observe `initial_claims`; loop layers —
/// layer 0 checks `initial_claims == final_claims[0].compute_gates()` (no
/// sumcheck); later layers sample a fresh `λ`, verify the batched layer sumcheck
/// against `claims.merged(λ)` and assert it equals `final_claims.compute_gates()
/// .merged(λ)`; each step observes `final_claims`, samples fold `r`, folds every
/// instance, and grows the single OOD point. The leftover claims are the
/// per-instance input-layer fractions.
pub fn verify<E: Engine>(
  proof: &LogupGkrProof<E>,
  _transcript: &mut E::TE,
) -> Result<LogupGkrOpeningClaim<E>, NovaError> {
  // Shape sanity the prover must satisfy (cheap, real):
  // one batched sumcheck per non-base layer, and every layer carries one split
  // claim per instance.
  if proof.final_claims.is_empty() || proof.sumchecks.len() + 1 != proof.final_claims.len() {
    return Err(NovaError::InvalidSumcheckProof);
  }
  let num_instances = proof.initial_claims.len();
  if num_instances == 0
    || proof
      .final_claims
      .iter()
      .any(|layer| layer.len() != num_instances)
  {
    return Err(NovaError::InvalidSumcheckProof);
  }

  let _ = (LAYER_SC_DEGREE, Fraction::<E::Scalar>::zero);
  unimplemented!("logup_gkr::verifier::verify is a stage-1 placeholder")
}
