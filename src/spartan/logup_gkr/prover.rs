//! Prover for the Logup-GKR fractional-sum argument.
//!
//! **Stage 1: placeholder only.** The signature is frozen against `proof.rs`
//! and `verifier.rs`; the body is filled in the implementation stage, reusing
//! Nova's batched-cubic sumcheck (`SumcheckProof::prove_batched_cubic`, which
//! already includes the Gruen eq-correction) for each layer, and absorbing
//! every leaf/root/claim into the transcript (the soundness red line).

use crate::errors::NovaError;
use crate::spartan::logup_gkr::layer::Layer;
use crate::spartan::logup_gkr::proof::{LogupGkrOpeningClaim, LogupGkrProof};
use crate::traits::Engine;

/// Proves the fractional-sum identity `Σ p/q = 0` for all logup instances in a
/// single batched tree (`_inputs` holds one input `Layer` per instance, e.g.
/// `[row, col]`), returning the proof and the shared opening claim to be batched
/// with the inner sumcheck.
///
/// # Stage 1
/// Not yet implemented — frozen signature with a placeholder body.
pub fn prove<E: Engine>(
  _inputs: Vec<Layer<E>>,
  _transcript: &mut E::TE,
) -> Result<(LogupGkrProof<E>, LogupGkrOpeningClaim<E>), NovaError> {
  unimplemented!("logup_gkr::prover::prove is a stage-1 placeholder")
}
