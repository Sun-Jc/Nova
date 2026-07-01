//! Logup-GKR fractional-sum memory-check argument.
//!
//! Replaces the inverse-logup memory-check in `ppsnark.rs` (the
//! `MemorySumcheckInstance` 6-route sumcheck + 4 inverse-polynomial
//! commitments) with a fractional-sum GKR tree per logup instance (`row`,
//! `col`). Projective fractions keep the circuit inversion-free, so the four
//! inverse commitments — measured at ~4.17s of a 15.7s prove at 1M
//! constraints (`perf-ppsnark-baseline.md`) — disappear.
//!
//! ## Module map
//! - [`fraction`]: projective fraction + 2-to-1 gate (pure).
//! - [`layer`]: input-layer types (leaves of each tree).
//! - [`proof`]: frozen proof/claim interface.
//! - [`prover`]: stage-1 placeholder; reuses `prove_batched_cubic` later.
//! - [`verifier`]: fold-down + root check, emits the shared opening claim.
//!
//! ## Boundary with ppSNARK (host reconcile contract)
//! The argument owns no commitment scheme. Its verifier returns a
//! [`proof::LogupGkrOpeningClaim`]: a single shared `eval_point` plus the
//! per-instance input-layer fractions `openings` (order `[row, col]`). The
//! **host** then:
//! 1. rerandomizes `L_row`/`L_col` at `eval_point` into a sumcheck batched with
//!    the inner sumcheck, and opens them (with the other columns) via HyperKZG
//!    at the shared point (see `rerandomize-batch-explained.md`);
//! 2. recomputes each instance's fraction from its opened `L`/`addr`/`ts`
//!    (`den = L·γ + addr + r`, `num = ts`) and checks it equals the matching
//!    entry of `openings` — the analogue of hp `reconcile_openings` /
//!    `eval_at_openings`;
//! 3. runs the `0/den` zero-sum balance check.
//! Steps 2–3 are the host's job, never the GKR verifier's.
//!
//! References: hyperplonk-logup-gkr (primary) and lambdaworks `gkr-logup`
//! (secondary); divergences are flagged at the use site.

pub mod fraction;
pub mod layer;
pub mod proof;
pub mod prover;
pub mod verifier;

pub use proof::{LogupGkrOpeningClaim, LogupGkrProof};
