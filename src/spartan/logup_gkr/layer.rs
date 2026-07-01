//! Layer type for the Logup-GKR fractional-sum tree.
//!
//! One [`Layer`] is a single level of a tree: two multilinear polynomials
//! (numerator, denominator). The same type serves the input layer, every
//! internal layer, and the output layer — they differ only in height (the
//! coefficient length halves each level, `N → N/2 → … → 1`). This mirrors
//! hyperplonk's single `Layer` type (`fractional_gkr/layer/mod.rs`), which has
//! no separate input/output variants.
//!
//! Each logup instance (`row` / `col`) starts from one input `Layer` whose
//! leaves are the projective cells `(num, den) = (ts, T+r)` on the table side
//! and `(1, W+r)` on the access side. Stage 1 fixes the type only; the tree
//! builder and its fold are written in the implementation stage.

use crate::spartan::polys::multilinear::MultilinearPolynomial;
use crate::traits::Engine;

/// One level of a fractional-sum tree: parallel numerator and denominator
/// multilinear polynomials over `{0,1}^{log len}`.
///
/// Storage is **two** MLEs. The "left/right" children that the fraction gate
/// consumes are *not* stored separately — they are the even/odd positions of
/// these two MLEs (`num[2i]`, `num[2i+1]`, `den[2i]`, `den[2i+1]`), read during
/// the fold. So a layer holds 2 polynomials; a per-layer *claim* exposes 4
/// values (`nL, nR, dL, dR`; see `proof::LayerFinalClaim`). Do not confuse the
/// two.
///
/// Invariant: `num.len() == den.len()` and both are a power of two; the padding
/// cell is `(0, 1)`, the additive identity of the fraction monoid.
///
/// # Numerator/denominator convention (read before constructing)
/// `num` is the **multiplicity** side (`ts` on the table side, `1` on the
/// access side); `den` is the **fingerprint** side (`T+r` / `W+r`, never
/// inverted). This is the SAME assignment as hyperplonk's
/// `new_input_layer(data, multiplicity)` — but note hp's constructor takes
/// `(data, multiplicity)` = `(den, num)`, the OPPOSITE positional order. To
/// make the swap impossible, this type has **no positional constructor**:
/// build it with field-init syntax so `num`/`den` are named explicitly, e.g.
/// `Layer { num: ts, den: t_plus_r }`.
pub struct Layer<E: Engine> {
  /// Numerator = multiplicity (`ts` on the table side, `1` on the access side).
  pub num: MultilinearPolynomial<E::Scalar>,
  /// Denominator = fingerprint (`T+r` / `W+r`), never inverted.
  pub den: MultilinearPolynomial<E::Scalar>,
}
