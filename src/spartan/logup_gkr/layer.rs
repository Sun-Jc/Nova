//! Input-layer construction for the Logup-GKR fractional-sum tree.
//!
//! Each logup instance (`row` / `col`) becomes one tree whose leaves are the
//! projective cells `(num, den) = (ts, T+r)` on the table side and
//! `(1, W+r)` on the access side. Stage 1 only fixes the type; the builder is
//! written in the implementation stage.

use crate::spartan::polys::multilinear::MultilinearPolynomial;
use crate::traits::Engine;

/// The input (leaf) layer of one fractional-sum tree: parallel numerator and
/// denominator multilinear polynomials over `{0,1}^{log N}`.
///
/// Invariant: `num.len() == den.len()` and both are a power of two; the padding
/// cell is `(0, 1)` so it is the additive identity of the fraction monoid.
///
/// # Numerator/denominator convention (read before constructing)
/// `num` is the **multiplicity** side (`ts` on the table side, `1` on the
/// access side); `den` is the **fingerprint** side (`T+r` / `W+r`, never
/// inverted). This is the SAME assignment as hyperplonk's
/// `new_input_layer(data, multiplicity)` — but note hp's constructor takes
/// `(data, multiplicity)` = `(den, num)`, the OPPOSITE positional order. To
/// make the swap impossible, this type has **no positional constructor**:
/// build it with field-init syntax so `num`/`den` are named explicitly, e.g.
/// `InputLayer { num: ts, den: t_plus_r }`.
pub struct InputLayer<E: Engine> {
  /// Numerator = multiplicity (`ts` on the table side, `1` on the access side).
  pub num: MultilinearPolynomial<E::Scalar>,
  /// Denominator = fingerprint (`T+r` / `W+r`), never inverted.
  pub den: MultilinearPolynomial<E::Scalar>,
}
