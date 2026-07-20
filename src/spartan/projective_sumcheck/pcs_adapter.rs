//! PCS boundary adapter (Stage 4): the zeta / subset-sum transform that lets a
//! standard **evaluation-basis** MLE commitment open a **coefficient-basis**
//! table at the projective sumcheck's reduced point.
//!
//! The projective sumcheck opens a witness as a coefficient-MLE:
//! ```text
//!   coeffMLE(vec, r) = Σ_b vec[b] · ∏_{i ∈ b} r_i.
//! ```
//! A standard MLE-PCS opens the evaluation-MLE:
//! ```text
//!   evalMLE(u, r) = Σ_b u[b] · eq̃(b, r).
//! ```
//! They coincide under the **zeta (subset-sum) transform** `u = ζ(vec)`:
//! ```text
//!   u[b] = Σ_{b' ⊆ b} vec[b'],   ⇒   evalMLE(ζ(vec), r) = coeffMLE(vec, r) ∀r.
//! ```
//! So a projective prover commits `ζ(vec)` for each witness and the existing
//! PCS opens it unchanged. The transform is `O(n · 2^n)` **at commit time
//! only** — outside the sumcheck hot path — so Route B stays transform-free in
//! the sumcheck itself (see the Stage-1 relation notes, §5).

use ff::Field;

/// In-place zeta (subset-sum) transform over the Boolean lattice:
/// `u[b] = Σ_{b' ⊆ b} vec[b']`, computed by the standard `n`-pass
/// dimension-wise prefix (the SOS/zeta transform dual to Möbius).
///
/// `vec.len()` must be a power of two (`2^n`). Runs in `O(n · 2^n)`.
pub fn zeta_transform<F: Field>(vec: &mut [F]) {
  let len = vec.len();
  debug_assert!(len.is_power_of_two(), "table length must be 2^n");
  let mut step = 1;
  while step < len {
    let mut base = 0;
    while base < len {
      for k in base..base + step {
        // Elements whose bit for this dimension is 1 accumulate the paired
        // element whose bit is 0.
        let lo = vec[k];
        vec[k + step] += lo;
      }
      base += step << 1;
    }
    step <<= 1;
  }
}

/// Returns `ζ(vec)` as a fresh vector, leaving the input untouched.
pub fn zeta<F: Field>(vec: &[F]) -> Vec<F> {
  let mut out = vec.to_vec();
  zeta_transform(&mut out);
  out
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::{provider::PallasEngine, traits::Engine};

  type E = PallasEngine;
  type Fr = <E as Engine>::Scalar;

  /// coeffMLE(vec, r) = Σ_b vec[b] ∏_{i∈b} r_i, with the low bit = X_0
  /// (matching the projective provers' binding order after reversal).
  fn coeff_mle(vec: &[Fr], r: &[Fr]) -> Fr {
    let n = vec.len();
    (0..n)
      .map(|b| {
        let mut acc = vec[b];
        for (i, ri) in r.iter().enumerate() {
          if (b >> i) & 1 == 1 {
            acc *= *ri;
          }
        }
        acc
      })
      .sum()
  }

  /// The core boundary identity: evalMLE(ζ(vec), r) == coeffMLE(vec, r) for all
  /// r. evalMLE uses the standard eq-weighted MLE (the PCS's opening).
  #[test]
  fn zeta_makes_eval_mle_equal_coeff_mle() {
    for num_vars in 1..=5usize {
      let n = 1usize << num_vars;
      let vec: Vec<Fr> = (0..n).map(|i| Fr::from((3 * i + 7) as u64)).collect();
      let r: Vec<Fr> = (0..num_vars)
        .map(|i| Fr::from((2 * i + 5) as u64))
        .collect();

      let coeff = coeff_mle(&vec, &r);

      // evalMLE(ζ(vec), r) = Σ_b u[b]·eq̃(b, r), with eq̃ in the SAME low-bit =
      // X_0 convention as coeff_mle (build it directly to avoid depending on the
      // library's internal bit order).
      let u = zeta(&vec);
      let eq_weight = |b: usize| -> Fr {
        let mut w = Fr::ONE;
        for (i, ri) in r.iter().enumerate() {
          w *= if (b >> i) & 1 == 1 {
            *ri
          } else {
            Fr::ONE - *ri
          };
        }
        w
      };
      let eval_mle: Fr = (0..n).map(|b| u[b] * eq_weight(b)).sum();

      assert_eq!(eval_mle, coeff, "num_vars={num_vars}");
    }
  }

  /// zeta is the identity's dual: ζ of a single monomial spreads to its
  /// supersets, and ζ(e_0) (constant term) fills every corner with 1.
  #[test]
  fn zeta_of_constant_term_is_all_ones() {
    let mut v = vec![Fr::ZERO; 8];
    v[0] = Fr::ONE;
    assert_eq!(zeta(&v), vec![Fr::ONE; 8]);
  }
}
