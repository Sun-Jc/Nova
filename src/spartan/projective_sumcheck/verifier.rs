//! The Projective SumCheck **verifier**: the round-by-round reduction.

use crate::{
  errors::NovaError,
  spartan::{polys::univariate::CompressedUniPoly, sumcheck::SumcheckProof},
  traits::{Engine, TranscriptEngineTrait},
};
use ff::Field;

/// Result of the SumCheck reduction: the sampled random point and the final
/// reduced claim `C_n`. The caller is responsible for checking
/// `C_n == G(r_0, ..., r_{n-1})` via its own oracle / PCS opening.
#[derive(Clone, Debug)]
pub struct ProjectiveSumcheckReduction<E: Engine> {
  /// The sampled evaluation point `r = (r_0, ..., r_{n-1})`.
  pub point: Vec<E::Scalar>,
  /// The final reduced claim `C_n = S_{n-1}(r_{n-1})`.
  pub final_claim: E::Scalar,
}

impl<E: Engine> ProjectiveSumcheckReduction<E> {
  /// Host-side final-oracle check: recompute the virtual polynomial `G(r)` from
  /// authenticated factor openings and compare it against the reduced claim
  /// `C_n`. This is the "second-order" kernel check that binds the SumCheck
  /// reduction to the committed polynomials — the analogue of `snark.rs`'s
  /// `claim_outer_final == expected` and the logup-GKR host's reconcile step.
  ///
  /// The virtual polynomial is a sum of scaled products of multilinear factors,
  /// normalized to a common projective degree `D` with the all-ones factor
  /// `U(X) = ∏_i (1 + X_i)` (see design §5):
  /// ```text
  ///   G(X) = Σ_t coeff_t · U(X)^{D - m_t} · ∏_{j ∈ term_t} F_j(X),
  ///   D = max_t m_t,   m_t = |term_t|.
  /// ```
  /// Evaluated at the sampled point `r = self.point`:
  /// ```text
  ///   G(r) = Σ_t coeff_t · U(r)^{D - m_t} · ∏_{j ∈ term_t} evals[j],
  ///   U(r) = ∏_i (1 + r_i).
  /// ```
  ///
  /// Arguments:
  /// - `evals`: the authenticated factor openings `F_j(r)`, one per factor, in
  ///   a fixed order (indexed by the product terms below). These come from PCS
  ///   openings at the point `r`; the caller must have verified them.
  /// - `terms`: the virtual polynomial as `(coeff_t, factor_indices)` pairs.
  ///   Each `factor_indices` lists which `evals` are multiplied in that term;
  ///   its length is the term's multiplicative degree `m_t`. `U(r)` is built
  ///   internally for the `D - m_t` homogenizing power — callers never pass `U`
  ///   as a factor.
  ///
  /// Returns `true` iff `G(r) == C_n`. Returns `false` (rather than panicking)
  /// if any factor index is out of range.
  ///
  /// Structured factors such as the projective equality polynomial
  /// `Eq^∞_ρ(r)` are *not* built here: fold them into `evals` as an ordinary
  /// opening (the caller evaluates them locally, `O(n)`, since they are public
  /// functions of `r`).
  pub fn verify_final_claim(&self, evals: &[E::Scalar], terms: &[(E::Scalar, Vec<usize>)]) -> bool {
    // U(r) = ∏_i (1 + r_i), evaluated locally from the sampled point.
    let u_at_r = self
      .point
      .iter()
      .fold(E::Scalar::ONE, |acc, r| acc * (E::Scalar::ONE + *r));

    // Common projective degree D = max multiplicative degree over all terms.
    let degree = terms.iter().map(|(_, idxs)| idxs.len()).max().unwrap_or(0);

    let mut g_at_r = E::Scalar::ZERO;
    for (coeff, idxs) in terms {
      // Product of the term's factor openings, scaled by its coefficient.
      let mut term_val = *coeff;
      for &j in idxs {
        let Some(&f_j) = evals.get(j) else {
          return false;
        };
        term_val *= f_j;
      }
      // Homogenize to degree D with U(r)^{D - m_t}.
      for _ in 0..(degree - idxs.len()) {
        term_val *= u_at_r;
      }
      g_at_r += term_val;
    }

    g_at_r == self.final_claim
  }
}

/// Verify a Projective SumCheck proof, running the round-by-round reduction.
///
/// Arguments:
/// - `initial_claim`: `C_0 = v`, the claimed projective sum.
/// - `degree_bounds`: per-round declared degree bound `D_i` (each `>= 1`);
///   `degree_bounds.len()` is the number of variables `n`.
/// - `proof`: the round messages, each a [`CompressedUniPoly`] compressed with
///   `UniPoly::compress_projective` (constant term omitted).
/// - `transcript`: Fiat-Shamir transcript; challenges are squeezed *after*
///   absorbing each round message.
///
/// On success returns the sampled point `r` and the reduced claim `C_n`. The
/// caller must still enforce `C_n == G(r)`.
pub fn verify<E: Engine>(
  initial_claim: E::Scalar,
  degree_bounds: &[usize],
  proof: &SumcheckProof<E>,
  transcript: &mut E::TE,
) -> Result<ProjectiveSumcheckReduction<E>, NovaError> {
  let num_vars = degree_bounds.len();
  let rounds: &[CompressedUniPoly<E::Scalar>] = proof.compressed_polys();
  if rounds.len() != num_vars {
    return Err(NovaError::InvalidSumcheckProof);
  }

  let mut claim = initial_claim;
  let mut point = Vec::with_capacity(num_vars);

  for (round_index, message) in rounds.iter().enumerate() {
    let degree = degree_bounds[round_index];
    // Projective compression requires D >= 1; the message stores exactly the
    // D non-constant coefficients [a_1, ..., a_D].
    if degree == 0 || message.stored_coeffs().len() != degree {
      return Err(NovaError::InvalidSumcheckProof);
    }

    // Reconstruct the full round polynomial via the projective identity
    // a_0 = C_i - a_D, reusing UniPoly as the dense container.
    let poly = message.decompress_projective(&claim);

    // Absorb the round message (all D+1 coeffs, incl. the linear term), then
    // squeeze the challenge (never before). See `absorb_round`.
    super::absorb_round::<E>(transcript, &poly);
    let r_i = transcript.squeeze(b"projective_sumcheck_challenge")?;

    // Reduce: C_{i+1} = S_i(r_i), reusing UniPoly's Horner evaluation.
    claim = poly.evaluate(&r_i);
    point.push(r_i);
  }

  Ok(ProjectiveSumcheckReduction {
    point,
    final_claim: claim,
  })
}
