//! Factorized Projective SumCheck prover (Stage 2): a target polynomial given
//! as a **sum of products of multilinear factors**, with an arbitrary-degree
//! product round kernel.
//!
//! ```text
//!   G(X) = Σ_t coeff_t · ∏_{j ∈ term_t} F_j(X)
//! ```
//! where each `F_j` is a multilinear factor stored as a dense coefficient
//! table (projective monomial basis) and shared by reference across terms. This
//! is the design doc's Phase-3 factorized prover, restricted to **Form B**: all
//! terms have the same multiplicative degree `D` (no structured factors, no `U`
//! homogenization yet — those arrive in Stage 3). Each round message therefore
//! has degree exactly `D`.
//!
//! It is written to satisfy the unmodified [`super::verifier::verify`] (same
//! transcript labels, absorb-then-squeeze order, `compress_projective` round
//! messages), and is differential-tested against
//! [`super::prover::prove_dense_multilinear`].

use crate::{
  spartan::{polys::univariate::UniPoly, sumcheck::SumcheckProof},
  traits::{Engine, TranscriptEngineTrait},
};
use ff::Field;

use super::prover::ProjectiveSumcheckProverOutput;

/// Multiplies the polynomial in `poly` (ascending coefficients, current degree
/// `current_degree`) by the linear factor `a0 + a1·T`, in place.
///
/// `poly[current_degree + 1]` must be present and zero on entry (the slot the
/// new leading term lands in). Design §9.3.
#[inline]
fn multiply_by_linear<F: Field>(poly: &mut [F], current_degree: usize, a0: F, a1: F) {
  for k in (0..=current_degree).rev() {
    let carry = poly[k] * a1;
    poly[k + 1] += carry;
    poly[k] *= a0;
  }
}

/// A factorized virtual polynomial `Σ_t coeff_t · ∏_{j ∈ term_t} F_j`.
///
/// `factors[j]` is factor `F_j`'s dense projective coefficient table (length
/// `2^num_vars`); `terms[t] = (coeff_t, factor_indices)` selects which factors
/// are multiplied in term `t`. All terms must have the same length `D` (Form B).
pub struct VirtualPolynomial<E: Engine> {
  num_vars: usize,
  /// Dense coefficient tables, one per factor; bound in place each round.
  factors: Vec<Vec<E::Scalar>>,
  /// Additive terms: `(coeff_t, indices into factors)`.
  terms: Vec<(E::Scalar, Vec<usize>)>,
}

impl<E: Engine> VirtualPolynomial<E> {
  /// Builds a virtual polynomial from factor tables and product terms.
  ///
  /// # Panics
  /// Panics if any factor table length is not `2^num_vars`, if `terms` is
  /// empty, if any factor index is out of range, or if the terms do not all
  /// share the same multiplicative degree (Form B requirement).
  pub fn new(
    num_vars: usize,
    factors: Vec<Vec<E::Scalar>>,
    terms: Vec<(E::Scalar, Vec<usize>)>,
  ) -> Self {
    let n = 1usize << num_vars;
    for f in &factors {
      assert_eq!(f.len(), n, "factor table length must be 2^num_vars");
    }
    assert!(
      !terms.is_empty(),
      "virtual polynomial needs at least one term"
    );
    let degree = terms[0].1.len();
    assert!(degree >= 1, "projective round degree D must be >= 1");
    for (_, idxs) in &terms {
      assert_eq!(idxs.len(), degree, "Form B: all terms must share degree D");
      for &j in idxs {
        assert!(j < factors.len(), "factor index out of range");
      }
    }
    Self {
      num_vars,
      factors,
      terms,
    }
  }

  /// The common multiplicative degree `D` of every term.
  pub fn degree(&self) -> usize {
    self.terms[0].1.len()
  }

  /// Runs the factorized Projective SumCheck prover.
  ///
  /// Mirrors the verifier's transcript discipline: each round builds the
  /// degree-`D` round polynomial `S_i` by summing, over every surviving suffix
  /// and every term, the product of the factors' linear slices; absorbs `S_i`;
  /// squeezes `r_i`; then binds every factor once at `r_i` in the monomial
  /// basis (`new = a_0 + r_i · a_1`).
  pub fn prove(mut self, transcript: &mut E::TE) -> ProjectiveSumcheckProverOutput<E> {
    let degree = self.degree();
    let num_vars = self.num_vars;

    // C_0 = Σ_b G[[b]]_D = Σ over all corners of Σ_t coeff_t ∏_j F_j[b].
    let initial_claim = self.initial_claim();

    let mut rounds = Vec::with_capacity(num_vars);
    let mut point = Vec::with_capacity(num_vars);

    let mut remaining = num_vars;
    for _ in 0..num_vars {
      let half = 1usize << (remaining - 1);

      // Round polynomial coefficients [a_0, ..., a_D].
      let mut round = vec![E::Scalar::ZERO; degree + 1];
      let mut scratch = vec![E::Scalar::ZERO; degree + 1];
      for suffix in 0..half {
        for (coeff, idxs) in &self.terms {
          // Build coeff · ∏_j (a0_j + a1_j T) into scratch.
          scratch.iter_mut().for_each(|c| *c = E::Scalar::ZERO);
          scratch[0] = *coeff;
          for (cur_deg, &j) in idxs.iter().enumerate() {
            let table = &self.factors[j];
            let a0 = table[suffix];
            let a1 = table[suffix + half];
            multiply_by_linear(&mut scratch, cur_deg, a0, a1);
          }
          for (r, s) in round.iter_mut().zip(scratch.iter()) {
            *r += *s;
          }
        }
      }

      let poly = UniPoly::<E::Scalar>::from_coeffs_no_trim(round)
        .expect("round polynomial has D+1 >= 2 coefficients");

      transcript.absorb(b"projective_sumcheck_round", &poly);
      let r_i = transcript
        .squeeze(b"projective_sumcheck_challenge")
        .expect("transcript squeeze failed");

      rounds.push(poly.compress_projective());
      point.push(r_i);

      // Bind every factor once at r_i: new[s] = table[s] + r_i · table[s+half].
      for table in &mut self.factors {
        let mut next = vec![E::Scalar::ZERO; half];
        for suffix in 0..half {
          next[suffix] = table[suffix] + r_i * table[suffix + half];
        }
        *table = next;
      }
      remaining -= 1;
    }

    // After n rounds each factor table holds a single value F_j(r); the reduced
    // claim is Σ_t coeff_t ∏_j F_j(r).
    let final_claim = self
      .terms
      .iter()
      .map(|(coeff, idxs)| idxs.iter().fold(*coeff, |acc, &j| acc * self.factors[j][0]))
      .sum();

    ProjectiveSumcheckProverOutput {
      proof: SumcheckProof::new(rounds),
      initial_claim,
      point,
      final_claim,
    }
  }

  /// `C_0 = Σ_b Σ_t coeff_t ∏_j F_j[b]` — the projective sum over all corners.
  fn initial_claim(&self) -> E::Scalar {
    let n = 1usize << self.num_vars;
    (0..n)
      .map(|b| {
        self
          .terms
          .iter()
          .map(|(coeff, idxs)| idxs.iter().fold(*coeff, |acc, &j| acc * self.factors[j][b]))
          .sum::<E::Scalar>()
      })
      .sum()
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::{
    provider::PallasEngine,
    spartan::projective_sumcheck::{prove_dense_multilinear, verify},
    traits::Engine,
  };
  use ff::Field;

  type E = PallasEngine;
  type Fr = <E as Engine>::Scalar;

  /// Form B, D=3: G = f·g·h (three multilinear factors). The factorized prover
  /// must (a) satisfy the verifier and (b) produce the same reduction as the
  /// dense oracle fed the corner table G[[b]] = f[b]·g[b]·h[b].
  #[test]
  fn form_b_product_matches_dense_oracle() {
    let num_vars = 3usize;
    let n = 1usize << num_vars;

    let f: Vec<Fr> = (0..n).map(|i| Fr::from((2 * i + 1) as u64)).collect();
    let g: Vec<Fr> = (0..n).map(|i| Fr::from((3 * i + 5) as u64)).collect();
    let h: Vec<Fr> = (0..n).map(|i| Fr::from((7 * i + 2) as u64)).collect();

    // Corner table for the dense oracle: coefficient-wise product (design §4).
    let corners: Vec<Fr> = (0..n).map(|i| f[i] * g[i] * h[i]).collect();

    // Factorized prover.
    let vp = VirtualPolynomial::<E>::new(num_vars, vec![f, g, h], vec![(Fr::ONE, vec![0, 1, 2])]);
    let mut ts_fac = <E as Engine>::TE::new(b"projsc_stage2");
    let out_fac = vp.prove(&mut ts_fac);

    // Dense oracle on the same corner table.
    let mut ts_den = <E as Engine>::TE::new(b"projsc_stage2");
    let out_den = prove_dense_multilinear::<E>(&corners, num_vars, &mut ts_den);

    // Only the initial claim is shared with the dense oracle: the dense oracle
    // runs a degree-1 sumcheck over the collapsed corner table while the
    // factorized prover runs degree-D, so their transcripts and sampled points
    // differ by construction.
    assert_eq!(out_fac.initial_claim, out_den.initial_claim);

    // And the factorized proof verifies.
    let degree_bounds = vec![3usize; num_vars];
    let mut ts_ver = <E as Engine>::TE::new(b"projsc_stage2");
    let reduction = verify::<E>(
      out_fac.initial_claim,
      &degree_bounds,
      &out_fac.proof,
      &mut ts_ver,
    )
    .unwrap();
    assert_eq!(reduction.point, out_fac.point);
    assert_eq!(reduction.final_claim, out_fac.final_claim);
  }

  /// Multiple additive terms of equal degree: G = f·g + h·k (D=2).
  #[test]
  fn form_b_sum_of_products() {
    let num_vars = 2usize;
    let n = 1usize << num_vars;
    let f: Vec<Fr> = (0..n).map(|i| Fr::from((i + 1) as u64)).collect();
    let g: Vec<Fr> = (0..n).map(|i| Fr::from((2 * i + 3) as u64)).collect();
    let h: Vec<Fr> = (0..n).map(|i| Fr::from((5 * i + 1) as u64)).collect();
    let k: Vec<Fr> = (0..n).map(|i| Fr::from((4 * i + 7) as u64)).collect();

    let corners: Vec<Fr> = (0..n).map(|i| f[i] * g[i] + h[i] * k[i]).collect();

    let vp = VirtualPolynomial::<E>::new(
      num_vars,
      vec![f, g, h, k],
      vec![(Fr::ONE, vec![0, 1]), (Fr::ONE, vec![2, 3])],
    );
    let mut ts_fac = <E as Engine>::TE::new(b"projsc_stage2b");
    let out_fac = vp.prove(&mut ts_fac);

    let mut ts_den = <E as Engine>::TE::new(b"projsc_stage2b");
    let out_den = prove_dense_multilinear::<E>(&corners, num_vars, &mut ts_den);

    assert_eq!(out_fac.initial_claim, out_den.initial_claim);

    let degree_bounds = vec![2usize; num_vars];
    let mut ts_ver = <E as Engine>::TE::new(b"projsc_stage2b");
    let reduction = verify::<E>(
      out_fac.initial_claim,
      &degree_bounds,
      &out_fac.proof,
      &mut ts_ver,
    )
    .unwrap();
    assert_eq!(reduction.point, out_fac.point);
    assert_eq!(reduction.final_claim, out_fac.final_claim);
  }

  /// Regression for the `from_coeffs` trailing-zero trimming bug (audit F1):
  /// a round whose leading coefficient `a_D` is zero must still produce a
  /// full-length `D`-coefficient message that the verifier accepts.
  ///
  /// Construct G = f·g (D=2) where `g`'s high half (the X_i^1 coefficients) is
  /// all zero in the first round, forcing `a_2 = Σ f_hi·g_hi = 0` there.
  #[test]
  fn form_b_zero_leading_coeff_round() {
    let num_vars = 2usize;
    let n = 1usize << num_vars;
    // f: arbitrary. g: zero on the high half of X_0 (indices with low bit 1),
    // i.e. g[suffix + half] = 0 for the first round → a_2 = 0 that round.
    let f: Vec<Fr> = (0..n).map(|i| Fr::from((i + 3) as u64)).collect();
    let half0 = n / 2;
    let g: Vec<Fr> = (0..n)
      .map(|i| {
        if i >= half0 {
          Fr::ZERO
        } else {
          Fr::from((2 * i + 1) as u64)
        }
      })
      .collect();

    let vp = VirtualPolynomial::<E>::new(num_vars, vec![f, g], vec![(Fr::ONE, vec![0, 1])]);
    let mut ts_fac = <E as Engine>::TE::new(b"projsc_f1");
    let out_fac = vp.prove(&mut ts_fac);

    // Every round message must carry exactly D = 2 stored coefficients.
    for m in out_fac.proof.compressed_polys() {
      assert_eq!(m.stored_coeffs().len(), 2);
    }

    let degree_bounds = vec![2usize; num_vars];
    let mut ts_ver = <E as Engine>::TE::new(b"projsc_f1");
    let reduction = verify::<E>(
      out_fac.initial_claim,
      &degree_bounds,
      &out_fac.proof,
      &mut ts_ver,
    )
    .unwrap();
    assert_eq!(reduction.point, out_fac.point);
    assert_eq!(reduction.final_claim, out_fac.final_claim);
  }
}
