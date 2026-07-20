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
use rayon::prelude::*;

use super::prover::ProjectiveSumcheckProverOutput;

/// Multiplies the polynomial in `poly` (ascending coefficients, current degree
/// `current_degree`) by the linear factor `a0 + a1·T`, in place.
///
/// `poly[current_degree + 1]` must be present and zero on entry (the slot the
/// new leading term lands in). Design §9.3.
#[inline]
pub(crate) fn multiply_by_linear<F: Field>(poly: &mut [F], current_degree: usize, a0: F, a1: F) {
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

  /// The number of variables `n`.
  pub fn num_vars(&self) -> usize {
    self.num_vars
  }

  /// Raises the common degree to `target` (`≥ current D`) by padding every term
  /// with the projective all-ones factor `U` (design §5). A no-op if already at
  /// `target`. Needed to bring instances to a shared degree before batching,
  /// where zero-padding would break the projective endpoint.
  ///
  /// # Panics
  /// Panics if `target < self.degree()`.
  pub fn lift_to_degree(&mut self, target: usize) {
    let cur = self.degree();
    assert!(target >= cur, "cannot lift to a lower degree");
    if target == cur {
      return;
    }
    // Reuse an existing all-ones factor if present, else append one.
    let n = 1usize << self.num_vars;
    let ones = vec![E::Scalar::ONE; n];
    let u_index = match self.factors.iter().position(|f| f == &ones) {
      Some(i) => i,
      None => {
        self.factors.push(ones);
        self.factors.len() - 1
      }
    };
    for (_, idxs) in &mut self.terms {
      for _ in 0..(target - cur) {
        idxs.push(u_index);
      }
    }
  }

  /// Builds a virtual polynomial from **mixed-degree** terms, homogenizing each
  /// to the common degree `D = max_t m_t` with the projective all-ones factor
  /// `U(X) = ∏_i (1 + X_i)` (design §5).
  ///
  /// `U` has corner coefficient `1` at every projective corner, so padding a
  /// term with `U^{D - m_t}` preserves every corner value while raising its
  /// declared degree to `D`. This is what lets a coefficient-basis pointwise
  /// relation mix factors of different multiplicative degree (e.g. `Az·Bz`
  /// against `u·Cz`) under one declared degree — see the Stage-1 relation
  /// ledger.
  ///
  /// Correctness-first: `U` is materialized as one dense `2^num_vars` factor and
  /// appended to `factors`; a structured (`O(n)`) `U` is a later perf item.
  ///
  /// # Panics
  /// Same guards as [`VirtualPolynomial::new`], except terms may have differing
  /// lengths (each `≥ 1`); the max length becomes `D`.
  pub fn new_homogenized(
    num_vars: usize,
    mut factors: Vec<Vec<E::Scalar>>,
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
    let degree = terms.iter().map(|(_, idxs)| idxs.len()).max().unwrap();
    assert!(degree >= 1, "projective round degree D must be >= 1");
    for (_, idxs) in &terms {
      assert!(!idxs.is_empty(), "each term needs at least one factor");
      for &j in idxs {
        assert!(j < factors.len(), "factor index out of range");
      }
    }

    // Append one dense all-ones factor U (every corner coefficient = 1).
    let u_index = factors.len();
    factors.push(vec![E::Scalar::ONE; n]);

    // Homogenize: pad each term with U^{D - m_t}.
    let terms = terms
      .into_iter()
      .map(|(coeff, mut idxs)| {
        let pad = degree - idxs.len();
        idxs.extend(std::iter::repeat(u_index).take(pad));
        (coeff, idxs)
      })
      .collect();

    Self {
      num_vars,
      factors,
      terms,
    }
  }

  /// Runs the factorized Projective SumCheck prover.
  ///
  /// Mirrors the verifier's transcript discipline: each round builds the
  /// degree-`D` round polynomial `S_i` by summing, over every surviving suffix
  /// and every term, the product of the factors' linear slices; absorbs `S_i`;
  /// squeezes `r_i`; then binds every factor once at `r_i` in the monomial
  /// basis (`new = a_0 + r_i · a_1`).
  pub fn prove(mut self, transcript: &mut E::TE) -> ProjectiveSumcheckProverOutput<E> {
    let num_vars = self.num_vars;

    // C_0 = Σ_b G[[b]]_D = Σ over all corners of Σ_t coeff_t ∏_j F_j[b].
    let initial_claim = self.initial_claim();

    let mut rounds = Vec::with_capacity(num_vars);
    let mut point = Vec::with_capacity(num_vars);

    let mut remaining = num_vars;
    for _ in 0..num_vars {
      // Full round polynomial coefficients [a_0, ..., a_D].
      let round = self.round_full_coeffs(remaining);

      let poly = UniPoly::<E::Scalar>::from_coeffs_no_trim(round)
        .expect("round polynomial has D+1 >= 2 coefficients");

      transcript.absorb(b"projective_sumcheck_round", &poly);
      let r_i = transcript
        .squeeze(b"projective_sumcheck_challenge")
        .expect("transcript squeeze failed");

      rounds.push(poly.compress_projective());
      point.push(r_i);

      self.bind(remaining, r_i);
      remaining -= 1;
    }

    // After n rounds each factor table holds a single value F_j(r); the reduced
    // claim is Σ_t coeff_t ∏_j F_j(r).
    let final_claim = self.reduced_claim();

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

  // --- Stepwise primitives (used by both `prove` and the batched prover) ---

  /// The public projective sum `C_0`.
  pub(crate) fn claim0(&self) -> E::Scalar {
    self.initial_claim()
  }

  /// Builds this round's full coefficient vector `[a_0, ..., a_D]` over the
  /// current (unbound) tables, with `remaining` variables left. Does not touch
  /// the transcript or bind anything. Parallelized over suffixes with a
  /// per-chunk scratch buffer (allocation-free in the inner loop).
  pub(crate) fn round_full_coeffs(&self, remaining: usize) -> Vec<E::Scalar> {
    let degree = self.degree();
    let half = 1usize << (remaining - 1);

    let fold_range = |lo: usize, hi: usize| -> Vec<E::Scalar> {
      let mut round = vec![E::Scalar::ZERO; degree + 1];
      let mut scratch = vec![E::Scalar::ZERO; degree + 1];
      for suffix in lo..hi {
        for (coeff, idxs) in &self.terms {
          scratch.iter_mut().for_each(|c| *c = E::Scalar::ZERO);
          scratch[0] = *coeff;
          for (cur_deg, &j) in idxs.iter().enumerate() {
            let table = &self.factors[j];
            multiply_by_linear(&mut scratch, cur_deg, table[suffix], table[suffix + half]);
          }
          for (r, s) in round.iter_mut().zip(scratch.iter()) {
            *r += *s;
          }
        }
      }
      round
    };

    // Parallel map-reduce over suffix chunks; small sizes stay serial.
    const PAR_THRESHOLD: usize = 1 << 10;
    if half < PAR_THRESHOLD {
      fold_range(0, half)
    } else {
      let chunk = (half / rayon::current_num_threads().max(1)).max(1);
      (0..half)
        .into_par_iter()
        .step_by(chunk)
        .map(|lo| fold_range(lo, (lo + chunk).min(half)))
        .reduce(
          || vec![E::Scalar::ZERO; degree + 1],
          |mut a, b| {
            for (x, y) in a.iter_mut().zip(b.iter()) {
              *x += *y;
            }
            a
          },
        )
    }
  }

  /// Binds every factor once at `r_i` in the monomial basis
  /// (`new = a_0 + r_i · a_1`), with `remaining` variables left before binding.
  /// In place: the lower half of each table is overwritten and truncated.
  pub(crate) fn bind(&mut self, remaining: usize, r_i: E::Scalar) {
    let half = 1usize << (remaining - 1);
    for table in &mut self.factors {
      let (lo, hi) = table.split_at_mut(half);
      lo.par_iter_mut()
        .zip(hi.par_iter())
        .for_each(|(a, b)| *a += r_i * *b);
      table.truncate(half);
    }
  }

  /// After all variables are bound, the reduced value `Σ_t coeff_t ∏_j F_j(r)`.
  pub(crate) fn reduced_claim(&self) -> E::Scalar {
    self
      .terms
      .iter()
      .map(|(coeff, idxs)| idxs.iter().fold(*coeff, |acc, &j| acc * self.factors[j][0]))
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

  /// Form A end-to-end: the design §13 coefficient-wise ZeroCheck
  /// `G = Eq^∞_ρ · (f·g − h·U)`, D=3, built with mixed-degree terms and
  /// `new_homogenized` auto-padding with U. `Eq^∞_ρ` is supplied as a factor
  /// table. The initial claim must be zero (the residual vanishes on every
  /// corner), the proof must verify, and `verify_final_claim` must reconstruct
  /// `G(r)` from the factor openings.
  #[test]
  fn form_a_zerocheck_eq_times_residual() {
    use crate::spartan::polys::eq_projective::EqPolynomialProjective;
    use crate::spartan::projective_sumcheck::ProjectiveSumcheckReduction;

    let num_vars = 2usize;
    let n = 1usize << num_vars;

    // Multilinear factors as coefficient tables (design §13: f=X+2Y etc.).
    // Corner index b: low bit = X_0. f[b], g[b], h[b] are the coefficients.
    let f = vec![Fr::ZERO, Fr::from(1), Fr::from(2), Fr::ZERO]; // X + 2Y
    let g = vec![Fr::ZERO, Fr::from(3), Fr::from(4), Fr::ZERO]; // 3X + 4Y
    let h = vec![Fr::ZERO, Fr::from(3), Fr::from(8), Fr::ZERO]; // 3X + 8Y  (= f∘g)

    let rho = vec![Fr::from(2), Fr::from(3)];
    let eq_table = EqPolynomialProjective::<Fr>::new(rho.clone()).evals();

    // Factors: [Eq, f, g, h]. Terms: +Eq·f·g and −Eq·h (mixed degree 3 vs 2 →
    // the second is homogenized by U to degree 3).
    let vp = VirtualPolynomial::<E>::new_homogenized(
      num_vars,
      vec![eq_table.clone(), f.clone(), g.clone(), h.clone()],
      vec![(Fr::ONE, vec![0, 1, 2]), (-Fr::ONE, vec![0, 3])],
    );
    // Residual f∘g − h = 0 at every corner ⇒ projective sum = 0.
    assert_eq!(vp.degree(), 3);

    let mut ts_fac = <E as Engine>::TE::new(b"projsc_forma");
    let out = vp.prove(&mut ts_fac);
    assert_eq!(out.initial_claim, Fr::ZERO);

    // Verify the reduction.
    let degree_bounds = vec![3usize; num_vars];
    let mut ts_ver = <E as Engine>::TE::new(b"projsc_forma");
    let reduction = verify::<E>(Fr::ZERO, &degree_bounds, &out.proof, &mut ts_ver).unwrap();
    assert_eq!(reduction.point, out.point);
    assert_eq!(reduction.final_claim, out.final_claim);

    // Final-oracle check: reconstruct G(r) from factor openings at r.
    // evals order [Eq, f, g, h]; U built internally by verify_final_claim.
    let r = &reduction.point;
    let eval = |t: &[Fr]| -> Fr {
      // coeff-basis MLE evaluate: Σ_b t[b] ∏_{i∈b} r_i. The prover binds the
      // high bit each round, so challenge r[k] pairs with bit (num_vars-1-k).
      (0..n)
        .map(|b| {
          let mut acc = t[b];
          for (k, rk) in r.iter().enumerate() {
            let bit = num_vars - 1 - k;
            if (b >> bit) & 1 == 1 {
              acc *= *rk;
            }
          }
          acc
        })
        .sum()
    };
    let evals = vec![eval(&eq_table), eval(&f), eval(&g), eval(&h)];
    let terms = vec![(Fr::ONE, vec![0usize, 1, 2]), (-Fr::ONE, vec![0usize, 3])];
    let red2 = ProjectiveSumcheckReduction::<E> {
      point: reduction.point.clone(),
      final_claim: reduction.final_claim,
    };
    assert!(red2.verify_final_claim(&evals, &terms));
  }
}
