//! Eq-factored projective sumcheck prover (Stage 5 optimization).
//!
//! ppSNARK's relations have the shape `Eq^∞_τ(X) · R(X)`, where `R` is a
//! sum-of-products of non-eq multilinear factors. The naive
//! [`VirtualPolynomial`](super::virtual_poly::VirtualPolynomial) binds `Eq^∞`
//! as a dense `2^n` factor and rebinds it every round. This module factors
//! `Eq^∞` out (Gruen-style, eprint 2024/108): eq is carried analytically as a
//! running scalar `eq_left` times a per-suffix weight, and contributes only its
//! linear factor `(1−τ_i)+τ_i·T` to each round polynomial — so the dense factor
//! set holds only `R`'s factors (degree `D−1` per term), and eq is never bound.
//!
//! Round `i` (current variable = high bit, `remaining` variables left):
//! ```text
//!   S_i(T) = eq_left · [(1−τ_cur)+τ_cur·T] · Σ_s eq_suffix[s] · R_slice(s,T),
//! ```
//! where `eq_left = ∏_{k bound} ((1−τ_k)+τ_k r_k)`, `eq_suffix` is the projective
//! eq corner table over the remaining suffix taus, and `R_slice(s,T)` is the
//! product of `R`'s factor linear slices (degree `D−1`). The eq linear factor
//! raises the message to degree `D`.
//!
//! Correctness is pinned by a differential test against the naive path (eq as
//! an explicit dense factor). The projective sumcheck verifier is untouched.

use crate::{
  spartan::{
    polys::univariate::UniPoly, projective_sumcheck::prover::ProjectiveSumcheckProverOutput,
    sumcheck::SumcheckProof,
  },
  traits::{Engine, TranscriptEngineTrait},
};
use ff::Field;
use rayon::prelude::*;

/// Natural-order projective eq corner table: `t[b] = ∏_i (bit_i(b) ? τ_i :
/// 1−τ_i)`, i.e. bit `i` corresponds to `τ_i` (low bit = `τ_0`), matching the
/// factored decomposition below. (`EqPolynomialProjective::evals` uses the
/// reverse tau order, so this dedicated builder keeps the current-var/suffix
/// split trivial: current var = `τ_cur`, suffix over `τ_0..τ_cur`.)
fn eq_corner_natural<F: Field>(taus: &[F]) -> Vec<F> {
  let n = taus.len();
  let mut evals = vec![F::ZERO; 1usize << n];
  evals[0] = F::ONE;
  for (i, &tau) in taus.iter().enumerate() {
    let block = 1usize << i;
    for b in 0..block {
      let lo = evals[b];
      evals[b] = lo * (F::ONE - tau);
      evals[b + block] = lo * tau;
    }
  }
  evals
}

/// A relation `Eq^∞_τ · R`, where `R = Σ_t coeff_t ∏_{j∈term_t} F_j` over the
/// non-eq factors, all terms sharing degree `Dr` (so the full message degree is
/// `D = Dr + 1`).
pub struct EqFactoredVirtualPolynomial<E: Engine> {
  num_vars: usize,
  taus: Vec<E::Scalar>,
  /// Non-eq factor tables (dense projective coefficients), bound each round.
  factors: Vec<Vec<E::Scalar>>,
  /// Terms over the non-eq factors; all share degree `Dr`.
  terms: Vec<(E::Scalar, Vec<usize>)>,
}

impl<E: Engine> EqFactoredVirtualPolynomial<E> {
  /// Builds an eq-factored relation. `taus.len()` must be `num_vars`, every
  /// factor table length `2^num_vars`, and all terms the same degree `Dr ≥ 1`.
  pub fn new(
    num_vars: usize,
    taus: Vec<E::Scalar>,
    factors: Vec<Vec<E::Scalar>>,
    terms: Vec<(E::Scalar, Vec<usize>)>,
  ) -> Self {
    assert_eq!(taus.len(), num_vars, "taus must have num_vars entries");
    let n = 1usize << num_vars;
    for f in &factors {
      assert_eq!(f.len(), n, "factor table length must be 2^num_vars");
    }
    assert!(!terms.is_empty(), "needs >= 1 term");
    let dr = terms[0].1.len();
    assert!(dr >= 1, "R term degree must be >= 1");
    for (_, idxs) in &terms {
      assert_eq!(idxs.len(), dr, "all R terms must share degree");
      for &j in idxs {
        assert!(j < factors.len(), "factor index out of range");
      }
    }
    Self {
      num_vars,
      taus,
      factors,
      terms,
    }
  }

  /// The full message degree `D = Dr + 1` (R's degree plus the eq linear factor).
  pub fn degree(&self) -> usize {
    self.terms[0].1.len() + 1
  }

  /// The projective sum `C_0 = Σ_b eq̃(τ,b) · R[b]`.
  fn initial_claim(&self, eq_full: &[E::Scalar]) -> E::Scalar {
    let n = 1usize << self.num_vars;
    (0..n)
      .into_par_iter()
      .map(|b| {
        let r: E::Scalar = self
          .terms
          .iter()
          .map(|(coeff, idxs)| idxs.iter().fold(*coeff, |acc, &j| acc * self.factors[j][b]))
          .sum();
        eq_full[b] * r
      })
      .sum()
  }

  /// Runs the eq-factored projective prover.
  pub fn prove(mut self, transcript: &mut E::TE) -> ProjectiveSumcheckProverOutput<E> {
    let num_vars = self.num_vars;
    let d = self.degree();
    let dr = d - 1;

    let eq_full = eq_corner_natural(&self.taus);
    let initial_claim = self.initial_claim(&eq_full);

    let mut rounds = Vec::with_capacity(num_vars);
    let mut point = Vec::with_capacity(num_vars);
    let mut eq_left = E::Scalar::ONE;

    let mut remaining = num_vars;
    for _ in 0..num_vars {
      let half = 1usize << (remaining - 1);
      let cur = remaining - 1; // current variable index (high bit)
      let (tau0, tau1) = (E::Scalar::ONE - self.taus[cur], self.taus[cur]);

      // Suffix eq weights over taus[0..cur] (low bit = tau_0), length `half`.
      let eq_suffix = eq_corner_natural(&self.taus[..cur]);

      // Inner = Σ_s eq_suffix[s] · R_slice(s, T), a degree-Dr polynomial.
      // Fast path: every term is exactly 2 factors (Dr = 2), the ppSNARK case.
      // Accumulate (c0, c1, c2) directly with no scratch alloc/zeroing.
      let all_deg2 = dr == 2 && self.terms.iter().all(|(_, idxs)| idxs.len() == 2);

      let fold_range = |lo: usize, hi: usize| -> Vec<E::Scalar> {
        if all_deg2 {
          let (mut c0, mut c1, mut c2) = (E::Scalar::ZERO, E::Scalar::ZERO, E::Scalar::ZERO);
          for suffix in lo..hi {
            let w = eq_suffix[suffix];
            if w == E::Scalar::ZERO {
              continue;
            }
            for (coeff, idxs) in &self.terms {
              let fa = &self.factors[idxs[0]];
              let fb = &self.factors[idxs[1]];
              let (a0, a1) = (fa[suffix], fa[suffix + half]);
              let (b0, b1) = (fb[suffix], fb[suffix + half]);
              let cw = *coeff * w;
              // (a0 + a1 T)(b0 + b1 T) scaled by cw.
              c0 += cw * a0 * b0;
              c1 += cw * (a0 * b1 + a1 * b0);
              c2 += cw * a1 * b1;
            }
          }
          return vec![c0, c1, c2];
        }
        // Generic path (any Dr): convolution via multiply_by_linear.
        let mut inner = vec![E::Scalar::ZERO; dr + 1];
        let mut scratch = vec![E::Scalar::ZERO; dr + 1];
        for suffix in lo..hi {
          let w = eq_suffix[suffix];
          if w == E::Scalar::ZERO {
            continue;
          }
          for (coeff, idxs) in &self.terms {
            scratch.iter_mut().for_each(|c| *c = E::Scalar::ZERO);
            scratch[0] = *coeff * w;
            for (cur_deg, &j) in idxs.iter().enumerate() {
              let table = &self.factors[j];
              super::virtual_poly::multiply_by_linear(
                &mut scratch,
                cur_deg,
                table[suffix],
                table[suffix + half],
              );
            }
            for (r, s) in inner.iter_mut().zip(scratch.iter()) {
              *r += *s;
            }
          }
        }
        inner
      };

      const PAR_THRESHOLD: usize = 1 << 10;
      let inner = if half < PAR_THRESHOLD {
        fold_range(0, half)
      } else {
        let chunk = (half / rayon::current_num_threads().max(1)).max(1);
        (0..half)
          .into_par_iter()
          .step_by(chunk)
          .map(|lo| fold_range(lo, (lo + chunk).min(half)))
          .reduce(
            || vec![E::Scalar::ZERO; dr + 1],
            |mut a, b| {
              for (x, y) in a.iter_mut().zip(b.iter()) {
                *x += *y;
              }
              a
            },
          )
      };

      // Multiply inner (degree Dr) by the eq linear factor (tau0 + tau1·T) and
      // by the running scalar eq_left → degree D message.
      let mut msg = vec![E::Scalar::ZERO; d + 1];
      for (k, c) in inner.iter().enumerate() {
        let cl = *c * eq_left;
        msg[k] += cl * tau0;
        msg[k + 1] += cl * tau1;
      }

      let poly =
        UniPoly::<E::Scalar>::from_coeffs_no_trim(msg).expect("degree D+1 >= 2 coefficients");

      transcript.absorb(b"projective_sumcheck_round", &poly);
      let r_i = transcript
        .squeeze(b"projective_sumcheck_challenge")
        .expect("transcript squeeze failed");

      rounds.push(poly.compress_projective());
      point.push(r_i);

      // Fold eq_left by the current variable's linear factor at r_i, and bind
      // the non-eq factors.
      eq_left *= tau0 + tau1 * r_i;
      for table in &mut self.factors {
        let (lo, hi) = table.split_at_mut(half);
        lo.par_iter_mut()
          .zip(hi.par_iter())
          .for_each(|(a, b)| *a += r_i * *b);
        table.truncate(half);
      }
      remaining -= 1;
    }

    // Final reduced claim: eq_left · R(r) = eq_left · Σ_t coeff_t ∏_j F_j(r).
    let r_of_r: E::Scalar = self
      .terms
      .iter()
      .map(|(coeff, idxs)| idxs.iter().fold(*coeff, |acc, &j| acc * self.factors[j][0]))
      .sum();
    let final_claim = eq_left * r_of_r;

    ProjectiveSumcheckProverOutput {
      proof: SumcheckProof::new(rounds),
      initial_claim,
      point,
      final_claim,
    }
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::{
    provider::PallasEngine,
    spartan::projective_sumcheck::{verify, VirtualPolynomial},
    traits::Engine,
  };
  use ff::Field;

  type E = PallasEngine;
  type Fr = <E as Engine>::Scalar;

  /// Differential test: the eq-factored prover must produce the SAME proof
  /// (initial claim, point, final claim) as the naive path that includes Eq^∞
  /// as an explicit dense factor. Relation: Eq^∞_τ · (f·g − h·U-less), i.e. two
  /// R-terms of equal degree 2.
  #[test]
  fn eq_factored_matches_naive() {
    let num_vars = 4usize;
    let n = 1usize << num_vars;
    let f: Vec<Fr> = (0..n).map(|i| Fr::from((2 * i + 1) as u64)).collect();
    let g: Vec<Fr> = (0..n).map(|i| Fr::from((3 * i + 5) as u64)).collect();
    let h: Vec<Fr> = (0..n).map(|i| Fr::from((7 * i + 2) as u64)).collect();
    let k: Vec<Fr> = (0..n).map(|i| Fr::from((i + 4) as u64)).collect();
    let taus: Vec<Fr> = (0..num_vars)
      .map(|i| Fr::from((5 * i + 3) as u64))
      .collect();

    // R = f·g + h·k  (Dr = 2), full degree D = 3.
    // Eq-factored path.
    let efvp = EqFactoredVirtualPolynomial::<E>::new(
      num_vars,
      taus.clone(),
      vec![f.clone(), g.clone(), h.clone(), k.clone()],
      vec![(Fr::ONE, vec![0, 1]), (Fr::ONE, vec![2, 3])],
    );
    let mut ts_e = <E as Engine>::TE::new(b"eqfac");
    let out_e = efvp.prove(&mut ts_e);

    // Naive path: eq as explicit factor 0, terms Eq·f·g and Eq·h·k (degree 3).
    // Use the SAME natural-order eq the factored prover uses.
    let eq = eq_corner_natural(&taus);
    let naive = VirtualPolynomial::<E>::new(
      num_vars,
      vec![eq, f, g, h, k],
      vec![(Fr::ONE, vec![0, 1, 2]), (Fr::ONE, vec![0, 3, 4])],
    );
    let mut ts_n = <E as Engine>::TE::new(b"eqfac");
    let out_n = naive.prove(&mut ts_n);

    assert_eq!(out_e.initial_claim, out_n.initial_claim);
    assert_eq!(out_e.point, out_n.point);
    assert_eq!(out_e.final_claim, out_n.final_claim);

    // And it verifies at D = 3.
    let degree_bounds = vec![3usize; num_vars];
    let mut ts_v = <E as Engine>::TE::new(b"eqfac");
    let reduction =
      verify::<E>(out_e.initial_claim, &degree_bounds, &out_e.proof, &mut ts_v).unwrap();
    assert_eq!(reduction.point, out_e.point);
    assert_eq!(reduction.final_claim, out_e.final_claim);
  }
}
