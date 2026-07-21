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

/// Builds the incremental prefix tables of the natural-order projective eq
/// polynomial: `out[k][b] = ∏_{i<k} (bit_i(b) ? τ_i : 1−τ_i)` for `k = 0..=len`,
/// with bit `i` corresponding to `τ_i` (low bit = `τ_0`). Each table is derived
/// from the previous by one doubling step, so building all of them costs
/// `Σ_k 2^k ≈ 2·2^len` field mults total.
fn eq_prefix_tables<F: Field>(taus: &[F]) -> Vec<Vec<F>> {
  let mut out = Vec::with_capacity(taus.len() + 1);
  out.push(vec![F::ONE]);
  for (i, &tau) in taus.iter().enumerate() {
    let prev = &out[i];
    let len = prev.len();
    let mut next = vec![F::ZERO; len * 2];
    for (b, &v) in prev.iter().enumerate() {
      next[b] = v * (F::ONE - tau);
      next[b + len] = v * tau;
    }
    out.push(next);
  }
  out
}

/// Split-half ("Gruen", eprint 2024/108 §3) representation of the natural-order
/// projective eq corner weights, in this module's own convention (bit `j ↔ τ_j`,
/// high index bound first). Stores the eq factor as two sides of `O(√N)` prefix
/// tables split at `mid = num_vars/2`, so the round kernel reads a suffix weight
/// by one multiplication of a low-side and a high-side entry — instead of
/// rebuilding a dense `O(2^cur)` suffix table (and reallocating) every round.
struct SplitEqNatural<F: Field> {
  mid: usize,
  /// `lo_prefix[k] = eq(τ_0..τ_{k-1})`, `k = 0..=mid`.
  lo_prefix: Vec<Vec<F>>,
  /// `hi_prefix[j] = eq(τ_mid..τ_{mid+j-1})`, `j = 0..=(num_vars − mid)`.
  hi_prefix: Vec<Vec<F>>,
}

impl<F: Field> SplitEqNatural<F> {
  fn new(taus: &[F]) -> Self {
    let mid = taus.len() / 2;
    let lo_prefix = eq_prefix_tables(&taus[..mid]);
    let hi_prefix = eq_prefix_tables(&taus[mid..]);
    Self {
      mid,
      lo_prefix,
      hi_prefix,
    }
  }

  /// Tables for the round whose surviving suffix spans the first `cur` variables:
  /// the weight of suffix index `s` is `lo[s & ((1<<bits)-1)] · hi[s >> bits]`.
  #[inline]
  fn round(&self, cur: usize) -> (&[F], &[F], usize) {
    if cur >= self.mid {
      (
        &self.lo_prefix[self.mid],
        &self.hi_prefix[cur - self.mid],
        self.mid,
      )
    } else {
      // hi_prefix[0] = [ONE]; the high factor collapses to 1.
      (&self.lo_prefix[cur], &self.hi_prefix[0], cur)
    }
  }

  /// Tables for the full eq over all `num_vars` variables (used for the initial
  /// claim): `eq(b) = lo[b & ((1<<mid)-1)] · hi[b >> mid]`.
  #[inline]
  fn full(&self) -> (&[F], &[F], usize) {
    (
      &self.lo_prefix[self.mid],
      self.hi_prefix.last().unwrap(),
      self.mid,
    )
  }
}

/// Natural-order projective eq corner table `t[b] = ∏_i (bit_i(b) ? τ_i :
/// 1−τ_i)`. Test-only reference builder mirroring [`SplitEqNatural`]; the prover
/// uses the split tables directly.
#[cfg(test)]
fn eq_corner_natural<F: Field>(taus: &[F]) -> Vec<F> {
  eq_prefix_tables(taus).pop().unwrap()
}

/// A relation `Eq^∞_τ · R`, where `R = Σ_t coeff_t ∏_{j∈term_t} F_j` over the
/// non-eq factors, all terms sharing degree `Dr` (so the full message degree is
/// `D = Dr + 1`).
pub struct EqFactoredVirtualPolynomial<E: Engine> {
  num_vars: usize,
  taus: Vec<E::Scalar>,
  /// Non-eq factor tables (dense projective coefficients), bound each round.
  factors: Vec<Vec<E::Scalar>>,
  /// Terms over the non-eq factors; term `t` carries `Dr − |term_t|` analytic U
  /// copies (mixed-degree; homogenized to `Dr = max_t |term_t|`).
  terms: Vec<(E::Scalar, Vec<usize>)>,
  // --- prover state (initialized at construction, advanced by `bind_round`) ---
  /// Split-half ("Gruen") eq tables, built once.
  split: SplitEqNatural<E::Scalar>,
  /// R-degree `Dr = max_t |term_t|`.
  dr: usize,
  /// Variables not yet bound (starts at `num_vars`).
  remaining: usize,
  /// Running eq prefix scalar `∏(τ0_k + τ1_k r_k)`.
  eq_left: E::Scalar,
  /// Running structured-U scalar `∏(1 + r_k)`.
  u_left: E::Scalar,
  /// Running claim `C_i = S_{i-1}(r_{i-1})` (BDDT hint).
  claim: E::Scalar,
  /// The projective sum `C_0` (cached for `claim0`).
  initial_claim: E::Scalar,
  two: E::Scalar,
  four: E::Scalar,
}

impl<E: Engine> EqFactoredVirtualPolynomial<E> {
  /// Builds an eq-factored relation with terms all of the same degree `Dr ≥ 1`
  /// (a special case of [`new_mixed`](Self::new_mixed) with no analytic-U
  /// homogenization). `taus.len()` must be `num_vars` and every factor table
  /// length `2^num_vars`.
  pub fn new(
    num_vars: usize,
    taus: Vec<E::Scalar>,
    factors: Vec<Vec<E::Scalar>>,
    terms: Vec<(E::Scalar, Vec<usize>)>,
  ) -> Self {
    let dr = terms[0].1.len();
    for (_, idxs) in &terms {
      assert_eq!(idxs.len(), dr, "all R terms must share degree");
    }
    Self::new_mixed(num_vars, taus, factors, terms)
  }

  /// Builds an eq-factored relation with **mixed-degree** terms, homogenized to
  /// the common R-degree `Dr = max_t |term_t|` with an **analytic** all-ones
  /// factor `U(X) = ∏_i (1 + X_i)`: a term of `m` real factors is treated as
  /// carrying `Dr − m` copies of `U`. `U` is never materialized as a dense
  /// `2^num_vars` table nor bound — the prover tracks the running scalar
  /// `u_left = ∏(1 + r_k)` and contributes `(1 + T)` per `U` copy each round
  /// (design §5; structured-`U` perf item). `taus.len()` must be `num_vars`,
  /// every factor table length `2^num_vars`, every real index in range, and
  /// `Dr ≥ 1`.
  pub fn new_mixed(
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
    let dr = terms.iter().map(|(_, idxs)| idxs.len()).max().unwrap();
    assert!(dr >= 1, "R term degree Dr must be >= 1");
    for (_, idxs) in &terms {
      for &j in idxs {
        assert!(j < factors.len(), "factor index out of range");
      }
    }

    let split = SplitEqNatural::<E::Scalar>::new(&taus);
    let initial_claim = Self::compute_initial_claim(num_vars, &factors, &terms, &split);
    let two = E::Scalar::ONE.double();
    let four = two.double();

    Self {
      num_vars,
      taus,
      factors,
      terms,
      split,
      dr,
      remaining: num_vars,
      eq_left: E::Scalar::ONE,
      u_left: E::Scalar::ONE,
      claim: initial_claim,
      initial_claim,
      two,
      four,
    }
  }

  /// The full message degree `D = Dr + 1` (R's degree `Dr = max_t |term_t|` plus
  /// the eq linear factor).
  pub fn degree(&self) -> usize {
    self.dr + 1
  }

  /// Raises the message degree to `target` (`≥ current`) by carrying more
  /// analytic `U` copies per term (`U(b) = 1` at every projective corner, so the
  /// projective sum — and thus the cached initial claim — is unchanged). Used to
  /// bring instances to a shared degree before batching.
  pub(crate) fn lift_to_degree(&mut self, target: usize) {
    assert!(target > self.dr, "cannot lift to a lower degree");
    self.dr = target - 1;
  }

  /// The number of variables `n`.
  pub fn num_vars(&self) -> usize {
    self.num_vars
  }

  /// The projective sum `C_0 = Σ_b eq̃(τ,b) · R[b]`.
  pub(crate) fn claim0(&self) -> E::Scalar {
    self.initial_claim
  }

  /// Computes `C_0 = Σ_b eq̃(τ,b) · R[b]` from the split-half eq tables (no dense
  /// `2^n` eq materialization). `R[b]` uses only real factors (`U[b] = 1`).
  fn compute_initial_claim(
    num_vars: usize,
    factors: &[Vec<E::Scalar>],
    terms: &[(E::Scalar, Vec<usize>)],
    split: &SplitEqNatural<E::Scalar>,
  ) -> E::Scalar {
    let n = 1usize << num_vars;
    let (lo, hi, bits) = split.full();
    let lomask = (1usize << bits) - 1;
    (0..n)
      .into_par_iter()
      .map(|b| {
        let eqb = lo[b & lomask] * hi[b >> bits];
        let r: E::Scalar = terms
          .iter()
          .map(|(coeff, idxs)| idxs.iter().fold(*coeff, |acc, &j| acc * factors[j][b]))
          .sum();
        eqb * r
      })
      .sum()
  }

  /// Builds the current round's degree-`D` message coefficients `[a_0, …, a_D]`
  /// from the unbound tables and the running eq/U scalars, deriving one inner
  /// coefficient from the running claim (BDDT). Pure: mutates no state and does
  /// not touch the transcript.
  pub(crate) fn round_message(&self) -> Vec<E::Scalar> {
    let dr = self.dr;
    let d = dr + 1;
    let remaining = self.remaining;
    let half = 1usize << (remaining - 1);
    let cur = remaining - 1; // current variable index (high bit)
    let (tau0, tau1) = (E::Scalar::ONE - self.taus[cur], self.taus[cur]);
    let eq_left = self.eq_left;
    let u_left = self.u_left;
    let claim = self.claim;
    let two = self.two;
    let four = self.four;

    // Suffix eq weights over taus[0..cur], read from the split-half tables:
    // eq_suffix[s] = eq_lo[s & lomask] · eq_hi[s >> split_bits].
    let (eq_lo, eq_hi, split_bits) = self.split.round(cur);
    let lomask = (1usize << split_bits) - 1;

    // Inner = Σ_s eq_suffix[s] · R_slice(s, T), a degree-Dr polynomial.
    // Fast path: Dr = 2 (the ppSNARK case); each term has 0, 1, or 2 real
    // factors, the rest being analytic U copies. Accumulate the endpoints
    // {inner(0)?, inner(1), inner(∞)} directly with no scratch alloc/zeroing.
    let all_deg2 = dr == 2;

    // BDDT (eprint 2025/1117 §6.2) in projective form. The round claim identity
    // `C = a_0 + a_D = eq_left·(tau0·c0 + tau1·c2)` gives one system-linear
    // constraint on the inner deg-2 coefficients, so we scan only the two
    // endpoints `{inner(1), inner(∞)=c2}` (2 N-sums) and recover the constant
    // `c0` from the running claim, then the cross term `c1 = inner(1) − c0 − c2`.
    // This drops the direct cross-term product (the most expensive per-cell
    // term). Falls back to also scanning `c0 = inner(0)` when `tau0` or
    // `eq_left` is zero (both challenge-measure-zero).
    let eq_left_inv: Option<E::Scalar> = eq_left.invert().into();
    let bddt = all_deg2 && tau0 != E::Scalar::ZERO && eq_left_inv.is_some();

    let fold_range = |lo: usize, hi: usize| -> Vec<E::Scalar> {
      if all_deg2 {
        // acc = [s_zero, s_one, s_inf]; s_zero stays 0 on the BDDT (2-sum) path.
        let (mut s_zero, mut s_one, mut s_inf) =
          (E::Scalar::ZERO, E::Scalar::ZERO, E::Scalar::ZERO);
        for suffix in lo..hi {
          let w = eq_lo[suffix & lomask] * eq_hi[suffix >> split_bits];
          if w == E::Scalar::ZERO {
            continue;
          }
          for (coeff, idxs) in &self.terms {
            let cw = *coeff * w;
            // R-term slice = (∏ real linear factors) · (u_left·(1+T))^{Dr−m}.
            match idxs.len() {
              2 => {
                let fa = &self.factors[idxs[0]];
                let fb = &self.factors[idxs[1]];
                let (a0, a1) = (fa[suffix], fa[suffix + half]);
                let (b0, b1) = (fb[suffix], fb[suffix + half]);
                s_inf += cw * a1 * b1;
                s_one += cw * (a0 + a1) * (b0 + b1);
                if !bddt {
                  s_zero += cw * a0 * b0;
                }
              }
              1 => {
                // real · U = u_left·(a0 + a1 T)·(1 + T).
                let fa = &self.factors[idxs[0]];
                let (a0, a1) = (fa[suffix], fa[suffix + half]);
                let cwu = cw * u_left;
                s_inf += cwu * a1;
                s_one += cwu * (a0 + a1) * two;
                if !bddt {
                  s_zero += cwu * a0;
                }
              }
              0 => {
                // U² = u_left²·(1 + T)².
                let cwu2 = cw * u_left * u_left;
                s_inf += cwu2;
                s_one += cwu2 * four;
                if !bddt {
                  s_zero += cwu2;
                }
              }
              _ => unreachable!("all_deg2 implies term has at most 2 real factors"),
            }
          }
        }
        return vec![s_zero, s_one, s_inf];
      }
      // Generic path (any Dr): convolution via multiply_by_linear.
      let mut inner = vec![E::Scalar::ZERO; dr + 1];
      let mut scratch = vec![E::Scalar::ZERO; dr + 1];
      for suffix in lo..hi {
        let w = eq_lo[suffix & lomask] * eq_hi[suffix >> split_bits];
        if w == E::Scalar::ZERO {
          continue;
        }
        for (coeff, idxs) in &self.terms {
          scratch.iter_mut().for_each(|c| *c = E::Scalar::ZERO);
          scratch[0] = *coeff * w;
          let mut cur_deg = 0;
          for &j in idxs {
            let table = &self.factors[j];
            super::virtual_poly::multiply_by_linear(
              &mut scratch,
              cur_deg,
              table[suffix],
              table[suffix + half],
            );
            cur_deg += 1;
          }
          // Analytic U homogenization: (u_left·(1+T)) per missing factor.
          for _ in 0..(dr - idxs.len()) {
            super::virtual_poly::multiply_by_linear(&mut scratch, cur_deg, u_left, u_left);
            cur_deg += 1;
          }
          for (r, s) in inner.iter_mut().zip(scratch.iter()) {
            *r += *s;
          }
        }
      }
      inner
    };

    const PAR_THRESHOLD: usize = 1 << 10;
    let agg = if half < PAR_THRESHOLD {
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

    // Reconstruct the inner deg-Dr coefficients. On the deg-2 fast path this
    // turns the scanned endpoints into [c0, c1, c2] (recovering c0 from the
    // claim under BDDT); the generic path already produced the coefficients.
    let inner = if all_deg2 {
      let s_one = agg[1];
      let c2 = agg[2];
      let c0 = if bddt {
        // c0 = (C·eq_left⁻¹ − tau1·c2)·tau0⁻¹.
        let tau0_inv = Option::<E::Scalar>::from(tau0.invert()).expect("tau0 != 0 under bddt");
        (claim * eq_left_inv.expect("eq_left invertible under bddt") - tau1 * c2) * tau0_inv
      } else {
        agg[0]
      };
      let c1 = s_one - c0 - c2;
      vec![c0, c1, c2]
    } else {
      agg
    };

    // Multiply inner (degree Dr) by the eq linear factor (tau0 + tau1·T) and
    // by the running scalar eq_left → the degree-D message coefficients.
    let mut msg = vec![E::Scalar::ZERO; d + 1];
    for (k, c) in inner.iter().enumerate() {
      let cl = *c * eq_left;
      msg[k] += cl * tau0;
      msg[k + 1] += cl * tau1;
    }
    msg
  }

  /// Advances the prover state by binding the current variable to `r`: updates
  /// the running claim from this instance's own message (`C_{i+1} = S_i(r)`),
  /// folds `eq_left` and `u_left`, binds the factor tables, and decrements
  /// `remaining`. `msg` must be this instance's [`round_message`](Self::round_message)
  /// output for the current round.
  pub(crate) fn bind_round(&mut self, r: E::Scalar, msg: &[E::Scalar]) {
    let cur = self.remaining - 1;
    let half = 1usize << cur;
    let (tau0, tau1) = (E::Scalar::ONE - self.taus[cur], self.taus[cur]);
    // C_{i+1} = S_i(r) via Horner over the ascending coefficients.
    let mut c_next = E::Scalar::ZERO;
    for coeff in msg.iter().rev() {
      c_next = c_next * r + *coeff;
    }
    self.claim = c_next;
    self.eq_left *= tau0 + tau1 * r;
    self.u_left *= E::Scalar::ONE + r;
    for table in &mut self.factors {
      let (lo, hi) = table.split_at_mut(half);
      lo.par_iter_mut()
        .zip(hi.par_iter())
        .for_each(|(a, b)| *a += r * *b);
      table.truncate(half);
    }
    self.remaining -= 1;
  }

  /// After all variables are bound, the reduced value `eq_left · R(r)`, where
  /// each term contributes `coeff · U(r)^{Dr−m} · ∏ real F_j(r)` and
  /// `U(r) = u_left = ∏(1 + r_k)`.
  pub(crate) fn reduced_claim(&self) -> E::Scalar {
    let dr = self.dr;
    let r_of_r: E::Scalar = self
      .terms
      .iter()
      .map(|(coeff, idxs)| {
        let base = idxs.iter().fold(*coeff, |acc, &j| acc * self.factors[j][0]);
        base * self.u_left.pow_vartime([(dr - idxs.len()) as u64])
      })
      .sum();
    self.eq_left * r_of_r
  }

  /// Runs the eq-factored projective prover (single instance).
  pub fn prove(mut self, transcript: &mut E::TE) -> ProjectiveSumcheckProverOutput<E> {
    let num_vars = self.num_vars;
    let initial_claim = self.initial_claim;
    let mut rounds = Vec::with_capacity(num_vars);
    let mut point = Vec::with_capacity(num_vars);
    for _ in 0..num_vars {
      let msg = self.round_message();
      let poly =
        UniPoly::<E::Scalar>::from_coeffs_no_trim(msg).expect("degree D+1 >= 2 coefficients");
      transcript.absorb(b"projective_sumcheck_round", &poly);
      let r_i = transcript
        .squeeze(b"projective_sumcheck_challenge")
        .expect("transcript squeeze failed");
      rounds.push(poly.compress_projective());
      point.push(r_i);
      self.bind_round(r_i, poly.coeffs());
    }
    let final_claim = self.reduced_claim();

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

  /// Differential check that the eq-factored prover (with BDDT claim-derivation)
  /// matches the naive dense-eq path over `R = f·g + h·k`, for a given tau
  /// vector. Shared by the general and the `tau = 1` fallback tests.
  fn check_matches_naive(num_vars: usize, taus: Vec<Fr>) {
    let n = 1usize << num_vars;
    let f: Vec<Fr> = (0..n).map(|i| Fr::from((2 * i + 1) as u64)).collect();
    let g: Vec<Fr> = (0..n).map(|i| Fr::from((3 * i + 5) as u64)).collect();
    let h: Vec<Fr> = (0..n).map(|i| Fr::from((7 * i + 2) as u64)).collect();
    let k: Vec<Fr> = (0..n).map(|i| Fr::from((i + 4) as u64)).collect();

    let efvp = EqFactoredVirtualPolynomial::<E>::new(
      num_vars,
      taus.clone(),
      vec![f.clone(), g.clone(), h.clone(), k.clone()],
      vec![(Fr::ONE, vec![0, 1]), (Fr::ONE, vec![2, 3])],
    );
    let mut ts_e = <E as Engine>::TE::new(b"eqfac_bddt");
    let out_e = efvp.prove(&mut ts_e);

    let eq = eq_corner_natural(&taus);
    let naive = VirtualPolynomial::<E>::new(
      num_vars,
      vec![eq, f, g, h, k],
      vec![(Fr::ONE, vec![0, 1, 2]), (Fr::ONE, vec![0, 3, 4])],
    );
    let mut ts_n = <E as Engine>::TE::new(b"eqfac_bddt");
    let out_n = naive.prove(&mut ts_n);

    assert_eq!(out_e.initial_claim, out_n.initial_claim);
    assert_eq!(out_e.point, out_n.point);
    assert_eq!(out_e.final_claim, out_n.final_claim);
    // Per-round messages must be byte-identical (proves BDDT reconstructs the
    // same round polynomial, not merely a verifying one).
    let polys_e = out_e.proof.compressed_polys();
    let polys_n = out_n.proof.compressed_polys();
    assert_eq!(polys_e.len(), polys_n.len());
    for (a, b) in polys_e.iter().zip(polys_n.iter()) {
      assert_eq!(a.stored_coeffs(), b.stored_coeffs());
    }
  }

  /// BDDT fast path (all `tau ≠ 1`) at a larger size must match the naive path
  /// byte-for-byte.
  #[test]
  fn eq_factored_bddt_matches_naive_large() {
    let num_vars = 8usize;
    let taus: Vec<Fr> = (0..num_vars)
      .map(|i| Fr::from((5 * i + 3) as u64))
      .collect();
    check_matches_naive(num_vars, taus);
  }

  /// The `tau = 1` fallback (`tau0 = 0`, BDDT can't invert) must still match the
  /// naive path. Places `tau = 1` at several positions, including the last
  /// variable bound (round 0) and the first.
  #[test]
  fn eq_factored_bddt_tau_one_fallback() {
    let num_vars = 6usize;
    for one_at in [0usize, 1, num_vars - 1] {
      let taus: Vec<Fr> = (0..num_vars)
        .map(|i| {
          if i == one_at {
            Fr::ONE
          } else {
            Fr::from((3 * i + 2) as u64)
          }
        })
        .collect();
      check_matches_naive(num_vars, taus);
    }
  }

  /// Structured-U path: an eq-factored relation with **mixed-degree** terms
  /// `Eq^∞_τ·(f·g − h·U)` (built via `new_mixed`, the second term carrying one
  /// analytic U) must match the naive path that materializes both eq and U as
  /// dense factors, byte-for-byte. Also covers the `tau = 1` fallback.
  #[test]
  fn eq_factored_structured_u_matches_naive() {
    for taus_pick in 0..2usize {
      let num_vars = 6usize;
      let n = 1usize << num_vars;
      let f: Vec<Fr> = (0..n).map(|i| Fr::from((2 * i + 1) as u64)).collect();
      let g: Vec<Fr> = (0..n).map(|i| Fr::from((3 * i + 5) as u64)).collect();
      let h: Vec<Fr> = (0..n).map(|i| Fr::from((7 * i + 2) as u64)).collect();
      let taus: Vec<Fr> = (0..num_vars)
        .map(|i| {
          if taus_pick == 1 && i == 2 {
            Fr::ONE // exercise the tau=1 fallback with structured U
          } else {
            Fr::from((5 * i + 3) as u64)
          }
        })
        .collect();

      // Eq-factored with analytic U: R = f·g − h·U (terms of degree 2 and 1).
      let efvp = EqFactoredVirtualPolynomial::<E>::new_mixed(
        num_vars,
        taus.clone(),
        vec![f.clone(), g.clone(), h.clone()],
        vec![(Fr::ONE, vec![0, 1]), (-Fr::ONE, vec![2])],
      );
      let mut ts_e = <E as Engine>::TE::new(b"eqfac_u");
      let out_e = efvp.prove(&mut ts_e);

      // Naive: eq (factor 0) and U (added by new_homogenized) are dense. Terms
      // Eq·f·g (degree 3) and −Eq·h (degree 2, homogenized by U to 3).
      let eq = eq_corner_natural(&taus);
      let naive = VirtualPolynomial::<E>::new_homogenized(
        num_vars,
        vec![eq, f, g, h],
        vec![(Fr::ONE, vec![0, 1, 2]), (-Fr::ONE, vec![0, 3])],
      );
      let mut ts_n = <E as Engine>::TE::new(b"eqfac_u");
      let out_n = naive.prove(&mut ts_n);

      assert_eq!(out_e.initial_claim, out_n.initial_claim);
      assert_eq!(out_e.point, out_n.point);
      assert_eq!(out_e.final_claim, out_n.final_claim);
      let polys_e = out_e.proof.compressed_polys();
      let polys_n = out_n.proof.compressed_polys();
      assert_eq!(polys_e.len(), polys_n.len());
      for (a, b) in polys_e.iter().zip(polys_n.iter()) {
        assert_eq!(a.stored_coeffs(), b.stored_coeffs());
      }
    }
  }
}
