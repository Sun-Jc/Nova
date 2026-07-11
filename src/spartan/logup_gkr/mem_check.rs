//! Host-layer memory-check verifier for ppSNARK, built on Logup-GKR.
//!
//! **Verifier-first.** This module is the *host* half of the Logup-GKR
//! memory-check: it wraps the frozen `logup_gkr::verify` (which owns no PCS)
//! and closes soundness by (1) recomputing each logup instance's input-layer
//! fraction from the columns the prover opens at the GKR evaluation point, (2)
//! checking those match what the GKR reduced to, and (3) running the fractional
//! balance check. It is written **before** any prover: the set of evaluations
//! [`MemCheckOpenings`] names here *is* the contract a prover must satisfy — the
//! prover must open exactly these columns at exactly [`eval_point`], and no
//! others.
//!
//! [`eval_point`]: LogupGkrOpeningClaim::eval_point
//!
//! ## The four sub-instances (why four, not two)
//! ppSNARK's memory-check is two logup relations (`row`, `col`), each a balance
//! `Σ ts/(T+r) = Σ 1/(W+r)`. We encode each relation as **two** height-N
//! fractional-sum sub-instances — a *table* side and an *access* side — so all
//! four share one GKR depth `log N` and the frozen single-batched-tree verifier
//! applies unchanged (it forbids uneven heights; here every side is exactly N
//! because ppSNARK pads every memory-check column to N in setup — see
//! `jcbase/ppsnark-pad-to-N.md`). A single 2N-leaf signed-multiplicity tree
//! (hyperplonk's shape) would instead emit a `log(2N)` point, which cannot be
//! rerandomized against the N-variable inner sumcheck; two N-leaf trees keep the
//! point at `log N`. Instance order is fixed:
//!
//! | idx | name        | num       | den                         |
//! |-----|-------------|-----------|-----------------------------|
//! | 0   | row_table   | `ts_row`  | `mem_row·γ + id + r`        |
//! | 1   | row_access  | `-1`      | `L_row·γ + addr_row + r`    |
//! | 2   | col_table   | `ts_col`  | `mem_col·γ + id + r`        |
//! | 3   | col_access  | `-1`      | `L_col·γ + addr_col + r`    |
//!
//! where `id = IdentityPolynomial` (the cell address `i`), `mem_row =
//! eq(r_outer_full, ·)`, `mem_col = z`, `addr_row = row`, `addr_col = col`, and
//! `(γ, r)` are the ppSNARK memory-check fingerprint challenges.
//!
//! ## Balance (the host's `0/den` check)
//! The GKR verifier reduces each sub-instance to a root fraction but does **not**
//! check the relation balances — that is this module's job. Balance is a
//! property of the whole relation, i.e. the *root* fractions (the proof's
//! `initial_claims`, bound to the reduction by the GKR verifier), not the
//! input-layer evaluations. The access side carries `num = -1`, so `root_table +
//! root_access = Σ ts/(T+r) − Σ 1/(W+r)`, which must vanish. In projective form
//! that means the *sum's numerator* is zero (the denominator, a product of
//! nonzero dens, cannot be), checked for the row pair `(0,1)` and the col pair
//! `(2,3)`. (This is looser than hyperplonk's exact `num == 0` on one merged
//! signed-multiplicity tree — audit item M1 — but equivalent for this padded
//! encoding.)

use crate::errors::NovaError;
use crate::spartan::logup_gkr::fraction::Fraction;
use crate::spartan::logup_gkr::proof::LogupGkrProof;
use crate::spartan::logup_gkr::verifier;
use crate::spartan::polys::eq::EqPolynomial;
use crate::spartan::polys::identity::IdentityPolynomial;
use crate::traits::Engine;
use ff::Field;

/// Fixed sub-instance count (`row_table, row_access, col_table, col_access`).
pub const NUM_SUB_INSTANCES: usize = 4;

/// The column evaluations a prover **must** open at the GKR
/// [`eval_point`](crate::spartan::logup_gkr::LogupGkrOpeningClaim::eval_point).
///
/// This is the prover's opening contract, defined by the verifier. Every field
/// is one column evaluated at the shared GKR point; the host reconstructs the
/// four sub-instance fractions from these and compares against the GKR's
/// reduced `openings`. The verifier recomputes the remaining pieces itself
/// (`id`, `mem_row = eq(r_outer_full, ·)`), so they are not opened here.
#[derive(Clone, Copy, Debug)]
pub struct MemCheckOpenings<E: Engine> {
  /// `L_row(eval_point)` — the row lookup column.
  pub eval_L_row: E::Scalar,
  /// `L_col(eval_point)` — the col lookup column.
  pub eval_L_col: E::Scalar,
  /// `row(eval_point)` — the row access-address column (`addr_row`).
  pub eval_row: E::Scalar,
  /// `col(eval_point)` — the col access-address column (`addr_col`).
  pub eval_col: E::Scalar,
  /// `ts_row(eval_point)` — the row multiplicity column.
  pub eval_ts_row: E::Scalar,
  /// `ts_col(eval_point)` — the col multiplicity column.
  pub eval_ts_col: E::Scalar,
  /// `mem_col(eval_point) = z(eval_point)` — the col table-value column.
  ///
  /// `mem_row` is `eq(r_outer_full, ·)`, which the verifier evaluates directly,
  /// so only `mem_col` needs opening.
  pub eval_mem_col: E::Scalar,
}

/// Verifies the ppSNARK memory-check via Logup-GKR.
///
/// Steps: (1) run the frozen GKR verifier to get the shared `eval_point` and the
/// four reduced input-layer fractions; (2) recompute those four fractions from
/// the prover's opened columns ([`MemCheckOpenings`]) and the fingerprint
/// `(gamma, r)`, and require they match; (3) check the two balances. Returns the
/// `eval_point` on success, so the caller can fold it into its batched PCS
/// opening set. The transcript must be positioned exactly as the prover left it
/// (GKR proof absorbed in the same slot).
///
/// `r_outer_full` is ppSNARK's extended outer challenge, defining `mem_row =
/// eq(r_outer_full, ·)`; the verifier evaluates it at `eval_point` itself.
pub fn verify<E: Engine>(
  proof: &LogupGkrProof<E>,
  gamma: E::Scalar,
  r: E::Scalar,
  r_outer_full: &[E::Scalar],
  openings: &MemCheckOpenings<E>,
  transcript: &mut E::TE,
) -> Result<Vec<E::Scalar>, NovaError> {
  // (1) Frozen GKR verifier: shape-check, root gates, per-layer sumchecks.
  let claim = verifier::verify::<E>(proof, transcript)?;
  let eval_point = claim.eval_point();
  let reduced = claim.openings();
  if reduced.len() != NUM_SUB_INSTANCES {
    return Err(NovaError::InvalidNumInstances);
  }

  // (2) Recompute the four input-layer fractions from the opened columns.
  // Pieces the verifier evaluates itself at eval_point:
  let eval_id = IdentityPolynomial::<E::Scalar>::new(eval_point.len()).evaluate(eval_point);
  let eval_mem_row = EqPolynomial::new(r_outer_full.to_vec()).evaluate(eval_point);

  let neg_one = -E::Scalar::ONE;

  // idx 0: row_table    num = ts_row    den = mem_row·γ + id + r
  let row_table = {
    let num = openings.eval_ts_row;
    let den = eval_mem_row * gamma + eval_id + r;
    Fraction::new(num, den)
  };

  // idx 1: row_access   num = -1        den = L_row·γ + addr_row + r
  let row_access = {
    let num = neg_one;
    let den = openings.eval_L_row * gamma + openings.eval_row + r;
    Fraction::new(num, den)
  };

  // idx 2: col_table    num = ts_col    den = mem_col·γ + id + r
  let col_table = {
    let num = openings.eval_ts_col;
    let den = openings.eval_mem_col * gamma + eval_id + r;
    Fraction::new(num, den)
  };

  // idx 3: col_access   num = -1        den = L_col·γ + addr_col + r
  let col_access = {
    let num = neg_one;
    let den = openings.eval_L_col * gamma + openings.eval_col + r;
    Fraction::new(num, den)
  };

  let recomputed = [row_table, row_access, col_table, col_access];

  // Each recomputed fraction must match the GKR-reduced input-layer fraction
  // (cross-multiplicative equality, since neither side is normalized).
  for (rc, red) in recomputed.iter().zip(reduced.iter()) {
    if rc.num * red.den != red.num * rc.den {
      return Err(NovaError::InvalidSumcheckProof);
    }
  }

  // (3) Balance: table side + access side = 0, for the row relation and the col
  // relation. Balance is a property of the whole relation, i.e. the *root*
  // fractions (`Σ ts/(T+r)` and `Σ -1/(W+r)`) — these are the GKR proof's
  // `initial_claims`, NOT the input-layer `recomputed` fractions from step (2)
  // (those are single-point evaluations at eval_point, a different quantity).
  // The roots are already bound to the reduction by the GKR verifier's root-gate
  // + per-layer sumcheck checks in step (1). The access side carries num = -1,
  // so `root_table + root_access = Σ ts/(T+r) − Σ 1/(W+r)`, which must vanish:
  // in projective form the sum's numerator is zero once its denominator (a
  // product of dens) is confirmed nonzero — checked explicitly just below.
  let [row_table_root, row_access_root, col_table_root, col_access_root] = proof.initial_claims[..]
  else {
    return Err(NovaError::InvalidNumInstances);
  };

  // Guard against a spurious `0/0` balance: `(t+a).num == 0` only certifies the
  // rational `t + a` is zero when its denominator `t.den · a.den` is nonzero.
  // Each root den is `Π_i (fingerprint_i + r)` with `r` a Fiat-Shamir challenge
  // drawn after the prover fixed its columns, so a zero factor happens with
  // probability ≤ N/|F| (negligible) for an honest prover and cannot be forced
  // by a malicious one — but the check is 4 comparisons and closes the path.
  let all_dens_nonzero = [
    row_table_root,
    row_access_root,
    col_table_root,
    col_access_root,
  ]
  .iter()
  .all(|f| f.den != E::Scalar::ZERO);
  if !all_dens_nonzero {
    return Err(NovaError::InvalidSumcheckProof);
  }

  let row_balanced = (row_table_root + row_access_root).num == E::Scalar::ZERO;
  let col_balanced = (col_table_root + col_access_root).num == E::Scalar::ZERO;
  if !row_balanced || !col_balanced {
    return Err(NovaError::InvalidSumcheckProof);
  }

  Ok(eval_point.to_vec())
}

#[cfg(test)]
mod tests {
  //! End-to-end host-verifier tests. Each builds four N-leaf sub-instances with
  //! a real (frozen) GKR prover, derives the [`MemCheckOpenings`] by evaluating
  //! the raw columns at the GKR `eval_point`, and checks `mem_check::verify`
  //! accepts a balanced witness and rejects a tampered one. The GKR prover is
  //! trusted here (it has its own round-trip tests); what is under test is the
  //! host reconcile + balance logic this module defines.
  use super::*;
  use crate::spartan::logup_gkr::layer::Layer;
  use crate::spartan::logup_gkr::prover;
  use crate::spartan::polys::multilinear::MultilinearPolynomial;
  use crate::traits::TranscriptEngineTrait;

  type E = crate::provider::Bn256EngineKZG;
  type Fr = <E as Engine>::Scalar;

  fn mle(v: Vec<Fr>) -> MultilinearPolynomial<Fr> {
    MultilinearPolynomial::new(v)
  }

  /// A hand-built N=4 memory-check witness whose row and col relations both
  /// balance. We choose the fingerprint pieces directly (γ, r and the per-cell
  /// columns) and derive `mem_row`/`mem_col`/dens so that `Σ ts/(T+r) =
  /// Σ 1/(W+r)` holds on each side.
  ///
  /// Construction: pick 4 distinct table dens `T[i]` freely; the access side is
  /// a multiset of reads into those cells with multiplicities `ts`, so the
  /// access dens are exactly the `T` values repeated per read. With N=4 reads
  /// over the 4 cells and `ts = [ts0..ts3]` summing to 4, the balance is
  /// `Σ ts[i]/(T[i]+r) = Σ_reads 1/(T[read]+r)` — identical multisets, so it
  /// holds by construction. Here we use `ts = [2,1,1,0]` and reads
  /// `[cell0, cell0, cell1, cell2]`.
  struct Witness {
    gamma: Fr,
    r: Fr,
    r_outer_full: Vec<Fr>,
    // raw N-length columns
    mem_row: Vec<Fr>,
    mem_col: Vec<Fr>,
    l_row: Vec<Fr>,
    l_col: Vec<Fr>,
    addr_row: Vec<Fr>,
    addr_col: Vec<Fr>,
    ts_row: Vec<Fr>,
    ts_col: Vec<Fr>,
  }

  impl Witness {
    // Build the four input layers [row_table, row_access, col_table, col_access].
    fn layers(&self) -> Vec<Layer<E>> {
      let n = self.mem_row.len();
      let id: Vec<Fr> = (0..n).map(|i| Fr::from(i as u64)).collect();
      let neg1 = -Fr::ONE;
      let den_table = |mem: &[Fr]| -> Vec<Fr> {
        (0..n)
          .map(|i| mem[i] * self.gamma + id[i] + self.r)
          .collect()
      };
      let den_access = |l: &[Fr], addr: &[Fr]| -> Vec<Fr> {
        (0..n)
          .map(|i| l[i] * self.gamma + addr[i] + self.r)
          .collect()
      };
      vec![
        Layer::<E> {
          num: mle(self.ts_row.clone()),
          den: mle(den_table(&self.mem_row)),
        },
        Layer::<E> {
          num: mle(vec![neg1; n]),
          den: mle(den_access(&self.l_row, &self.addr_row)),
        },
        Layer::<E> {
          num: mle(self.ts_col.clone()),
          den: mle(den_table(&self.mem_col)),
        },
        Layer::<E> {
          num: mle(vec![neg1; n]),
          den: mle(den_access(&self.l_col, &self.addr_col)),
        },
      ]
    }

    fn openings(&self, pt: &[Fr]) -> MemCheckOpenings<E> {
      let ev = |v: &[Fr]| mle(v.to_vec()).evaluate(pt);
      MemCheckOpenings {
        eval_L_row: ev(&self.l_row),
        eval_L_col: ev(&self.l_col),
        eval_row: ev(&self.addr_row),
        eval_col: ev(&self.addr_col),
        eval_ts_row: ev(&self.ts_row),
        eval_ts_col: ev(&self.ts_col),
        eval_mem_col: ev(&self.mem_col),
      }
    }
  }

  // A balanced N=4 witness. mem_row is eq(r_outer_full, ·) so the verifier can
  // recompute it; we set r_outer_full = [0,0] giving mem_row = [1,0,0,0].
  fn balanced_witness() -> Witness {
    let r_outer_full = vec![Fr::ZERO, Fr::ZERO];
    let mem_row = EqPolynomial::new(r_outer_full.clone()).evals(); // [1,0,0,0]
    let mem_col = vec![Fr::from(5), Fr::from(6), Fr::from(7), Fr::from(8)];
    // reads = [cell0, cell0, cell1, cell2]; ts = [2,1,1,0].
    let reads = [0usize, 0, 1, 2];
    let ts_row = vec![Fr::from(2), Fr::from(1), Fr::from(1), Fr::ZERO];
    let ts_col = ts_row.clone();
    let addr_row: Vec<Fr> = reads.iter().map(|&i| Fr::from(i as u64)).collect();
    let addr_col = addr_row.clone();
    // access lookup value = the table value at the read cell.
    let l_row: Vec<Fr> = reads.iter().map(|&i| mem_row[i]).collect();
    let l_col: Vec<Fr> = reads.iter().map(|&i| mem_col[i]).collect();
    Witness {
      gamma: Fr::from(3),
      r: Fr::from(9),
      r_outer_full,
      mem_row,
      mem_col,
      l_row,
      l_col,
      addr_row,
      addr_col,
      ts_row,
      ts_col,
    }
  }

  fn run(w: &Witness) -> Result<Vec<Fr>, NovaError> {
    let mut tr_p = <E as Engine>::TE::new(b"memcheck-test");
    let (proof, claim) = prover::prove::<E>(w.layers(), &mut tr_p).expect("prove");
    let openings = w.openings(claim.eval_point());
    let mut tr_v = <E as Engine>::TE::new(b"memcheck-test");
    verify::<E>(&proof, w.gamma, w.r, &w.r_outer_full, &openings, &mut tr_v)
  }

  #[test]
  fn accepts_balanced_witness() {
    let w = balanced_witness();
    let pt = run(&w).expect("host verifier must accept a balanced witness");
    assert_eq!(pt.len(), 2, "eval_point has log N = 2 variables");
  }

  #[test]
  fn rejects_tampered_ts() {
    // Break the row balance by bumping a multiplicity: Σ ts/(T+r) no longer
    // equals Σ 1/(W+r). The GKR proof is still built from the tampered layers,
    // so the balance check (not reconcile) is what fails.
    let mut w = balanced_witness();
    w.ts_row[0] += Fr::ONE;
    assert!(run(&w).is_err(), "must reject an unbalanced multiplicity");
  }

  #[test]
  fn rejects_mismatched_opening() {
    // Keep the layers balanced but feed the host a wrong opened column, so the
    // reconcile step (recomputed fraction vs GKR opening) fails.
    let w = balanced_witness();
    let mut tr_p = <E as Engine>::TE::new(b"memcheck-test");
    let (proof, claim) = prover::prove::<E>(w.layers(), &mut tr_p).expect("prove");
    let mut openings = w.openings(claim.eval_point());
    openings.eval_L_row += Fr::ONE; // inconsistent with the committed layer
    let mut tr_v = <E as Engine>::TE::new(b"memcheck-test");
    assert!(
      verify::<E>(&proof, w.gamma, w.r, &w.r_outer_full, &openings, &mut tr_v).is_err(),
      "must reject an opening that disagrees with the GKR reduction"
    );
  }
}
