//! Bridge layer wiring Logup-GKR to ppSNARK's memory-check.
//!
//! This module is the *host* half of the Logup-GKR memory-check — the link
//! between the pure, ppSNARK-agnostic `logup_gkr` argument (which owns no PCS)
//! and `ppsnark` (which owns the commitments, the inner sumcheck, and the
//! rerandomization of `L_row`/`L_col`). It lives outside `logup_gkr/` on
//! purpose: that module stays a standalone fractional-sum argument, while this
//! one depends on both sides.
//!
//! **Verifier-first.** The verifier half here wraps the frozen
//! `logup_gkr::verify` and closes soundness by (1) recomputing each logup
//! instance's input-layer fraction from the columns the prover opens at the GKR
//! evaluation point, (2) checking those match what the GKR reduced to, and (3)
//! running the fractional balance check. It is written **before** any prover:
//! the set of evaluations [`MemCheckOpenings`] names here *is* the contract a
//! prover must satisfy — the prover must open exactly these columns at exactly
//! the GKR eval_point, and no others.
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
use crate::spartan::logup_gkr::layer::Layer;
use crate::spartan::logup_gkr::proof::LogupGkrProof;
use crate::spartan::logup_gkr::verifier;
use crate::spartan::polys::eq::EqPolynomial;
use crate::spartan::polys::identity::IdentityPolynomial;
use crate::spartan::polys::multilinear::MultilinearPolynomial;
use crate::spartan::sumcheck::eq_sumcheck::EqSumCheckInstance;
use crate::spartan::sumcheck::{SumcheckEngine, SumcheckProof};
use crate::traits::Engine;
use ff::Field;
use rayon::prelude::*;

/// Fixed sub-instance count (`row_table, row_access, col_table, col_access`).
pub const NUM_SUB_INSTANCES: usize = 4;

/// Number of columns the rerandomize instance carries from the GKR `eval_point`
/// into the inner sumcheck. These are every column reconcile needs at the GKR
/// point that the verifier cannot self-compute (it computes only `mem_row = eq`
/// and the identity `id`): `L_row, L_col, addr_row, addr_col, ts_row, ts_col,
/// mem_col`. The order is fixed by [`MemCheckOpenings::rerand_claims`]
/// [`MemCheckOpenings::rerand_claims`] and must match between prover and
/// verifier.
pub const NUM_RERAND_COLUMNS: usize = 7;

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

impl<E: Engine> MemCheckOpenings<E> {
  /// The claimed column values at the GKR `eval_point`, in the fixed
  /// [`NUM_RERAND_COLUMNS`] order the rerandomize instance and the verifier both
  /// use: `[L_row, L_col, addr_row, addr_col, ts_row, ts_col, mem_col]`.
  pub fn rerand_claims(&self) -> [E::Scalar; NUM_RERAND_COLUMNS] {
    [
      self.eval_L_row,
      self.eval_L_col,
      self.eval_row,
      self.eval_col,
      self.eval_ts_row,
      self.eval_ts_col,
      self.eval_mem_col,
    ]
  }
}

/// The prover's memory-check witness: the raw length-N columns fed into
/// Logup-GKR.
///
/// These are exactly ppSNARK's padded memory-check columns (all length
/// `N = 2^{log N}`; see `jcbase/ppsnark-pad-to-N.md`). [`build_input_layers`]
/// turns them into the four GKR input layers, and [`prove`] consumes the
/// witness, opening the subset named by [`MemCheckOpenings`] at the shared
/// point. Field meanings match the four-sub-instance table in the module docs:
/// - `mem_row = eq(r_outer_full, ·)`, `mem_col = z` (table values);
/// - `L_row`/`L_col` the lookup columns, `addr_row = row`/`addr_col = col` the
///   access addresses, `ts_row`/`ts_col` the multiplicities.
#[derive(Clone)]
pub struct MemCheckWitness<E: Engine> {
  /// Row table values `mem_row = eq(r_outer_full, ·)`.
  pub mem_row: Vec<E::Scalar>,
  /// Col table values `mem_col = z`.
  pub mem_col: Vec<E::Scalar>,
  /// Row lookup column `L_row`.
  pub L_row: Vec<E::Scalar>,
  /// Col lookup column `L_col`.
  pub L_col: Vec<E::Scalar>,
  /// Row access addresses `addr_row = row`.
  pub addr_row: Vec<E::Scalar>,
  /// Col access addresses `addr_col = col`.
  pub addr_col: Vec<E::Scalar>,
  /// Row multiplicities `ts_row`.
  pub ts_row: Vec<E::Scalar>,
  /// Col multiplicities `ts_col`.
  pub ts_col: Vec<E::Scalar>,
}

/// Builds the four GKR input layers `[row_table, row_access, col_table,
/// col_access]` from the raw columns and the fingerprint `(gamma, r)`.
///
/// This is the prover-side dual of the verifier's per-instance reconstruction
/// in [`verify`]: it must produce, for every leaf `i`, exactly the fractions the
/// verifier recomputes at `eval_point`. The layers, in order (matching the
/// module-doc table and [`NUM_SUB_INSTANCES`]):
///
/// | idx | name        | num        | den                         |
/// |-----|-------------|------------|-----------------------------|
/// | 0   | row_table   | `ts_row`   | `mem_row·γ + id + r`        |
/// | 1   | row_access  | `-1`       | `L_row·γ + addr_row + r`    |
/// | 2   | col_table   | `ts_col`   | `mem_col·γ + id + r`        |
/// | 3   | col_access  | `-1`       | `L_col·γ + addr_col + r`    |
///
/// `id[i] = i` is the cell address. All columns must share one length `N`, a
/// power of two (`debug_assert`ed); the returned layers each have `log N`
/// variables, the shared GKR depth.
pub fn build_input_layers<E: Engine>(
  cols: &MemCheckWitness<E>,
  gamma: E::Scalar,
  r: E::Scalar,
) -> Vec<Layer<E>> {
  let n = cols.mem_row.len();
  debug_assert!(
    n.is_power_of_two() && n >= 2,
    "N must be a power of two >= 2"
  );
  for col in [
    &cols.mem_col,
    &cols.L_row,
    &cols.L_col,
    &cols.addr_row,
    &cols.addr_col,
    &cols.ts_row,
    &cols.ts_col,
  ] {
    debug_assert_eq!(col.len(), n, "all memory-check columns must share length N");
  }

  let neg_one = -E::Scalar::ONE;
  // Table-side den: mem·γ + id + r, where id[i] = i (the cell address).
  //
  // `id[i]` is built by per-chunk accumulation instead of `Scalar::from(i)` per
  // element: each chunk pays ONE `from` for its base index, then walks its cells
  // with a field `+ ONE`, which is far cheaper than a u64→Montgomery conversion.
  // The chunks run in parallel (`par_chunks_mut`), so this keeps full width
  // while dropping N `from`s to `N / chunk_size`.
  let one = E::Scalar::ONE;
  let chunk_size = 1 + n / rayon::current_num_threads().max(1);
  let den_table = |mem: &[E::Scalar]| -> Vec<E::Scalar> {
    let mut out = vec![E::Scalar::ZERO; n];
    out
      .par_chunks_mut(chunk_size)
      .enumerate()
      .for_each(|(c, chunk)| {
        let mut id = E::Scalar::from((c * chunk_size) as u64); // base index of this chunk
        for (out_i, mem_i) in chunk.iter_mut().zip(&mem[c * chunk_size..]) {
          *out_i = *mem_i * gamma + id + r;
          id += one;
        }
      });
    out
  };
  // Access-side den: L·γ + addr + r.
  let den_access = |l: &[E::Scalar], addr: &[E::Scalar]| -> Vec<E::Scalar> {
    (0..n)
      .into_par_iter()
      .map(|i| l[i] * gamma + addr[i] + r)
      .collect()
  };
  let mle = |v: Vec<E::Scalar>| MultilinearPolynomial::new(v);

  vec![
    // idx 0: row_table
    Layer::<E> {
      num: mle(cols.ts_row.clone()),
      den: mle(den_table(&cols.mem_row)),
    },
    // idx 1: row_access
    Layer::<E> {
      num: mle(vec![neg_one; n]),
      den: mle(den_access(&cols.L_row, &cols.addr_row)),
    },
    // idx 2: col_table
    Layer::<E> {
      num: mle(cols.ts_col.clone()),
      den: mle(den_table(&cols.mem_col)),
    },
    // idx 3: col_access
    Layer::<E> {
      num: mle(vec![neg_one; n]),
      den: mle(den_access(&cols.L_col, &cols.addr_col)),
    },
  ]
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

/// Prover output for the Logup-GKR memory-check, before ppSNARK integration.
///
/// Bundles the three things the three downstream consumers need:
/// - `proof`: the GKR proof, absorbed into the SNARK and replayed by [`verify`];
/// - `openings`: the column evaluations at the GKR `eval_point`, the
///   [`MemCheckOpenings`] the host reconcile step checks. In a standalone run
///   these are opened directly at `eval_point`; once wired into ppSNARK the same
///   values arrive at the shared inner point via `rerandomize` instead;
/// - `rerandomize`: the [`RerandomizeSumcheckInstance`] that moves `L_row`/
///   `L_col` from `eval_point` into the inner sumcheck bundle (unused by the
///   standalone verifier, produced here so the full plumbing is exercised).
pub struct MemCheckProverOutput<E: Engine> {
  /// The GKR fractional-sum proof.
  pub proof: LogupGkrProof<E>,
  /// Column evaluations at the GKR `eval_point` (host reconcile input).
  pub openings: MemCheckOpenings<E>,
  /// L_row/L_col opening-point reduction into the inner sumcheck.
  pub rerandomize: RerandomizeSumcheckInstance<E>,
  /// The shared GKR evaluation point (length `log N`).
  pub eval_point: Vec<E::Scalar>,
}

/// Proves the ppSNARK memory-check via Logup-GKR from the raw columns.
///
/// This is the prover-side entry point mirroring [`verify`]. It:
/// 1. builds the four GKR input layers ([`build_input_layers`]);
/// 2. runs the frozen GKR prover to fold them and emit the proof plus the shared
///    `eval_point`;
/// 3. evaluates the opened columns at `eval_point` to form [`MemCheckOpenings`];
/// 4. builds the [`RerandomizeSumcheckInstance`] that will later carry
///    `L_row`/`L_col` into the inner sumcheck.
///
/// The transcript must be in the same state the verifier expects at the GKR
/// slot (the GKR prover absorbs exactly what [`verify`] replays). `(gamma, r)`
/// are ppSNARK's memory-check fingerprint challenges.
pub fn prove<E: Engine>(
  witness: MemCheckWitness<E>,
  gamma: E::Scalar,
  r: E::Scalar,
  transcript: &mut E::TE,
) -> Result<MemCheckProverOutput<E>, NovaError> {
  // (1)+(2) Build the four input layers and fold them through the GKR trees.
  let layers = build_input_layers(&witness, gamma, r);
  let (proof, claim) = crate::spartan::logup_gkr::prover::prove::<E>(layers, transcript)?;
  let eval_point = claim.eval_point().to_vec();

  // (3) Assemble the opened columns at the shared point. The prover is honest
  // and owns every column, so each opening is taken by the cheapest route that
  // yields the same value:
  // - `ts_row`/`ts_col` ARE the table-side numerators the GKR reduction already
  //   produced (openings 0, 2), so reuse them directly;
  // - `addr_row`/`addr_col`/`mem_col` are fused into the GKR dens
  //   (`L·γ + addr + r`, `mem·γ + id + r`), so invert those closed forms instead
  //   of re-evaluating the MLEs — one field op vs an N-wide evaluation;
  // - `L_row`/`L_col` must be evaluated directly (they seed the rerandomize
  //   claims and are the other unknown in the access dens), so they stay `ev`.
  // This is a pure prover-side shortcut: soundness lives in the verifier, which
  // opens each column against its own commitment (see `verify`).
  let ev = |v: &[E::Scalar]| MultilinearPolynomial::evaluate_with(v, &eval_point);
  let [row_table, row_access, col_table, col_access] = claim.openings()[..] else {
    return Err(NovaError::InvalidNumInstances);
  };
  let eval_L_row = ev(&witness.L_row);
  let eval_L_col = ev(&witness.L_col);
  let eval_id = IdentityPolynomial::<E::Scalar>::new(eval_point.len()).evaluate(&eval_point);
  let gamma_inv = gamma.invert().expect("fingerprint gamma is nonzero");
  let openings = MemCheckOpenings {
    eval_L_row,
    eval_L_col,
    // row_access.den = L_row·γ + addr_row + r  ⇒  addr_row = den − L_row·γ − r
    eval_row: row_access.den - eval_L_row * gamma - r,
    // col_access.den = L_col·γ + addr_col + r  ⇒  addr_col = den − L_col·γ − r
    eval_col: col_access.den - eval_L_col * gamma - r,
    eval_ts_row: row_table.num, // row_table numerator
    eval_ts_col: col_table.num, // col_table numerator
    // col_table.den = mem_col·γ + id + r  ⇒  mem_col = (den − id − r)·γ⁻¹
    eval_mem_col: (col_table.den - eval_id - r) * gamma_inv,
  };

  // (4) Rerandomize instance: carry every column reconcile needs from the GKR
  // eval_point into the inner sumcheck. Columns and claims share the fixed
  // RERAND order [L_row, L_col, addr_row, addr_col, ts_row, ts_col, mem_col].
  // The witness is consumed here, so its columns are moved in (no clone).
  let claims = openings.rerand_claims().to_vec();
  let columns = vec![
    witness.L_row,
    witness.L_col,
    witness.addr_row,
    witness.addr_col,
    witness.ts_row,
    witness.ts_col,
    witness.mem_col,
  ];
  let rerandomize = RerandomizeSumcheckInstance::new(eval_point.clone(), columns, claims);

  Ok(MemCheckProverOutput {
    proof,
    openings,
    rerandomize,
    eval_point,
  })
}

/// Rerandomizes the GKR verifier's per-column evaluation requests at the GKR
/// `eval_point` into the inner sumcheck, so every column reconcile needs is
/// opened at the shared inner point `r_inner_batched` instead of at
/// `eval_point`.
///
/// # Why several columns, not just L
/// The host reconcile ([`verify`]) checks the four GKR-reduced fractions at
/// `eval_point`. Each fraction's num/den is a fingerprint of several columns
/// (`ts`, `L`, `addr`, `mem_col`); the verifier can self-compute only `mem_row =
/// eq(r_outer_full, ·)` and the identity `id`. Every other column it needs at
/// `eval_point` must be carried there. So this instance rerandomizes **all** of
/// them (order fixed by [`MemCheckOpenings::rerand_claims`]): each column `X` becomes one
/// sumcheck `Σ_y eq(eval_point, y) · X(y)` over the same `y ∈ {0,1}^{log N}`
/// domain as the inner ABC/E sumcheck, folding into the shared `prove_helper`
/// bundle and landing at `r_inner_batched`. This is exactly ppSNARK's E-claim
/// mechanism (`Σ eq(r_outer, y) · E(y)`) applied to each column, all sharing one
/// `eq_sumcheck` because they use the same `eval_point`.
///
/// # Degree
/// Each summand `eq(eval_point, ·) · X(·)` is a product of two multilinears, so
/// the true round-polynomial degree is **2**. [`prove_helper`] hardcodes degree
/// 3 (it interpolates every instance with `from_evals_deg3` and asserts all
/// bundled instances share a degree), so [`Self::degree`] reports 3 and the
/// quadratic evals carry a zero cubic coefficient — identical to how the E-claim
/// rides in the degree-3 inner instance. Running outside `prove_helper` could
/// reclaim the extra sample point (see HANDOFF O1).
///
/// [`prove_helper`]: super::ppsnark
pub struct RerandomizeSumcheckInstance<E: Engine> {
  /// Transparent `eq(eval_point, ·)` factor, shared by all columns.
  eq_sumcheck: EqSumCheckInstance<E>,
  /// The columns being rerandomized, order [`MemCheckOpenings::rerand_claims`].
  polys: Vec<MultilinearPolynomial<E::Scalar>>,
  /// Running claim per column (BDDT, eprint 2025/1117 §6.2).
  running_claims: Vec<E::Scalar>,
  /// Saved `[p(0), 0, p(-1)]` per column, used by [`Self::bound`].
  saved_evals: Vec<[E::Scalar; 3]>,
}

impl<E: Engine> RerandomizeSumcheckInstance<E> {
  /// Builds the instance from the GKR `eval_point`, the columns (order
  /// [`MemCheckOpenings::rerand_claims`]), and their claimed values `X(eval_point)` (the GKR
  /// verifier's requested values, which seed the running claims and are the
  /// instance's initial sumcheck claims). Every column must have length
  /// `N = 2^{eval_point.len()}`.
  pub fn new(
    eval_point: Vec<E::Scalar>,
    columns: Vec<Vec<E::Scalar>>,
    claims: Vec<E::Scalar>,
  ) -> Self {
    assert_eq!(columns.len(), claims.len());
    let saved_evals = vec![[E::Scalar::ZERO; 3]; columns.len()];
    Self {
      eq_sumcheck: EqSumCheckInstance::new(eval_point),
      polys: columns
        .into_iter()
        .map(MultilinearPolynomial::new)
        .collect(),
      running_claims: claims,
      saved_evals,
    }
  }
}

impl<E: Engine> SumcheckEngine<E> for RerandomizeSumcheckInstance<E> {
  fn initial_claims(&self) -> Vec<E::Scalar> {
    self.running_claims.clone()
  }

  fn degree(&self) -> usize {
    // True degree is 2 (eq · X); reported as 3 to ride in the degree-3
    // prove_helper bundle. See the type docs and HANDOFF O1.
    3
  }

  fn size(&self) -> usize {
    let n = self.polys[0].len();
    debug_assert!(self.polys.iter().all(|p| p.len() == n));
    n
  }

  fn evaluation_points(&mut self) -> Vec<Vec<E::Scalar>> {
    // Each column is one quadratic `eq(eval_point, ·) · X(·)`, sampled the same
    // way as the E-claim. The cubic coefficient is zero (degree 2).
    let evals: Vec<[E::Scalar; 3]> = self
      .polys
      .par_iter()
      .zip(self.running_claims.par_iter())
      .map(|(poly, &claim)| {
        let (e0, _, einf) = self
          .eq_sumcheck
          .evaluation_points_quadratic_with_one_input(poly, claim);
        [e0, E::Scalar::ZERO, einf]
      })
      .collect();

    self.saved_evals = evals.clone();
    evals.into_iter().map(|e| e.to_vec()).collect()
  }

  fn bound(&mut self, r: &E::Scalar) {
    self.running_claims = self
      .running_claims
      .iter()
      .zip(self.saved_evals.iter())
      .map(|(&claim, saved)| SumcheckProof::<E>::update_claim(claim, saved, r))
      .collect();

    self
      .polys
      .par_iter_mut()
      .for_each(|poly| poly.bind_poly_var_top(r));

    self.eq_sumcheck.bound(r);
  }

  fn final_claims(&self) -> Vec<Vec<E::Scalar>> {
    self.polys.iter().map(|p| vec![p[0]]).collect()
  }
}

#[cfg(test)]
mod tests {
  //! End-to-end host-verifier tests. Each builds four N-leaf sub-instances with
  //! a real (frozen) GKR prover, derives the [`MemCheckOpenings`] by evaluating
  //! the raw columns at the GKR `eval_point`, and checks `mem_check::verify`
  //! accepts a balanced witness and rejects a tampered one. The GKR prover is
  //! trusted here (it has its own round-trip tests); what is under test is the
  //! End-to-end tests through the top-level [`prove`]/[`verify`] pair. Each
  //! builds a balanced N=4 witness, proves it (four sub-instances folded by the
  //! frozen GKR prover), and checks the host verifier accepts it and rejects
  //! tampered multiplicities or mismatched openings. What is under test is this
  //! module's own logic — `build_input_layers`, reconcile, balance, and the
  //! rerandomize claims — with the GKR prover/verifier trusted (own tests).
  use super::*;
  use crate::traits::TranscriptEngineTrait;

  type E = crate::provider::Bn256EngineKZG;
  type Fr = <E as Engine>::Scalar;

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
    cols: MemCheckWitness<E>,
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
    let L_row: Vec<Fr> = reads.iter().map(|&i| mem_row[i]).collect();
    let L_col: Vec<Fr> = reads.iter().map(|&i| mem_col[i]).collect();
    Witness {
      gamma: Fr::from(3),
      r: Fr::from(9),
      r_outer_full,
      cols: MemCheckWitness {
        mem_row,
        mem_col,
        L_row,
        L_col,
        addr_row,
        addr_col,
        ts_row,
        ts_col,
      },
    }
  }

  fn run(w: &Witness) -> Result<Vec<Fr>, NovaError> {
    let mut tr_p = <E as Engine>::TE::new(b"memcheck-test");
    let out = prove::<E>(w.cols.clone(), w.gamma, w.r, &mut tr_p).expect("prove");
    let mut tr_v = <E as Engine>::TE::new(b"memcheck-test");
    verify::<E>(
      &out.proof,
      w.gamma,
      w.r,
      &w.r_outer_full,
      &out.openings,
      &mut tr_v,
    )
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
    w.cols.ts_row[0] += Fr::ONE;
    assert!(run(&w).is_err(), "must reject an unbalanced multiplicity");
  }

  #[test]
  fn rejects_mismatched_opening() {
    // Keep the layers balanced but feed the host a wrong opened column, so the
    // reconcile step (recomputed fraction vs GKR opening) fails.
    let w = balanced_witness();
    let mut tr_p = <E as Engine>::TE::new(b"memcheck-test");
    let mut out = prove::<E>(w.cols.clone(), w.gamma, w.r, &mut tr_p).expect("prove");
    out.openings.eval_L_row += Fr::ONE; // inconsistent with the committed layer
    let mut tr_v = <E as Engine>::TE::new(b"memcheck-test");
    assert!(
      verify::<E>(
        &out.proof,
        w.gamma,
        w.r,
        &w.r_outer_full,
        &out.openings,
        &mut tr_v
      )
      .is_err(),
      "must reject an opening that disagrees with the GKR reduction"
    );
  }

  #[test]
  fn rerandomize_claims_match_openings() {
    // The rerandomize instance's initial claims must be exactly the claimed
    // column values at eval_point, in the fixed RERAND order.
    let w = balanced_witness();
    let mut tr_p = <E as Engine>::TE::new(b"memcheck-test");
    let out = prove::<E>(w.cols.clone(), w.gamma, w.r, &mut tr_p).expect("prove");
    let claims = out.rerandomize.initial_claims();
    assert_eq!(claims, out.openings.rerand_claims().to_vec());
    assert_eq!(claims.len(), NUM_RERAND_COLUMNS);
  }
}
