//! This module implements `RelaxedR1CSSNARK` traits using a spark-based approach to prove evaluations of
//! sparse multilinear polynomials involved in Spartan's sum-check protocol, thereby providing a preprocessing SNARK
//! The verifier in this preprocessing SNARK maintains a commitment to R1CS matrices. This is beneficial when using a
//! polynomial commitment scheme in which the verifier's costs is succinct.
//! The SNARK implemented here is described in the MicroNova paper.
use crate::traits::evm_serde::EvmCompatSerde;
use crate::{
  digest::{DigestComputer, SimpleDigestible},
  errors::NovaError,
  r1cs::{R1CSShape, RelaxedR1CSInstance, RelaxedR1CSWitness},
  spartan::{
    math::Math,
    mem_check_logup,
    mem_check_logup_gkr::{self, MemCheckWitness},
    polys::{
      eq::EqPolynomial,
      masked_eq::MaskedEqPolynomial,
      multilinear::{MultilinearPolynomial, SparsePolynomial},
      univariate::{CompressedUniPoly, UniPoly},
    },
    powers,
    sumcheck::{eq_sumcheck::EqSumCheckInstance, SumcheckEngine, SumcheckProof},
    PolyEvalInstance, PolyEvalWitness,
  },
  traits::{
    commitment::{CommitmentEngineTrait, Len},
    evaluation::EvaluationEngineTrait,
    snark::{DigestHelperTrait, RelaxedR1CSSNARKTrait},
    Engine, TranscriptEngineTrait, TranscriptReprTrait,
  },
  zip_with, Commitment, CommitmentKey,
};
use core::cmp::max;
use ff::Field;
use itertools::Itertools as _;
use once_cell::sync::OnceCell;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use serde_with::serde_as;

fn padded<E: Engine>(v: &[E::Scalar], n: usize, e: &E::Scalar) -> Vec<E::Scalar> {
  let mut v_padded = vec![*e; n];
  for (i, v_i) in v.iter().enumerate() {
    v_padded[i] = *v_i;
  }
  v_padded
}

/// A type that holds `R1CSShape` in a form amenable to memory checking
#[derive(Clone, Serialize, Deserialize)]
#[serde(bound = "")]
pub struct R1CSShapeSparkRepr<E: Engine> {
  /// size of the vectors
  pub N: usize,

  /// dense representation of row indices
  pub row: Vec<E::Scalar>,
  /// dense representation of column indices
  pub col: Vec<E::Scalar>,
  /// values from A matrix
  pub val_A: Vec<E::Scalar>,
  /// values from B matrix
  pub val_B: Vec<E::Scalar>,
  /// values from C matrix
  pub val_C: Vec<E::Scalar>,

  /// timestamp polynomial for rows
  pub ts_row: Vec<E::Scalar>,
  /// timestamp polynomial for columns
  pub ts_col: Vec<E::Scalar>,
}

/// A type that holds a commitment to a sparse polynomial
#[derive(Clone, Serialize, Deserialize)]
#[serde(bound = "")]
pub struct R1CSShapeSparkCommitment<E: Engine> {
  /// size of each vector
  pub N: usize,

  // commitments to the dense representation
  /// commitment to row indices
  pub comm_row: Commitment<E>,
  /// commitment to column indices
  pub comm_col: Commitment<E>,
  /// commitment to A matrix values
  pub comm_val_A: Commitment<E>,
  /// commitment to B matrix values
  pub comm_val_B: Commitment<E>,
  /// commitment to C matrix values
  pub comm_val_C: Commitment<E>,

  // commitments to the timestamp polynomials
  /// commitment to row timestamps
  pub comm_ts_row: Commitment<E>,
  /// commitment to column timestamps
  pub comm_ts_col: Commitment<E>,
}

impl<E: Engine> TranscriptReprTrait<E::GE> for R1CSShapeSparkCommitment<E> {
  fn to_transcript_bytes(&self) -> Vec<u8> {
    [
      self.comm_row,
      self.comm_col,
      self.comm_val_A,
      self.comm_val_B,
      self.comm_val_C,
      self.comm_ts_row,
      self.comm_ts_col,
    ]
    .as_slice()
    .to_transcript_bytes()
  }
}

impl<E: Engine> R1CSShapeSparkRepr<E> {
  /// represents `R1CSShape` in a Spark-friendly format amenable to memory checking
  pub fn new(S: &R1CSShape<E>) -> R1CSShapeSparkRepr<E> {
    let N = {
      let total_nz = S.A.len() + S.B.len() + S.C.len();
      max(total_nz, max(2 * S.num_vars, S.num_cons)).next_power_of_two()
    };

    let (mut row, mut col) = (vec![0; N], vec![N - 1; N]); // we make col lookup into the last entry of z, so we commit to zeros

    for (i, (r, c, _)) in S.A.iter().chain(S.B.iter()).chain(S.C.iter()).enumerate() {
      row[i] = r;
      col[i] = c;
    }
    let val_A = {
      let mut val = vec![E::Scalar::ZERO; N];
      for (i, (_, _, v)) in S.A.iter().enumerate() {
        val[i] = v;
      }
      val
    };

    let val_B = {
      let mut val = vec![E::Scalar::ZERO; N];
      for (i, (_, _, v)) in S.B.iter().enumerate() {
        val[S.A.len() + i] = v;
      }
      val
    };

    let val_C = {
      let mut val = vec![E::Scalar::ZERO; N];
      for (i, (_, _, v)) in S.C.iter().enumerate() {
        val[S.A.len() + S.B.len() + i] = v;
      }
      val
    };

    // timestamp calculation routine
    let timestamp_calc = |num_ops: usize, num_cells: usize, addr_trace: &[usize]| -> Vec<usize> {
      let mut ts = vec![0usize; num_cells];

      assert!(num_ops >= addr_trace.len());
      for addr in addr_trace {
        assert!(*addr < num_cells);
        ts[*addr] += 1;
      }
      ts
    };

    // timestamp polynomials for row
    let (ts_row, ts_col) =
      rayon::join(|| timestamp_calc(N, N, &row), || timestamp_calc(N, N, &col));

    // a routine to turn a vector of usize into a vector scalars
    let to_vec_scalar = |v: &[usize]| -> Vec<E::Scalar> {
      (0..v.len())
        .map(|i| E::Scalar::from(v[i] as u64))
        .collect::<Vec<E::Scalar>>()
    };

    R1CSShapeSparkRepr {
      N,

      // dense representation
      row: to_vec_scalar(&row),
      col: to_vec_scalar(&col),
      val_A,
      val_B,
      val_C,

      // timestamp polynomials
      ts_row: to_vec_scalar(&ts_row),
      ts_col: to_vec_scalar(&ts_col),
    }
  }

  /// Commits to the R1CSShapeSparkRepr using the provided commitment key
  pub fn commit(&self, ck: &CommitmentKey<E>) -> R1CSShapeSparkCommitment<E> {
    let comm_vec: Vec<Commitment<E>> = [
      &self.row,
      &self.col,
      &self.val_A,
      &self.val_B,
      &self.val_C,
      &self.ts_row,
      &self.ts_col,
    ]
    .par_iter()
    .map(|v| E::CE::commit(ck, v, &E::Scalar::ZERO))
    .collect();

    R1CSShapeSparkCommitment {
      N: self.row.len(),
      comm_row: comm_vec[0],
      comm_col: comm_vec[1],
      comm_val_A: comm_vec[2],
      comm_val_B: comm_vec[3],
      comm_val_C: comm_vec[4],
      comm_ts_row: comm_vec[5],
      comm_ts_col: comm_vec[6],
    }
  }

  // computes evaluation oracles
  fn evaluation_oracles(
    &self,
    S: &R1CSShape<E>,
    r_outer: &[E::Scalar],
    z: &[E::Scalar],
  ) -> (
    Vec<E::Scalar>,
    Vec<E::Scalar>,
    Vec<E::Scalar>,
    Vec<E::Scalar>,
  ) {
    let mem_row = EqPolynomial::new(r_outer.to_vec()).evals();
    let mem_col = padded::<E>(z, self.N, &E::Scalar::ZERO);

    let (L_row, L_col) = {
      let mut L_row = vec![mem_row[0]; self.N]; // we place mem_row[0] since resized row is appended with 0s
      let mut L_col = vec![mem_col[self.N - 1]; self.N]; // we place mem_col[N-1] since resized col is appended with N-1

      for (i, (val_r, val_c)) in S
        .A
        .iter()
        .chain(S.B.iter())
        .chain(S.C.iter())
        .map(|(r, c, _)| (mem_row[r], mem_col[c]))
        .enumerate()
      {
        L_row[i] = val_r;
        L_col[i] = val_c;
      }
      (L_row, L_col)
    };

    (mem_row, mem_col, L_row, L_col)
  }
}

/// The [WitnessBoundSumcheck] ensures that the witness polynomial W defined over n = log(N) variables,
/// is zero outside of the first `num_vars = 2^m` entries.
///
/// # Details
///
/// The `W` polynomial is padded with zeros to size N = 2^n.
/// The `masked_eq` polynomials is defined as with regards to a random challenge `tau` as
/// the eq(tau) polynomial, where the first 2^m evaluations to 0.
///
/// The instance is given by
///  `0 = ∑_{0≤i<2^n} masked_eq[i] * W[i]`.
/// It is equivalent to the expression
///  `0 = ∑_{2^m≤i<2^n} eq[i] * W[i]`
/// Since `eq` is random, the instance is only satisfied if `W[2^{m}..] = 0`.
pub struct WitnessBoundSumcheck<E: Engine> {
  poly_W: MultilinearPolynomial<E::Scalar>,
  poly_masked_eq: MultilinearPolynomial<E::Scalar>,
}

impl<E: Engine> WitnessBoundSumcheck<E> {
  /// Creates a new WitnessBoundSumcheck instance
  pub fn new(tau: Vec<E::Scalar>, poly_W_padded: Vec<E::Scalar>, num_vars: usize) -> Self {
    let num_vars_log = num_vars.log_2();
    // When num_vars = num_rounds, we shouldn't have to prove anything
    // but we still want this instance to compute the evaluation of W
    let num_rounds = poly_W_padded.len().log_2();
    assert!(num_vars_log < num_rounds);

    let poly_masked_eq_evals =
      MaskedEqPolynomial::new(&EqPolynomial::new(tau), num_vars_log).evals();

    Self {
      poly_W: MultilinearPolynomial::new(poly_W_padded),
      poly_masked_eq: MultilinearPolynomial::new(poly_masked_eq_evals),
    }
  }
}
impl<E: Engine> SumcheckEngine<E> for WitnessBoundSumcheck<E> {
  fn initial_claims(&self) -> Vec<E::Scalar> {
    vec![E::Scalar::ZERO]
  }

  fn degree(&self) -> usize {
    3
  }

  fn size(&self) -> usize {
    assert_eq!(self.poly_W.len(), self.poly_masked_eq.len());
    self.poly_W.len()
  }

  fn evaluation_points(&mut self) -> Vec<Vec<E::Scalar>> {
    // masked_eq * W is a quadratic polynomial (A * B)
    let (eval_point_0, eval_point_inf) =
      SumcheckProof::<E>::compute_eval_points_quadratic(&self.poly_masked_eq, &self.poly_W);

    // bound_coeff is always 0 for quadratic polynomials
    vec![vec![eval_point_0, E::Scalar::ZERO, eval_point_inf]]
  }

  fn bound(&mut self, r: &E::Scalar) {
    [&mut self.poly_W, &mut self.poly_masked_eq]
      .par_iter_mut()
      .for_each(|poly| poly.bind_poly_var_top(r));
  }

  fn final_claims(&self) -> Vec<Vec<E::Scalar>> {
    vec![vec![self.poly_W[0], self.poly_masked_eq[0]]]
  }
}

/// Inner batched sumcheck instance for PPSNARK
///
/// Proves two claims:
/// 1. ABC claim: `claim = Σ L_row(y) * L_col(y) * val(y)` where val = val_A + c*val_B + c²*val_C
/// 2. E claim: `claim_E = Σ eq(r_outer, y) * E(y)` to rerandomize the E opening point
pub struct InnerBatchedSumcheckInstance<E: Engine> {
  /// The claim value for ABC: eval_Az + c * eval_Bz + c² * eval_Cz
  pub claim: E::Scalar,
  /// Row lookup polynomial
  pub poly_L_row: MultilinearPolynomial<E::Scalar>,
  /// Column lookup polynomial
  pub poly_L_col: MultilinearPolynomial<E::Scalar>,
  /// Value polynomial (val_A + c*val_B + c²*val_C)
  pub poly_val: MultilinearPolynomial<E::Scalar>,

  /// The claim value for E: eval_E_at_r_outer
  pub claim_E: E::Scalar,
  /// Factored eq polynomial eq(r_outer, .) (Gruen, eprint 2024/108)
  eq_sumcheck: EqSumCheckInstance<E>,
  /// Error polynomial E
  pub poly_E: MultilinearPolynomial<E::Scalar>,

  /// Running claim for E sub-claim (BDDT, eprint 2025/1117 Section 6.2)
  running_claim_E: E::Scalar,
  /// Saved eval points for E sub-claim
  saved_evals_E: [E::Scalar; 3],
}

impl<E: Engine> InnerBatchedSumcheckInstance<E> {
  /// Create a new inner batched sumcheck instance
  pub fn new(
    claim: E::Scalar,
    L_row: Vec<E::Scalar>,
    L_col: Vec<E::Scalar>,
    val: Vec<E::Scalar>,
    claim_E: E::Scalar,
    r_outer: Vec<E::Scalar>,
    E: Vec<E::Scalar>,
  ) -> Self {
    Self {
      claim,
      poly_L_row: MultilinearPolynomial::new(L_row),
      poly_L_col: MultilinearPolynomial::new(L_col),
      poly_val: MultilinearPolynomial::new(val),
      claim_E,
      eq_sumcheck: EqSumCheckInstance::new(r_outer),
      poly_E: MultilinearPolynomial::new(E),
      running_claim_E: claim_E,
      saved_evals_E: [E::Scalar::ZERO; 3],
    }
  }
}

impl<E: Engine> SumcheckEngine<E> for InnerBatchedSumcheckInstance<E> {
  fn initial_claims(&self) -> Vec<E::Scalar> {
    vec![self.claim, self.claim_E]
  }

  fn degree(&self) -> usize {
    3
  }

  fn size(&self) -> usize {
    assert_eq!(self.poly_L_row.len(), self.poly_val.len());
    assert_eq!(self.poly_L_row.len(), self.poly_L_col.len());
    assert_eq!(self.poly_L_row.len(), self.poly_E.len());
    self.poly_L_row.len()
  }

  fn evaluation_points(&mut self) -> Vec<Vec<E::Scalar>> {
    let (poly_A, poly_B, poly_C) = (&self.poly_L_row, &self.poly_L_col, &self.poly_val);

    // E claim: 1 N-scaling sum + claim-derived (BDDT, eprint 2025/1117 Section 6.2)
    let ((eval_point_0, bound_coeff, eval_point_inf), (eval_E_0, eval_E_bound_coeff, eval_E_inf)) =
      rayon::join(
        || SumcheckProof::<E>::compute_eval_points_cubic(poly_A, poly_B, poly_C),
        || {
          self
            .eq_sumcheck
            .evaluation_points_quadratic_with_one_input(&self.poly_E, self.running_claim_E)
        },
      );

    self.saved_evals_E = [eval_E_0, eval_E_bound_coeff, eval_E_inf];

    vec![
      vec![eval_point_0, bound_coeff, eval_point_inf],
      // E claim is degree-2 (quadratic), so cubic coefficient is zero
      vec![eval_E_0, E::Scalar::ZERO, eval_E_inf],
    ]
  }

  fn bound(&mut self, r: &E::Scalar) {
    self.running_claim_E =
      SumcheckProof::<E>::update_claim(self.running_claim_E, &self.saved_evals_E, r);

    [
      &mut self.poly_L_row,
      &mut self.poly_L_col,
      &mut self.poly_val,
      &mut self.poly_E,
    ]
    .par_iter_mut()
    .for_each(|poly| poly.bind_poly_var_top(r));

    self.eq_sumcheck.bound(r);
  }

  fn final_claims(&self) -> Vec<Vec<E::Scalar>> {
    vec![
      vec![self.poly_L_row[0], self.poly_L_col[0]],
      vec![self.poly_E[0]],
    ]
  }
}

/// A type that represents the prover's key
#[derive(Clone, Serialize, Deserialize)]
#[serde(bound = "")]
pub struct ProverKey<E: Engine, EE: EvaluationEngineTrait<E>> {
  pk_ee: EE::ProverKey,
  S_repr: R1CSShapeSparkRepr<E>,
  S_comm: R1CSShapeSparkCommitment<E>,
  vk_digest: E::Scalar, // digest of verifier's key
}

/// A type that represents the verifier's key
#[derive(Clone, Serialize, Deserialize)]
#[serde(bound = "")]
pub struct VerifierKey<E: Engine, EE: EvaluationEngineTrait<E>> {
  num_cons: usize,
  num_vars: usize,
  vk_ee: EE::VerifierKey,
  S_comm: R1CSShapeSparkCommitment<E>,
  #[serde(skip, default = "OnceCell::new")]
  digest: OnceCell<E::Scalar>,
}

impl<E: Engine, EE: EvaluationEngineTrait<E>> SimpleDigestible for VerifierKey<E, EE> {}

/// A succinct proof of knowledge of a witness to a relaxed R1CS instance
/// The proof is produced using Spartan's combination of the sum-check and
/// the commitment to a vector viewed as a polynomial commitment
#[serde_as]
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(bound = "")]
pub struct RelaxedR1CSSNARK<E: Engine, EE: EvaluationEngineTrait<E>> {
  // commitment to oracles for memory reads
  comm_L_row: Commitment<E>,
  comm_L_col: Commitment<E>,

  // Memory-check proof data. Bundled per implementation so the SNARK gates the
  // whole memory-check behind one field: inverse-logup (feature `logup-no-gkr`)
  // carries the four inverse-oracle commitments + evals; Logup-GKR (default)
  // carries the fractional-sum proof + the rerandomize claims.
  mem_check: mem_check_logup::LogupProofData<E>,
  mem_check_gkr: mem_check_logup_gkr::GkrProofData<E>,

  // outer sum-check proof
  sc_outer: SumcheckProof<E>,

  // claims about Az, Bz, and Cz virtual polynomials at the end of sum-check
  #[serde_as(as = "EvmCompatSerde")]
  eval_Az_at_r_outer: E::Scalar,
  #[serde_as(as = "EvmCompatSerde")]
  eval_Bz_at_r_outer: E::Scalar,
  #[serde_as(as = "EvmCompatSerde")]
  eval_Cz_at_r_outer: E::Scalar,
  #[serde_as(as = "EvmCompatSerde")]
  eval_E_at_r_outer: E::Scalar,

  // inner batched sum-check proof
  sc_inner_batched: SumcheckProof<E>,

  // claims from the end of sum-check
  #[serde_as(as = "EvmCompatSerde")]
  eval_E: E::Scalar,
  #[serde_as(as = "EvmCompatSerde")]
  eval_L_row: E::Scalar,
  #[serde_as(as = "EvmCompatSerde")]
  eval_L_col: E::Scalar,
  #[serde_as(as = "EvmCompatSerde")]
  eval_val_A: E::Scalar,
  #[serde_as(as = "EvmCompatSerde")]
  eval_val_B: E::Scalar,
  #[serde_as(as = "EvmCompatSerde")]
  eval_val_C: E::Scalar,

  #[serde_as(as = "EvmCompatSerde")]
  eval_W: E::Scalar,

  // Shared address/multiplicity evaluations (both memory-check implementations
  // open these at the inner point).
  #[serde_as(as = "EvmCompatSerde")]
  eval_row: E::Scalar, // address
  #[serde_as(as = "EvmCompatSerde")]
  eval_ts_row: E::Scalar,
  #[serde_as(as = "EvmCompatSerde")]
  eval_col: E::Scalar, // address
  #[serde_as(as = "EvmCompatSerde")]
  eval_ts_col: E::Scalar,

  // a PCS evaluation argument
  eval_arg: EE::EvaluationArgument,
}

impl<E: Engine, EE: EvaluationEngineTrait<E>> RelaxedR1CSSNARK<E, EE> {
  /// Batched inner sum-check prover for 4 instances: memory, inner_batched,
  /// witness, and the Logup-GKR rerandomize instance (L_row/L_col opening-point
  /// reduction). All four share the same size and (reported) degree, so they
  /// fold into one RLC'd sumcheck landing at the shared point `r_inner_batched`.
  fn prove_helper<T1, T2, T3, T4>(
    mem: &mut T1,
    inner: &mut T2,
    witness: &mut T3,
    rerand: &mut T4,
    transcript: &mut E::TE,
  ) -> Result<
    (
      SumcheckProof<E>,
      Vec<E::Scalar>,
      Vec<Vec<E::Scalar>>,
      Vec<Vec<E::Scalar>>,
      Vec<Vec<E::Scalar>>,
      Vec<Vec<E::Scalar>>,
    ),
    NovaError,
  >
  where
    T1: SumcheckEngine<E>,
    T2: SumcheckEngine<E>,
    T3: SumcheckEngine<E>,
    T4: SumcheckEngine<E>,
  {
    // sanity checks
    assert_eq!(mem.size(), inner.size());
    assert_eq!(mem.size(), witness.size());
    assert_eq!(mem.size(), rerand.size());
    assert_eq!(mem.degree(), inner.degree());
    assert_eq!(mem.degree(), witness.degree());
    assert_eq!(mem.degree(), rerand.degree());

    // these claims are already added to the transcript, so we do not need to add
    let claims = mem
      .initial_claims()
      .into_iter()
      .chain(inner.initial_claims())
      .chain(witness.initial_claims())
      .chain(rerand.initial_claims())
      .collect::<Vec<E::Scalar>>();

    let s = transcript.squeeze(b"r")?;
    let coeffs = powers::<E>(&s, claims.len());

    // compute the joint claim
    let claim = zip_with!((claims.iter(), coeffs.iter()), |c_1, c_2| *c_1 * c_2).sum();

    let mut e = claim;
    let mut r: Vec<E::Scalar> = Vec::new();
    let mut cubic_polys: Vec<CompressedUniPoly<E::Scalar>> = Vec::new();
    let num_rounds = mem.size().log_2();
    for _ in 0..num_rounds {
      let (evals_mem, (evals_inner, (evals_witness, evals_rerand))) = rayon::join(
        || mem.evaluation_points(),
        || {
          rayon::join(
            || inner.evaluation_points(),
            || {
              rayon::join(
                || witness.evaluation_points(),
                || rerand.evaluation_points(),
              )
            },
          )
        },
      );

      let evals: Vec<Vec<E::Scalar>> = evals_mem
        .into_iter()
        .chain(evals_inner)
        .chain(evals_witness)
        .chain(evals_rerand)
        .collect::<Vec<Vec<E::Scalar>>>();
      assert_eq!(evals.len(), claims.len());

      let evals_combined_0 = (0..evals.len()).map(|i| evals[i][0] * coeffs[i]).sum();
      let evals_combined_bound_coeff = (0..evals.len()).map(|i| evals[i][1] * coeffs[i]).sum();
      let evals_combined_inf = (0..evals.len()).map(|i| evals[i][2] * coeffs[i]).sum();

      let evals = vec![
        evals_combined_0,
        e - evals_combined_0,
        evals_combined_bound_coeff,
        evals_combined_inf,
      ];

      let poly = UniPoly::from_evals_deg3(&evals);

      // append the prover's message to the transcript
      transcript.absorb(b"p", &poly);

      // derive the verifier's challenge for the next round
      let r_i = transcript.squeeze(b"c")?;
      r.push(r_i);

      let _ = rayon::join(
        || mem.bound(&r_i),
        || {
          rayon::join(
            || inner.bound(&r_i),
            || rayon::join(|| witness.bound(&r_i), || rerand.bound(&r_i)),
          )
        },
      );

      e = poly.evaluate(&r_i);
      cubic_polys.push(poly.compress());
    }

    let mem_claims = mem.final_claims();
    let inner_claims = inner.final_claims();
    let witness_claims = witness.final_claims();
    let rerand_claims = rerand.final_claims();

    Ok((
      SumcheckProof::new(cubic_polys),
      r,
      mem_claims,
      inner_claims,
      witness_claims,
      rerand_claims,
    ))
  }
}

impl<E: Engine, EE: EvaluationEngineTrait<E>> VerifierKey<E, EE> {
  fn new(
    num_cons: usize,
    num_vars: usize,
    S_comm: R1CSShapeSparkCommitment<E>,
    vk_ee: EE::VerifierKey,
  ) -> Self {
    VerifierKey {
      num_cons,
      num_vars,
      S_comm,
      vk_ee,
      digest: Default::default(),
    }
  }
}
impl<E: Engine, EE: EvaluationEngineTrait<E>> DigestHelperTrait<E> for VerifierKey<E, EE> {
  /// Returns the digest of the verifier's key
  fn digest(&self) -> E::Scalar {
    self
      .digest
      .get_or_try_init(|| {
        let dc = DigestComputer::new(self);
        dc.digest()
      })
      .cloned()
      .expect("Failure to retrieve digest!")
  }
}

impl<E: Engine, EE: EvaluationEngineTrait<E>> RelaxedR1CSSNARKTrait<E> for RelaxedR1CSSNARK<E, EE> {
  type ProverKey = ProverKey<E, EE>;
  type VerifierKey = VerifierKey<E, EE>;

  fn ck_floor() -> Box<dyn for<'a> Fn(&'a R1CSShape<E>) -> usize> {
    Box::new(|shape: &R1CSShape<E>| -> usize {
      // the commitment key should be large enough to commit to the R1CS matrices
      shape.A.len() + shape.B.len() + shape.C.len()
    })
  }

  fn setup(
    ck: &CommitmentKey<E>,
    S: &R1CSShape<E>,
  ) -> Result<(Self::ProverKey, Self::VerifierKey), NovaError> {
    // check the provided commitment key meets minimal requirements
    if ck.length() < Self::ck_floor()(S) {
      return Err(NovaError::InvalidCommitmentKeyLength);
    }
    let (pk_ee, vk_ee) = EE::setup(ck)?;

    // pad the R1CS matrices
    let S = S.pad();

    let S_repr = R1CSShapeSparkRepr::new(&S);
    let S_comm = S_repr.commit(ck);

    let vk = VerifierKey::new(S.num_cons, S.num_vars, S_comm.clone(), vk_ee);

    let pk = ProverKey {
      pk_ee,
      S_repr,
      S_comm,
      vk_digest: vk.digest(),
    };

    Ok((pk, vk))
  }

  /// produces a succinct proof of satisfiability of a `RelaxedR1CS` instance
  fn prove(
    ck: &CommitmentKey<E>,
    pk: &Self::ProverKey,
    S: &R1CSShape<E>,
    U: &RelaxedR1CSInstance<E>,
    W: &RelaxedR1CSWitness<E>,
  ) -> Result<Self, NovaError> {
    // pad the R1CSShape
    let S = S.pad();
    // sanity check that R1CSShape has all required size characteristics
    assert!(S.is_regular_shape());

    let W = W.pad(&S); // pad the witness
    let mut transcript = E::TE::new(b"RelaxedR1CSSNARK");

    // append the verifier key (which includes commitment to R1CS matrices) and the RelaxedR1CSInstance to the transcript
    transcript.absorb(b"vk", &pk.vk_digest);
    transcript.absorb(b"U", U);

    // compute the full satisfying assignment by concatenating W.W, U.u, and U.X
    let z = [W.W.clone(), vec![U.u], U.X.clone()].concat();

    // compute Az, Bz, Cz
    let (Az, Bz, Cz) = S.multiply_vec(&z)?;

    // Shortened outer sum-check: run for log(num_cons) rounds instead of log(N).
    // Az, Bz, Cz, E are naturally size num_cons and zero-padded to N, so the
    // sum over {0,1}^log(N) collapses to {0,1}^log(num_cons).
    let num_rounds_outer = S.num_cons.log_2();
    let num_rounds_inner = pk.S_repr.N.log_2();
    let tau = (0..num_rounds_outer)
      .map(|_| transcript.squeeze(b"t"))
      .collect::<Result<Vec<_>, NovaError>>()?;

    // Step 1: Outer sum-check (standalone, degree-3) on size-m polynomials
    // Proves: 0 = Σ_{x∈{0,1}^log(m)} eq(τ,x) * (Az(x) * Bz(x) - (u·Cz(x) + E(x)))
    let uCz_E: Vec<E::Scalar> = Cz
      .iter()
      .zip(W.E.iter())
      .map(|(cz, e)| U.u * *cz + *e)
      .collect();
    let mut poly_Az = MultilinearPolynomial::new(Az);
    let mut poly_Bz = MultilinearPolynomial::new(Bz);
    let mut poly_uCz_E = MultilinearPolynomial::new(uCz_E);

    let (sc_outer, r_outer, claims_outer) = SumcheckProof::<E>::prove_cubic_with_three_inputs(
      &E::Scalar::ZERO,
      tau,
      &mut poly_Az,
      &mut poly_Bz,
      &mut poly_uCz_E,
      &mut transcript,
    )?;

    // Claims from the shortened outer sum-check (evaluations at the log(m)-length point)
    let eval_Az_at_r_outer = claims_outer[0];
    let eval_Bz_at_r_outer = claims_outer[1];
    let eval_Cz_at_r_outer = MultilinearPolynomial::evaluate_with(&Cz, &r_outer);
    let eval_E_at_r_outer = claims_outer[2] - U.u * eval_Cz_at_r_outer;

    // Absorb outer sum-check claims into transcript
    transcript.absorb(
      b"e",
      &[
        eval_Az_at_r_outer,
        eval_Bz_at_r_outer,
        eval_Cz_at_r_outer,
        eval_E_at_r_outer,
      ]
      .as_slice(),
    );

    // Squeeze random padding challenges and extend r_outer to length log(N).
    // For zero-padded polys P of size m within N: P(r_outer_full) = factor · P(r_outer)
    // where factor = Π(1 - r_pad_j) and r_outer_full = (r_pad, r_outer).
    // r_pad occupies the top (MSB) positions since padding variables are the high-order bits.
    let num_pad_rounds = num_rounds_inner
      .checked_sub(num_rounds_outer)
      .ok_or(NovaError::InvalidSumcheckProof)?;
    let r_pad = (0..num_pad_rounds)
      .map(|_| transcript.squeeze(b"p"))
      .collect::<Result<Vec<E::Scalar>, NovaError>>()?;
    let r_outer_full: Vec<E::Scalar> = r_pad.iter().chain(r_outer.iter()).cloned().collect();
    let factor: E::Scalar = r_pad
      .iter()
      .fold(E::Scalar::ONE, |acc, r| acc * (E::Scalar::ONE - r));

    // Pad E and W to size N for inner sum-check and PCS
    let E = padded::<E>(&W.E, pk.S_repr.N, &E::Scalar::ZERO);
    let W = padded::<E>(&W.W, pk.S_repr.N, &E::Scalar::ZERO);

    // -----------------------------------------------------------------------
    // Step 2: Prepare the batched inner batched sum-check (memory + inner_batched + witness)
    // -----------------------------------------------------------------------

    // Compute evaluation oracles at r_outer_full (the extended outer challenge)
    // L_row(i) = eq(r_outer_full, row(i)) for all i
    // L_col(i) = z(col(i)) for all i, where z is the full satisfying assignment
    let (mem_row, mem_col, L_row, L_col) = pk.S_repr.evaluation_oracles(&S, &r_outer_full, &z);
    let (comm_L_row, comm_L_col) = rayon::join(
      || E::CE::commit(ck, &L_row, &E::Scalar::ZERO),
      || E::CE::commit(ck, &L_col, &E::Scalar::ZERO),
    );

    // Absorb commitments to L_row and L_col
    transcript.absorb(b"e", &vec![comm_L_row, comm_L_col].as_slice());

    // Squeeze challenge for batching inner batched ABC claims
    let c = transcript.squeeze(b"c")?;

    // Squeeze challenges for memory sum-check
    let gamma = transcript.squeeze(b"g")?;
    let r = transcript.squeeze(b"r")?;

    let (mut inner_batched_sc_inst, mem_res) = rayon::join(
      || {
        // Inner batched sum-check instance for:
        // (a) ABC claim: factor·(v_A + c·v_B + c²·v_C) = Σ L_row(y) * (val_A + c·val_B + c²·val_C)(y) * L_col(y)
        // (b) E claim: factor·eval_E = Σ eq(r_outer_full, y) * E(y)
        // The claims are scaled by factor because the inner sum-check uses r_outer_full
        // and eq(r_outer_full, j) = factor · eq(r_outer, j) for j < m.
        let val = zip_with!(
          par_iter,
          (pk.S_repr.val_A, pk.S_repr.val_B, pk.S_repr.val_C),
          |v_a, v_b, v_c| *v_a + c * *v_b + c * c * *v_c
        )
        .collect::<Vec<E::Scalar>>();

        InnerBatchedSumcheckInstance::new(
          factor * (eval_Az_at_r_outer + c * eval_Bz_at_r_outer + c * c * eval_Cz_at_r_outer),
          L_row.clone(),
          L_col.clone(),
          val,
          factor * eval_E_at_r_outer,
          r_outer_full.clone(), // eq challenges for factored eq (Gruen)
          E.clone(),
        )
      },
      || {
        // Inverse-logup memory sum-check instance to prove L_row/L_col are
        // well-formed.
        mem_check_logup::prove_step::<E>(
          ck,
          r,
          gamma,
          &mem_row,
          &pk.S_repr.row,
          &L_row,
          &pk.S_repr.ts_row,
          &mem_col,
          &pk.S_repr.col,
          &L_col,
          &pk.S_repr.ts_col,
          num_rounds_inner,
          &mut transcript,
        )
      },
    );

    let (mut mem_sc_inst, comm_mem_oracles, mem_oracles) = mem_res?;

    // Witness bound sum-check using r_outer_full as the random evaluation point
    let mut witness_sc_inst =
      WitnessBoundSumcheck::new(r_outer_full.clone(), W.clone(), S.num_vars);

    // Logup-GKR memory-check tree fold: folds the four sub-instances into GKR
    // proofs and emits the rerandomize instance carrying the reconcile columns
    // into the inner sumcheck below. Reuses the same fingerprint (gamma, r); the
    // tree fold absorbs into the transcript here, BEFORE the inner sumcheck.
    let logup_gkr_witness = MemCheckWitness::<E> {
      mem_row: mem_row.clone(),
      mem_col: mem_col.clone(),
      L_row: L_row.clone(),
      L_col: L_col.clone(),
      addr_row: pk.S_repr.row.clone(),
      addr_col: pk.S_repr.col.clone(),
      ts_row: pk.S_repr.ts_row.clone(),
      ts_col: pk.S_repr.ts_col.clone(),
    };
    let (mut rerandomize_sc_inst, gkr_proof_data, _gkr_rerand_claims) =
      mem_check_logup_gkr::prove_step::<E>(logup_gkr_witness, gamma, r, &mut transcript)?;

    // -----------------------------------------------------------------------
    // Step 3: Run the batched inner batched sum-check (4 instances)
    // -----------------------------------------------------------------------
    let (
      sc_inner_batched,
      r_inner_batched,
      claims_mem,
      claims_inner_batched,
      claims_witness,
      claims_rerand,
    ) = Self::prove_helper(
      &mut mem_sc_inst,
      &mut inner_batched_sc_inst,
      &mut witness_sc_inst,
      &mut rerandomize_sc_inst,
      &mut transcript,
    )?;

    // Claims from the inner batched sum-check
    let eval_L_row = claims_inner_batched[0][0];
    let eval_L_col = claims_inner_batched[0][1];
    let eval_E = claims_inner_batched[1][0]; // E(r_inner_batched) — rerandomized to open at same point

    // The rerandomize instance's L_row/L_col at r_inner_batched must equal the
    // inner ABC path's — they are the same polynomial at the same point.
    debug_assert_eq!(claims_rerand[0][0], eval_L_row);
    debug_assert_eq!(claims_rerand[1][0], eval_L_col);

    let eval_t_plus_r_inv_row = claims_mem[0][0];
    let eval_w_plus_r_inv_row = claims_mem[0][1];
    let eval_ts_row = claims_mem[0][2];

    let eval_t_plus_r_inv_col = claims_mem[1][0];
    let eval_w_plus_r_inv_col = claims_mem[1][1];
    let eval_ts_col = claims_mem[1][2];
    let eval_W = claims_witness[0][0];

    // Compute evaluations at r_inner_batched that did not come for free from the sum-check
    let (eval_val_A, eval_val_B, eval_val_C, eval_row, eval_col) = {
      let e = MultilinearPolynomial::multi_evaluate_with(
        &[
          &pk.S_repr.val_A,
          &pk.S_repr.val_B,
          &pk.S_repr.val_C,
          &pk.S_repr.row,
          &pk.S_repr.col,
        ],
        &r_inner_batched,
      );
      (e[0], e[1], e[2], e[3], e[4])
    };

    // All evaluations are at r_inner_batched — fold into one claim for batch PCS opening
    let eval_vec = vec![
      eval_W,
      eval_E,
      eval_L_row,
      eval_L_col,
      eval_val_A,
      eval_val_B,
      eval_val_C,
      eval_t_plus_r_inv_row,
      eval_row,
      eval_w_plus_r_inv_row,
      eval_ts_row,
      eval_t_plus_r_inv_col,
      eval_col,
      eval_w_plus_r_inv_col,
      eval_ts_col,
    ];

    let comm_vec = [
      U.comm_W,
      U.comm_E,
      comm_L_row,
      comm_L_col,
      pk.S_comm.comm_val_A,
      pk.S_comm.comm_val_B,
      pk.S_comm.comm_val_C,
      comm_mem_oracles[0],
      pk.S_comm.comm_row,
      comm_mem_oracles[1],
      pk.S_comm.comm_ts_row,
      comm_mem_oracles[2],
      pk.S_comm.comm_col,
      comm_mem_oracles[3],
      pk.S_comm.comm_ts_col,
    ];
    let poly_vec = [
      &W,
      &E,
      &L_row,
      &L_col,
      &pk.S_repr.val_A,
      &pk.S_repr.val_B,
      &pk.S_repr.val_C,
      mem_oracles[0].as_ref(),
      &pk.S_repr.row,
      mem_oracles[1].as_ref(),
      &pk.S_repr.ts_row,
      mem_oracles[2].as_ref(),
      &pk.S_repr.col,
      mem_oracles[3].as_ref(),
      &pk.S_repr.ts_col,
    ];
    transcript.absorb(b"e", &eval_vec.as_slice()); // comm_vec is already in the transcript
    let c = transcript.squeeze(b"c")?;
    let w: PolyEvalWitness<E> = PolyEvalWitness::batch(&poly_vec, &c);
    let u: PolyEvalInstance<E> =
      PolyEvalInstance::batch(&comm_vec, &r_inner_batched, &eval_vec, &c);

    let eval_arg = EE::prove(
      ck,
      &pk.pk_ee,
      &mut transcript,
      &u.c,
      &w.p,
      &r_inner_batched,
      &u.e,
    )?;

    Ok(RelaxedR1CSSNARK {
      comm_L_row,
      comm_L_col,

      mem_check: mem_check_logup::LogupProofData {
        comm_t_plus_r_inv_row: comm_mem_oracles[0],
        comm_w_plus_r_inv_row: comm_mem_oracles[1],
        comm_t_plus_r_inv_col: comm_mem_oracles[2],
        comm_w_plus_r_inv_col: comm_mem_oracles[3],
        eval_t_plus_r_inv_row,
        eval_w_plus_r_inv_row,
        eval_t_plus_r_inv_col,
        eval_w_plus_r_inv_col,
      },
      mem_check_gkr: gkr_proof_data,

      sc_outer,

      eval_Az_at_r_outer,
      eval_Bz_at_r_outer,
      eval_Cz_at_r_outer,
      eval_E_at_r_outer,

      sc_inner_batched,

      eval_E,
      eval_L_row,
      eval_L_col,
      eval_val_A,
      eval_val_B,
      eval_val_C,

      eval_W,

      eval_row,
      eval_ts_row,
      eval_col,
      eval_ts_col,

      eval_arg,
    })
  }

  /// verifies a proof of satisfiability of a `RelaxedR1CS` instance
  fn verify(&self, vk: &Self::VerifierKey, U: &RelaxedR1CSInstance<E>) -> Result<(), NovaError> {
    let mut transcript = E::TE::new(b"RelaxedR1CSSNARK");

    // append the verifier key (including commitment to R1CS matrices) and the RelaxedR1CSInstance to the transcript
    transcript.absorb(b"vk", &vk.digest());
    transcript.absorb(b"U", U);

    // Shortened outer sum-check runs for log(num_cons) rounds
    let num_rounds_outer = vk.num_cons.log_2();
    let num_rounds_inner = vk.S_comm.N.log_2();
    let tau = (0..num_rounds_outer)
      .map(|_| transcript.squeeze(b"t"))
      .collect::<Result<Vec<_>, NovaError>>()?;

    // -----------------------------------------------------------------------
    // Step 1: Verify the shortened outer sum-check
    // Claim: 0 = Σ_{x∈{0,1}^log(m)} eq(τ,x) * (Az(x) * Bz(x) - (u·Cz(x) + E(x)))
    // -----------------------------------------------------------------------
    let (claim_sc_outer_final, r_outer) =
      self
        .sc_outer
        .verify(E::Scalar::ZERO, num_rounds_outer, 3, &mut transcript)?;

    // Check outer sum-check final claim
    let eq_tau_at_r_outer = EqPolynomial::new(tau).evaluate(&r_outer);
    let claim_sc_outer_expected = eq_tau_at_r_outer
      * (self.eval_Az_at_r_outer * self.eval_Bz_at_r_outer
        - U.u * self.eval_Cz_at_r_outer
        - self.eval_E_at_r_outer);
    if claim_sc_outer_expected != claim_sc_outer_final {
      return Err(NovaError::InvalidSumcheckProof);
    }

    // Absorb outer sum-check claims
    transcript.absorb(
      b"e",
      &[
        self.eval_Az_at_r_outer,
        self.eval_Bz_at_r_outer,
        self.eval_Cz_at_r_outer,
        self.eval_E_at_r_outer,
      ]
      .as_slice(),
    );

    // Squeeze random padding challenges and extend r_outer to length log(N)
    // r_pad occupies the top (MSB) positions since padding variables are the high-order bits.
    let num_pad_rounds = num_rounds_inner
      .checked_sub(num_rounds_outer)
      .ok_or(NovaError::InvalidSumcheckProof)?;
    let r_pad = (0..num_pad_rounds)
      .map(|_| transcript.squeeze(b"p"))
      .collect::<Result<Vec<E::Scalar>, NovaError>>()?;
    let r_outer_full: Vec<E::Scalar> = r_pad.iter().chain(r_outer.iter()).cloned().collect();
    let factor: E::Scalar = r_pad
      .iter()
      .fold(E::Scalar::ONE, |acc, r| acc * (E::Scalar::ONE - r));

    // -----------------------------------------------------------------------
    // Step 2: Verify the batched inner batched sum-check (memory + inner_batched + witness)
    // -----------------------------------------------------------------------

    // Absorb commitments to L_row and L_col
    transcript.absorb(b"e", &vec![self.comm_L_row, self.comm_L_col].as_slice());

    // Squeeze challenge for batching inner batched ABC claims
    let c = transcript.squeeze(b"c")?;

    let gamma = transcript.squeeze(b"g")?;
    let r = transcript.squeeze(b"r")?;

    // Inverse-logup memory-check, transcript phase: absorb the inverse-oracle
    // commitments and squeeze the memory sumcheck's eq randomness `rho`.
    let rho =
      mem_check_logup::verify_pre_inner::<E>(&self.mem_check, num_rounds_inner, &mut transcript)?;

    // Logup-GKR memory-check, transcript phase: replay the GKR proof, reconcile
    // the prover-claimed columns against the reduction, run the balance check,
    // and absorb the claimed values — all before the inner sumcheck's `s`. The
    // reconcile/balance are transcript-free algebra; the claimed values are
    // still bound to the committed columns by the rerandomize sub-claims below.
    let gkr_eval_point = mem_check_logup_gkr::verify_pre_inner::<E>(
      &self.mem_check_gkr,
      gamma,
      r,
      &r_outer_full,
      &mut transcript,
    )?;

    // 16 claims: 6 memory + 2 inner (ABC + E) + 1 witness + 7 rerandomize.
    let num_claims = mem_check_logup_gkr::RERAND_BASE + mem_check_logup_gkr::NUM_RERAND_COLUMNS;
    let s = transcript.squeeze(b"r")?;
    let coeffs = powers::<E>(&s, num_claims);

    // Compute the combined initial claim. The memory-check slot contributes one
    // term (`mem_initial`); the ABC/E terms are shared. Under Logup-GKR the slot
    // term is the rerandomize columns' claimed values (coeffs 9-15); under
    // inverse-logup it is zero (the six memory routes prove `0 = Σ...`).
    // Claim 6: inner ABC = factor * (eval_Az + c * eval_Bz + c² * eval_Cz)
    // Claim 7: inner E = factor * eval_E
    // Claim 8: witness (zero)
    // The factor accounts for zero-padding: eval_P(r_outer_full) = factor * eval_P(r_outer)
    let claim_inner_batched_ABC = factor
      * (self.eval_Az_at_r_outer + c * self.eval_Bz_at_r_outer + c * c * self.eval_Cz_at_r_outer);
    let mem_initial = mem_check_logup_gkr::verify_initial_claim::<E>(&self.mem_check_gkr, &coeffs);
    let claim = coeffs[6] * claim_inner_batched_ABC
      + coeffs[7] * factor * self.eval_E_at_r_outer
      + mem_initial;

    // Verify inner batched sum-check
    let (claim_sc_inner_batched_final, r_inner_batched) =
      self
        .sc_inner_batched
        .verify(claim, num_rounds_inner, 3, &mut transcript)?;

    // Verify inner batched sum-check final claim
    let claim_sc_inner_batched_expected = {
      let eq_r_outer = EqPolynomial::new(r_outer_full.clone());
      let eq_r_outer_at_r_inner_batched = eq_r_outer.evaluate(&r_inner_batched);
      let taus_masked_bound_r_inner_batched =
        MaskedEqPolynomial::new(&eq_r_outer, vk.num_vars.log_2()).evaluate(&r_inner_batched);

      // mem_col = z at r_inner_batched, reconstructed from eval_W and public IO.
      let eval_mem_col_at_r_inner = {
        let (factor, r_inner_batched_unpad) = {
          let l = vk.S_comm.N.log_2() - (2 * vk.num_vars).log_2();
          let mut factor = E::Scalar::ONE;
          for r_p in r_inner_batched.iter().take(l) {
            factor *= E::Scalar::ONE - r_p
          }
          (factor, r_inner_batched[l..].to_vec())
        };
        let eval_X = {
          let X = vec![U.u]
            .into_iter()
            .chain(U.X.iter().cloned())
            .collect::<Vec<E::Scalar>>();
          let poly_X = SparsePolynomial::new(r_inner_batched_unpad.len() - 1, X);
          poly_X.evaluate(&r_inner_batched_unpad[1..])
        };
        self.eval_W + factor * r_inner_batched_unpad[0] * eval_X
      };

      // Inverse-logup memory claims (coeffs 0-5).
      let public_io: Vec<E::Scalar> = vec![U.u].into_iter().chain(U.X.iter().cloned()).collect();
      // Memory-check slot's final-claim contribution. Under Logup-GKR it is the
      // rerandomize columns' claims (coeffs 9-15); under inverse-logup it is the
      // six memory routes (coeffs 0-5). Both are computed here for now (the two
      // implementations run alongside each other); a feature gate later selects
      // one by replacing this single `mem_final` term.
      let mem_final_logup = mem_check_logup::verify_final_claim::<E>(
        &self.mem_check,
        &coeffs,
        rho,
        gamma,
        r,
        &r_outer_full,
        &r_inner_batched,
        num_rounds_inner,
        vk.num_vars,
        vk.S_comm.N,
        self.eval_W,
        self.eval_L_row,
        self.eval_L_col,
        self.eval_row,
        self.eval_col,
        self.eval_ts_row,
        self.eval_ts_col,
        &public_io,
      );

      // Inner batched ABC claim (coeff 6): L_row * L_col * (val_A + c·val_B + c²·val_C)
      let claim_inner_batched_ABC_final = coeffs[6]
        * self.eval_L_row
        * self.eval_L_col
        * (self.eval_val_A + c * self.eval_val_B + c * c * self.eval_val_C);

      // Inner batched E claim (coeff 7): eq(r_outer_full, r_inner_batched) * E(r_inner_batched)
      let claim_inner_batched_E_final = coeffs[7] * eq_r_outer_at_r_inner_batched * self.eval_E;

      // Witness claim (coeff 8)
      let claim_witness_final = coeffs[8] * taus_masked_bound_r_inner_batched * self.eval_W;

      // Logup-GKR rerandomize claims (coeffs 9-15).
      let rerand_col_evals = [
        self.eval_L_row,
        self.eval_L_col,
        self.eval_row,
        self.eval_col,
        self.eval_ts_row,
        self.eval_ts_col,
        eval_mem_col_at_r_inner,
      ];
      let mem_final_gkr = mem_check_logup_gkr::verify_final_claim::<E>(
        &coeffs,
        &gkr_eval_point,
        &r_inner_batched,
        &rerand_col_evals,
      );

      // The memory-check slot's final contribution (both implementations for now).
      let mem_final = mem_final_logup + mem_final_gkr;

      mem_final + claim_inner_batched_ABC_final + claim_inner_batched_E_final + claim_witness_final
    };

    if claim_sc_inner_batched_expected != claim_sc_inner_batched_final {
      return Err(NovaError::InvalidSumcheckProof);
    }

    // Verify polynomial openings — all at r_inner_batched
    let eval_vec = vec![
      self.eval_W,
      self.eval_E,
      self.eval_L_row,
      self.eval_L_col,
      self.eval_val_A,
      self.eval_val_B,
      self.eval_val_C,
      self.mem_check.eval_t_plus_r_inv_row,
      self.eval_row,
      self.mem_check.eval_w_plus_r_inv_row,
      self.eval_ts_row,
      self.mem_check.eval_t_plus_r_inv_col,
      self.eval_col,
      self.mem_check.eval_w_plus_r_inv_col,
      self.eval_ts_col,
    ];
    let comm_vec = [
      U.comm_W,
      U.comm_E,
      self.comm_L_row,
      self.comm_L_col,
      vk.S_comm.comm_val_A,
      vk.S_comm.comm_val_B,
      vk.S_comm.comm_val_C,
      self.mem_check.comm_t_plus_r_inv_row,
      vk.S_comm.comm_row,
      self.mem_check.comm_w_plus_r_inv_row,
      vk.S_comm.comm_ts_row,
      self.mem_check.comm_t_plus_r_inv_col,
      vk.S_comm.comm_col,
      self.mem_check.comm_w_plus_r_inv_col,
      vk.S_comm.comm_ts_col,
    ];
    transcript.absorb(b"e", &eval_vec.as_slice()); // comm_vec is already in the transcript
    let c = transcript.squeeze(b"c")?;
    let u: PolyEvalInstance<E> =
      PolyEvalInstance::batch(&comm_vec, &r_inner_batched, &eval_vec, &c);

    // verify
    EE::verify(
      &vk.vk_ee,
      &mut transcript,
      &u.c,
      &r_inner_batched,
      &u.e,
      &self.eval_arg,
    )?;

    Ok(())
  }
}
