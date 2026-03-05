//! This module defines a collection of traits that define the behavior of a commitment engine
//! We require the commitment engine to provide a commitment to vectors with a single group element
#[cfg(feature = "io")]
use crate::provider::ptau::PtauFileError;
use crate::traits::{AbsorbInRO2Trait, AbsorbInROTrait, Engine, TranscriptReprTrait};
use core::{
  fmt::Debug,
  ops::{Add, Mul, MulAssign, Range},
};
use num_integer::Integer;
use num_traits::ToPrimitive;
use rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator, ParallelIterator};
use serde::{Deserialize, Serialize};

/// A helper trait for types implementing scalar multiplication.
pub trait ScalarMul<Rhs, Output = Self>: Mul<Rhs, Output = Output> + MulAssign<Rhs> {}

impl<T, Rhs, Output> ScalarMul<Rhs, Output> for T where T: Mul<Rhs, Output = Output> + MulAssign<Rhs>
{}

/// This trait defines the behavior of the commitment
pub trait CommitmentTrait<E: Engine>:
  Clone
  + Copy
  + Debug
  + Default
  + PartialEq
  + Eq
  + Send
  + Sync
  + TranscriptReprTrait<E::GE>
  + Serialize
  + for<'de> Deserialize<'de>
  + AbsorbInROTrait<E>
  + AbsorbInRO2Trait<E>
  + Add<Self, Output = Self>
  + ScalarMul<E::Scalar>
{
  /// Returns the coordinate representation of the commitment
  fn to_coordinates(&self) -> (E::Base, E::Base, bool);
}

/// A trait that helps determine the length of a structure.
/// Note this does not impose any memory representation constraints on the structure.
pub trait Len {
  /// Returns the length of the structure.
  fn length(&self) -> usize;
}

/// A trait that ties different pieces of the commitment generation together
pub trait CommitmentEngineTrait<E: Engine>: Clone + Send + Sync {
  /// Holds the type of the commitment key
  /// The key should quantify its length in terms of group generators.
  type CommitmentKey: Len + Clone + Debug + Send + Sync + Serialize + for<'de> Deserialize<'de>;

  /// Holds the type of the derandomization key
  type DerandKey: Clone + Debug + Send + Sync + Serialize + for<'de> Deserialize<'de>;

  /// Holds the type of the commitment
  type Commitment: CommitmentTrait<E>;

  /// Load keys
  #[cfg(feature = "io")]
  fn load_setup(
    reader: &mut (impl std::io::Read + std::io::Seek),
    label: &'static [u8],
    n: usize,
  ) -> Result<Self::CommitmentKey, PtauFileError>;

  /// Saves the key to the provided writer.
  #[cfg(feature = "io")]
  fn save_setup(
    ck: &Self::CommitmentKey,
    writer: &mut (impl std::io::Write + std::io::Seek),
  ) -> Result<(), PtauFileError>;

  /// Samples a new commitment key of a specified size.
  ///
  /// # Errors
  ///
  /// Returns an error if the setup cannot be performed (e.g., HyperKZG in production
  /// builds without the `test-utils` feature).
  fn setup(label: &'static [u8], n: usize)
    -> Result<Self::CommitmentKey, crate::errors::NovaError>;

  /// Extracts the blinding generator
  fn derand_key(ck: &Self::CommitmentKey) -> Self::DerandKey;

  /// Commits to the provided vector using the provided generators and random blind
  fn commit(ck: &Self::CommitmentKey, v: &[E::Scalar], r: &E::Scalar) -> Self::Commitment;

  /// Batch commits to the provided vectors using the provided generators and random blind
  fn batch_commit(
    ck: &Self::CommitmentKey,
    v: &[Vec<E::Scalar>],
    r: &[E::Scalar],
  ) -> Vec<Self::Commitment> {
    assert!(v.len() == r.len());
    v.par_iter()
      .zip(r.par_iter())
      .map(|(v_i, r_i)| Self::commit(ck, v_i, r_i))
      .collect()
  }

  /// Commits to the provided vector of sparse binary scalars using the provided generators and random blind
  fn commit_sparse_binary(
    ck: &Self::CommitmentKey,
    non_zero_indices: &[usize],
    r: &E::Scalar,
  ) -> Self::Commitment;

  /// Holds the type of the precomputed one-hot table.
  ///
  /// The table precomputes pairwise sums of generators for each consecutive pair of
  /// blocks (t=2). This halves the number of point additions during batch commits.
  type OneHotTable: Send + Sync;

  /// Commits to a one-hot structured vector where each block of `block_size` elements
  /// has exactly one non-zero (=1) entry at the specified offset.
  ///
  /// `block_offsets[i]` gives the offset within block `i` (`0 <= offset < block_size`).
  /// The global index for block `i` is `i * block_size + block_offsets[i]`.
  fn commit_one_hot(
    ck: &Self::CommitmentKey,
    block_size: usize,
    block_offsets: &[usize],
    r: &E::Scalar,
  ) -> Self::Commitment;

  /// Batch commits to multiple one-hot structured vectors.
  ///
  /// Default implementation uses `par_iter` over individual `commit_one_hot` calls.
  fn batch_commit_one_hot(
    ck: &Self::CommitmentKey,
    block_size: usize,
    all_offsets: &[Vec<usize>],
    r: &[E::Scalar],
  ) -> Vec<Self::Commitment> {
    assert_eq!(all_offsets.len(), r.len());
    all_offsets
      .par_iter()
      .zip(r.par_iter())
      .map(|(offsets, r_i)| Self::commit_one_hot(ck, block_size, offsets, r_i))
      .collect()
  }

  /// Builds a precomputed one-hot table for the given commitment key.
  ///
  /// The table precomputes K² pairwise sums for each consecutive pair of blocks,
  /// enabling `batch_commit_one_hot_with_table` to halve point additions.
  ///
  /// One-time cost: `(num_blocks / 2) × K²` point additions.
  /// Memory: `(num_blocks / 2) × K²` affine points (for K=16, ~8 GB at 2^20 blocks).
  fn build_one_hot_table(
    ck: &Self::CommitmentKey,
    block_size: usize,
    num_blocks: usize,
  ) -> Self::OneHotTable;

  /// Batch commits using a precomputed one-hot table with transposed iteration.
  ///
  /// Combines precomputed pair tables (t=2) with block-major iteration order.
  /// Each pair of blocks requires only 1 table lookup instead of 2 point additions,
  /// and the transposed order keeps per-commit accumulators hot in cache.
  ///
  /// Benchmark (BN256, K=16, 2^20 blocks, 600 commits):
  /// - `batch_commit_one_hot` (default): 24.8s
  /// - `batch_commit_one_hot_with_table`: 9.3s (2.66× faster)
  fn batch_commit_one_hot_with_table(
    ck: &Self::CommitmentKey,
    table: &Self::OneHotTable,
    all_offsets: &[Vec<usize>],
    r: &[E::Scalar],
  ) -> Vec<Self::Commitment>;

  /// Commits to the provided vector of "small" scalars (at most 64 bits) using the provided generators and random blind
  fn commit_small<T: Integer + Into<u64> + Copy + Sync + ToPrimitive>(
    ck: &Self::CommitmentKey,
    v: &[T],
    r: &E::Scalar,
  ) -> Self::Commitment;

  /// Commits to the provided vector of "small" scalars (at most 64 bits) using the provided generators and random blind (range)
  fn commit_small_range<T: Integer + Into<u64> + Copy + Sync + ToPrimitive>(
    ck: &Self::CommitmentKey,
    v: &[T],
    r: &E::Scalar,
    range: Range<usize>,
    max_num_bits: usize,
  ) -> Self::Commitment;

  /// Batch commits to the provided vectors of "small" scalars (at most 64 bits) using the provided generators and random blind
  fn batch_commit_small<T: Integer + Into<u64> + Copy + Sync + ToPrimitive>(
    ck: &Self::CommitmentKey,
    v: &[Vec<T>],
    r: &[E::Scalar],
  ) -> Vec<Self::Commitment> {
    assert!(v.len() == r.len());
    v.par_iter()
      .zip(r.par_iter())
      .map(|(v_i, r_i)| Self::commit_small(ck, v_i, r_i))
      .collect()
  }

  /// Remove given blind from commitment
  fn derandomize(
    dk: &Self::DerandKey,
    commit: &Self::Commitment,
    r: &E::Scalar,
  ) -> Self::Commitment;

  /// Returns the coordinates of each generator in the commitment key.
  ///
  /// This method extracts the (x, y) coordinates of each generator point
  /// in the commitment key. This is useful for operations that need direct
  /// access to the underlying elliptic curve points, such as in-circuit
  /// verification of polynomial evaluations.
  ///
  /// # Panics
  ///
  /// Panics if any generator point is the point at infinity.
  fn ck_to_coordinates(ck: &Self::CommitmentKey) -> Vec<(E::Base, E::Base)>;

  /// Returns the generators as projective group elements.
  ///
  /// This provides full group-operation access (add, double, scalar mul)
  /// on the commitment key generators, useful for precomputing correction
  /// points in circuit optimizations.
  fn ck_to_group_elements(ck: &Self::CommitmentKey) -> Vec<E::GE>;
}
