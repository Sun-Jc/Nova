//! Benchmarking sparse vs dense commitment for boolean vectors using HyperKZG over BN254.
//!
//! This benchmark compares:
//! - Sparse commitment: commits directly from non-zero indices
//! - Dense commitment: constructs full vector from indices, then commits using optimized 1-bit path
//! - Vector construction only: measures just the time to build the dense vector from indices
//!
//! Parameters:
//! - Vector length: 64 to 2^25
//! - Density: percentage of non-zero entries (10%, 25%, 50%)
use core::time::Duration;
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use nova_snark::{
  provider::Bn256EngineKZG,
  traits::{commitment::CommitmentEngineTrait, Engine},
};
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha20Rng;

criterion_group! {
name = sparse_commit;
config = Criterion::default().warm_up_time(Duration::from_millis(3000)).sample_size(10);
targets = bench_sparse_commit
}

criterion_main!(sparse_commit);

/// Fixed seed for reproducibility
const SEED: u64 = 0xDEADBEEF_CAFEBABE;

/// Densities to test (as percentages)
const DENSITIES: &[u32] = &[5, 10, 25, 50];

/// Generate non-zero indices for a boolean vector with specified density.
/// Returns only the indices where the value is 1.
fn generate_non_zero_indices(size: usize, density_percent: u32, seed: u64) -> Vec<usize> {
  let mut rng = ChaCha20Rng::seed_from_u64(seed);
  let mut non_zero_indices = Vec::new();

  for i in 0..size {
    // Generate a number 0-99, if < density_percent, it's a 1
    let val: u32 = rng.gen_range(0..100);
    if val < density_percent {
      non_zero_indices.push(i);
    }
  }

  non_zero_indices
}

/// Construct a dense vector from non-zero indices.
/// This is the operation we want to time separately.
/// Returns Vec<u8> with 0s and 1s for use with commit_small_range.
#[inline(never)]
fn make_dense_vector(size: usize, non_zero_indices: &[usize]) -> Vec<u8> {
  let mut dense = vec![0u8; size];
  for &idx in non_zero_indices {
    dense[idx] = 1;
  }
  dense
}

fn bench_sparse_commit(c: &mut Criterion) {
  type E = Bn256EngineKZG;

  let max_size = 1 << 25;

  // Setup commitment key for the maximum size
  println!("Setting up commitment key for size {}...", max_size);
  let ck = <E as Engine>::CE::setup(b"sparse_commit_bench", max_size);
  println!("Commitment key setup complete.");

  let zero = <E as Engine>::Scalar::zero();

  // Pre-generate non-zero indices for all (size, density) combinations
  println!("Pre-generating non-zero indices for all configurations...");

  // Benchmark across different sizes: 64, 128, 256, ..., up to 1<<25
  let mut size = 64usize;
  while size <= max_size {
    for &density in DENSITIES {
      // Use a combined seed for reproducibility across (size, density) pairs
      let combined_seed = SEED
        .wrapping_add(size as u64)
        .wrapping_mul(density as u64 + 1);
      let non_zero_indices = generate_non_zero_indices(size, density, combined_seed);

      let actual_density = (non_zero_indices.len() as f64 / size as f64) * 100.0;
      println!(
        "Size {}, target density {}%: {} non-zero elements ({:.2}% actual)",
        size,
        density,
        non_zero_indices.len(),
        actual_density
      );

      // Benchmark 1: Just vector construction (make_vec)
      c.bench_function(&format!("make_vec_{size}_d{density}"), |b| {
        b.iter(|| black_box(make_dense_vector(size, &non_zero_indices)))
      });

      // Benchmark 2: Sparse commitment (commit only, no vector construction)
      c.bench_function(&format!("sparse_commit_{size}_d{density}"), |b| {
        b.iter(|| {
          black_box(<E as Engine>::CE::commit_sparse_binary(
            &ck,
            &non_zero_indices,
            &zero,
          ))
        })
      });

      // Benchmark 3: Dense commitment (vector construction + commit using 1-bit optimized path)
      c.bench_function(&format!("dense_commit_{size}_d{density}"), |b| {
        b.iter(|| {
          let dense_vec = make_dense_vector(size, &non_zero_indices);
          black_box(<E as Engine>::CE::commit_small_range(
            &ck,
            &dense_vec,
            &zero,
            0..size,
            1, // max_num_bits = 1 for boolean values
          ))
        })
      });
    }

    // Double twice the size for next iteration
    size *= 4;
  }
}
