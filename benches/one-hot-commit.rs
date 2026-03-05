//! Benchmarks for one-hot sparse binary commitment.
//!
//! Compares three approaches:
//! 1. `commit_sparse_binary` — existing generic sparse binary path
//! 2. `commit_one_hot` — new one-hot-aware path with prefetch
//! 3. `OneHotPrecompTable` — precomputed pair-table lookup (amortized over many commits)
//!
//! All benchmarks use K=16 (block size) on BN254 and Pallas curves.
use core::time::Duration;
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use nova_snark::{
  provider::{msm::OneHotPrecompTable, Bn256EngineKZG, PallasEngine},
  traits::{commitment::CommitmentEngineTrait, Engine},
};
use rand::Rng;

criterion_group! {
  name = one_hot;
  config = Criterion::default().warm_up_time(Duration::from_millis(2000)).sample_size(10);
  targets = bench_one_hot_bn256, bench_one_hot_pallas
}

criterion_main!(one_hot);

/// Generate random one-hot offsets for the given number of blocks and block size.
fn random_offsets(num_blocks: usize, block_size: usize) -> Vec<usize> {
  let mut rng = rand::thread_rng();
  (0..num_blocks)
    .map(|_| rng.gen::<usize>() % block_size)
    .collect()
}

/// Convert one-hot offsets to global indices (for use with commit_sparse_binary).
fn offsets_to_global_indices(offsets: &[usize], block_size: usize) -> Vec<usize> {
  offsets
    .iter()
    .enumerate()
    .map(|(i, &offset)| i * block_size + offset)
    .collect()
}

fn bench_one_hot_bn256(c: &mut Criterion) {
  type E = Bn256EngineKZG;
  let k: usize = 16;
  let label = "bn256";

  let sizes: Vec<usize> = vec![1 << 12, 1 << 14, 1 << 16, 1 << 18, 1 << 20];
  let max_size = *sizes.last().unwrap();

  let ck = <E as Engine>::CE::setup(b"bench_one_hot", max_size).unwrap();
  let zero = <E as Engine>::Scalar::default();

  let max_blocks = max_size / k;
  let offsets = random_offsets(max_blocks, k);
  let global_indices = offsets_to_global_indices(&offsets, k);

  let mut group = c.benchmark_group(format!("one_hot_{label}_K{k}"));

  for &size in &sizes {
    let num_blocks = size / k;

    group.bench_with_input(BenchmarkId::new("sparse_binary", size), &size, |b, _| {
      let indices: Vec<usize> = global_indices[..num_blocks].to_vec();
      b.iter(|| {
        black_box(<E as Engine>::CE::commit_sparse_binary(
          &ck, &indices, &zero,
        ))
      });
    });

    group.bench_with_input(BenchmarkId::new("one_hot", size), &size, |b, _| {
      let offs = &offsets[..num_blocks];
      b.iter(|| black_box(<E as Engine>::CE::commit_one_hot(&ck, k, offs, &zero)));
    });

    group.bench_with_input(BenchmarkId::new("precomp_table", size), &size, |b, _| {
      let table = OneHotPrecompTable::new(&ck.ck()[..size], k, num_blocks);
      let offs = &offsets[..num_blocks];
      b.iter(|| black_box(table.batch_add(offs)));
    });
  }

  group.finish();

  // Table build time benchmark
  let mut build_group = c.benchmark_group(format!("one_hot_{label}_K{k}_table_build"));
  for &size in &sizes {
    let num_blocks = size / k;
    build_group.bench_with_input(BenchmarkId::new("build_table", size), &size, |b, _| {
      b.iter(|| black_box(OneHotPrecompTable::new(&ck.ck()[..size], k, num_blocks)));
    });
  }
  build_group.finish();
}

fn bench_one_hot_pallas(c: &mut Criterion) {
  type E = PallasEngine;
  let k: usize = 16;
  let label = "pallas";

  let sizes: Vec<usize> = vec![1 << 12, 1 << 14, 1 << 16, 1 << 18, 1 << 20];
  let max_size = *sizes.last().unwrap();

  let ck = <E as Engine>::CE::setup(b"bench_one_hot", max_size).unwrap();
  let zero = <E as Engine>::Scalar::default();

  let max_blocks = max_size / k;
  let offsets = random_offsets(max_blocks, k);
  let global_indices = offsets_to_global_indices(&offsets, k);

  let mut group = c.benchmark_group(format!("one_hot_{label}_K{k}"));

  for &size in &sizes {
    let num_blocks = size / k;

    group.bench_with_input(BenchmarkId::new("sparse_binary", size), &size, |b, _| {
      let indices: Vec<usize> = global_indices[..num_blocks].to_vec();
      b.iter(|| {
        black_box(<E as Engine>::CE::commit_sparse_binary(
          &ck, &indices, &zero,
        ))
      });
    });

    group.bench_with_input(BenchmarkId::new("one_hot", size), &size, |b, _| {
      let offs = &offsets[..num_blocks];
      b.iter(|| black_box(<E as Engine>::CE::commit_one_hot(&ck, k, offs, &zero)));
    });

    group.bench_with_input(BenchmarkId::new("precomp_table", size), &size, |b, _| {
      let table = OneHotPrecompTable::new(&ck.ck()[..size], k, num_blocks);
      let offs = &offsets[..num_blocks];
      b.iter(|| black_box(table.batch_add(offs)));
    });
  }

  group.finish();

  // Table build time benchmark
  let mut build_group = c.benchmark_group(format!("one_hot_{label}_K{k}_table_build"));
  for &size in &sizes {
    let num_blocks = size / k;
    build_group.bench_with_input(BenchmarkId::new("build_table", size), &size, |b, _| {
      b.iter(|| black_box(OneHotPrecompTable::new(&ck.ck()[..size], k, num_blocks)));
    });
  }
  build_group.finish();
}
