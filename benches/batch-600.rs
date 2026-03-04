//! Benchmark: 600 one-hot commits at real scale.
//!
//! K=16, 2^20 blocks, vector length = 2^24.
//! Compares:
//! 1. Sequential sparse_binary (baseline)
//! 2. Parallel sparse_binary (par_iter over 600 commits)
//! 3. Parallel transpose (block-major order)
//! 4. Parallel transpose + precomp table (t=2)
use core::time::Duration;
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use nova_snark::{
  provider::{
    msm::{batch_add_one_hot_precomp_transpose, batch_add_one_hot_transpose, OneHotPrecompTable},
    Bn256EngineKZG,
  },
  traits::{commitment::CommitmentEngineTrait, Engine},
};
use rand::Rng;
use rayon::prelude::*;

criterion_group! {
  name = batch600;
  config = Criterion::default().warm_up_time(Duration::from_millis(3000)).sample_size(10);
  targets = bench_batch_600
}

criterion_main!(batch600);

fn random_offsets(num_blocks: usize, block_size: usize) -> Vec<usize> {
  let mut rng = rand::thread_rng();
  (0..num_blocks)
    .map(|_| rng.gen::<usize>() % block_size)
    .collect()
}

fn offsets_to_global_indices(offsets: &[usize], block_size: usize) -> Vec<usize> {
  offsets
    .iter()
    .enumerate()
    .map(|(i, &offset)| i * block_size + offset)
    .collect()
}

fn bench_batch_600(c: &mut Criterion) {
  type E = Bn256EngineKZG;
  let k: usize = 16;
  let num_blocks: usize = 1 << 20;
  let total_size = num_blocks * k;
  let num_commits = 600;

  println!("Setting up {total_size} generators...");
  let ck = <E as Engine>::CE::setup(b"bench_batch", total_size).unwrap();
  let zero = <E as Engine>::Scalar::default();
  println!("Setup done. Generating {num_commits} random offset vectors...");

  let all_offsets: Vec<Vec<usize>> = (0..num_commits)
    .map(|_| random_offsets(num_blocks, k))
    .collect();

  let all_global_indices: Vec<Vec<usize>> = all_offsets
    .iter()
    .map(|offsets| offsets_to_global_indices(offsets, k))
    .collect();

  println!("Building precomp table t=2...");
  let table = OneHotPrecompTable::new(&ck.ck()[..total_size], k, num_blocks);
  println!("Table built. Starting benchmarks...");

  let mut group = c.benchmark_group(format!("batch_{num_commits}_commits_K{k}"));

  // 1. Sequential: 600 sparse_binary commits one after another
  group.bench_function("sequential_sparse", |b| {
    b.iter(|| {
      for indices in &all_global_indices {
        black_box(<E as Engine>::CE::commit_sparse_binary(&ck, indices, &zero));
      }
    });
  });

  // 2. Parallel: par_iter over 600 commits
  group.bench_function("parallel_sparse", |b| {
    b.iter(|| {
      let results: Vec<_> = all_global_indices
        .par_iter()
        .map(|indices| <E as Engine>::CE::commit_sparse_binary(&ck, indices, &zero))
        .collect();
      black_box(results);
    });
  });

  // 3. Parallel transpose (block-major, no precomp)
  group.bench_function("parallel_transpose", |b| {
    b.iter(|| {
      let refs: Vec<&[usize]> = all_offsets.iter().map(|v| v.as_slice()).collect();
      let results = batch_add_one_hot_transpose(&ck.ck()[..total_size], k, &refs);
      black_box(results);
    });
  });

  // 4. Parallel transpose + precomp table (t=2)
  group.bench_function("precomp_transpose", |b| {
    b.iter(|| {
      let refs: Vec<&[usize]> = all_offsets.iter().map(|v| v.as_slice()).collect();
      let results = batch_add_one_hot_precomp_transpose(&table, &refs);
      black_box(results);
    });
  });

  group.finish();
}
