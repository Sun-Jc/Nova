//! Quick benchmark: one-hot commit at small scale (N≈300, K=16).
//!
//! Tests whether precomputed tables help when there are only ~18 blocks.
use core::time::Duration;
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use nova_snark::{
  provider::{msm::OneHotPrecompTable, Bn256EngineKZG, PallasEngine},
  traits::{commitment::CommitmentEngineTrait, Engine},
};
use rand::Rng;

criterion_group! {
  name = one_hot_small;
  config = Criterion::default().warm_up_time(Duration::from_millis(1000)).sample_size(50);
  targets = bench_small_bn256, bench_small_pallas
}

criterion_main!(one_hot_small);

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

fn bench_small_bn256(c: &mut Criterion) {
  type E = Bn256EngineKZG;
  let k: usize = 16;
  let label = "bn256";

  // Small sizes around N=300
  let sizes: Vec<usize> = vec![
    16 * 4,   // N=64,   4 blocks
    16 * 8,   // N=128,  8 blocks
    16 * 18,  // N=288,  18 blocks (≈300)
    16 * 20,  // N=320,  20 blocks
    16 * 32,  // N=512,  32 blocks
    16 * 64,  // N=1024, 64 blocks
    16 * 128, // N=2048, 128 blocks
  ];

  let max_size = *sizes.last().unwrap();
  let ck = <E as Engine>::CE::setup(b"bench_small", max_size).unwrap();
  let zero = <E as Engine>::Scalar::default();

  let max_blocks = max_size / k;
  let offsets = random_offsets(max_blocks, k);
  let global_indices = offsets_to_global_indices(&offsets, k);

  let mut group = c.benchmark_group(format!("small_{label}_K{k}"));

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
}

fn bench_small_pallas(c: &mut Criterion) {
  type E = PallasEngine;
  let k: usize = 16;
  let label = "pallas";

  let sizes: Vec<usize> = vec![
    16 * 4,   // N=64,   4 blocks
    16 * 8,   // N=128,  8 blocks
    16 * 18,  // N=288,  18 blocks (≈300)
    16 * 20,  // N=320,  20 blocks
    16 * 32,  // N=512,  32 blocks
    16 * 64,  // N=1024, 64 blocks
    16 * 128, // N=2048, 128 blocks
  ];

  let max_size = *sizes.last().unwrap();
  let ck = <E as Engine>::CE::setup(b"bench_small", max_size).unwrap();
  let zero = <E as Engine>::Scalar::default();

  let max_blocks = max_size / k;
  let offsets = random_offsets(max_blocks, k);
  let global_indices = offsets_to_global_indices(&offsets, k);

  let mut group = c.benchmark_group(format!("small_{label}_K{k}"));

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
}
