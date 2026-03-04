//! Benchmark: one-hot commit at the real target scale.
//!
//! K=16, 2^20 blocks (non-zero count), vector length = 2^24.
//! Compares sparse_binary, one_hot, precomp_table (t=2), and precomp_table3 (t=3).
use core::time::Duration;
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use nova_snark::{
  provider::{
    msm::{OneHotPrecompTable, OneHotPrecompTable3},
    Bn256EngineKZG,
  },
  traits::{commitment::CommitmentEngineTrait, Engine},
};
use rand::Rng;

criterion_group! {
  name = one_hot_real;
  config = Criterion::default().warm_up_time(Duration::from_millis(3000)).sample_size(10);
  targets = bench_real_bn256
}

criterion_main!(one_hot_real);

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

fn bench_real_bn256(c: &mut Criterion) {
  type E = Bn256EngineKZG;
  let k: usize = 16;

  // 2^20 blocks × K=16 = 2^24 total elements
  let num_blocks: usize = 1 << 20;
  let total_size = num_blocks * k;

  println!(
    "Setting up commitment key with {total_size} generators ({num_blocks} blocks, K={k})..."
  );
  let ck = <E as Engine>::CE::setup(b"bench_real", total_size).unwrap();
  let zero = <E as Engine>::Scalar::default();
  println!("Setup complete.");

  let offsets = random_offsets(num_blocks, k);
  let global_indices = offsets_to_global_indices(&offsets, k);

  let mut group = c.benchmark_group(format!("real_bn256_K{k}_blocks{num_blocks}"));

  // 1. Baseline: sparse_binary
  group.bench_function("sparse_binary", |b| {
    b.iter(|| {
      black_box(<E as Engine>::CE::commit_sparse_binary(
        &ck,
        &global_indices,
        &zero,
      ))
    });
  });

  // 2. One-hot with prefetch
  group.bench_function("one_hot", |b| {
    b.iter(|| black_box(<E as Engine>::CE::commit_one_hot(&ck, k, &offsets, &zero)));
  });

  // 3. Precomputed table t=2 (build outside benchmark loop)
  println!("Building precomputed table t=2 for {num_blocks} blocks...");
  let table2 = OneHotPrecompTable::new(&ck.ck()[..total_size], k, num_blocks);
  println!("Table t=2 built.");

  group.bench_function("precomp_t2", |b| {
    b.iter(|| black_box(table2.batch_add(&offsets)));
  });

  // 4. Precomputed table t=3 (build outside benchmark loop)
  println!("Building precomputed table t=3 for {num_blocks} blocks...");
  let table3 = OneHotPrecompTable3::new(&ck.ck()[..total_size], k, num_blocks);
  println!("Table t=3 built.");

  group.bench_function("precomp_t3", |b| {
    b.iter(|| black_box(table3.batch_add(&offsets)));
  });

  group.finish();
}
