//! Benchmark: one-hot batch commitment.
//!
//! Block size K=16, 2^20 blocks (16M elements).
//! Compares sparse_binary baseline vs batch_commit_one_hot at batch sizes 1, 16, 64, 256.
use core::time::Duration;
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use nova_snark::{
  provider::Bn256EngineKZG,
  traits::{commitment::CommitmentEngineTrait, Engine},
};
use rand::Rng;

criterion_group! {
  name = one_hot_commit;
  config = Criterion::default().warm_up_time(Duration::from_millis(3000)).sample_size(10);
  targets = bench_one_hot
}

criterion_main!(one_hot_commit);

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

fn bench_one_hot(c: &mut Criterion) {
  type E = Bn256EngineKZG;
  let k: usize = 16;
  let num_blocks: usize = 1 << 20;
  let total_size = num_blocks * k;
  let zero = <E as Engine>::Scalar::default();

  let ck = <E as Engine>::CE::setup(b"bench_one_hot", total_size).unwrap();

  let mut group = c.benchmark_group(format!("one_hot_K{k}_2^20"));

  // Baseline: single commit via commit_sparse_binary
  let offsets = random_offsets(num_blocks, k);
  let global_indices = offsets_to_global_indices(&offsets, k);
  group.bench_function("sparse_binary_x1", |b| {
    b.iter(|| {
      black_box(<E as Engine>::CE::commit_sparse_binary(
        &ck,
        &global_indices,
        &zero,
      ))
    });
  });

  // batch_commit_one_hot at four batch sizes
  for batch_size in [1, 16, 64, 256] {
    let all_offsets: Vec<Vec<usize>> = (0..batch_size)
      .map(|_| random_offsets(num_blocks, k))
      .collect();
    let all_r = vec![zero; batch_size];

    group.bench_function(&format!("batch_one_hot_x{batch_size}"), |b| {
      b.iter(|| {
        black_box(<E as Engine>::CE::batch_commit_one_hot(
          &ck,
          k,
          &all_offsets,
          &all_r,
        ))
      });
    });
  }

  group.finish();
}
