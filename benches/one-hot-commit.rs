//! Benchmark: one-hot batch commitment.
//!
//! Block size K=16, varying scales (2^16, 2^18, 2^20 blocks).
//! Compares commit-major, block-major, and precomp+block-major
//! at batch sizes 1, 16, 64, 256.
use core::time::Duration;
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use nova_snark::{
  provider::{
    msm::{batch_add_one_hot_block, batch_add_precomp_block, OneHotPrecompTable},
    Bn256EngineKZG,
  },
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

fn bench_one_hot(c: &mut Criterion) {
  type E = Bn256EngineKZG;
  let k: usize = 16;

  for log_blocks in [16, 18, 20] {
    let num_blocks: usize = 1 << log_blocks;
    let total_size = num_blocks * k;
    let zero = <E as Engine>::Scalar::default();

    let ck = <E as Engine>::CE::setup(b"bench_one_hot", total_size).unwrap();

    // Build precomp table (outside bench loop)
    let table = OneHotPrecompTable::new(ck.ck(), k, num_blocks);

    let mut group = c.benchmark_group(format!("one_hot_K{k}_2^{log_blocks}"));

    for batch_size in [1, 16, 64, 256] {
      let all_offsets: Vec<Vec<usize>> = (0..batch_size)
        .map(|_| random_offsets(num_blocks, k))
        .collect();
      let all_r = vec![zero; batch_size];

      // Commit-major (baseline)
      group.bench_function(format!("commit_major_x{batch_size}"), |b| {
        b.iter(|| {
          black_box(<E as Engine>::CE::batch_commit_one_hot(
            &ck,
            k,
            &all_offsets,
            &all_r,
          ))
        });
      });

      // Block-major
      group.bench_function(format!("block_major_x{batch_size}"), |b| {
        b.iter(|| {
          black_box(batch_add_one_hot_block(ck.ck(), k, &all_offsets))
        });
      });

      // Precomp + block-major
      group.bench_function(format!("precomp_block_x{batch_size}"), |b| {
        b.iter(|| {
          black_box(batch_add_precomp_block(&table, &all_offsets))
        });
      });
    }

    group.finish();
  }
}
