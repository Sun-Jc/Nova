//! Benchmark: batch_sparse_binary_msm vs N independent batch_add calls
//!
//! Simulates the one-hot-per-slot structure: 2^20 bases, slot_size=16,
//! N instances each selecting one base per slot.
use core::time::Duration;
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use halo2curves::{
  bn256::{G1Affine as Affine, G1 as Point},
  group::Group,
};
use nova_snark::provider::msm::{batch_add, batch_sparse_binary_msm};
use rand::Rng;
use rayon::prelude::*;

criterion_group! {
  name = batch_binary;
  config = Criterion::default().warm_up_time(Duration::from_millis(3000)).sample_size(10);
  targets = bench_batch_binary_msm
}

criterion_main!(batch_binary);

fn bench_batch_binary_msm(c: &mut Criterion) {
  let slot_size: usize = 16;

  for log_bases in [16, 20] {
    let num_bases: usize = 1 << log_bases;
    let num_slots = num_bases / slot_size;

    // Generate random bases
    let bases: Vec<Affine> = (0..num_bases)
      .into_par_iter()
      .map(|_| {
        let mut rng = rand::thread_rng();
        Affine::from(Point::random(&mut rng))
      })
      .collect();

    for num_instances in [4, 16, 64, 256] {
      // Generate random one-hot-per-slot index sets
      let index_sets: Vec<Vec<usize>> = (0..num_instances)
        .into_par_iter()
        .map(|_| {
          let mut rng = rand::thread_rng();
          (0..num_slots)
            .map(|slot| slot * slot_size + rng.gen_range(0..slot_size))
            .collect()
        })
        .collect();

      let index_refs: Vec<&[usize]> = index_sets.iter().map(|v| v.as_slice()).collect();

      let label = format!("2^{log_bases}_N{num_instances}");

      // Benchmark: N independent batch_add calls
      c.bench_with_input(
        BenchmarkId::new("independent_batch_add", &label),
        &(&bases, &index_sets),
        |b, &(bases, index_sets)| {
          b.iter(|| {
            let results: Vec<_> = index_sets
              .iter()
              .map(|indices| black_box(batch_add(bases, indices)))
              .collect();
            black_box(results)
          })
        },
      );

      // Benchmark: batch_sparse_binary_msm
      c.bench_with_input(
        BenchmarkId::new("batch_sparse_binary_msm", &label),
        &(&bases, &index_refs),
        |b, &(bases, index_refs)| {
          b.iter(|| black_box(batch_sparse_binary_msm(bases, index_refs, slot_size)))
        },
      );
    }
  }
}
