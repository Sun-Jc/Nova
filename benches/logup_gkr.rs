//! Benchmark for the Logup-GKR tree construction (`Layer::build_tree`, which
//! repeatedly calls the parallelized `fold_up`) over BN254 scalars.
//!
//! Run: `cargo bench --bench logup_gkr`.
use core::time::Duration;
use criterion::{criterion_group, criterion_main, BatchSize, BenchmarkId, Criterion};
use ff::Field;
use nova_snark::{
  provider::Bn256EngineKZG,
  spartan::{logup_gkr::layer::Layer, polys::multilinear::MultilinearPolynomial},
  traits::Engine,
};
use rand::rngs::StdRng;
use rand::SeedableRng;
use rayon::prelude::*;

criterion_group! {
  name = logup_gkr;
  config = Criterion::default().warm_up_time(Duration::from_millis(1000)).sample_size(20);
  targets = bench_build_tree
}
criterion_main!(logup_gkr);

fn bench_build_tree(c: &mut Criterion) {
  type E = Bn256EngineKZG;
  type Fr = <E as Engine>::Scalar;

  let mut group = c.benchmark_group("logup-gkr-build-tree");

  for &log_n in &[16usize, 18, 20, 22] {
    let n = 1usize << log_n;

    // Random, invertible-safe input layer: num arbitrary, den nonzero.
    let (num, den): (Vec<Fr>, Vec<Fr>) = (0..n)
      .into_par_iter()
      .map(|i| {
        let mut rng = StdRng::seed_from_u64(i as u64);
        (Fr::random(&mut rng), Fr::random(&mut rng) + Fr::ONE)
      })
      .unzip();

    group.throughput(criterion::Throughput::Elements(n as u64));
    group.bench_with_input(BenchmarkId::from_parameter(log_n), &log_n, |b, _| {
      b.iter_batched(
        || Layer::<E> {
          num: MultilinearPolynomial::new(num.clone()),
          den: MultilinearPolynomial::new(den.clone()),
        },
        |layer| layer.build_tree(),
        BatchSize::SmallInput,
      )
    });
  }
  group.finish();
}
