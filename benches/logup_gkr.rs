//! Benchmark for the Logup-GKR tree construction (`Layer::build_tree`, which
//! repeatedly calls the parallelized `fold_up`) over BN254 scalars.
//!
//! Run: `cargo bench --bench logup_gkr`.
use core::time::Duration;
use criterion::{black_box, criterion_group, criterion_main, BatchSize, BenchmarkId, Criterion};
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
  targets = bench_build_tree, bench_fold_crossover
}
criterion_main!(logup_gkr);

// One fraction-add fold of a length-`2n` layer into length-`n`, serial.
fn fold_serial<Fr: Field>(num: &[Fr], den: &[Fr], n: usize) -> (Vec<Fr>, Vec<Fr>) {
  (0..n)
    .map(|i| {
      let d = den[i] * den[i + n];
      let p = num[i] * den[i + n] + num[i + n] * den[i];
      (p, d)
    })
    .unzip()
}

// Same fold, parallel.
fn fold_parallel<Fr: Field>(num: &[Fr], den: &[Fr], n: usize) -> (Vec<Fr>, Vec<Fr>) {
  (0..n)
    .into_par_iter()
    .map(|i| {
      let d = den[i] * den[i + n];
      let p = num[i] * den[i + n] + num[i + n] * den[i];
      (p, d)
    })
    .unzip()
}

// Sweep a single fold across sizes bracketing PARALLEL_THRESHOLD (=4096),
// serial vs parallel, so the empirical crossover = the optimal threshold.
fn bench_fold_crossover(c: &mut Criterion) {
  type Fr = <Bn256EngineKZG as Engine>::Scalar;
  let mut group = c.benchmark_group("logup-gkr-fold-crossover");

  for &n in &[8192usize, 32768, 65536, 131072, 262144, 524288] {
    let len = 2 * n;
    let (num, den): (Vec<Fr>, Vec<Fr>) = (0..len)
      .into_par_iter()
      .map(|i| {
        let mut rng = StdRng::seed_from_u64(i as u64);
        (Fr::random(&mut rng), Fr::random(&mut rng) + Fr::ONE)
      })
      .unzip();

    group.bench_with_input(BenchmarkId::new("serial", n), &n, |b, &n| {
      b.iter(|| fold_serial(black_box(&num), black_box(&den), n))
    });
    group.bench_with_input(BenchmarkId::new("parallel", n), &n, |b, &n| {
      b.iter(|| fold_parallel(black_box(&num), black_box(&den), n))
    });
  }
  group.finish();
}

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
