//! Measure exact table build time for 2^20 blocks, K=16.
use core::time::Duration;
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use nova_snark::{
  provider::{msm::OneHotPrecompTable, Bn256EngineKZG},
  traits::{commitment::CommitmentEngineTrait, Engine},
};

criterion_group! {
  name = table_build;
  config = Criterion::default().warm_up_time(Duration::from_millis(1000)).sample_size(10);
  targets = bench_table_build
}

criterion_main!(table_build);

fn bench_table_build(c: &mut Criterion) {
  type E = Bn256EngineKZG;
  let k: usize = 16;
  let num_blocks: usize = 1 << 20;
  let total_size = num_blocks * k;

  println!("Setting up {total_size} generators...");
  let ck = <E as Engine>::CE::setup(b"bench_build", total_size).unwrap();
  println!("Setup done. Benchmarking table build for {num_blocks} blocks...");

  let mut group = c.benchmark_group("table_build_full");
  group.bench_function(&format!("build_{num_blocks}_blocks"), |b| {
    b.iter(|| {
      black_box(OneHotPrecompTable::new(
        &ck.ck()[..total_size],
        k,
        num_blocks,
      ))
    });
  });
  group.finish();
}
