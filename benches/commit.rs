//! Benchmarking the commit times for Pallas (Pasta) curve using
//! halo2curves library and the Nova-provided commitment engine, on a range of scalar bit-widths
use core::time::Duration;
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use halo2curves::{
  ff::Field,
  group::{Curve, Group},
  msm::msm_best,
  pasta::{Fq as Scalar, Pallas},
};
use nova_snark::{
  provider::PallasEngine,
  traits::{commitment::CommitmentEngineTrait, Engine},
};
use rand::Rng;
use rayon::prelude::*;

criterion_group! {
name = commit;
config = Criterion::default().warm_up_time(Duration::from_millis(3000)).sample_size(10);
targets = bench_commit
}

criterion_main!(commit);

fn bench_commit(c: &mut Criterion) {
  type E = PallasEngine;

  let sizes: Vec<usize> = vec![8, 1024, 1 << 15, 1 << 20, 1 << 24];
  let max = *sizes.last().unwrap();

  // sample bases for the purpose of testing
  let ck = <E as Engine>::CE::setup(b"test_from_label", max).unwrap();

  let zero = <E as Engine>::Scalar::zero();

  // generate random affine bases for msm_best baseline
  // (Pedersen CommitmentKey does not expose bases publicly)
  let bases_affine = (0..max)
    .into_par_iter()
    .map(|_| {
      let mut rng = rand::thread_rng();
      Pallas::random(&mut rng).to_affine()
    })
    .collect::<Vec<_>>();

  // random scalars that are in the set {0, 1}
  let scalars_u1 = (0..max)
    .into_par_iter()
    .map(|_| {
      let mut rng = rand::thread_rng();
      rng.gen::<u16>() % 2
    })
    .collect::<Vec<_>>();

  let scalars_u1_field = scalars_u1
    .iter()
    .map(|&x| Scalar::from(x as u64))
    .collect::<Vec<_>>();

  // pre-compute non-zero indices for commit_sparse_binary benchmark
  let non_zero_indices_u1: Vec<usize> = scalars_u1
    .iter()
    .enumerate()
    .filter(|(_, &v)| v != 0)
    .map(|(i, _)| i)
    .collect();

  // 10-bit scalars that are in the set {0, ..., 2^10-1}
  let scalars_u10 = (0..max)
    .into_par_iter()
    .map(|_| {
      let mut rng = rand::thread_rng();
      rng.gen::<u16>() % (1 << 10)
    })
    .collect::<Vec<_>>();

  let scalars_u10_field = scalars_u10
    .iter()
    .map(|&x| Scalar::from(x as u64))
    .collect::<Vec<_>>();

  // random scalars that are in the set {0, ..., 2^16-1}
  let scalars_u16 = (0..max)
    .into_par_iter()
    .map(|_| {
      let mut rng = rand::thread_rng();
      rng.gen::<u16>()
    })
    .collect::<Vec<_>>();

  let scalars_u16_field = scalars_u16
    .iter()
    .map(|&x| Scalar::from(x as u64))
    .collect::<Vec<_>>();

  // random scalars that are in the set {0, ..., 2^32-1}
  let scalars_u32 = (0..max)
    .into_par_iter()
    .map(|_| {
      let mut rng = rand::thread_rng();
      rng.gen::<u32>()
    })
    .collect::<Vec<_>>();

  let scalars_u32_field = scalars_u32
    .iter()
    .map(|&x| Scalar::from(x as u64))
    .collect::<Vec<_>>();

  // random scalars that are in the set {0, ..., 2^64-1}
  let scalars_u64 = (0..max)
    .into_par_iter()
    .map(|_| {
      let mut rng = rand::thread_rng();
      rng.gen::<u64>()
    })
    .collect::<Vec<_>>();

  let scalars_u64_field = scalars_u64
    .iter()
    .map(|&x| Scalar::from(x))
    .collect::<Vec<_>>();

  // random scalars in the set {0, ..., p-1}, where p is the modulus for the
  // scalar field of Pallas
  let scalars_random_field = (0..max)
    .into_par_iter()
    .map(|_| {
      let mut rng = rand::thread_rng();
      Scalar::random(&mut rng)
    })
    .collect::<Vec<_>>();

  for &size in &sizes {
    c.bench_function(&format!("halo2curves_commit_u1_{size}"), |b| {
      b.iter(|| black_box(msm_best(&scalars_u1_field[..size], &bases_affine[..size])))
    });

    c.bench_function(&format!("nova_generic_commit_u1_{size}"), |b| {
      b.iter(|| {
        black_box(<E as Engine>::CE::commit(
          &ck,
          &scalars_u1_field[..size],
          &zero,
        ))
      })
    });

    c.bench_function(&format!("nova_specialized_commit_u1_{size}"), |b| {
      b.iter(|| {
        black_box(<E as Engine>::CE::commit_small(
          &ck,
          &scalars_u1[..size],
          &zero,
        ))
      })
    });

    let nz_count = non_zero_indices_u1.partition_point(|&i| i < size);
    c.bench_function(&format!("nova_batch_add_commit_u1_{size}"), |b| {
      b.iter(|| {
        black_box(<E as Engine>::CE::commit_sparse_binary(
          &ck,
          &non_zero_indices_u1[..nz_count],
          &zero,
        ))
      })
    });

    c.bench_function(&format!("halo2curves_commit_u10_{size}"), |b| {
      b.iter(|| black_box(msm_best(&scalars_u10_field[..size], &bases_affine[..size])))
    });

    c.bench_function(&format!("nova_generic_commit_u10_{size}"), |b| {
      b.iter(|| {
        black_box(<E as Engine>::CE::commit(
          &ck,
          &scalars_u10_field[..size],
          &zero,
        ))
      })
    });

    c.bench_function(&format!("nova_specialized_commit_u10_{size}"), |b| {
      b.iter(|| {
        black_box(<E as Engine>::CE::commit_small(
          &ck,
          &scalars_u10[..size],
          &zero,
        ))
      })
    });

    c.bench_function(&format!("halo2curves_commit_u16_{size}"), |b| {
      b.iter(|| black_box(msm_best(&scalars_u16_field[..size], &bases_affine[..size])))
    });

    c.bench_function(&format!("nova_generic_commit_u16_{size}"), |b| {
      b.iter(|| {
        black_box(<E as Engine>::CE::commit(
          &ck,
          &scalars_u16_field[..size],
          &zero,
        ))
      })
    });

    c.bench_function(&format!("nova_specialized_commit_u16_{size}"), |b| {
      b.iter(|| {
        black_box(<E as Engine>::CE::commit_small(
          &ck,
          &scalars_u16[..size],
          &zero,
        ))
      })
    });

    c.bench_function(&format!("halo2curves_commit_u32_{size}"), |b| {
      b.iter(|| black_box(msm_best(&scalars_u32_field[..size], &bases_affine[..size])))
    });

    c.bench_function(&format!("nova_generic_commit_u32_{size}"), |b| {
      b.iter(|| {
        black_box(<E as Engine>::CE::commit(
          &ck,
          &scalars_u32_field[..size],
          &zero,
        ))
      })
    });

    c.bench_function(&format!("nova_specialized_commit_u32_{size}"), |b| {
      b.iter(|| {
        black_box(<E as Engine>::CE::commit_small(
          &ck,
          &scalars_u32[..size],
          &zero,
        ))
      })
    });

    c.bench_function(&format!("halo2curves_commit_u64_{size}"), |b| {
      b.iter(|| black_box(msm_best(&scalars_u64_field[..size], &bases_affine[..size])))
    });

    c.bench_function(&format!("nova_generic_commit_u64_{size}"), |b| {
      b.iter(|| {
        black_box(<E as Engine>::CE::commit(
          &ck,
          &scalars_u64_field[..size],
          &zero,
        ))
      })
    });

    c.bench_function(&format!("nova_specialized_commit_u64_{size}"), |b| {
      b.iter(|| {
        black_box(<E as Engine>::CE::commit_small(
          &ck,
          &scalars_u64[..size],
          &zero,
        ))
      })
    });

    c.bench_function(&format!("halo2curves_commit_random_{size}"), |b| {
      b.iter(|| {
        black_box(msm_best(
          &scalars_random_field[..size],
          &bases_affine[..size],
        ))
      })
    });

    c.bench_function(&format!("nova_generic_commit_random_{size}"), |b| {
      b.iter(|| {
        black_box(<E as Engine>::CE::commit(
          &ck,
          &scalars_random_field[..size],
          &zero,
        ))
      })
    });
  }
}
