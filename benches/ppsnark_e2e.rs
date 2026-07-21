//! Level-2 end-to-end benchmark: **eval-basis ppSNARK vs projective ppSNARK**
//! on identical workloads, via `DirectSNARK` over a `NonTrivialCircuit`.
//!
//! - eval: `spartan::ppsnark::RelaxedR1CSSNARK` with HyperKZG.
//! - projective: `spartan::ppsnark_projective::RelaxedR1CSSNARKProjective` with
//!   the coefficient-form HyperKZG adapter (`HyperKZGCoeffAdapter`).
//!
//! Both prove/verify the same statement; this measures the whole-SNARK Prove
//! and Verify cost side by side.
//!
//! Requires the `test-utils` feature (HyperKZG `setup` is disabled in
//! production builds; `test-utils` enables the insecure random-tau setup):
//!   `cargo bench --features test-utils --bench ppsnark_e2e`
#![allow(non_snake_case)]

use criterion::*;
use nova_snark::{
  provider::{coeff_eval_adapter::HyperKZGCoeffAdapter, Bn256EngineKZG},
  spartan::{
    direct::DirectSNARK, ppsnark::RelaxedR1CSSNARK, ppsnark_projective::RelaxedR1CSSNARKProjective,
  },
  traits::{circuit::NonTrivialCircuit, Engine},
};
use std::hint::black_box;
use std::time::Duration;

type E = Bn256EngineKZG;
type Fr = <E as Engine>::Scalar;

// eval-basis ppSNARK with plain HyperKZG.
type EvalEE = nova_snark::provider::hyperkzg::EvaluationEngine<E>;
type EvalSNARK = RelaxedR1CSSNARK<E, EvalEE>;

// projective ppSNARK with the coefficient-form HyperKZG adapter.
type ProjEE = HyperKZGCoeffAdapter<E>;
type ProjSNARK = RelaxedR1CSSNARKProjective<E, ProjEE>;

criterion_group! {
  name = ppsnark_e2e;
  config = Criterion::default().warm_up_time(Duration::from_millis(3000));
  targets = bench_e2e
}
criterion_main!(ppsnark_e2e);

const NUM_SAMPLES: usize = 10;

/// NonTrivialCircuit(num_cons) computes z -> z^(2^num_cons) by repeated squaring.
fn circuit_output(num_cons: usize, z0: Fr) -> Fr {
  let mut x = z0;
  for _ in 0..num_cons {
    x = x * x;
  }
  x
}

fn bench_e2e(c: &mut Criterion) {
  for &num_cons in [8192usize, 16384, 32768, 65536, 131072, 262144].iter() {
    let mut group = c.benchmark_group(format!("ppsnark-e2e-CircuitSize-{num_cons}"));
    group.sample_size(NUM_SAMPLES);

    let circuit = NonTrivialCircuit::<Fr>::new(num_cons);
    let z0 = Fr::from(42);
    let input = vec![z0];
    let io = vec![z0, circuit_output(num_cons, z0)];

    // --- eval-basis ppSNARK ---
    let (eval_pk, eval_vk) =
      DirectSNARK::<E, EvalSNARK, NonTrivialCircuit<Fr>>::setup(circuit.clone()).unwrap();
    group.bench_function("eval/Prove", |b| {
      b.iter(|| {
        let res = DirectSNARK::<E, EvalSNARK, _>::prove(
          black_box(&eval_pk),
          black_box(circuit.clone()),
          black_box(&input),
        );
        assert!(res.is_ok());
      })
    });
    let eval_snark =
      DirectSNARK::<E, EvalSNARK, _>::prove(&eval_pk, circuit.clone(), &input).unwrap();
    group.bench_function("eval/Verify", |b| {
      b.iter(|| {
        assert!(eval_snark
          .verify(black_box(&eval_vk), black_box(&io))
          .is_ok())
      });
    });

    // --- projective ppSNARK ---
    let (proj_pk, proj_vk) =
      DirectSNARK::<E, ProjSNARK, NonTrivialCircuit<Fr>>::setup(circuit.clone()).unwrap();
    group.bench_function("projective/Prove", |b| {
      b.iter(|| {
        let res = DirectSNARK::<E, ProjSNARK, _>::prove(
          black_box(&proj_pk),
          black_box(circuit.clone()),
          black_box(&input),
        );
        assert!(res.is_ok());
      })
    });
    let proj_snark =
      DirectSNARK::<E, ProjSNARK, _>::prove(&proj_pk, circuit.clone(), &input).unwrap();
    group.bench_function("projective/Verify", |b| {
      b.iter(|| {
        assert!(proj_snark
          .verify(black_box(&proj_vk), black_box(&io))
          .is_ok())
      });
    });

    group.finish();
  }
}
