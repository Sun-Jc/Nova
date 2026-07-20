//! Level-1 benchmark: projective vs eval-basis prover cost for the ppSNARK
//! **outer sumcheck**, on identical data at matched sizes.
//!
//! The eval-basis side runs `SumcheckProof::prove_cubic_with_three_inputs`
//! (the actual ppSNARK outer prover, Gruen + BDDT optimized). The projective
//! side builds the coefficient-basis outer virtual polynomial
//! `Eq^∞_τ·(Az·Bz − u·Cz·U − E·U)` and runs the factorized projective prover.
//! Both prove the same relation over `{0,1}^n` vs `{0,∞}^n`; this measures the
//! per-sumcheck prover speedup before end-to-end wiring.
//!
//! Run: `cargo criterion --bench ppsnark_projective`.
#![allow(non_snake_case)]

use criterion::*;
use ff::Field;
use nova_snark::{
  provider::Bn256EngineKZG,
  spartan::{
    polys::multilinear::MultilinearPolynomial, ppsnark_projective::build_outer,
    sumcheck::SumcheckProof,
  },
  traits::{Engine, TranscriptEngineTrait},
};
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::hint::black_box;
use std::time::Duration;

type E = Bn256EngineKZG;
type Fr = <E as Engine>::Scalar;

criterion_group! {
  name = ppsnark_projective;
  config = Criterion::default().warm_up_time(Duration::from_millis(2000));
  targets = bench_outer_prover
}
criterion_main!(ppsnark_projective);

fn bench_outer_prover(c: &mut Criterion) {
  for &num_vars in [12usize, 14, 16, 18, 20].iter() {
    let n = 1usize << num_vars;

    // Shared data.
    let az: Vec<Fr> = (0..n)
      .into_par_iter()
      .map(|i| Fr::from((2 * i + 1) as u64))
      .collect();
    let bz: Vec<Fr> = (0..n)
      .into_par_iter()
      .map(|i| Fr::from((3 * i + 4) as u64))
      .collect();
    let cz: Vec<Fr> = (0..n)
      .into_par_iter()
      .map(|i| Fr::from((i + 2) as u64))
      .collect();
    let u = Fr::from(7);
    // Satisfying E = Az∘Bz − u·Cz.
    let e: Vec<Fr> = (0..n)
      .into_par_iter()
      .map(|i| az[i] * bz[i] - u * cz[i])
      .collect();
    let tau: Vec<Fr> = (0..num_vars)
      .map(|i| Fr::from((5 * i + 3) as u64))
      .collect();

    let mut group = c.benchmark_group(format!("outer-sumcheck-vars-{num_vars}"));
    group.sample_size(10);

    // Eval-basis: prove_cubic_with_three_inputs over eq(τ)·(Az·Bz − uCz_E).
    group.bench_function("eval_basis", |b| {
      b.iter(|| {
        let uCz_E: Vec<Fr> = (0..n).map(|i| u * cz[i] + e[i]).collect();
        let mut pa = MultilinearPolynomial::new(az.clone());
        let mut pb = MultilinearPolynomial::new(bz.clone());
        let mut pc = MultilinearPolynomial::new(uCz_E);
        let mut ts = <E as Engine>::TE::new(b"bench");
        let res = SumcheckProof::<E>::prove_cubic_with_three_inputs(
          black_box(&Fr::ZERO),
          black_box(tau.clone()),
          &mut pa,
          &mut pb,
          &mut pc,
          &mut ts,
        );
        assert!(res.is_ok());
      })
    });

    // Projective: factorized prover on the coeff-basis outer virtual polynomial.
    group.bench_function("projective", |b| {
      b.iter(|| {
        let vp = build_outer::<E>(
          num_vars,
          az.clone(),
          bz.clone(),
          cz.clone(),
          e.clone(),
          u,
          &tau,
        );
        let mut ts = <E as Engine>::TE::new(b"bench");
        let out = vp.prove(&mut ts);
        black_box(out.final_claim);
      })
    });

    group.finish();
  }
}
