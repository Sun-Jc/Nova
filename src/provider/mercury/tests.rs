use std::time::Instant;

use ff::Field;
use rand_core::OsRng;
use rayon::iter::{IntoParallelIterator, ParallelIterator};

use crate::provider::{hyperkzg, mercury};
use crate::spartan::polys::multilinear::MultilinearPolynomial;
use crate::traits::commitment::CommitmentEngineTrait;
use crate::traits::evaluation::EvaluationEngineTrait;
use crate::traits::{Engine, TranscriptEngineTrait};
use crate::{provider::Bn256EngineKZG, spartan::polys::univariate::UniPoly};

type F = halo2curves::bn256::Fr;
type E = Bn256EngineKZG;
type EE = mercury::EvaluationEngine<E>;

#[test]
fn test_mercury_ee() {
  let log_n = 15;
  let n = 1 << log_n;
  let poly = UniPoly {
    coeffs: (0..n)
      .into_par_iter()
      .map(|_| F::random(OsRng))
      .collect::<Vec<_>>(),
  };
  let point = (0..log_n).map(|_| F::random(OsRng)).collect::<Vec<_>>();

  let ck = <<E as Engine>::CE as CommitmentEngineTrait<E>>::CommitmentKey::setup_from_rng(
    b"test", n, OsRng,
  );

  let (pk, vk) = EE::setup(&ck);

  let eval = MultilinearPolynomial::new(poly.coeffs.clone()).evaluate(&point);

  let mut transcript = <E as Engine>::TE::new(b"test");

  let comm = <E as Engine>::CE::commit(&ck, &poly.coeffs, &F::ZERO);

  let start = Instant::now();
  let arg = EE::prove(
    &ck,
    &pk,
    &mut transcript,
    &comm,
    &poly.coeffs,
    &point,
    &eval,
  )
  .unwrap();
  let dur = start.elapsed();

  println!("Mercury: {:?}", &dur);

  {
    let (pk, vk) = hyperkzg::EvaluationEngine::setup(&ck);

    let mut transcript = <E as Engine>::TE::new(b"test");

    let start = Instant::now();
    let arg = hyperkzg::EvaluationEngine::<E>::prove(
      &ck,
      &pk,
      &mut transcript,
      &comm,
      &poly.coeffs,
      &point,
      &eval,
    )
    .unwrap();
    let dur = start.elapsed();

    println!("HyperKZG: {:?}", &dur);

    let mut transcript = <E as Engine>::TE::new(b"test");

    hyperkzg::EvaluationEngine::<E>::verify(&vk, &mut transcript, &comm, &point, &eval, &arg)
      .unwrap();
  }

  let mut transcript = <E as Engine>::TE::new(b"test");

  EE::verify(&vk, &mut transcript, &comm, &point, &eval, &arg).unwrap();
}
