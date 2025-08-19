use crate::{
  errors::NovaError,
  spartan::polys::{
    multilinear::MultilinearPolynomial,
    univariate::{CompressedUniPoly, UniPoly},
  },
  traits::{Engine, TranscriptEngineTrait},
};
use ff::Field;
use itertools::Itertools as _;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

/// Defines a trait for implementing sum-check in a generic manner
pub trait SumcheckEngine<E: Engine>: Send + Sync {
  /// returns the initial claims
  fn initial_claims(&self) -> Vec<E::Scalar>;

  /// degree of the sum-check polynomial
  fn degree(&self) -> usize;

  /// the size of the polynomials
  fn size(&self) -> usize;

  /// returns evaluation points at 0, 2, d-1 (where d is the degree of the sum-check polynomial)
  fn evaluation_points(&self) -> Vec<Vec<E::Scalar>>;

  /// bounds a variable in the constituent polynomials
  fn bound(&mut self, r: &E::Scalar);

  /// returns the final claims
  fn final_claims(&self) -> Vec<Vec<E::Scalar>>;
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(bound = "")]
pub(crate) struct SumcheckProof<E: Engine> {
  compressed_polys: Vec<CompressedUniPoly<E::Scalar>>,
}

impl<E: Engine> SumcheckProof<E> {
  pub fn new(compressed_polys: Vec<CompressedUniPoly<E::Scalar>>) -> Self {
    Self { compressed_polys }
  }

  pub fn verify(
    &self,
    claim: E::Scalar,
    num_rounds: usize,
    degree_bound: usize,
    transcript: &mut E::TE,
  ) -> Result<(E::Scalar, Vec<E::Scalar>), NovaError> {
    let mut e = claim;
    let mut r: Vec<E::Scalar> = Vec::new();

    // verify that there is a univariate polynomial for each round
    if self.compressed_polys.len() != num_rounds {
      return Err(NovaError::InvalidSumcheckProof);
    }

    for i in 0..self.compressed_polys.len() {
      let poly = self.compressed_polys[i].decompress(&e);

      // verify degree bound
      if poly.degree() != degree_bound {
        return Err(NovaError::InvalidSumcheckProof);
      }

      // we do not need to check if poly(0) + poly(1) = e, as
      // decompress() call above already ensures that holds
      debug_assert_eq!(poly.eval_at_zero() + poly.eval_at_one(), e);

      // append the prover's message to the transcript
      transcript.absorb(b"p", &poly);

      //derive the verifier's challenge for the next round
      let r_i = transcript.squeeze(b"c")?;

      r.push(r_i);

      // evaluate the claimed degree-ell polynomial at r_i
      e = poly.evaluate(&r_i);
    }

    Ok((e, r))
  }

  pub fn verify_batch(
    &self,
    claims: &[E::Scalar],
    num_rounds: &[usize],
    coeffs: &[E::Scalar],
    degree_bound: usize,
    transcript: &mut E::TE,
  ) -> Result<(E::Scalar, Vec<E::Scalar>), NovaError> {
    let num_instances = claims.len();
    assert_eq!(num_rounds.len(), num_instances);
    assert_eq!(coeffs.len(), num_instances);

    // n = maxᵢ{nᵢ}
    let num_rounds_max = *num_rounds.iter().max().unwrap();

    // Random linear combination of claims,
    // where each claim is scaled by 2^{n-nᵢ} to account for the padding.
    //
    // claim = ∑ᵢ coeffᵢ⋅2^{n-nᵢ}⋅cᵢ
    let claim = zip_with!(
      (
        zip_with!(iter, (claims, num_rounds), |claim, num_rounds| {
          let scaling_factor = 1 << (num_rounds_max - num_rounds);
          E::Scalar::from(scaling_factor as u64) * claim
        }),
        coeffs.iter()
      ),
      |scaled_claim, coeff| scaled_claim * coeff
    )
    .sum();

    self.verify(claim, num_rounds_max, degree_bound, transcript)
  }

  #[inline]
  fn compute_eval_points_quad<F>(
    poly_A: &MultilinearPolynomial<E::Scalar>,
    poly_B: &MultilinearPolynomial<E::Scalar>,
    comb_func: &F,
  ) -> (E::Scalar, E::Scalar)
  where
    F: Fn(&E::Scalar, &E::Scalar) -> E::Scalar + Sync,
  {
    let len = poly_A.len() / 2;
    (0..len)
      .into_par_iter()
      .map(|i| {
        // eval 0: bound_func is A(low)
        let eval_point_0 = comb_func(&poly_A[i], &poly_B[i]);

        // eval 2: bound_func is -A(low) + 2*A(high)
        let poly_A_bound_point = poly_A[len + i] + poly_A[len + i] - poly_A[i];
        let poly_B_bound_point = poly_B[len + i] + poly_B[len + i] - poly_B[i];
        let eval_point_2 = comb_func(&poly_A_bound_point, &poly_B_bound_point);
        (eval_point_0, eval_point_2)
      })
      .reduce(
        || (E::Scalar::ZERO, E::Scalar::ZERO),
        |a, b| (a.0 + b.0, a.1 + b.1),
      )
  }

  pub fn prove_quad<F>(
    claim: &E::Scalar,
    num_rounds: usize,
    poly_A: &mut MultilinearPolynomial<E::Scalar>,
    poly_B: &mut MultilinearPolynomial<E::Scalar>,
    comb_func: F,
    transcript: &mut E::TE,
  ) -> Result<(Self, Vec<E::Scalar>, Vec<E::Scalar>), NovaError>
  where
    F: Fn(&E::Scalar, &E::Scalar) -> E::Scalar + Sync,
  {
    let mut r: Vec<E::Scalar> = Vec::new();
    let mut polys: Vec<CompressedUniPoly<E::Scalar>> = Vec::new();
    let mut claim_per_round = *claim;
    for _ in 0..num_rounds {
      let poly = {
        let (eval_point_0, eval_point_2) =
          Self::compute_eval_points_quad(poly_A, poly_B, &comb_func);

        let evals = vec![eval_point_0, claim_per_round - eval_point_0, eval_point_2];
        UniPoly::from_evals(&evals)
      };

      // append the prover's message to the transcript
      transcript.absorb(b"p", &poly);

      //derive the verifier's challenge for the next round
      let r_i = transcript.squeeze(b"c")?;
      r.push(r_i);
      polys.push(poly.compress());

      // Set up next round
      claim_per_round = poly.evaluate(&r_i);

      // bind all tables to the verifier's challenge
      rayon::join(
        || poly_A.bind_poly_var_top(&r_i),
        || poly_B.bind_poly_var_top(&r_i),
      );
    }

    Ok((
      SumcheckProof {
        compressed_polys: polys,
      },
      r,
      vec![poly_A[0], poly_B[0]],
    ))
  }

  pub fn prove_quad_batch<F>(
    claims: &[E::Scalar],
    num_rounds: &[usize],
    mut poly_A_vec: Vec<MultilinearPolynomial<E::Scalar>>,
    mut poly_B_vec: Vec<MultilinearPolynomial<E::Scalar>>,
    coeffs: &[E::Scalar],
    comb_func: F,
    transcript: &mut E::TE,
  ) -> Result<(Self, Vec<E::Scalar>, (Vec<E::Scalar>, Vec<E::Scalar>)), NovaError>
  where
    F: Fn(&E::Scalar, &E::Scalar) -> E::Scalar + Sync,
  {
    let num_claims = claims.len();

    assert_eq!(num_rounds.len(), num_claims);
    assert_eq!(poly_A_vec.len(), num_claims);
    assert_eq!(poly_B_vec.len(), num_claims);
    assert_eq!(coeffs.len(), num_claims);

    for (i, &num_rounds) in num_rounds.iter().enumerate() {
      let expected_size = 1 << num_rounds;

      // Direct indexing with the assumption that the index will always be in bounds
      let a = &poly_A_vec[i];
      let b = &poly_B_vec[i];

      for (l, polyname) in [(a.len(), "poly_A_vec"), (b.len(), "poly_B_vec")].iter() {
        assert_eq!(
          *l, expected_size,
          "Mismatch in size for {polyname} at index {i}"
        );
      }
    }

    let num_rounds_max = *num_rounds.iter().max().unwrap();
    let mut e = zip_with!(
      iter,
      (claims, num_rounds, coeffs),
      |claim, num_rounds, coeff| {
        let scaled_claim = E::Scalar::from((1 << (num_rounds_max - num_rounds)) as u64) * claim;
        scaled_claim * coeff
      }
    )
    .sum();
    let mut r: Vec<E::Scalar> = Vec::new();
    let mut quad_polys: Vec<CompressedUniPoly<E::Scalar>> = Vec::new();

    for current_round in 0..num_rounds_max {
      let remaining_rounds = num_rounds_max - current_round;
      let evals: Vec<(E::Scalar, E::Scalar)> = zip_with!(
        par_iter,
        (num_rounds, claims, poly_A_vec, poly_B_vec),
        |num_rounds, claim, poly_A, poly_B| {
          if remaining_rounds <= *num_rounds {
            Self::compute_eval_points_quad(poly_A, poly_B, &comb_func)
          } else {
            let remaining_variables = remaining_rounds - num_rounds - 1;
            let scaled_claim = E::Scalar::from((1 << remaining_variables) as u64) * claim;
            (scaled_claim, scaled_claim)
          }
        }
      )
      .collect();

      let evals_combined_0 = (0..evals.len()).map(|i| evals[i].0 * coeffs[i]).sum();
      let evals_combined_2 = (0..evals.len()).map(|i| evals[i].1 * coeffs[i]).sum();

      let evals = vec![evals_combined_0, e - evals_combined_0, evals_combined_2];
      let poly = UniPoly::from_evals(&evals);

      // append the prover's message to the transcript
      transcript.absorb(b"p", &poly);

      // derive the verifier's challenge for the next round
      let r_i = transcript.squeeze(b"c")?;
      r.push(r_i);

      // bound all tables to the verifier's challenge
      zip_with_for_each!(
        (
          num_rounds.par_iter(),
          poly_A_vec.par_iter_mut(),
          poly_B_vec.par_iter_mut()
        ),
        |num_rounds, poly_A, poly_B| {
          if remaining_rounds <= *num_rounds {
            let _ = rayon::join(
              || poly_A.bind_poly_var_top(&r_i),
              || poly_B.bind_poly_var_top(&r_i),
            );
          }
        }
      );

      e = poly.evaluate(&r_i);
      quad_polys.push(poly.compress());
    }
    poly_A_vec.iter().for_each(|p| assert_eq!(p.len(), 1));
    poly_B_vec.iter().for_each(|p| assert_eq!(p.len(), 1));

    let poly_A_final = poly_A_vec
      .into_iter()
      .map(|poly| poly[0])
      .collect::<Vec<_>>();
    let poly_B_final = poly_B_vec
      .into_iter()
      .map(|poly| poly[0])
      .collect::<Vec<_>>();

    let eval_expected = zip_with!(
      iter,
      (poly_A_final, poly_B_final, coeffs),
      |eA, eB, coeff| comb_func(eA, eB) * coeff
    )
    .sum::<E::Scalar>();
    assert_eq!(e, eval_expected);

    let claims_prod = (poly_A_final, poly_B_final);

    Ok((SumcheckProof::new(quad_polys), r, claims_prod))
  }

  #[inline]
  pub fn compute_eval_points_cubic<F>(
    poly_A: &MultilinearPolynomial<E::Scalar>,
    poly_B: &MultilinearPolynomial<E::Scalar>,
    poly_C: &MultilinearPolynomial<E::Scalar>,
    comb_func: &F,
  ) -> (E::Scalar, E::Scalar, E::Scalar)
  where
    F: Fn(&E::Scalar, &E::Scalar, &E::Scalar) -> E::Scalar + Sync,
  {
    let len = poly_A.len() / 2;
    (0..len)
      .into_par_iter()
      .map(|i| {
        // eval 0: bound_func is A(low)
        let eval_point_0 = comb_func(&poly_A[i], &poly_B[i], &poly_C[i]);

        // eval 2: bound_func is -A(low) + 2*A(high)
        let poly_A_bound_point = poly_A[len + i] + poly_A[len + i] - poly_A[i];
        let poly_B_bound_point = poly_B[len + i] + poly_B[len + i] - poly_B[i];
        let poly_C_bound_point = poly_C[len + i] + poly_C[len + i] - poly_C[i];
        let eval_point_2 = comb_func(
          &poly_A_bound_point,
          &poly_B_bound_point,
          &poly_C_bound_point,
        );

        // eval 3: bound_func is -2A(low) + 3A(high); computed incrementally with bound_func applied to eval(2)
        let poly_A_bound_point = poly_A_bound_point + poly_A[len + i] - poly_A[i];
        let poly_B_bound_point = poly_B_bound_point + poly_B[len + i] - poly_B[i];
        let poly_C_bound_point = poly_C_bound_point + poly_C[len + i] - poly_C[i];
        let eval_point_3 = comb_func(
          &poly_A_bound_point,
          &poly_B_bound_point,
          &poly_C_bound_point,
        );
        (eval_point_0, eval_point_2, eval_point_3)
      })
      .reduce(
        || (E::Scalar::ZERO, E::Scalar::ZERO, E::Scalar::ZERO),
        |a, b| (a.0 + b.0, a.1 + b.1, a.2 + b.2),
      )
  }

  #[inline]
  pub fn compute_eval_points_cubic_with_additive_term<F>(
    poly_A: &MultilinearPolynomial<E::Scalar>,
    poly_B: &MultilinearPolynomial<E::Scalar>,
    poly_C: &MultilinearPolynomial<E::Scalar>,
    poly_D: &MultilinearPolynomial<E::Scalar>,
    comb_func: &F,
  ) -> (E::Scalar, E::Scalar, E::Scalar)
  where
    F: Fn(&E::Scalar, &E::Scalar, &E::Scalar, &E::Scalar) -> E::Scalar + Sync,
  {
    let len = poly_A.len() / 2;
    (0..len)
      .into_par_iter()
      .map(|i| {
        // eval 0: bound_func is A(low)
        let eval_point_0 = comb_func(&poly_A[i], &poly_B[i], &poly_C[i], &poly_D[i]);

        // eval 2: bound_func is -A(low) + 2*A(high)
        let poly_A_bound_point = poly_A[len + i] + poly_A[len + i] - poly_A[i];
        let poly_B_bound_point = poly_B[len + i] + poly_B[len + i] - poly_B[i];
        let poly_C_bound_point = poly_C[len + i] + poly_C[len + i] - poly_C[i];
        let poly_D_bound_point = poly_D[len + i] + poly_D[len + i] - poly_D[i];
        let eval_point_2 = comb_func(
          &poly_A_bound_point,
          &poly_B_bound_point,
          &poly_C_bound_point,
          &poly_D_bound_point,
        );

        // eval 3: bound_func is -2A(low) + 3A(high); computed incrementally with bound_func applied to eval(2)
        let poly_A_bound_point = poly_A_bound_point + poly_A[len + i] - poly_A[i];
        let poly_B_bound_point = poly_B_bound_point + poly_B[len + i] - poly_B[i];
        let poly_C_bound_point = poly_C_bound_point + poly_C[len + i] - poly_C[i];
        let poly_D_bound_point = poly_D_bound_point + poly_D[len + i] - poly_D[i];
        let eval_point_3 = comb_func(
          &poly_A_bound_point,
          &poly_B_bound_point,
          &poly_C_bound_point,
          &poly_D_bound_point,
        );
        (eval_point_0, eval_point_2, eval_point_3)
      })
      .reduce(
        || (E::Scalar::ZERO, E::Scalar::ZERO, E::Scalar::ZERO),
        |a, b| (a.0 + b.0, a.1 + b.1, a.2 + b.2),
      )
  }

  #[inline]
  #[allow(unused_variables, unused)]
  pub fn compute_eval_points_cubic_with_eq_partial<F>(
    poly_eq_left: &MultilinearPolynomial<E::Scalar>,
    poly_eq_right: &MultilinearPolynomial<E::Scalar>,
    poly_A: &MultilinearPolynomial<E::Scalar>,
    poly_B: &MultilinearPolynomial<E::Scalar>,
    poly_C: &MultilinearPolynomial<E::Scalar>,
    comb_func: F,
  ) -> (E::Scalar, E::Scalar, E::Scalar)
  where
    F: Fn(&E::Scalar, &E::Scalar, &E::Scalar) -> E::Scalar + Sync,
  {
    let half_l = poly_eq_right.len() as usize;
    let half_p = half_l * poly_eq_left.len();

    assert_eq!(half_p * 2, poly_A.len());

    let (zero_A, one_A) = poly_A.Z.split_at(half_p);
    let (zero_B, one_B) = poly_B.Z.split_at(half_p);
    let (zero_C, one_C) = poly_C.Z.split_at(half_p);

    let zip_A = zero_A.par_iter().zip(one_A);
    let zip_B = zero_B.par_iter().zip(one_B);
    let zip_C = zero_C.par_iter().zip(one_C);

    dbg!(poly_A.len());
    dbg!(zero_A.len());

    // dbg!(poly_eq_left);
    // dbg!(poly_eq_right);

    zip_A
      .zip(zip_B)
      .zip(zip_C)
      .enumerate()
      .map(
        |(i, (((zero_a, one_a), (zero_b, one_b)), (zero_c, one_c)))| {
          let one_a_double = one_a.double();
          let one_b_double = one_b.double();
          let one_c_double = one_c.double();
          let one_a_triple = one_a_double + one_a;
          let one_b_triple = one_b_double + one_b;
          let one_c_triple = one_c_double + one_c;
          let zero_a_double = zero_a.double();
          let zero_b_double = zero_b.double();
          let zero_c_double = zero_c.double();

          let e_a_0 = zero_a;
          let e_b_0 = zero_b;
          let e_c_0 = zero_c;

          let e_a_2 = one_a_double - zero_a;
          let e_b_2 = one_b_double - zero_b;
          let e_c_2 = one_c_double - zero_c;

          let e_a_3 = one_a_triple - zero_a_double;
          let e_b_3 = one_b_triple - zero_b_double;
          let e_c_3 = one_c_triple - zero_c_double;

          let eval_0 = comb_func(&e_a_0, &e_b_0, &e_c_0);
          let eval_2 = comb_func(&e_a_2, &e_b_2, &e_c_2);
          let eval_3 = comb_func(&e_a_3, &e_b_3, &e_c_3);

          let right_id = i % half_l;
          let left_id = i >> (half_l.ilog2() as usize);

          assert_eq!(i, right_id + left_id * half_l);
          assert!(half_l.is_power_of_two());

          let left = poly_eq_left[left_id];
          let right = poly_eq_right[right_id];

          let eval_0 = eval_0 * left * right;
          let eval_2 = eval_2 * left * right;
          let eval_3 = eval_3 * left * right;

          (eval_0, eval_2, eval_3)
        },
      )
      .reduce(
        || (E::Scalar::ZERO, E::Scalar::ZERO, E::Scalar::ZERO),
        |a, b| (a.0 + b.0, a.1 + b.1, a.2 + b.2),
      )
  }

  #[inline]
  #[allow(unused_variables, unused)]
  pub fn compute_eval_points_cubic_with_eq_partial_2<F>(
    poly_eq_right: &MultilinearPolynomial<E::Scalar>,
    poly_A: &MultilinearPolynomial<E::Scalar>,
    poly_B: &MultilinearPolynomial<E::Scalar>,
    poly_C: &MultilinearPolynomial<E::Scalar>,
    comb_func: F,
  ) -> (E::Scalar, E::Scalar, E::Scalar)
  where
    F: Fn(&E::Scalar, &E::Scalar, &E::Scalar) -> E::Scalar + Sync,
  {
    let half_l = poly_eq_right.len() as usize;
    let half_p = half_l / 2;

    assert_eq!(poly_eq_right.len() * 2, poly_A.len());

    let (zero_A, one_A) = poly_A.Z.split_at(half_l);
    let (zero_B, one_B) = poly_B.Z.split_at(half_l);
    let (zero_C, one_C) = poly_C.Z.split_at(half_l);

    let zip_A = zero_A.par_iter().zip(one_A);
    let zip_B = zero_B.par_iter().zip(one_B);
    let zip_C = zero_C.par_iter().zip(one_C);

    dbg!(poly_A.len());
    dbg!(zero_A.len());
    assert_eq!(zero_A.len() * 2, poly_A.len());

    // dbg!(poly_eq_left);
    // dbg!(poly_eq_right);

    zip_A
      .zip(zip_B)
      .zip(zip_C)
      .enumerate()
      .map(
        |(i, (((zero_a, one_a), (zero_b, one_b)), (zero_c, one_c)))| {
          let one_a_double = one_a.double();
          let one_b_double = one_b.double();
          let one_c_double = one_c.double();
          let one_a_triple = one_a_double + one_a;
          let one_b_triple = one_b_double + one_b;
          let one_c_triple = one_c_double + one_c;
          let zero_a_double = zero_a.double();
          let zero_b_double = zero_b.double();
          let zero_c_double = zero_c.double();

          let e_a_0 = zero_a;
          let e_b_0 = zero_b;
          let e_c_0 = zero_c;

          let e_a_2 = one_a_double - zero_a;
          let e_b_2 = one_b_double - zero_b;
          let e_c_2 = one_c_double - zero_c;

          let e_a_3 = one_a_triple - zero_a_double;
          let e_b_3 = one_b_triple - zero_b_double;
          let e_c_3 = one_c_triple - zero_c_double;

          let eval_0 = comb_func(&e_a_0, &e_b_0, &e_c_0);
          let eval_2 = comb_func(&e_a_2, &e_b_2, &e_c_2);
          let eval_3 = comb_func(&e_a_3, &e_b_3, &e_c_3);

          let right_id = i;
          assert!(half_l.is_power_of_two());

          let right = poly_eq_right[right_id];

          let eval_0 = eval_0 * right;
          let eval_2 = eval_2 * right;
          let eval_3 = eval_3 * right;

          (eval_0, eval_2, eval_3)
        },
      )
      .reduce(
        || (E::Scalar::ZERO, E::Scalar::ZERO, E::Scalar::ZERO),
        |a, b| (a.0 + b.0, a.1 + b.1, a.2 + b.2),
      )
  }

  pub fn prove_cubic_with_additive_term<F>(
    claim: &E::Scalar,
    num_rounds: usize,
    poly_A: &mut MultilinearPolynomial<E::Scalar>,
    poly_B: &mut MultilinearPolynomial<E::Scalar>,
    poly_C: &mut MultilinearPolynomial<E::Scalar>,
    poly_D: &mut MultilinearPolynomial<E::Scalar>,
    comb_func: F,
    transcript: &mut E::TE,
  ) -> Result<(Self, Vec<E::Scalar>, Vec<E::Scalar>), NovaError>
  where
    F: Fn(&E::Scalar, &E::Scalar, &E::Scalar, &E::Scalar) -> E::Scalar + Sync,
  {
    let mut r: Vec<E::Scalar> = Vec::new();
    let mut polys: Vec<CompressedUniPoly<E::Scalar>> = Vec::new();
    let mut claim_per_round = *claim;

    for _ in 0..num_rounds {
      let poly = {
        // Make an iterator returning the contributions to the evaluations
        let (eval_point_0, eval_point_2, eval_point_3) =
          Self::compute_eval_points_cubic_with_additive_term(
            poly_A, poly_B, poly_C, poly_D, &comb_func,
          );

        let evals = vec![
          eval_point_0,
          claim_per_round - eval_point_0,
          eval_point_2,
          eval_point_3,
        ];
        UniPoly::from_evals(&evals)
      };

      // append the prover's message to the transcript
      transcript.absorb(b"p", &poly);

      //derive the verifier's challenge for the next round
      let r_i = transcript.squeeze(b"c")?;
      r.push(r_i);
      polys.push(poly.compress());

      // Set up next round
      claim_per_round = poly.evaluate(&r_i);

      // bound all tables to the verifier's challenge
      rayon::join(
        || {
          rayon::join(
            || poly_A.bind_poly_var_top(&r_i),
            || poly_B.bind_poly_var_top(&r_i),
          )
        },
        || {
          rayon::join(
            || poly_C.bind_poly_var_top(&r_i),
            || poly_D.bind_poly_var_top(&r_i),
          )
        },
      );
    }

    Ok((
      SumcheckProof {
        compressed_polys: polys,
      },
      r,
      vec![poly_A[0], poly_B[0], poly_C[0], poly_D[0]],
    ))
  }
}
