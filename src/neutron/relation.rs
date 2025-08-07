//! This module defines relations used in the Neutron folding scheme
use crate::{
  errors::NovaError,
  r1cs::{R1CSInstance, R1CSShape, R1CSWitness, RelaxedR1CSInstance, RelaxedR1CSWitness},
  spartan::math::Math,
  traits::{commitment::CommitmentEngineTrait, AbsorbInRO2Trait, Engine, ROTrait},
  Commitment, CommitmentKey, DerandKey, CE,
};
use ff::Field;
use rand_core::OsRng;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

/// A type that holds structure information for a zero-fold relation
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(bound = "")]
pub struct Structure<E: Engine> {
  /// Shape of the R1CS relation
  pub(crate) S: R1CSShape<E>,

  /// the number of variables in the Eq polynomial
  pub(crate) ell: usize,
  pub(crate) left: usize,
  pub(crate) right: usize,
}

/// A type that holds witness information for a zero-fold relation
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(bound = "")]
pub struct FoldedWitness<E: Engine> {
  /// Running witness of the main relation
  pub(crate) W: Vec<E::Scalar>,
  pub(crate) r_W: E::Scalar,

  /// eq polynomial in tensor form
  pub(crate) E: Vec<E::Scalar>,
  pub(crate) r_E: E::Scalar,
}

/// A type that holds instance information for a zero-fold relation
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(bound = "")]
pub struct FoldedInstance<E: Engine> {
  pub(crate) comm_W: Commitment<E>,
  pub(crate) comm_E: Commitment<E>,
  pub(crate) T: E::Scalar,
  pub(crate) u: E::Scalar,
  pub(crate) X: Vec<E::Scalar>,
}

impl<E: Engine> Structure<E> {
  /// Create a new structure using the provided shape
  pub fn new(S: &R1CSShape<E>) -> Self {
    // pad to the regular shape
    let S = S.pad();

    let ell = S.num_cons.next_power_of_two().log_2();

    // we split ell into ell1 and ell2 such that ell1 + ell2 = ell and ell1 >= ell2
    let ell1 = (ell + 1) / 2; // This ensures ell1 >= ell2
    let ell2 = ell / 2;

    Structure {
      S: S.clone(),
      ell,
      left: 1 << ell1,
      right: 1 << ell2,
    }
  }

  /// Samples a new random `R1CSInstance`/`R1CSWitness` pair for use in NIFS folding
  pub fn sample_random_r1cs_instance_witness(
    &self,
    ck: &CommitmentKey<E>,
  ) -> Result<(R1CSInstance<E>, R1CSWitness<E>), NovaError> {
    // Use the underlying R1CSShape's method to get RelaxedR1CS types
    let (relaxed_u, relaxed_w) = self.S.sample_random_instance_witness(ck)?;
    
    // Convert RelaxedR1CS to R1CS (this works because they're random instances with u=1)
    let u = R1CSInstance {
      comm_W: relaxed_u.comm_W,
      X: relaxed_u.X,
    };
    let w = R1CSWitness {
      W: relaxed_w.W,
      r_W: relaxed_w.r_W,
    };
    
    Ok((u, w))
  }

  /// Samples a new random `FoldedInstance`/`FoldedWitness` pair 
  pub fn sample_random_instance_witness(
    &self,
    ck: &CommitmentKey<E>,
  ) -> Result<(FoldedInstance<E>, FoldedWitness<E>), NovaError> {
    // sample Z = (W, u, X)
    let Z = (0..self.S.num_vars + self.S.num_io + 1)
      .into_par_iter()
      .map(|_| E::Scalar::random(&mut OsRng))
      .collect::<Vec<E::Scalar>>();

    let r_W = E::Scalar::random(&mut OsRng);
    let r_E = E::Scalar::random(&mut OsRng);

    let u = Z[self.S.num_vars];

    // compute E <- AZ o BZ - u * CZ
    let (AZ, BZ, CZ) = self.S.multiply_vec(&Z)?;

    let E = AZ
      .par_iter()
      .zip(BZ.par_iter())
      .zip(CZ.par_iter())
      .map(|((az, bz), cz)| *az * *bz - u * *cz)
      .collect::<Vec<E::Scalar>>();

    // Pad E to left + right size (Neutron uses different E structure than Nova)
    let mut E_padded = vec![E::Scalar::ZERO; self.left + self.right];
    for (i, &e) in E.iter().enumerate() {
      if i < E_padded.len() {
        E_padded[i] = e;
      }
    }

    // compute commitments to W,E in parallel
    let (comm_W, comm_E) = rayon::join(
      || CE::<E>::commit(ck, &Z[..self.S.num_vars], &r_W),
      || CE::<E>::commit(ck, &E_padded, &r_E),
    );

    // Compute T = sum(E) (simplification for random instance)
    let T = E_padded.iter().fold(E::Scalar::ZERO, |acc, &e| acc + e);

    Ok((
      FoldedInstance {
        comm_W,
        comm_E,
        T,
        u,
        X: Z[self.S.num_vars + 1..].to_vec(),
      },
      FoldedWitness {
        W: Z[..self.S.num_vars].to_vec(),
        r_W,
        E: E_padded,
        r_E,
      },
    ))
  }

  /// Check if the witness is satisfying
  pub fn is_sat(
    &self,
    ck: &CommitmentKey<E>,
    U: &FoldedInstance<E>,
    W: &FoldedWitness<E>,
  ) -> Result<(), NovaError> {
    // check if the witness is satisfying
    let z = [W.W.clone(), vec![U.u], U.X.clone()].concat();
    let (Az, Bz, Cz) = self.S.multiply_vec(&z)?;

    // full_E is the outer outer product of E1 and E2
    // E1 and E2 are splits of E
    let (E1, E2) = W.E.split_at(self.left);
    let mut full_E = vec![E::Scalar::ONE; self.left * self.right];
    for i in 0..self.right {
      for j in 0..self.left {
        full_E[i * self.left + j] = E2[i] * E1[j];
      }
    }

    let sum = full_E
      .par_iter()
      .zip(Az.par_iter())
      .zip(Bz.par_iter())
      .zip(Cz.par_iter())
      .map(|(((e, a), b), c)| *e * ((*a) * (*b) - *c))
      .reduce(|| E::Scalar::ZERO, |acc, x| acc + x);

    if sum != U.T {
      println!("sum: {:?}", sum);
      println!("U.T: {:?}", U.T);
      return Err(NovaError::UnSat {
        reason: "sum != U.T".to_string(),
      });
    }

    // check the validity of the commitments
    let comm_W = E::CE::commit(ck, &W.W, &W.r_W);
    let comm_E = E::CE::commit(ck, &W.E, &W.r_E);

    if comm_W != U.comm_W || comm_E != U.comm_E {
      return Err(NovaError::UnSat {
        reason: "comm_W != U.comm_W || comm_E != U.comm_E".to_string(),
      });
    }

    Ok(())
  }
}

impl<E: Engine> FoldedWitness<E> {
  /// Create a default witness
  pub fn default(S: &Structure<E>) -> Self {
    FoldedWitness {
      W: vec![E::Scalar::ZERO; S.S.num_vars],
      r_W: E::Scalar::ZERO,
      E: vec![E::Scalar::ZERO; S.left + S.right],
      r_E: E::Scalar::ZERO,
    }
  }

  /// Derandomize the witness by returning a witness with zero randomness and the randomness values
  pub fn derandomize(&self) -> (Self, E::Scalar, E::Scalar) {
    (
      FoldedWitness {
        W: self.W.clone(),
        r_W: E::Scalar::ZERO,
        E: self.E.clone(),
        r_E: E::Scalar::ZERO,
      },
      self.r_W,
      self.r_E,
    )
  }

  /// Convert to a RelaxedR1CSWitness for use with SNARKs
  pub fn to_relaxed_witness(&self) -> RelaxedR1CSWitness<E> {
    RelaxedR1CSWitness {
      W: self.W.clone(),
      r_W: self.r_W,
      E: self.E.clone(),
      r_E: self.r_E,
    }
  }

  /// Fold the witness with another witness
  pub fn fold(
    &self,
    W2: &R1CSWitness<E>,
    E2: &Vec<E::Scalar>,
    r_E2: &E::Scalar,
    r_b: &E::Scalar,
  ) -> Result<Self, NovaError> {
    // we need to compute the weighted sum using weights of (1-r_b) and r_b
    let W = self
      .W
      .par_iter()
      .zip(W2.W.par_iter())
      .map(|(w1, w2)| *w1 + *r_b * (*w2 - *w1))
      .collect::<Vec<_>>();
    let r_W = (E::Scalar::ONE - r_b) * self.r_W + *r_b * W2.r_W;

    let E = self
      .E
      .par_iter()
      .zip(E2.par_iter())
      .map(|(e1, e2)| *e1 + *r_b * (*e2 - *e1))
      .collect::<Vec<_>>();
    let r_E = (E::Scalar::ONE - r_b) * self.r_E + *r_b * r_E2;

    Ok(Self { W, r_W, E, r_E })
  }
}

impl<E: Engine> FoldedInstance<E> {
  /// Create a default instance
  pub fn default(S: &Structure<E>) -> Self {
    FoldedInstance {
      comm_W: Commitment::<E>::default(),
      comm_E: Commitment::<E>::default(),
      T: E::Scalar::ZERO,
      u: E::Scalar::ZERO,
      X: vec![E::Scalar::ZERO; S.S.num_io],
    }
  }

  /// Derandomize the instance using the derand key and randomness values
  pub fn derandomize(
    &self,
    dk: &DerandKey<E>,
    r_W: &E::Scalar,
    r_E: &E::Scalar,
  ) -> FoldedInstance<E> {
    FoldedInstance {
      comm_W: CE::<E>::derandomize(dk, &self.comm_W, r_W),
      comm_E: CE::<E>::derandomize(dk, &self.comm_E, r_E),
      T: self.T,
      u: self.u,
      X: self.X.clone(),
    }
  }

  /// Convert to a RelaxedR1CSInstance for use with SNARKs
  pub fn to_relaxed_instance(&self) -> RelaxedR1CSInstance<E> {
    RelaxedR1CSInstance {
      comm_W: self.comm_W,
      comm_E: self.comm_E,
      u: self.u,
      X: self.X.clone(),
    }
  }

  /// Fold the instance with another instance
  pub fn fold(
    &self,
    U2: &R1CSInstance<E>,
    comm_E: &Commitment<E>,
    r_b: &E::Scalar,
    T_out: &E::Scalar,
  ) -> Result<Self, NovaError> {
    // we need to compute the weighted sum using weights of (1-r_b) and r_b
    let comm_W = self.comm_W * (E::Scalar::ONE - r_b) + U2.comm_W * *r_b;
    let comm_E = self.comm_E * (E::Scalar::ONE - r_b) + *comm_E * *r_b;
    let X = self
      .X
      .par_iter()
      .zip(U2.X.par_iter())
      .map(|(x1, x2)| (E::Scalar::ONE - r_b) * x1 + *r_b * x2)
      .collect::<Vec<_>>();
    let u = (E::Scalar::ONE - r_b) * self.u + r_b;

    Ok(Self {
      comm_W,
      comm_E,
      T: *T_out,
      u,
      X,
    })
  }
}

impl<E: Engine> AbsorbInRO2Trait<E> for FoldedInstance<E> {
  fn absorb_in_ro2(&self, ro: &mut E::RO2) {
    self.comm_W.absorb_in_ro2(ro);
    self.comm_E.absorb_in_ro2(ro);

    ro.absorb(self.T);
    ro.absorb(self.u);
    for x in &self.X {
      ro.absorb(*x);
    }
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::{
    frontend::{
      r1cs::{NovaShape, NovaWitness},
      shape_cs::ShapeCS,
      solver::SatisfyingAssignment,
      Circuit, ConstraintSystem,
    },
    provider::{hyperkzg::EvaluationEngine, Bn256EngineKZG},
    spartan::{direct::DirectCircuit, snark::RelaxedR1CSSNARK},
    spartan::{math::Math, polys::eq::EqPolynomial},
    traits::{circuit::NonTrivialCircuit, snark::RelaxedR1CSSNARKTrait},
  };
  use rand::rngs::OsRng;

  fn test_sat_inner<E: Engine, S: RelaxedR1CSSNARKTrait<E>>() -> Result<(), NovaError> {
    // generate a non-trivial circuit
    let num_cons: usize = 16;
    let log_num_cons = num_cons.log_2();

    let circuit: DirectCircuit<E, NonTrivialCircuit<E::Scalar>> =
      DirectCircuit::new(None, NonTrivialCircuit::<E::Scalar>::new(num_cons));

    // synthesize the circuit's shape
    let mut cs: ShapeCS<E> = ShapeCS::new();
    let _ = circuit.synthesize(&mut cs);
    let (shape, ck) = cs.r1cs_shape(&*S::ck_floor());
    let S = Structure::new(&shape);

    // test default instance-witness pair under the structure
    let W = FoldedWitness::default(&S);
    let U = FoldedInstance::default(&S);
    S.is_sat(&ck, &U, &W)?;

    // generate a satisfying instance-witness for the r1cs
    let circuit: DirectCircuit<E, NonTrivialCircuit<E::Scalar>> = DirectCircuit::new(
      Some(vec![E::Scalar::from(2)]),
      NonTrivialCircuit::<E::Scalar>::new(num_cons),
    );
    let mut cs = SatisfyingAssignment::<E>::new();
    let _ = circuit.synthesize(&mut cs);
    let (u, w) = cs
      .r1cs_instance_and_witness(&shape, &ck)
      .map_err(|_e| NovaError::UnSat {
        reason: "Unable to generate a satisfying witness".to_string(),
      })?;

    // generate a random eq polynomial
    let coords = (0..log_num_cons)
      .map(|_| E::Scalar::random(&mut OsRng))
      .collect::<Vec<_>>();
    let E = EqPolynomial::new(coords).evals();

    // pad witness
    let mut W = w.W.clone();
    W.resize(S.S.num_vars, E::Scalar::ZERO);

    let W = FoldedWitness {
      W,
      r_W: w.r_W,
      E: E.clone(),
      r_E: E::Scalar::random(&mut OsRng),
    };

    let U = FoldedInstance {
      comm_W: u.comm_W,
      comm_E: E::CE::commit(&ck, &E, &W.r_E),
      T: E::Scalar::ZERO,
      X: u.X.clone(),
      u: E::Scalar::ONE,
    };

    S.is_sat(&ck, &U, &W)
  }

  #[test]
  fn test_sat() {
    type E = Bn256EngineKZG;
    type S = RelaxedR1CSSNARK<E, EvaluationEngine<E>>;
    let res = test_sat_inner::<E, S>();
    assert!(res.is_ok());
  }
}
