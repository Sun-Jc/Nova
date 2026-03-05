//! This module provides a multi-scalar multiplication routine.
//!
//! For full-field-element MSM, we use a signed-decomposition + bit-width-partitioning strategy
//!
//! 1. **Signed scalar decomposition**: For each scalar `s`, we compare `num_bits(s)` vs
//!    `num_bits(p - s)` and use whichever representation is smaller, negating the base point
//!    if we use `p - s`. This halves the effective scalar range.
//! 2. **Bit-width partitioning**: Scalars are routed to the optimal algorithm based on their
//!    actual bit-width after signed reduction: binary accumulation for 0/1, single-window
//!    bucket sort for ≤10 bits, multi-window Pippenger for ≤64 bits, and wNAF for the rest.
//! 3. **wNAF (windowed non-adjacent form)**: Large scalars use signed digit decomposition
//!    with XYZZ bucket coordinates for ~18% faster accumulation.
//! 4. **XYZZ bucket coordinates**: Extended Jacobian `(X, Y, ZZ, ZZZ)` provides cheaper
//!    mixed addition (7M + 2S) compared to standard Jacobian (~11M + 5S for proj+affine).
//!
//! The MSM implementations (for integer types and field types) are adapted from halo2/jolt.
#![allow(unsafe_code)]
use ff::{Field, PrimeField};
use halo2curves::{group::Group, CurveAffine};
use num_integer::Integer;
use num_traits::{ToPrimitive, Zero};
use rayon::{current_num_threads, prelude::*};

// ==================================================================================
// XYZZ (Extended Jacobian) Bucket coordinates
// ==================================================================================

/// Extended Jacobian (XYZZ) coordinates for efficient MSM bucket accumulation.
///
/// Stores `(X, Y, ZZ, ZZZ)` where `ZZ = Z²` and `ZZZ = Z³` for a Jacobian point
/// with coordinates `(X/ZZ, Y/ZZZ)` in affine.
///
/// Mixed addition (affine + XYZZ) costs 7M + 2S vs ~11M + 5S for standard projective+affine.
/// Formula source: <https://www.hyperelliptic.org/EFD/g1p/auto-shortw-xyzz.html>
///
/// **Assumes `a = 0`** in the curve equation `y² = x³ + ax + b`, which holds for
/// all curves used in Nova (BN254, Grumpkin, Pallas, Vesta, secp256k1, secq256k1).
#[derive(Copy, Clone)]
struct BucketXYZZ<F: Field> {
  x: F,
  y: F,
  zz: F,
  zzz: F,
}

impl<F: Field> BucketXYZZ<F> {
  /// The point at infinity (identity).
  #[inline]
  fn zero() -> Self {
    Self {
      x: F::ONE,
      y: F::ONE,
      zz: F::ZERO,
      zzz: F::ZERO,
    }
  }

  /// Check if this is the identity.
  #[inline]
  fn is_zero(&self) -> bool {
    self.zz == F::ZERO
  }

  /// Double in place (dbl-2008-s-1, assumes a=0).
  /// Cost: 2M + 5S + 7add
  fn double_in_place(&mut self) {
    if self.is_zero() {
      return;
    }
    // U = 2*Y1
    let u = self.y.double();
    // V = U^2
    let v = u.square();
    // W = U*V
    let w = u * v;
    // S = X1*V
    let s = self.x * v;
    // M = 3*X1^2 (a=0, so no a*ZZ^2 term)
    let x_sq = self.x.square();
    let m = x_sq.double() + x_sq;
    // X3 = M^2 - 2*S
    self.x = m.square() - s.double();
    // Y3 = M*(S - X3) - W*Y1
    self.y = m * (s - self.x) - w * self.y;
    // ZZ3 = V*ZZ1
    self.zz *= v;
    // ZZZ3 = W*ZZZ1
    self.zzz *= w;
  }

  /// XYZZ += XYZZ (full addition, add-2008-s).
  fn add_assign_bucket(&mut self, other: &Self) {
    if other.is_zero() {
      return;
    }
    if self.is_zero() {
      *self = *other;
      return;
    }
    // U1 = X1*ZZ2, U2 = X2*ZZ1
    let u1 = self.x * other.zz;
    let u2 = other.x * self.zz;
    // S1 = Y1*ZZZ2, S2 = Y2*ZZZ1
    let s1 = self.y * other.zzz;
    let s2 = other.y * self.zzz;

    if u1 == u2 {
      if s1 == s2 {
        self.double_in_place();
      } else {
        *self = Self::zero();
      }
      return;
    }
    let p = u2 - u1;
    let r = s2 - s1;
    let pp = p.square();
    let ppp = p * pp;
    let q = u1 * pp;
    self.x = r.square() - ppp - q.double();
    self.y = r * (q - self.x) - s1 * ppp;
    self.zz = self.zz * other.zz * pp;
    self.zzz = self.zzz * other.zzz * ppp;
  }
}

/// Mixed addition: BucketXYZZ += CurveAffine point (madd-2008-s).
/// Cost: 7M + 2S
#[inline]
fn bucket_add_affine<C: CurveAffine>(bucket: &mut BucketXYZZ<C::Base>, p: &C) {
  if bool::from(p.is_identity()) {
    return;
  }
  let coords = p.coordinates().unwrap();
  let px = *coords.x();
  let py = *coords.y();

  if bucket.is_zero() {
    bucket.x = px;
    bucket.y = py;
    bucket.zz = C::Base::ONE;
    bucket.zzz = C::Base::ONE;
    return;
  }
  // U2 = X2*ZZ1, S2 = Y2*ZZZ1
  let u2 = px * bucket.zz;
  let s2 = py * bucket.zzz;

  if bucket.x == u2 {
    if bucket.y == s2 {
      bucket.double_in_place();
    } else {
      *bucket = BucketXYZZ::zero();
    }
    return;
  }
  let p_val = u2 - bucket.x;
  let r = s2 - bucket.y;
  let pp = p_val.square();
  let ppp = p_val * pp;
  let q = bucket.x * pp;
  bucket.x = r.square() - ppp - q.double();
  bucket.y = r * (q - bucket.x) - bucket.y * ppp;
  bucket.zz *= pp;
  bucket.zzz *= ppp;
}

/// Convert XYZZ bucket to projective curve point.
///
/// Computes affine coordinates `(X/ZZ, Y/ZZZ)` then converts to projective.
/// Only called O(windows) times per thread, so the field inversion cost is negligible.
#[inline]
fn bucket_to_curve<C: CurveAffine>(bucket: &BucketXYZZ<C::Base>) -> C::CurveExt {
  if bucket.is_zero() {
    return C::CurveExt::identity();
  }
  let zz_inv = bucket.zz.invert().unwrap();
  let zzz_inv = bucket.zzz.invert().unwrap();
  let x = bucket.x * zz_inv;
  let y = bucket.y * zzz_inv;
  C::from_xy(x, y)
    .expect("XYZZ bucket should produce a valid curve point")
    .into()
}

// ==================================================================================
// Scalar utilities
// ==================================================================================

/// Count significant bits in a field element (from its little-endian repr).
#[inline]
fn scalar_num_bits<F: PrimeField>(s: &F) -> u32 {
  let repr = s.to_repr();
  let bytes = repr.as_ref();
  for i in (0..bytes.len()).rev() {
    if bytes[i] != 0 {
      return i as u32 * 8 + (8 - bytes[i].leading_zeros());
    }
  }
  0
}

/// Extract the low 64 bits from a field element's little-endian representation.
#[inline]
fn repr_low_u64<F: PrimeField>(s: &F) -> u64 {
  let repr = s.to_repr();
  let bytes = repr.as_ref();
  let mut buf = [0u8; 8];
  let len = bytes.len().min(8);
  buf[..len].copy_from_slice(&bytes[..len]);
  u64::from_le_bytes(buf)
}

// ==================================================================================
// Main MSM with signed decomposition + bit-width partitioning
// ==================================================================================

/// Performs an optimized multi-scalar multiplication for full field element scalars.
///
/// This function uses signed scalar decomposition to halve the effective scalar range,
/// then partitions scalars by bit-width to route each group to the optimal MSM algorithm.
/// Small scalars (≤64 bits) use Nova's bucket-sort MSMs; large scalars delegate to
/// halo2curves' `msm_best`.
///
/// Adapted from the Jolt MSM implementation, with optimizations for Nova's use cases.
pub fn msm<C: CurveAffine>(coeffs: &[C::Scalar], bases: &[C]) -> C::Curve {
  assert_eq!(coeffs.len(), bases.len());
  let n = coeffs.len();
  if n == 0 {
    return C::Curve::identity();
  }

  // For very small inputs, use the simple fallback
  if n <= 16 {
    return msm_simple(coeffs, bases);
  }

  // Group indices: 0=unit_pos, 1=unit_neg, 2=pos≤8, 3=neg≤8,
  // 4=pos≤16, 5=neg≤16, 6=pos≤32, 7=neg≤32, 8=pos≤64, 9=neg≤64, 10=large
  const NUM_GROUPS: usize = 11;

  // Phase 1: Classify each scalar in parallel
  // Encode as u64: top 4 bits = group, bottom 60 bits = original index
  let classified: Vec<u64> = coeffs
    .par_iter()
    .enumerate()
    .filter_map(|(i, s)| {
      if bool::from(s.is_zero()) {
        return None;
      }
      let neg_s = -(*s);
      let bits_s = scalar_num_bits(s);
      let bits_neg = scalar_num_bits(&neg_s);

      let group = if bits_s <= 1 {
        0u8 // unit positive
      } else if bits_neg <= 1 {
        1u8 // unit negative
      } else if bits_s <= 8 {
        2u8
      } else if bits_neg <= 8 {
        3u8
      } else if bits_s <= 16 {
        4u8
      } else if bits_neg <= 16 {
        5u8
      } else if bits_s <= 32 {
        6u8
      } else if bits_neg <= 32 {
        7u8
      } else if bits_s <= 64 {
        8u8
      } else if bits_neg <= 64 {
        9u8
      } else {
        10u8 // large
      };
      Some(((i as u64) & 0x0FFF_FFFF_FFFF_FFFF) | ((group as u64) << 60))
    })
    .collect();

  if classified.is_empty() {
    return C::Curve::identity();
  }

  // Phase 2: Sort by group for efficient partitioning
  let mut classified = classified;
  classified.par_sort_unstable_by_key(|v| (v >> 60) as u8);

  let extract_group = |v: u64| (v >> 60) as u8;
  let extract_index = |v: u64| (v & 0x0FFF_FFFF_FFFF_FFFF) as usize;

  // Find partition boundaries
  let mut boundaries = [0usize; NUM_GROUPS + 1];
  {
    let mut pos = 0;
    for g in 0..NUM_GROUPS as u8 {
      boundaries[g as usize] = pos;
      pos += classified[pos..].partition_point(|v| extract_group(*v) <= g);
    }
    boundaries[NUM_GROUPS] = classified.len();
  }

  // Helper to extract (bases, u64_scalars) for a group range
  let extract_u64_group = |start: usize, end: usize, negate: bool| -> (Vec<C>, Vec<u64>) {
    classified[start..end]
      .iter()
      .map(|&v| {
        let idx = extract_index(v);
        let b = bases[idx];
        let s = if negate { -coeffs[idx] } else { coeffs[idx] };
        (b, repr_low_u64(&s))
      })
      .unzip()
  };

  // Helper to extract (bases, bool_scalars) for unit groups
  let extract_binary_group = |start: usize, end: usize| -> Vec<C> {
    classified[start..end]
      .iter()
      .map(|&v| bases[extract_index(v)])
      .collect()
  };

  // Phase 3: Compute MSM for each group in parallel
  // Binary groups (unit scalars): just accumulate bases
  let (g0_start, g0_end) = (boundaries[0], boundaries[1]);
  let (g1_start, g1_end) = (boundaries[1], boundaries[2]);

  // Small-scalar groups
  let (g2_start, g2_end) = (boundaries[2], boundaries[3]);
  let (g3_start, g3_end) = (boundaries[3], boundaries[4]);
  let (g4_start, g4_end) = (boundaries[4], boundaries[5]);
  let (g5_start, g5_end) = (boundaries[5], boundaries[6]);
  let (g6_start, g6_end) = (boundaries[6], boundaries[7]);
  let (g7_start, g7_end) = (boundaries[7], boundaries[8]);
  let (g8_start, g8_end) = (boundaries[8], boundaries[9]);
  let (g9_start, g9_end) = (boundaries[9], boundaries[10]);

  // Large scalar group
  let (g10_start, g10_end) = (boundaries[10], boundaries[11]);

  // Execute all groups in parallel using nested rayon joins
  let (binary_result, small_and_large_result) = rayon::join(
    || {
      // Binary: pos - neg
      let (pos, neg) = rayon::join(
        || {
          let bases_pos = extract_binary_group(g0_start, g0_end);
          accumulate_bases::<C>(&bases_pos)
        },
        || {
          let bases_neg = extract_binary_group(g1_start, g1_end);
          accumulate_bases::<C>(&bases_neg)
        },
      );
      pos - neg
    },
    || {
      let (small_result, large_result) = rayon::join(
        || {
          // Small scalar groups: compute each pair (positive - negative)
          let ((r8, r16), (r32, r64)) = rayon::join(
            || {
              rayon::join(
                || {
                  let (pos_b, pos_s) = extract_u64_group(g2_start, g2_end, false);
                  let (neg_b, neg_s) = extract_u64_group(g3_start, g3_end, true);
                  msm_small_with_max_num_bits(&pos_s, &pos_b, 8)
                    - msm_small_with_max_num_bits(&neg_s, &neg_b, 8)
                },
                || {
                  let (pos_b, pos_s) = extract_u64_group(g4_start, g4_end, false);
                  let (neg_b, neg_s) = extract_u64_group(g5_start, g5_end, true);
                  msm_small_with_max_num_bits(&pos_s, &pos_b, 16)
                    - msm_small_with_max_num_bits(&neg_s, &neg_b, 16)
                },
              )
            },
            || {
              rayon::join(
                || {
                  let (pos_b, pos_s) = extract_u64_group(g6_start, g6_end, false);
                  let (neg_b, neg_s) = extract_u64_group(g7_start, g7_end, true);
                  msm_small_with_max_num_bits(&pos_s, &pos_b, 32)
                    - msm_small_with_max_num_bits(&neg_s, &neg_b, 32)
                },
                || {
                  let (pos_b, pos_s) = extract_u64_group(g8_start, g8_end, false);
                  let (neg_b, neg_s) = extract_u64_group(g9_start, g9_end, true);
                  msm_small_with_max_num_bits(&pos_s, &pos_b, 64)
                    - msm_small_with_max_num_bits(&neg_s, &neg_b, 64)
                },
              )
            },
          );
          r8 + r16 + r32 + r64
        },
        || {
          // Large scalars: delegate to halo2curves' optimized MSM
          if g10_start >= g10_end {
            return C::Curve::identity();
          }
          let (large_bases, large_coeffs): (Vec<C>, Vec<C::Scalar>) = classified
            [g10_start..g10_end]
            .iter()
            .map(|&v| {
              let idx = extract_index(v);
              (bases[idx], coeffs[idx])
            })
            .unzip();
          halo2curves::msm::msm_best(&large_coeffs, &large_bases)
        },
      );
      small_result + large_result
    },
  );

  binary_result + small_and_large_result
}

/// Simple MSM fallback for very small inputs.
fn msm_simple<C: CurveAffine>(coeffs: &[C::Scalar], bases: &[C]) -> C::Curve {
  coeffs
    .iter()
    .zip(bases.iter())
    .fold(C::Curve::identity(), |acc, (coeff, base)| {
      acc + *base * coeff
    })
}

/// Accumulate bases (sum of affine points, for binary MSM).
fn accumulate_bases<C: CurveAffine>(bases: &[C]) -> C::Curve {
  let num_threads = current_num_threads();
  if bases.is_empty() {
    return C::Curve::identity();
  }
  if bases.len() > num_threads {
    let chunk = bases.len().div_ceil(num_threads);
    bases
      .par_chunks(chunk)
      .map(|chunk| {
        chunk.iter().fold(C::Curve::identity(), |mut acc, b| {
          acc += *b;
          acc
        })
      })
      .reduce(C::Curve::identity, |a, b| a + b)
  } else {
    bases.iter().fold(C::Curve::identity(), |mut acc, b| {
      acc += *b;
      acc
    })
  }
}

fn num_bits(n: usize) -> usize {
  if n == 0 {
    0
  } else {
    (n.ilog2() + 1) as usize
  }
}

// ==================================================================================
// Small-scalar MSM with XYZZ buckets
// ==================================================================================

/// Multi-scalar multiplication using the best algorithm for the given scalars.
pub fn msm_small<C: CurveAffine, T: Integer + Into<u64> + Copy + Sync + ToPrimitive>(
  scalars: &[T],
  bases: &[C],
) -> C::Curve {
  let max_num_bits = num_bits(scalars.iter().max().unwrap().to_usize().unwrap());
  msm_small_with_max_num_bits(scalars, bases, max_num_bits)
}

/// Multi-scalar multiplication using the best algorithm for the given scalars.
pub fn msm_small_with_max_num_bits<
  C: CurveAffine,
  T: Integer + Into<u64> + Copy + Sync + ToPrimitive,
>(
  scalars: &[T],
  bases: &[C],
  max_num_bits: usize,
) -> C::Curve {
  assert_eq!(bases.len(), scalars.len());

  match max_num_bits {
    0 => C::identity().into(),
    1 => msm_binary(scalars, bases),
    2..=10 => msm_10(scalars, bases, max_num_bits),
    11..=32 => msm_small_rest(scalars, bases, max_num_bits),
    _ => {
      // For >32-bit scalars, halo2curves' msm_best is faster than our
      // bucket-sort Pippenger (e.g., 192ms vs 244ms at u64, 2^20 points).
      let field_scalars: Vec<C::ScalarExt> = scalars
        .iter()
        .map(|s| C::ScalarExt::from((*s).into()))
        .collect();
      halo2curves::msm::msm_best(&field_scalars, bases)
    }
  }
}

fn msm_binary<C: CurveAffine, T: Integer + Sync>(scalars: &[T], bases: &[C]) -> C::Curve {
  assert_eq!(scalars.len(), bases.len());
  let num_threads = current_num_threads();
  let process_chunk = |scalars: &[T], bases: &[C]| {
    let mut acc = C::Curve::identity();
    scalars
      .iter()
      .zip(bases.iter())
      .filter(|(scalar, _)| !scalar.is_zero())
      .for_each(|(_, base)| {
        acc += *base;
      });
    acc
  };

  if scalars.len() > num_threads {
    let chunk = scalars.len() / num_threads;
    scalars
      .par_chunks(chunk)
      .zip(bases.par_chunks(chunk))
      .map(|(scalars, bases)| process_chunk(scalars, bases))
      .reduce(C::Curve::identity, |sum, evl| sum + evl)
  } else {
    process_chunk(scalars, bases)
  }
}

/// MSM optimized for up to 10-bit scalars, using XYZZ bucket coordinates.
fn msm_10<C: CurveAffine, T: Into<u64> + Zero + Copy + Sync>(
  scalars: &[T],
  bases: &[C],
  max_num_bits: usize,
) -> C::Curve {
  fn msm_10_serial<C: CurveAffine, T: Into<u64> + Zero + Copy>(
    scalars: &[T],
    bases: &[C],
    max_num_bits: usize,
  ) -> C::Curve {
    let num_buckets: usize = 1 << max_num_bits;
    let mut buckets: Vec<BucketXYZZ<C::Base>> = vec![BucketXYZZ::zero(); num_buckets];

    scalars
      .iter()
      .zip(bases.iter())
      .filter(|(scalar, _base)| !scalar.is_zero())
      .for_each(|(scalar, base)| {
        let bucket_index: u64 = (*scalar).into();
        bucket_add_affine::<C>(&mut buckets[bucket_index as usize], base);
      });

    let mut result: BucketXYZZ<C::Base> = BucketXYZZ::zero();
    let mut running_sum: BucketXYZZ<C::Base> = BucketXYZZ::zero();
    for b in buckets.into_iter().skip(1).rev() {
      running_sum.add_assign_bucket(&b);
      result.add_assign_bucket(&running_sum);
    }
    bucket_to_curve::<C>(&result)
  }

  let num_threads = current_num_threads();
  if scalars.len() > num_threads {
    let chunk_size = scalars.len() / num_threads;
    scalars
      .par_chunks(chunk_size)
      .zip(bases.par_chunks(chunk_size))
      .map(|(scalars_chunk, bases_chunk)| msm_10_serial(scalars_chunk, bases_chunk, max_num_bits))
      .reduce(C::Curve::identity, |sum, evl| sum + evl)
  } else {
    msm_10_serial(scalars, bases, max_num_bits)
  }
}

fn msm_small_rest<C: CurveAffine, T: Into<u64> + Zero + Copy + Sync>(
  scalars: &[T],
  bases: &[C],
  max_num_bits: usize,
) -> C::Curve {
  fn msm_small_rest_serial<C: CurveAffine, T: Into<u64> + Zero + Copy>(
    scalars: &[T],
    bases: &[C],
    max_num_bits: usize,
  ) -> C::Curve {
    let mut c = if bases.len() < 32 {
      3
    } else {
      compute_ln(bases.len()) + 2
    };

    if max_num_bits == 32 || max_num_bits == 64 {
      c = 8;
    }

    let scalars_and_bases_iter = scalars.iter().zip(bases).filter(|(s, _base)| !s.is_zero());
    let window_starts: Vec<usize> = (0..max_num_bits).step_by(c).collect();

    // Each window is of size `c`.
    // We divide up the bits 0..num_bits into windows of size `c`, and
    // process each such window.
    let window_sums: Vec<C::CurveExt> = window_starts
      .iter()
      .map(|&w_start| {
        let mut res: BucketXYZZ<C::Base> = BucketXYZZ::zero();
        // We don't need the "zero" bucket, so we only have 2^c - 1 buckets.
        let mut buckets: Vec<BucketXYZZ<C::Base>> = vec![BucketXYZZ::zero(); (1 << c) - 1];
        // This clone is cheap, because the iterator contains just a
        // pointer and an index into the original vectors.
        scalars_and_bases_iter.clone().for_each(|(&scalar, base)| {
          let scalar: u64 = scalar.into();
          if scalar == 1 {
            // We only process unit scalars once in the first window.
            if w_start == 0 {
              bucket_add_affine::<C>(&mut res, base);
            }
          } else {
            let mut scalar = scalar;

            // We right-shift by w_start, thus getting rid of the
            // lower bits.
            scalar >>= w_start;

            // We mod the remaining bits by 2^{window size}, thus taking `c` bits.
            scalar %= 1 << c;

            // If the scalar is non-zero, we update the corresponding
            // bucket.
            // (Recall that `buckets` doesn't have a zero bucket.)
            if scalar != 0 {
              bucket_add_affine::<C>(&mut buckets[(scalar - 1) as usize], base);
            }
          }
        });

        // Prefix sum using XYZZ coordinates
        let mut running_sum: BucketXYZZ<C::Base> = BucketXYZZ::zero();
        for b in buckets.into_iter().rev() {
          running_sum.add_assign_bucket(&b);
          res.add_assign_bucket(&running_sum);
        }
        bucket_to_curve::<C>(&res)
      })
      .collect();

    // We store the sum for the lowest window.
    let lowest = *window_sums.first().unwrap();

    // We're traversing windows from high to low.
    lowest
      + window_sums[1..]
        .iter()
        .rev()
        .fold(C::CurveExt::identity(), |mut total, sum_i| {
          total += sum_i;
          for _ in 0..c {
            total = total.double();
          }
          total
        })
  }

  let num_threads = current_num_threads();
  if scalars.len() > num_threads {
    let chunk_size = scalars.len() / num_threads;
    scalars
      .par_chunks(chunk_size)
      .zip(bases.par_chunks(chunk_size))
      .map(|(scalars_chunk, bases_chunk)| {
        msm_small_rest_serial(scalars_chunk, bases_chunk, max_num_bits)
      })
      .reduce(C::Curve::identity, |sum, evl| sum + evl)
  } else {
    msm_small_rest_serial(scalars, bases, max_num_bits)
  }
}

fn compute_ln(a: usize) -> usize {
  // log2(a) * ln(2)
  if a == 0 {
    0 // Handle edge case where log2 is undefined
  } else {
    a.ilog2() as usize * 69 / 100
  }
}

#[inline(always)]
pub(crate) fn batch_add<C: CurveAffine>(bases: &[C], one_indices: &[usize]) -> C::Curve {
  fn add_chunk<C: CurveAffine>(bases: impl Iterator<Item = C>) -> C::Curve {
    let mut acc = C::Curve::identity();
    for base in bases {
      acc += base;
    }
    acc
  }

  let num_chunks = rayon::current_num_threads();
  let chunk_size = (one_indices.len() + num_chunks).div_ceil(num_chunks);

  let comm = one_indices
    .par_chunks(chunk_size)
    .into_par_iter()
    .map(|chunk| add_chunk(chunk.iter().map(|index| bases[*index])))
    .reduce(C::Curve::identity, |sum, evl| sum + evl);

  comm
}

// ==================================================================================
// One-hot sparse binary commitment: batch addition with structural optimizations
// ==================================================================================

/// Software prefetch hint for read access.
///
/// On x86_64, issues a `prefetcht0` instruction to bring data into all cache levels.
/// On other architectures, this is a no-op.
#[inline(always)]
fn prefetch_read<T>(ptr: *const T) {
  #[cfg(target_arch = "x86_64")]
  unsafe {
    core::arch::x86_64::_mm_prefetch(ptr as *const i8, core::arch::x86_64::_MM_HINT_T0);
  }
  #[cfg(not(target_arch = "x86_64"))]
  {
    let _ = ptr; // suppress unused warning
  }
}

/// Batch addition of affine points for one-hot structured vectors.
///
/// A one-hot vector of length `num_blocks * block_size` has exactly one non-zero (=1) entry
/// per block of size `block_size`. `block_offsets[i]` gives the offset within block `i`
/// (must satisfy `0 <= block_offsets[i] < block_size`).
///
/// The global index for block `i` is `i * block_size + block_offsets[i]`.
///
/// This function uses software prefetch hints to hide memory latency from strided access
/// into the (potentially large) generator array.
pub(crate) fn batch_add_one_hot<C: CurveAffine>(
  bases: &[C],
  block_size: usize,
  block_offsets: &[usize],
) -> C::Curve {
  assert!(block_size > 0, "block_size must be positive");
  let num_blocks = block_offsets.len();
  if num_blocks == 0 {
    return C::Curve::identity();
  }

  let num_threads = current_num_threads();

  if num_blocks > num_threads {
    let chunk_size = num_blocks.div_ceil(num_threads);

    let bucket = (0..num_blocks)
      .into_par_iter()
      .chunks(chunk_size)
      .map(|block_indices| {
        let mut bucket = BucketXYZZ::<C::Base>::zero();
        for (local_idx, block_idx) in block_indices.iter().enumerate() {
          // Prefetch the point for the next iteration to hide memory latency.
          // The stride between accessed generators is block_size * sizeof(C),
          // which is typically 64-128 bytes * block_size — well beyond a cache line.
          if local_idx + 1 < block_indices.len() {
            let next_block = block_indices[local_idx + 1];
            let next_global = next_block * block_size + block_offsets[next_block];
            prefetch_read(bases.as_ptr().wrapping_add(next_global));
          }
          let global_idx = block_idx * block_size + block_offsets[*block_idx];
          debug_assert!(global_idx < bases.len());
          bucket_add_affine::<C>(&mut bucket, &bases[global_idx]);
        }
        bucket
      })
      .reduce(BucketXYZZ::zero, |mut a, b| {
        a.add_assign_bucket(&b);
        a
      });

    bucket_to_curve::<C>(&bucket)
  } else {
    // Small number of blocks: sequential with prefetch
    let mut bucket = BucketXYZZ::<C::Base>::zero();
    for i in 0..num_blocks {
      if i + 1 < num_blocks {
        let next_global = (i + 1) * block_size + block_offsets[i + 1];
        prefetch_read(bases.as_ptr().wrapping_add(next_global));
      }
      let global_idx = i * block_size + block_offsets[i];
      bucket_add_affine::<C>(&mut bucket, &bases[global_idx]);
    }
    if bucket.is_zero() {
      C::Curve::identity()
    } else {
      bucket_to_curve::<C>(&bucket)
    }
  }
}

/// Precomputed table for one-hot block pairs.
///
/// For each consecutive pair of blocks `(2i, 2i+1)`, the table stores all `K * K`
/// possible sums `bases[2i*K + a] + bases[(2i+1)*K + b]` for `a, b in 0..K`.
/// This turns two point additions into a single table lookup.
///
/// Memory layout: `tables[pair_idx][a * K + b]` = `bases[2*pair_idx*K + a] + bases[(2*pair_idx+1)*K + b]`
///
/// For K=16, each pair needs 256 precomputed points (~16 KB per pair).
pub struct OneHotPrecompTable<C: CurveAffine> {
  /// Precomputed sums for each block pair: tables[pair][a * K + b]
  tables: Vec<Vec<C>>,
  /// Block size K
  block_size: usize,
  /// Total number of blocks
  num_blocks: usize,
  /// Leftover bases for the odd-one-out block (if num_blocks is odd)
  leftover_bases: Option<Vec<C>>,
}

impl<C: CurveAffine> OneHotPrecompTable<C> {
  /// Build precomputation tables for the given generator bases and block size.
  ///
  /// `bases` is the full generator array of length `num_blocks * block_size`.
  /// For each pair of consecutive blocks, precomputes all K² possible sums.
  ///
  /// Cost: `num_pairs * K² * (1 point addition)` precomputation.
  /// For K=16, num_blocks=1024: 512 pairs × 256 entries = 131072 additions.
  pub fn new(bases: &[C], block_size: usize, num_blocks: usize) -> Self {
    assert!(block_size > 0);
    assert!(bases.len() >= num_blocks * block_size);

    let num_pairs = num_blocks / 2;
    let k = block_size;

    // Precompute tables in parallel, one per pair
    let tables: Vec<Vec<C>> = (0..num_pairs)
      .into_par_iter()
      .map(|pair_idx| {
        let block_a_start = 2 * pair_idx * k;
        let block_b_start = (2 * pair_idx + 1) * k;

        let mut table = Vec::with_capacity(k * k);
        for a in 0..k {
          let pa = bases[block_a_start + a];
          for b in 0..k {
            let pb = bases[block_b_start + b];
            // Store the affine sum: pa + pb
            let sum: C = (C::Curve::identity() + pa + pb).into();
            table.push(sum);
          }
        }
        table
      })
      .collect();

    // Handle the odd block (if num_blocks is odd)
    let leftover_bases = if num_blocks % 2 == 1 {
      let last_block_start = (num_blocks - 1) * k;
      Some(bases[last_block_start..last_block_start + k].to_vec())
    } else {
      None
    };

    Self {
      tables,
      block_size: k,
      num_blocks,
      leftover_bases,
    }
  }

  /// Commit using precomputed tables.
  ///
  /// `block_offsets[i]` is the offset within block `i` (0 <= offset < K).
  ///
  /// For each pair of blocks, performs a single table lookup instead of two point additions.
  /// The final accumulation of pair results uses XYZZ coordinates with prefetch.
  pub fn batch_add(&self, block_offsets: &[usize]) -> C::Curve {
    assert_eq!(block_offsets.len(), self.num_blocks);
    let k = self.block_size;
    let num_pairs = self.num_blocks / 2;

    if num_pairs == 0 {
      // Only 0 or 1 blocks
      if self.num_blocks == 1 {
        if let Some(ref leftover) = self.leftover_bases {
          return C::Curve::identity() + leftover[block_offsets[0]];
        }
      }
      return C::Curve::identity();
    }

    let num_threads = current_num_threads();

    let pair_result = if num_pairs > num_threads {
      let chunk_size = num_pairs.div_ceil(num_threads);

      let bucket = (0..num_pairs)
        .into_par_iter()
        .chunks(chunk_size)
        .map(|pair_indices| {
          let mut bucket = BucketXYZZ::<C::Base>::zero();
          for (local_idx, &pair_idx) in pair_indices.iter().enumerate() {
            // Prefetch the next table entry
            if local_idx + 1 < pair_indices.len() {
              let next_pair = pair_indices[local_idx + 1];
              let next_a = block_offsets[2 * next_pair];
              let next_b = block_offsets[2 * next_pair + 1];
              let next_table_idx = next_a * k + next_b;
              prefetch_read(self.tables[next_pair].as_ptr().wrapping_add(next_table_idx));
            }

            let a = block_offsets[2 * pair_idx];
            let b = block_offsets[2 * pair_idx + 1];
            debug_assert!(a < k && b < k);
            let table_idx = a * k + b;
            bucket_add_affine::<C>(&mut bucket, &self.tables[pair_idx][table_idx]);
          }
          bucket
        })
        .reduce(BucketXYZZ::zero, |mut a, b| {
          a.add_assign_bucket(&b);
          a
        });

      bucket_to_curve::<C>(&bucket)
    } else {
      // Sequential with prefetch
      let mut bucket = BucketXYZZ::<C::Base>::zero();
      for pair_idx in 0..num_pairs {
        if pair_idx + 1 < num_pairs {
          let next_a = block_offsets[2 * (pair_idx + 1)];
          let next_b = block_offsets[2 * (pair_idx + 1) + 1];
          let next_table_idx = next_a * k + next_b;
          prefetch_read(
            self.tables[pair_idx + 1]
              .as_ptr()
              .wrapping_add(next_table_idx),
          );
        }

        let a = block_offsets[2 * pair_idx];
        let b = block_offsets[2 * pair_idx + 1];
        debug_assert!(a < k && b < k);
        let table_idx = a * k + b;
        bucket_add_affine::<C>(&mut bucket, &self.tables[pair_idx][table_idx]);
      }
      if bucket.is_zero() {
        C::Curve::identity()
      } else {
        bucket_to_curve::<C>(&bucket)
      }
    };

    // Add the odd block if present
    if let Some(ref leftover) = self.leftover_bases {
      let last_offset = block_offsets[self.num_blocks - 1];
      debug_assert!(last_offset < k);
      pair_result + leftover[last_offset]
    } else {
      pair_result
    }
  }
}

/// Precomputed table for one-hot block triples (t=3).
///
/// For each consecutive triple of blocks `(3i, 3i+1, 3i+2)`, the table stores all `K³`
/// possible sums `bases[3i*K + a] + bases[(3i+1)*K + b] + bases[(3i+2)*K + c]`.
/// This turns three point additions into a single table lookup.
///
/// For K=16, each triple needs 4096 precomputed points (~256 KB per triple).
pub struct OneHotPrecompTable3<C: CurveAffine> {
  /// Precomputed sums for each block triple: tables[triple][a * K² + b * K + c]
  tables: Vec<Vec<C>>,
  /// Block size K
  block_size: usize,
  /// Total number of blocks
  num_blocks: usize,
  /// Number of full triples
  num_triples: usize,
  /// Leftover blocks (0, 1, or 2 blocks that don't form a full triple)
  leftover_bases: Vec<Vec<C>>,
}

impl<C: CurveAffine> OneHotPrecompTable3<C> {
  /// Build t=3 precomputation tables.
  ///
  /// Cost: `num_triples * K³` point additions.
  /// For K=16, 2^20 blocks: 349,525 triples × 4096 entries ≈ 1.43 billion additions.
  pub fn new(bases: &[C], block_size: usize, num_blocks: usize) -> Self {
    assert!(block_size > 0);
    assert!(bases.len() >= num_blocks * block_size);

    let k = block_size;
    let num_triples = num_blocks / 3;
    let leftover_count = num_blocks % 3;

    // Precompute tables in parallel, one per triple
    let tables: Vec<Vec<C>> = (0..num_triples)
      .into_par_iter()
      .map(|triple_idx| {
        let ba = 3 * triple_idx * k;
        let bb = (3 * triple_idx + 1) * k;
        let bc = (3 * triple_idx + 2) * k;

        let mut table = Vec::with_capacity(k * k * k);
        for a in 0..k {
          let pa = bases[ba + a];
          for b in 0..k {
            let pab = C::Curve::identity() + pa + bases[bb + b];
            for c in 0..k {
              let sum: C = (pab + bases[bc + c]).into();
              table.push(sum);
            }
          }
        }
        table
      })
      .collect();

    // Collect leftover blocks (those not forming a full triple)
    let leftover_bases: Vec<Vec<C>> = (0..leftover_count)
      .map(|i| {
        let block_start = (num_triples * 3 + i) * k;
        bases[block_start..block_start + k].to_vec()
      })
      .collect();

    Self {
      tables,
      block_size: k,
      num_blocks,
      num_triples,
      leftover_bases,
    }
  }

  /// Commit using t=3 precomputed tables.
  pub fn batch_add(&self, block_offsets: &[usize]) -> C::Curve {
    assert_eq!(block_offsets.len(), self.num_blocks);
    let k = self.block_size;
    let k2 = k * k;

    if self.num_triples == 0 {
      // No full triples, just add leftovers
      let mut acc = C::Curve::identity();
      for (i, bases) in self.leftover_bases.iter().enumerate() {
        let block_idx = i;
        acc += bases[block_offsets[block_idx]];
      }
      return acc;
    }

    let num_threads = current_num_threads();

    let triple_result = if self.num_triples > num_threads {
      let chunk_size = self.num_triples.div_ceil(num_threads);

      let bucket = (0..self.num_triples)
        .into_par_iter()
        .chunks(chunk_size)
        .map(|triple_indices| {
          let mut bucket = BucketXYZZ::<C::Base>::zero();
          for (local_idx, &triple_idx) in triple_indices.iter().enumerate() {
            if local_idx + 1 < triple_indices.len() {
              let next = triple_indices[local_idx + 1];
              let na = block_offsets[3 * next];
              let nb = block_offsets[3 * next + 1];
              let nc = block_offsets[3 * next + 2];
              prefetch_read(
                self.tables[next]
                  .as_ptr()
                  .wrapping_add(na * k2 + nb * k + nc),
              );
            }

            let a = block_offsets[3 * triple_idx];
            let b = block_offsets[3 * triple_idx + 1];
            let c = block_offsets[3 * triple_idx + 2];
            debug_assert!(a < k && b < k && c < k);
            let table_idx = a * k2 + b * k + c;
            bucket_add_affine::<C>(&mut bucket, &self.tables[triple_idx][table_idx]);
          }
          bucket
        })
        .reduce(BucketXYZZ::zero, |mut a, b| {
          a.add_assign_bucket(&b);
          a
        });

      bucket_to_curve::<C>(&bucket)
    } else {
      let mut bucket = BucketXYZZ::<C::Base>::zero();
      for triple_idx in 0..self.num_triples {
        if triple_idx + 1 < self.num_triples {
          let next = triple_idx + 1;
          let na = block_offsets[3 * next];
          let nb = block_offsets[3 * next + 1];
          let nc = block_offsets[3 * next + 2];
          prefetch_read(
            self.tables[next]
              .as_ptr()
              .wrapping_add(na * k2 + nb * k + nc),
          );
        }

        let a = block_offsets[3 * triple_idx];
        let b = block_offsets[3 * triple_idx + 1];
        let c = block_offsets[3 * triple_idx + 2];
        debug_assert!(a < k && b < k && c < k);
        let table_idx = a * k2 + b * k + c;
        bucket_add_affine::<C>(&mut bucket, &self.tables[triple_idx][table_idx]);
      }
      if bucket.is_zero() {
        C::Curve::identity()
      } else {
        bucket_to_curve::<C>(&bucket)
      }
    };

    // Add leftover blocks
    let mut result = triple_result;
    for (i, bases) in self.leftover_bases.iter().enumerate() {
      let block_idx = self.num_triples * 3 + i;
      result += bases[block_offsets[block_idx]];
    }
    result
  }
}

// ==================================================================================
// Batch one-hot commit: parallel transpose (block-major order)
// ==================================================================================

/// Batch addition for multiple one-hot vectors simultaneously, using transposed iteration order.
///
/// Instead of computing each commit independently (commit-major order),
/// this iterates in block-major order: for each block, load the generators once
/// and distribute to all commits that select them. This reduces generator memory
/// reads from `num_commits × num_blocks` (random access) to `num_blocks × K` (sequential).
///
/// # Arguments
/// * `bases` — full generator array of length `num_blocks × block_size`
/// * `block_size` — K (elements per block)
/// * `all_offsets` — `all_offsets[commit_idx][block_idx]` = offset within that block
///
/// # Returns
/// One projective point per commit: `results[i]` = sum of selected generators for commit `i`.
pub fn batch_add_one_hot_transpose<C: CurveAffine>(
  bases: &[C],
  block_size: usize,
  all_offsets: &[&[usize]],
) -> Vec<C::Curve> {
  let num_commits = all_offsets.len();
  if num_commits == 0 {
    return vec![];
  }

  let num_blocks = all_offsets[0].len();
  let k = block_size;
  assert!(bases.len() >= num_blocks * k);
  for offsets in all_offsets.iter() {
    assert_eq!(offsets.len(), num_blocks);
  }

  let num_threads = current_num_threads();
  let blocks_per_thread = num_blocks.div_ceil(num_threads);

  // Each thread processes a range of blocks and maintains its own set of accumulators.
  // Accumulator per commit = one BucketXYZZ (~128 bytes).
  // 600 commits × 128 bytes = ~75 KB per thread → fits in L1/L2 cache.
  let partial_results: Vec<Vec<BucketXYZZ<C::Base>>> = (0..num_threads)
    .into_par_iter()
    .map(|thread_idx| {
      let block_start = thread_idx * blocks_per_thread;
      let block_end = (block_start + blocks_per_thread).min(num_blocks);
      if block_start >= block_end {
        return vec![BucketXYZZ::zero(); num_commits];
      }

      let mut accumulators: Vec<BucketXYZZ<C::Base>> = vec![BucketXYZZ::zero(); num_commits];

      for block_b in block_start..block_end {
        let base_idx = block_b * k;

        // Prefetch next block's generators
        if block_b + 1 < block_end {
          prefetch_read(bases.as_ptr().wrapping_add((block_b + 1) * k));
        }

        // For each commit, look up its offset for this block,
        // fetch the corresponding generator, and add to its accumulator.
        // No branching — every commit does exactly one addition per block.
        for (i, offsets) in all_offsets.iter().enumerate() {
          let v = offsets[block_b];
          debug_assert!(v < k);
          bucket_add_affine::<C>(&mut accumulators[i], &bases[base_idx + v]);
        }
      }

      accumulators
    })
    .collect();

  // Merge partial results from all threads: for each commit, sum across threads
  (0..num_commits)
    .into_par_iter()
    .map(|commit_idx| {
      let mut merged = BucketXYZZ::<C::Base>::zero();
      for thread_result in &partial_results {
        merged.add_assign_bucket(&thread_result[commit_idx]);
      }
      if merged.is_zero() {
        C::Curve::identity()
      } else {
        bucket_to_curve::<C>(&merged)
      }
    })
    .collect()
}

// ==================================================================================
// Batch one-hot commit: parallel transpose with precomputed tables
// ==================================================================================

/// Batch addition for multiple one-hot vectors using precomputed pair tables + transposed iteration.
///
/// Combines two optimizations:
/// 1. **Precomputed tables (t=2)**: each pair of blocks has K² precomputed sums,
///    so each commit needs only 1 table lookup per pair instead of 2 point additions.
/// 2. **Transposed iteration**: iterates in pair-major order across all commits,
///    keeping accumulators hot in cache.
///
/// This halves the number of point additions compared to plain transpose:
/// `num_commits × (num_blocks / 2)` additions instead of `num_commits × num_blocks`.
///
/// # Arguments
/// * `table` — precomputed `OneHotPrecompTable` (t=2) for the generators
/// * `all_offsets` — `all_offsets[commit_idx][block_idx]` = offset within that block
///
/// # Returns
/// One projective point per commit.
pub fn batch_add_one_hot_precomp_transpose<C: CurveAffine>(
  table: &OneHotPrecompTable<C>,
  all_offsets: &[&[usize]],
) -> Vec<C::Curve> {
  let num_commits = all_offsets.len();
  if num_commits == 0 {
    return vec![];
  }

  let k = table.block_size;
  let num_blocks = table.num_blocks;
  let num_pairs = num_blocks / 2;

  for offsets in all_offsets.iter() {
    assert_eq!(offsets.len(), num_blocks);
  }

  if num_pairs == 0 {
    // 0 or 1 blocks: handle trivially
    return all_offsets
      .iter()
      .map(|offsets| {
        if num_blocks == 1 {
          if let Some(ref leftover) = table.leftover_bases {
            return C::Curve::identity() + leftover[offsets[0]];
          }
        }
        C::Curve::identity()
      })
      .collect();
  }

  let num_threads = current_num_threads();
  let pairs_per_thread = num_pairs.div_ceil(num_threads);

  // Each thread processes a range of pairs, maintaining its own 600 accumulators.
  let partial_results: Vec<Vec<BucketXYZZ<C::Base>>> = (0..num_threads)
    .into_par_iter()
    .map(|thread_idx| {
      let pair_start = thread_idx * pairs_per_thread;
      let pair_end = (pair_start + pairs_per_thread).min(num_pairs);
      if pair_start >= pair_end {
        return vec![BucketXYZZ::zero(); num_commits];
      }

      let mut accumulators: Vec<BucketXYZZ<C::Base>> = vec![BucketXYZZ::zero(); num_commits];

      for pair_idx in pair_start..pair_end {
        // Prefetch next pair's table
        if pair_idx + 1 < pair_end {
          prefetch_read(table.tables[pair_idx + 1].as_ptr());
        }

        let pair_table = &table.tables[pair_idx];

        // For each commit, compute the table index from its two block offsets
        // and add the precomputed sum to the accumulator. One lookup per pair.
        for (i, offsets) in all_offsets.iter().enumerate() {
          let a = offsets[2 * pair_idx];
          let b = offsets[2 * pair_idx + 1];
          debug_assert!(a < k && b < k);
          let table_idx = a * k + b;
          bucket_add_affine::<C>(&mut accumulators[i], &pair_table[table_idx]);
        }
      }

      accumulators
    })
    .collect();

  // Merge partial results from all threads
  let mut results: Vec<C::Curve> = (0..num_commits)
    .into_par_iter()
    .map(|commit_idx| {
      let mut merged = BucketXYZZ::<C::Base>::zero();
      for thread_result in &partial_results {
        merged.add_assign_bucket(&thread_result[commit_idx]);
      }
      if merged.is_zero() {
        C::Curve::identity()
      } else {
        bucket_to_curve::<C>(&merged)
      }
    })
    .collect();

  // Handle the leftover odd block (if num_blocks is odd)
  if let Some(ref leftover) = table.leftover_bases {
    let last_block = num_blocks - 1;
    for (i, offsets) in all_offsets.iter().enumerate() {
      let offset = offsets[last_block];
      debug_assert!(offset < k);
      results[i] += leftover[offset];
    }
  }

  results
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::provider::{
    bn256_grumpkin::{bn256, grumpkin},
    pasta::{pallas, vesta},
    secp_secq::{secp256k1, secq256k1},
  };
  use ff::Field;
  use halo2curves::{group::Group, CurveAffine};
  use rand_core::OsRng;

  fn test_general_msm_with<F: Field, A: CurveAffine<ScalarExt = F>>() {
    let n = 8;
    let coeffs = (0..n).map(|_| F::random(OsRng)).collect::<Vec<_>>();
    let bases = (0..n)
      .map(|_| A::from(A::generator() * F::random(OsRng)))
      .collect::<Vec<_>>();

    assert_eq!(coeffs.len(), bases.len());
    let naive = coeffs
      .iter()
      .zip(bases.iter())
      .fold(A::CurveExt::identity(), |acc, (coeff, base)| {
        acc + *base * coeff
      });
    let msm = msm(&coeffs, &bases);

    assert_eq!(naive, msm)
  }

  #[test]
  fn test_general_msm() {
    test_general_msm_with::<pallas::Scalar, pallas::Affine>();
    test_general_msm_with::<vesta::Scalar, vesta::Affine>();
    test_general_msm_with::<bn256::Scalar, bn256::Affine>();
    test_general_msm_with::<grumpkin::Scalar, grumpkin::Affine>();
    test_general_msm_with::<secp256k1::Scalar, secp256k1::Affine>();
    test_general_msm_with::<secq256k1::Scalar, secq256k1::Affine>();
  }

  fn test_msm_ux_with<F: PrimeField, A: CurveAffine<ScalarExt = F>>() {
    let n = 8;
    let bases = (0..n)
      .map(|_| A::from(A::generator() * F::random(OsRng)))
      .collect::<Vec<_>>();

    for bit_width in [1, 4, 8, 10, 16, 20, 32, 40, 64] {
      println!("bit_width: {bit_width}");
      assert!(bit_width <= 64); // Ensure we don't overflow F::from
      let mask = if bit_width == 64 {
        u64::MAX
      } else {
        (1u64 << bit_width) - 1
      };
      let coeffs: Vec<u64> = (0..n)
        .map(|_| rand::random::<u64>() & mask)
        .collect::<Vec<_>>();
      let coeffs_scalar: Vec<F> = coeffs.iter().map(|b| F::from(*b)).collect::<Vec<_>>();
      let general = msm(&coeffs_scalar, &bases);
      let integer = msm_small(&coeffs, &bases);

      assert_eq!(general, integer);
    }
  }

  #[test]
  fn test_msm_ux() {
    test_msm_ux_with::<pallas::Scalar, pallas::Affine>();
    test_msm_ux_with::<vesta::Scalar, vesta::Affine>();
    test_msm_ux_with::<bn256::Scalar, bn256::Affine>();
    test_msm_ux_with::<grumpkin::Scalar, grumpkin::Affine>();
    test_msm_ux_with::<secp256k1::Scalar, secp256k1::Affine>();
    test_msm_ux_with::<secq256k1::Scalar, secq256k1::Affine>();
  }

  // ==================================================================================
  // One-hot batch add tests
  // ==================================================================================

  /// Test that batch_add_one_hot matches naive computation (sum of selected generators).
  fn test_one_hot_correctness_with<A: CurveAffine>() {
    let k: usize = 16;
    for &num_blocks in &[1, 2, 3, 7, 16, 64, 128] {
      let n = num_blocks * k;
      let bases: Vec<A> = (0..n)
        .map(|_| A::from(A::generator() * A::Scalar::random(OsRng)))
        .collect();

      // Generate random one-hot offsets
      let offsets: Vec<usize> = (0..num_blocks)
        .map(|_| (rand::random::<usize>()) % k)
        .collect();

      // Compute with one-hot function
      let result_one_hot = batch_add_one_hot::<A>(&bases, k, &offsets);

      // Compute naively: sum of bases[i * k + offsets[i]] for each block
      let naive = offsets
        .iter()
        .enumerate()
        .fold(A::CurveExt::identity(), |acc, (i, &offset)| {
          acc + bases[i * k + offset]
        });

      assert_eq!(
        result_one_hot, naive,
        "one-hot mismatch at num_blocks={num_blocks}"
      );

      // Also verify it matches batch_add with equivalent indices
      let global_indices: Vec<usize> = offsets
        .iter()
        .enumerate()
        .map(|(i, &offset)| i * k + offset)
        .collect();
      let result_sparse = batch_add::<A>(&bases, &global_indices);
      assert_eq!(
        result_one_hot, result_sparse,
        "one-hot vs sparse mismatch at num_blocks={num_blocks}"
      );
    }
  }

  #[test]
  fn test_one_hot_correctness() {
    test_one_hot_correctness_with::<pallas::Affine>();
    test_one_hot_correctness_with::<bn256::Affine>();
    test_one_hot_correctness_with::<secp256k1::Affine>();
  }

  /// Test that all offsets pointing to the same position still works correctly.
  fn test_one_hot_same_offset_with<A: CurveAffine>() {
    let k: usize = 16;
    let num_blocks = 32;
    let n = num_blocks * k;
    let bases: Vec<A> = (0..n)
      .map(|_| A::from(A::generator() * A::Scalar::random(OsRng)))
      .collect();

    // All offsets = 0 (first element of each block)
    let offsets = vec![0usize; num_blocks];
    let result = batch_add_one_hot::<A>(&bases, k, &offsets);
    let naive = (0..num_blocks).fold(A::CurveExt::identity(), |acc, i| acc + bases[i * k]);
    assert_eq!(result, naive);

    // All offsets = k-1 (last element of each block)
    let offsets = vec![k - 1; num_blocks];
    let result = batch_add_one_hot::<A>(&bases, k, &offsets);
    let naive = (0..num_blocks).fold(A::CurveExt::identity(), |acc, i| acc + bases[i * k + k - 1]);
    assert_eq!(result, naive);
  }

  #[test]
  fn test_one_hot_same_offset() {
    test_one_hot_same_offset_with::<pallas::Affine>();
    test_one_hot_same_offset_with::<bn256::Affine>();
  }

  /// Test precomputed table correctness against naive and batch_add_one_hot.
  fn test_one_hot_precomp_with<A: CurveAffine>() {
    let k: usize = 16;
    for &num_blocks in &[1, 2, 3, 4, 7, 15, 16, 33, 64] {
      let n = num_blocks * k;
      let bases: Vec<A> = (0..n)
        .map(|_| A::from(A::generator() * A::Scalar::random(OsRng)))
        .collect();

      let table = OneHotPrecompTable::new(&bases, k, num_blocks);

      let offsets: Vec<usize> = (0..num_blocks)
        .map(|_| (rand::random::<usize>()) % k)
        .collect();

      let result_precomp = table.batch_add(&offsets);
      let result_one_hot = batch_add_one_hot::<A>(&bases, k, &offsets);

      assert_eq!(
        result_precomp, result_one_hot,
        "precomp mismatch at num_blocks={num_blocks}"
      );
    }
  }

  #[test]
  fn test_one_hot_precomp() {
    test_one_hot_precomp_with::<pallas::Affine>();
    test_one_hot_precomp_with::<bn256::Affine>();
    test_one_hot_precomp_with::<secp256k1::Affine>();
  }

  /// Test single-block edge cases for precomputed tables.
  #[test]
  fn test_one_hot_precomp_single_block() {
    type A = pallas::Affine;
    let k: usize = 16;
    let bases: Vec<A> = (0..k)
      .map(|_| A::from(A::generator() * pallas::Scalar::random(OsRng)))
      .collect();

    let table = OneHotPrecompTable::new(&bases, k, 1);
    for offset in 0..k {
      let result = table.batch_add(&[offset]);
      // Verify against batch_add_one_hot which we already tested
      let expected = batch_add_one_hot::<A>(&bases, k, &[offset]);
      assert_eq!(result, expected);
    }
  }

  /// Test t=3 precomputed table correctness.
  fn test_one_hot_precomp3_with<A: CurveAffine>() {
    let k: usize = 16;
    for &num_blocks in &[3, 4, 6, 7, 9, 15, 16, 33, 64] {
      let n = num_blocks * k;
      let bases: Vec<A> = (0..n)
        .map(|_| A::from(A::generator() * A::Scalar::random(OsRng)))
        .collect();

      let table = OneHotPrecompTable3::new(&bases, k, num_blocks);

      let offsets: Vec<usize> = (0..num_blocks)
        .map(|_| (rand::random::<usize>()) % k)
        .collect();

      let result_precomp = table.batch_add(&offsets);
      let result_one_hot = batch_add_one_hot::<A>(&bases, k, &offsets);

      assert_eq!(
        result_precomp, result_one_hot,
        "precomp3 mismatch at num_blocks={num_blocks}"
      );
    }
  }

  #[test]
  fn test_one_hot_precomp3() {
    test_one_hot_precomp3_with::<pallas::Affine>();
    test_one_hot_precomp3_with::<bn256::Affine>();
    test_one_hot_precomp3_with::<secp256k1::Affine>();
  }

  /// Test batch transpose correctness: results must match individual one-hot commits.
  fn test_batch_transpose_with<A: CurveAffine>() {
    let k: usize = 16;
    let num_blocks = 64;
    let n = num_blocks * k;
    let bases: Vec<A> = (0..n)
      .map(|_| A::from(A::generator() * A::Scalar::random(OsRng)))
      .collect();

    for &num_commits in &[1, 2, 10, 50, 100] {
      // Generate random offsets for each commit
      let all_offsets: Vec<Vec<usize>> = (0..num_commits)
        .map(|_| {
          (0..num_blocks)
            .map(|_| rand::random::<usize>() % k)
            .collect()
        })
        .collect();

      let all_offsets_refs: Vec<&[usize]> = all_offsets.iter().map(|v| v.as_slice()).collect();

      // Compute with transpose batch
      let batch_results = batch_add_one_hot_transpose::<A>(&bases, k, &all_offsets_refs);

      // Compute each individually
      for (i, offsets) in all_offsets.iter().enumerate() {
        let individual = batch_add_one_hot::<A>(&bases, k, offsets);
        assert_eq!(
          batch_results[i], individual,
          "transpose mismatch at commit={i}, num_commits={num_commits}"
        );
      }
    }
  }

  #[test]
  fn test_batch_transpose() {
    test_batch_transpose_with::<pallas::Affine>();
    test_batch_transpose_with::<bn256::Affine>();
  }

  /// Test precomp+transpose correctness.
  fn test_batch_precomp_transpose_with<A: CurveAffine>() {
    let k: usize = 16;
    let num_blocks = 64;
    let n = num_blocks * k;
    let bases: Vec<A> = (0..n)
      .map(|_| A::from(A::generator() * A::Scalar::random(OsRng)))
      .collect();

    let table = OneHotPrecompTable::new(&bases, k, num_blocks);

    for &num_commits in &[1, 2, 10, 50, 100] {
      let all_offsets: Vec<Vec<usize>> = (0..num_commits)
        .map(|_| {
          (0..num_blocks)
            .map(|_| rand::random::<usize>() % k)
            .collect()
        })
        .collect();

      let all_offsets_refs: Vec<&[usize]> = all_offsets.iter().map(|v| v.as_slice()).collect();

      // Compute with precomp+transpose
      let batch_results = batch_add_one_hot_precomp_transpose::<A>(&table, &all_offsets_refs);

      // Compute each individually
      for (i, offsets) in all_offsets.iter().enumerate() {
        let individual = batch_add_one_hot::<A>(&bases, k, offsets);
        assert_eq!(
          batch_results[i], individual,
          "precomp_transpose mismatch at commit={i}, num_commits={num_commits}"
        );
      }
    }
  }

  #[test]
  fn test_batch_precomp_transpose() {
    test_batch_precomp_transpose_with::<pallas::Affine>();
    test_batch_precomp_transpose_with::<bn256::Affine>();
  }
}
