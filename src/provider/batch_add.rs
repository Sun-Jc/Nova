use ff::Field;
use halo2curves::CurveAffine;
use halo2curves::CurveExt;

#[inline(always)]
pub(crate) fn batch_add_custom<C: CurveAffine>(n: usize, points: &mut [(C::Base, C::Base)]) -> C {
  if n == 1 {
    return C::from_xy(points[0].0, points[0].1).unwrap();
  } else if n == 2 {
    let p1 = C::from_xy(points[0].0, points[0].1).unwrap();
    let p2 = C::from_xy(points[1].0, points[1].1).unwrap();
    return (p1 + p2).into();
  }
  assert!(n > 2);

  let mut n = n;
  let (points, rest) = if !n.is_power_of_two() {
    let next_power = n.next_power_of_two();
    n = next_power / 2;
    let (points, rest) = points.split_at_mut(n);
    (points, Some(rest))
  } else {
    (points, None)
  };

  let mut buffer = vec![C::Base::ZERO; n];
  let round = n.next_power_of_two().trailing_zeros();

  for r in (1..=round).rev() {
    let mid = 1 << (r - 1);

    let mut acc = C::Base::ONE;

    for i in 0..mid {
      let (x1, y1) = points[i];
      let (x2, y2) = points[i + mid];

      let denom = match x1 == x2 && y1 == y2 {
        true => y1 + y1,
        false => x2 - x1,
      };

      let i2 = i * 2;
      buffer[i2] = denom;
      buffer[i2 + 1] = acc;

      acc *= denom;
    }

    acc = acc.invert().unwrap();

    for i in (0..mid).rev() {
      let i2: usize = i * 2;

      let inv = buffer[i2 + 1] * acc;

      acc *= buffer[i2];

      let (x1, y1) = points[i];
      let (x2, y2) = points[i + mid];

      let nom = match x1 == x2 && y1 == y2 {
        true => {
          let x1_square = x1 * x1;
          x1_square + x1_square + x1_square + C::a()
        }
        false => y2 - y1,
      };

      let m = nom * inv;
      let x3 = m * m - x1 - x2;
      let y3 = m * (x1 - x3) - y1;

      points[i] = (x3, y3);
    }
  }

  let res = points[0];
  let res = C::from_xy(res.0, res.1).unwrap();

  if let Some(rest) = rest {
    (res + batch_add_custom::<C>(rest.len(), rest)).into()
  } else {
    res
  }
}
