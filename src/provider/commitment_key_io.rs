use std::{
  io::{self, Read, Write},
  path::Path, sync::mpsc::{channel, Receiver},
};

use ff::PrimeField;
use halo2curves::{bn256::Bn256, group::Group, CurveAffine};
use rayon::iter::{IntoParallelIterator, IntoParallelRefIterator, ParallelIterator};

use crate::traits::Engine;

use super::{
  bn256_grumpkin::{bn256, grumpkin},
  pasta::pallas,
  secp_secq::{secp256k1, secq256k1},
  traits::DlogGroup,
};

fn encode<G: DlogGroup>() -> Vec<u8> {
  let one = G::gen().affine();
  let zero = G::zero().affine();
  let mut encoded = G::encode(&one);
  encoded.extend(G::encode(&zero));
  encoded
}

/// Returns the ID of the DLog group
pub fn id_of<G: DlogGroup>() -> Option<u8> {
  let bn256 = encode::<bn256::Point>();
  let grumpkin = encode::<grumpkin::Point>();
  let palas = encode::<pallas::Point>();
  let vesta = encode::<pallas::Point>();
  let secp = encode::<secp256k1::Point>();
  let secq = encode::<secq256k1::Point>();

  let self_encode = encode::<G>();

  if bn256 == self_encode {
    Some(1)
  } else if grumpkin == self_encode {
    Some(2)
  } else if palas == self_encode {
    Some(3)
  } else if vesta == self_encode {
    Some(4)
  } else if secp == self_encode {
    Some(5)
  } else if secq == self_encode {
    Some(6)
  } else {
    None
  }
}

/// A trait for saving and loading commitment keys
pub trait CommitmentKeyIO: Sized {
  /// Saves the commitment key to the provided writer
  fn save_to(&self, writer: &mut impl Write) -> Result<(), io::Error>;
  /// Loads the commitment key from the provided reader
  fn load_from(reader: &mut impl Read, needed_num_ck: usize) -> Result<Self, io::Error>;
  /// Checks the sanity of the key file
  fn check_sanity_of_file(path: impl AsRef<Path>, required_len_ck_list: usize) -> bool;
}

use rayon::spawn;

fn read_one<R: Read>(
  mut reader: R,
) -> io::Result<(R, Receiver<bn256::Affine>)>
{
  let mut xy = [0u8; 64];
  reader.read_exact(&mut xy)?;

  let (tx, rx) = channel();
  spawn(move || {
    let mut x = [0u8; 32];
    let mut y = [0u8; 32];

    x.copy_from_slice(&xy[..32]);
    y.copy_from_slice(&xy[32..]);

    let x = bn256::Base::from_repr(x.into()).unwrap();
    let y = bn256::Base::from_repr(y.into()).unwrap();
    let p = bn256::Affine::from_xy(x, y).unwrap();

    tx.send(p).unwrap();
  });

  Ok((reader, rx))
}

fn read_one3<R: Read>(
  mut reader: R,
) -> io::Result<(R, [u8; 64])>
{
  let mut xy = [0u8; 64];
  reader.read_exact(&mut xy)?;

  Ok((reader, xy))
}

pub fn load_ck_vec2(
  reader: &mut impl Read,
  ck_len: usize,
) -> Vec<bn256::Affine>
{
  let mut rxs = Vec::with_capacity(ck_len);
  let mut reader = reader;
  for _ in 0..ck_len {
    let (reader2, rx) = read_one(reader).unwrap();
    rxs.push(rx);
    reader = reader2;
  }
  rxs.into_par_iter().map(|rx| rx.recv().unwrap()).collect()
  // rxs.into_iter().map(|rx| rx.recv().unwrap()).collect()
}

pub fn load_ck_vec3(
  reader: &mut impl Read,
  ck_len: usize,
) -> Vec<bn256::Affine>
{
  let mut data = Vec::with_capacity(ck_len);
  let mut reader = reader;
  for _ in 0..ck_len {
    let (reader2, d) = read_one3(reader).unwrap();
    data.push(d);
    reader = reader2;
  }
  data.into_par_iter().map(|rx| {
    let x = bn256::Base::from_repr(rx[..32].into()).unwrap();
    let y = bn256::Base::from_repr(rx[32..].into()).unwrap();
    bn256::Affine::from_xy(x, y).unwrap()
  }).collect()
  // rxs.into_iter().map(|rx| rx.recv().unwrap()).collect()
}


pub fn load_ck_vec<E: Engine>(
  reader: &mut impl Read,
  ck_len: usize,
) -> Vec<<E::GE as DlogGroup>::AffineGroupElement>
where
  E::GE: DlogGroup,
{
  let mut ck = Vec::with_capacity(ck_len);
  for _ in 0..ck_len {
    let p = <E::GE as DlogGroup>::decode_from(reader).unwrap();
    ck.push(p);
  }
  ck
}

pub fn write_ck_vec<E: Engine>(
  writer: &mut impl Write,
  ck: &[<E::GE as DlogGroup>::AffineGroupElement],
) -> Result<(), io::Error>
where
  E::GE: DlogGroup,
{
  let chunk_size = rayon::current_num_threads() * 1000;
  ck.chunks(chunk_size).for_each(|cks| {
    let bins = cks
      .par_iter()
      .map(<E::GE as DlogGroup>::encode)
      .flatten()
      .collect::<Vec<_>>();
    writer.write_all(&bins).unwrap();
  });
  Ok(())
}
