#![allow(non_snake_case)]
//! Per-phase timing + peak-RSS profile of `ppsnark::prove`, comparing the
//! default Logup-GKR memory-check against the `logup-no-gkr` inverse-logup path.
//!
//! Custom harness (not criterion): drives `DirectSNARK::prove` at several
//! circuit sizes, repeating each `REPEATS` times to average out noise, and uses
//! the library's `phase-profile` hook to record wall-clock time per phase and
//! peak RSS (`ru_maxrss`) at each phase boundary.
//!
//! Which memory-check path is measured depends on the build feature:
//!   cargo bench --bench ppsnark_phases --features phase-profile
//!   cargo bench --bench ppsnark_phases --features "phase-profile,logup-no-gkr"
//!
//! Emit machine-readable JSON lines (one object per size) by setting
//! `PHASE_JSON=1`.

use nova_snark::{
  provider::Bn256EngineKZG,
  spartan::direct::DirectSNARK,
  spartan::phase_profile,
  traits::{circuit::NonTrivialCircuit, Engine},
};
use std::sync::{Arc, Mutex};
use std::time::Instant;

type E = Bn256EngineKZG;
type EE = nova_snark::provider::hyperkzg::EvaluationEngine<E>;
type S = nova_snark::spartan::ppsnark::RelaxedR1CSSNARK<E, EE>;

/// Phase boundary labels emitted by `ppsnark::prove`, in order. The interval
/// *ending* at label `i+1` is attributed to the phase named by label `i+1`
/// (except "start", which just resets the clock).
const PHASES: &[&str] = &[
  "outer_sc",
  "eval_oracles",
  "commit_L",
  "mem_check",
  "witness_bound",
  "inner_batch",
  "final_evals",
  "open",
];

/// Peak RSS in bytes, via `getrusage(RUSAGE_SELF).ru_maxrss`. On Linux the
/// field is KiB; on macOS it is bytes. Monotonic non-decreasing peak.
fn peak_rss_bytes() -> u64 {
  let mut usage: libc::rusage = unsafe { std::mem::zeroed() };
  let ret = unsafe { libc::getrusage(libc::RUSAGE_SELF, &mut usage) };
  if ret != 0 {
    return 0;
  }
  let maxrss = usage.ru_maxrss as u64;
  if cfg!(target_os = "macos") {
    maxrss
  } else {
    maxrss * 1024
  }
}

/// Accumulated per-phase samples across repeats.
#[derive(Default, Clone)]
struct Acc {
  /// nanoseconds per phase, one Vec per phase
  times_ns: Vec<Vec<u128>>,
  /// peak RSS (bytes) observed at each phase boundary, one Vec per phase
  rss_peak: Vec<Vec<u64>>,
}

fn median(xs: &mut [u128]) -> u128 {
  xs.sort_unstable();
  let n = xs.len();
  if n == 0 {
    0
  } else if n % 2 == 1 {
    xs[n / 2]
  } else {
    (xs[n / 2 - 1] + xs[n / 2]) / 2
  }
}

fn max_u64(xs: &[u64]) -> u64 {
  xs.iter().copied().max().unwrap_or(0)
}

fn label() -> &'static str {
  if cfg!(feature = "logup-no-gkr") {
    "no-gkr"
  } else {
    "gkr"
  }
}

fn main() {
  // Sizes: mirror benches/ppsnark.rs but keep the run tractable. Override with
  // PHASE_SIZES="16384,65536".
  let sizes: Vec<usize> = std::env::var("PHASE_SIZES")
    .ok()
    .map(|s| {
      s.split(',')
        .filter_map(|t| t.trim().parse().ok())
        .collect::<Vec<_>>()
    })
    .filter(|v: &Vec<usize>| !v.is_empty())
    .unwrap_or_else(|| vec![16384, 65536, 262144, 1048576]);

  let repeats: usize = std::env::var("PHASE_REPEATS")
    .ok()
    .and_then(|s| s.parse().ok())
    .unwrap_or(10);

  let json = std::env::var("PHASE_JSON").is_ok();

  eprintln!(
    "# ppsnark prove phase profile — path={}, repeats={}, sizes={:?}",
    label(),
    repeats,
    sizes
  );

  for &num_cons in &sizes {
    let circuit = NonTrivialCircuit::<<E as Engine>::Scalar>::new(num_cons);
    let input = vec![<E as Engine>::Scalar::from(42)];
    let (pk, _vk) =
      DirectSNARK::<E, S, NonTrivialCircuit<<E as Engine>::Scalar>>::setup(circuit.clone())
        .unwrap();

    // Warm up once (allocator, caches) — not recorded.
    let _ = DirectSNARK::prove(&pk, circuit.clone(), &input).unwrap();

    let acc = Arc::new(Mutex::new(Acc {
      times_ns: vec![Vec::with_capacity(repeats); PHASES.len()],
      rss_peak: vec![Vec::with_capacity(repeats); PHASES.len()],
    }));

    // Shared clock, reset at "start", advanced at each boundary.
    let last = Arc::new(Mutex::new(Instant::now()));
    let idx = Arc::new(Mutex::new(0usize));

    {
      let acc = acc.clone();
      let last = last.clone();
      let idx = idx.clone();
      phase_profile::set_hook(Some(Box::new(move |lbl: &str| {
        let now = Instant::now();
        if lbl == "start" {
          *last.lock().unwrap() = now;
          *idx.lock().unwrap() = 0;
          return;
        }
        let mut l = last.lock().unwrap();
        let dt = now.duration_since(*l).as_nanos();
        *l = now;
        let mut i = idx.lock().unwrap();
        let phase_i = *i;
        *i += 1;
        if phase_i < PHASES.len() {
          debug_assert_eq!(PHASES[phase_i], lbl, "phase order mismatch");
          let mut a = acc.lock().unwrap();
          a.times_ns[phase_i].push(dt);
          a.rss_peak[phase_i].push(peak_rss_bytes());
        }
      })));
    }

    for _ in 0..repeats {
      let _ = DirectSNARK::prove(&pk, circuit.clone(), &input).unwrap();
    }

    phase_profile::set_hook(None);

    let a = acc.lock().unwrap();
    let mut total_ns = 0u128;
    // (name, median_ns, min_ns, max_ns, rss_peak_bytes)
    let mut rows: Vec<(String, u128, u128, u128, u64)> = Vec::new();
    for (i, name) in PHASES.iter().enumerate() {
      let mut t = a.times_ns[i].clone();
      let med = median(&mut t);
      let mn = t.first().copied().unwrap_or(0); // sorted by median()
      let mx = t.last().copied().unwrap_or(0);
      let rss = max_u64(&a.rss_peak[i]);
      total_ns += med;
      rows.push((name.to_string(), med, mn, mx, rss));
    }
    let peak = rows.iter().map(|r| r.4).max().unwrap_or(0);

    if json {
      let phases_json: Vec<String> = rows
        .iter()
        .map(|(n, t, mn, mx, r)| {
          format!(
            "{{\"phase\":\"{}\",\"median_ms\":{:.3},\"min_ms\":{:.3},\"max_ms\":{:.3},\"rss_peak_mb\":{:.1}}}",
            n,
            *t as f64 / 1e6,
            *mn as f64 / 1e6,
            *mx as f64 / 1e6,
            *r as f64 / 1048576.0
          )
        })
        .collect();
      println!(
        "{{\"path\":\"{}\",\"num_cons\":{},\"total_ms\":{:.3},\"peak_rss_mb\":{:.1},\"phases\":[{}]}}",
        label(),
        num_cons,
        total_ns as f64 / 1e6,
        peak as f64 / 1048576.0,
        phases_json.join(",")
      );
    } else {
      println!(
        "\n=== path={} num_cons={} (2^{}) repeats={} ===",
        label(),
        num_cons,
        (num_cons as f64).log2().round() as u32,
        repeats
      );
      println!(
        "  {:<14} {:>11} {:>11} {:>11} {:>7} {:>12}",
        "phase", "median_ms", "min_ms", "max_ms", "pct", "rss_peak_mb"
      );
      for (n, t, mn, mx, r) in &rows {
        let pct = if total_ns > 0 {
          *t as f64 / total_ns as f64 * 100.0
        } else {
          0.0
        };
        println!(
          "  {:<14} {:>11.3} {:>11.3} {:>11.3} {:>6.1}% {:>12.1}",
          n,
          *t as f64 / 1e6,
          *mn as f64 / 1e6,
          *mx as f64 / 1e6,
          pct,
          *r as f64 / 1048576.0
        );
      }
      println!(
        "  {:<14} {:>11.3} {:>11} {:>11} {:>7} {:>12.1}",
        "TOTAL",
        total_ns as f64 / 1e6,
        "",
        "",
        "",
        peak as f64 / 1048576.0
      );
    }
  }
}
