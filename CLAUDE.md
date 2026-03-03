# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Nova is a high-speed recursive SNARK (Succinct Non-Interactive Argument of Knowledge) library implementing Incrementally Verifiable Computation (IVC) via folding schemes. Developed by Microsoft Research, written in Rust (edition 2021, MSRV 1.79.0).

## Build & Test Commands

```bash
# Build
cargo build --release

# Run all tests (use --release, crypto tests are slow in debug)
cargo test --release

# Run a single test
cargo test --release <test_name>
# e.g. cargo test --release test_general_msm

# Run tests for experimental NeutronNova
cargo test --release --features experimental neutron

# Lint
cargo clippy --all-targets -- -D warnings

# Format check
cargo fmt --all -- --check

# Build docs
cargo doc --no-deps --document-private-items

# Benchmarks (requires cargo-criterion)
cargo criterion --bench recursive-snark
cargo criterion --bench compressed-snark
cargo criterion --bench commit

# Run examples
cargo run --release --example minroot
cargo run --release --example hashchain
```

## Feature Flags

| Feature | Purpose |
|---|---|
| `io` (default) | Load/save commitment keys from Powers of Tau files |
| `evm` | EVM-compatible big-endian serialization |
| `experimental` | NeutronNova variant |
| `test-utils` | Insecure random tau generation (testing only) |
| `flamegraph` | pprof2 criterion flamegraph profiling |
| `blitzar` | GPU-accelerated MSM for BN256 via blitzar crate |

WASM target: `cargo build --no-default-features --target wasm32-unknown-unknown`

## Code Style

- **rustfmt**: 2-space indentation, Unix newlines, `use_try_shorthand = true`
- **`#![forbid(unsafe_code)]`** — no unsafe code anywhere
- **`#![allow(non_snake_case)]`** — uppercase variable names are used extensively for mathematical notation (e.g., `A`, `B`, `L`, `R`, `T`, `W`)
- **`#![deny(missing_docs)]`** — all public items must be documented
- Clippy: `type-complexity-threshold = 9999`, `too-many-arguments-threshold = 20`
- No `print_stdout`/`print_stderr` outside tests

## Architecture

### Generics & Curve Cycles

Everything is generic over `Engine` types. Nova requires **curve cycles** where `E1::Base == E2::Scalar` and vice versa, enabling efficient in-circuit verification. The primary generic pattern is:

```rust
PublicParams<E1, E2, C>       // E1 = primary engine, E2 = secondary, C = step circuit
RecursiveSNARK<E1, E2, C>     // Incremental proof object
CompressedSNARK<E1, E2, S1, S2>  // S1, S2 = SNARK traits for compression
```

### Key Modules

- **`src/nova/`** — Core IVC: `PublicParams`, `RecursiveSNARK`, `CompressedSNARK`, NIFS (Non-Interactive Folding Scheme). The augmented circuit lives in `nova/circuit/`.
- **`src/spartan/`** — Spartan SNARK backend. `snark.rs` (non-preprocessing) and `ppsnark.rs` (MicroSpartan preprocessing variant). Contains sumcheck protocol and multilinear polynomial types.
- **`src/r1cs/`** — R1CS constraint system: `R1CSShape` (matrices A, B, C in CSR sparse format), `R1CSInstance`, `R1CSWitness`, and their relaxed variants for folding.
- **`src/frontend/`** — Circuit DSL: `ConstraintSystem` trait, `LinearCombination`, `Variable`. Gadgets for boolean, arithmetic, SHA256, and Poseidon hash (ported from Neptune, ~15 files).
- **`src/provider/`** — Concrete implementations of engine traits for 6 curve pairs. Contains commitment engines, MSM, hash functions, and transcript implementations.
- **`src/traits/`** — Core trait definitions: `Engine`, `Group`, `ROTrait`, `TranscriptEngineTrait`, `CommitmentEngineTrait`, `StepCircuit`, `RelaxedR1CSSNARKTrait`.
- **`src/gadgets/`** — In-circuit gadgets for elliptic curve operations (`ecc.rs`, 45KB) and non-native field arithmetic.
- **`src/neutron/`** — Experimental NeutronNova variant (feature-gated behind `experimental`).

### Provider Architecture (Curves & Commitment Schemes)

Curve implementations use two macros in `provider/traits.rs`:
- **`impl_traits!`** — Full implementation including `DlogGroupExt` (MSM dispatch to `provider/msm.rs`)
- **`impl_traits_no_dlog_ext!`** — Everything except MSM, allowing custom MSM (e.g., BN256 with blitzar GPU)

Available engines and their commitment schemes:

| Engine | Curves | Commitment |
|---|---|---|
| `Bn256EngineKZG` | BN256/Grumpkin | HyperKZG |
| `Bn256EngineIPA` | BN256/Grumpkin | Pedersen+IPA |
| `PallasEngine` | Pallas/Vesta | Pedersen+IPA |
| `Secp256k1Engine` | Secp256k1/Secq256k1 | Pedersen+IPA |

### MSM (Multi-Scalar Multiplication) — Performance Critical Path

`src/provider/msm.rs` is the performance bottleneck. It implements a signed-decomposition + bit-width-partitioning strategy:

1. Signed scalar decomposition (halves effective scalar range)
2. Classify scalars by bit-width into 11 groups (binary, ≤8, ≤16, ≤32, ≤64, large)
3. Route each group to optimal algorithm: binary accumulation, single-window bucket sort (≤10 bits), multi-window Pippenger (≤64 bits), wNAF with XYZZ bucket coordinates (large)

The MSM is consumed via the `DlogGroupExt` trait's `vartime_multiscalar_mul` method, primarily called from Pedersen commitments (`provider/pedersen.rs`) and the IPA/HyperKZG evaluation arguments.

BN256 has a special case in `bn256_grumpkin.rs`: it uses `impl_traits_no_dlog_ext!` and manually implements `DlogGroupExt`, allowing conditional dispatch to blitzar GPU (`#[cfg(feature = "blitzar")]`).

### Commitment Engines

- **Pedersen** (`pedersen.rs`): Vector commitment using MSM. Used by all IPA-based engines.
- **HyperKZG** (`hyperkzg.rs`): Pairing-based polynomial commitment. BN256 only.
- **Mercury** (`mercury.rs`): Optimized KZG variant. BN256 only.

### Parallelism

Uses `rayon` throughout. MSM parallelizes over windows and chunks. Batch operations (`batch_vartime_multiscalar_mul`) parallelize over instances.

### x86_64 ASM Optimization

On x86_64, `halo2curves` compiles with the `asm` feature for ~50% faster field arithmetic. Non-x86_64 targets use pure Rust. Controlled via `[target.'cfg(...)'.dependencies]` in Cargo.toml.
