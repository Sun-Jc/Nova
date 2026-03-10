# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Nova is a high-speed recursive SNARK library implementing Incrementally Verifiable Computation (IVC) via folding schemes. It is developed by Microsoft Research. The prover's work to update proofs is constant regardless of steps executed, and the verifier circuit is ~10,000 multiplication gates.

**Crate name:** `nova-snark` (v0.59.0)
**Rust edition:** 2021, **MSRV:** 1.79.0
**Base commit:** `3ca2f0636f8ea2be0c1b46e005b1f7fd8042a258` — all subsequent changes are based on this commit.
**Reference branch:** `opt/one-hot-precomp-transpose` — reference for patterns and approaches, not the working branch.

## Optimization Goal

**Scope:** One-Hot Commitment Batch — optimize batch commitment of one-hot structured vectors.

Three methods under exploration:
1. **XYZZ coordinates** — use extended Jacobian `(X, Y, ZZ, ZZZ)` for cheaper bucket accumulation (7M+2S vs ~11M+5S)
2. **Batch block** — transpose iteration order from per-vector to per-block, improving cache locality across batch
3. **Precompute** — precompute pairwise sums of generators for adjacent blocks, halving point additions

**Benchmark baseline:** commit `6a3f01f` on branch `feat/batch-commit-one-hot`. ~38ms/vector (K=16, 2^20 blocks, BN254). Run: `cargo bench --bench one-hot-commit --features test-utils`

## Build & Development Commands

```bash
# Build (release recommended — crypto code is extremely slow in debug)
cargo build --release
cargo build --release --examples --benches

# Run all tests (always use --release)
cargo test --release

# Run a specific test by name
cargo test --release <test_name>
# Filter to tests containing "ivc" in a specific module
cargo test --release ivc -- --nocapture

# Run experimental NeutronNova tests
cargo test --release --features experimental neutron

# Lint and format
cargo fmt --all -- --check
cargo clippy --all-targets -- -D warnings

# Generate docs (CI checks this)
cargo doc --no-deps --document-private-items

# Run examples
cargo run --release --example minroot
cargo run --release --example hashchain

# Run benchmarks (Criterion)
cargo bench --bench recursive-snark
cargo bench --bench compressed-snark
cargo bench --bench sha256
cargo bench --bench ppsnark
cargo bench --bench commit

# WASM build (no default features, no ASM)
cargo build --no-default-features --target wasm32-unknown-unknown

# Verify all CI-checked feature combos locally
cargo check --features evm
cargo check --no-default-features
cargo check --features "io,evm,experimental,test-utils"
```

**CI runs 9 jobs**: build, build-wasm, test, test-experimental, fmt, clippy, docs, spelling (via `typos` crate), and check-features. All must pass before merge.

## Cargo Features

| Feature | Description |
|---------|-------------|
| `io` (default) | Powers of Tau file I/O for HyperKZG/Mercury setup |
| `evm` | EVM-friendly serialization support |
| `experimental` | NeutronNova (unstable, may break between releases) |
| `test-utils` | Insecure random tau generation — **never use in production** |
| `flamegraph` | Criterion flamegraph profiling |
| `blitzar` | GPU-accelerated MSM via Blitzar (implicit feature from optional dep) |

CI validates all feature combinations (see Build & Development Commands above).

## Architecture

### Curve Cycles (enabling recursive proofs)
Nova operates over a **cycle of elliptic curves** where curve A's base field = curve B's scalar field and vice versa. The library parameterizes proofs over two `Engine` types (`E1`, `E2`) with this relationship enforced at the type level:
```
E1: Engine<Base = <E2 as Engine>::Scalar>
E2: Engine<Base = <E1 as Engine>::Scalar>
```
Three curve cycles are supported: **Pallas/Vesta**, **BN254/Grumpkin** (pairing-capable), **secp256k1/secq256k1**.

### Core Module Map

- **`src/nova/`** — Main IVC scheme: `PublicParams` (setup), `RecursiveSNARK` (incremental proving), `CompressedSNARK` (Spartan compression). The `nifs.rs` submodule implements Non-Interactive Folding.
- **`src/spartan/`** — Spartan SNARK for compressing IVC proofs. Two variants: `snark.rs` (non-preprocessing) and `ppsnark.rs` (preprocessing/MicroSpartan for on-chain verification). Uses sumcheck protocol (`sumcheck.rs`) and multilinear polynomial operations (`polys/`).
- **`src/traits/`** — Core trait hierarchy. `Engine` bundles curve group, fields, RO (random oracle), transcript, and commitment engine. `StepCircuit` is what users implement to define their computation. `CommitmentEngineTrait` and `RelaxedR1CSSNARKTrait` abstract commitment schemes and SNARK backends. `EvaluationEngineTrait` (in `evaluation.rs`) defines polynomial evaluation arguments — this is how commitment schemes (Pedersen/IPA, HyperKZG, Mercury) prove polynomial evaluations.
- **`src/provider/`** — Concrete implementations of traits:
  - Curves: `pasta.rs`, `bn256_grumpkin.rs`, `secp_secq.rs`
  - Commitments: `pedersen.rs` (+ IPA in `ipa_pc.rs`), `hyperkzg.rs` (pairing-based), `mercury.rs` (optimized KZG)
  - Hash/RO: `poseidon.rs` (circuit-friendly), `keccak.rs` (transcript)
  - `msm.rs` — Multi-scalar multiplication
  - `ptau.rs` — Powers of Tau file loading
  - `blitzar.rs` — Optional GPU acceleration
- **`src/frontend/`** — Bellman-style circuit API. Users build circuits with `ConstraintSystem`, `AllocatedNum`, `Boolean`. Includes gadgets for SHA256, Poseidon sponge, and test utilities (`test_shape_cs`, `util_cs/`).
- **`src/gadgets/`** — Low-level in-circuit primitives: `ecc.rs` (elliptic curve operations), `nonnative/` (cross-field arithmetic — critical for recursive verification).
- **`src/r1cs/`** — R1CS representation: `R1CSShape` (circuit structure), `R1CSInstance`/`R1CSWitness` (concrete assignments), `RelaxedR1CSInstance`/`RelaxedR1CSWitness` (folded form). Sparse matrix representation in `sparse.rs`.
- **`src/neutron/`** — Experimental NeutronNova IVC variant (feature-gated behind `experimental`).

### Implementing a StepCircuit

Users define their IVC computation by implementing `StepCircuit<F>` (in `src/traits/circuit.rs`). The two required methods:
- `arity()` — number of input/output elements per step (input and output vectors must be the same size)
- `synthesize(cs, z)` — given a constraint system and input `z_i`, produce output `z_{i+1}` using bellman-style gadgets (`AllocatedNum`, `Boolean`, etc.)

**Important**: `StepCircuit::synthesize` must NOT call `inputize`/`alloc_io` — pass outputs via the return value only, or you'll get `InvalidStepCircuitIO`.

See `examples/minroot.rs` for a complete working example. `TrivialCircuit` (identity) and `NonTrivialCircuit` (repeated squaring) in `src/traits/circuit.rs` are useful for testing.

### Key Data Flow
```
User StepCircuit → PublicParams::setup() → RecursiveSNARK::new() / .prove_step()
  → [fold with random instance for ZK] → CompressedSNARK::prove() → .verify()
```

Nova uses a **primary** circuit (the user's `StepCircuit`) running on E1 and a **secondary** (augmented) circuit on E2 that verifies the primary's folding proof. The secondary circuit is always a `TrivialCircuit`. Both circuits are augmented internally with verification gadgets — this is where the nonnative arithmetic in `gadgets/nonnative/` is critical, since each circuit must verify operations over the other curve's field.

### Crate-Level Type Aliases

`src/lib.rs` defines convenience aliases used throughout the codebase — when you see bare `CommitmentKey<E>` or `Commitment<E>`, they resolve through the `Engine`'s associated commitment engine:
```rust
type CommitmentKey<E> = <<E as Engine>::CE as CommitmentEngineTrait<E>>::CommitmentKey;
type Commitment<E>    = <<E as Engine>::CE as CommitmentEngineTrait<E>>::Commitment;
type DerandKey<E>     = <<E as Engine>::CE as CommitmentEngineTrait<E>>::DerandKey;
type CE<E>            = <E as Engine>::CE;
```

### Commitment Scheme Selection
- **Pedersen + IPA**: Works on all curves, no trusted setup needed
- **HyperKZG**: Requires pairing curves (BN254) + Powers of Tau setup, faster verification
- **Mercury**: Optimized KZG variant, even faster verification than HyperKZG

## Code Conventions

- **English only in code** — All source code, comments, doc comments, error messages, string literals, and commit messages must be in English. No Chinese, Japanese, Korean, or other non-ASCII scripts in code. (Conversations with the developer may use Chinese or English.)
- **Maximize reuse** — Reuse existing functions, types, and patterns before introducing new ones. Prefer extending or composing existing abstractions.
- **Minimize diff** — Keep changes as small and focused as possible. Do not refactor unrelated code alongside feature work.
- **No anti-patterns** — Avoid code smells such as unnecessary clones, redundant allocations, copy-paste duplication, or overly complex generics where simpler solutions exist.
- **2-space indentation**, Unix newlines, `use_try_shorthand = true` (see `rustfmt.toml`)
- **`#![deny(unsafe_code)]`** — No unsafe Rust in this crate (note: `deny` not `forbid`, but the convention is to never use `unsafe`)
- **`#![deny(missing_docs)]`** — All public items must have doc comments
- **`#![deny(warnings, unused, future_incompatible, nonstandard_style, rust_2018_idioms)]`** — strict lint baseline
- **`#![allow(non_snake_case)]`** — Mathematical notation (uppercase variables like `A`, `B`, `W`, `U`) is intentional and expected
- **No `print!`/`println!` outside tests** — enforced by clippy (`clippy::print_stdout`, `clippy::print_stderr`) via `cfg_attr(not(test), warn(...))`
- Clippy is configured permissively for crypto code: type complexity threshold 9999, max function args 20 (`.clippy.toml`)
- **`NovaError` is `#[non_exhaustive]`** — always use wildcard arms when matching against it
- **`#[serde(bound = "")]`** is used pervasively on generic structs — this overrides Serde's default trait bound inference, which would otherwise require bounds that don't hold for the complex generic types in this codebase
- **EVM serialization**: When the `evm` feature is enabled, `CustomSerdeTrait` (in `traits/evm_serde.rs`) overrides default serde for curve/field types via `serde_with::serde_as`. New types exposed for EVM must implement `CustomSerdeTrait`.
- **Code simplicity over minor performance gains** — this is a core project value. Performance PRs require reproducible benchmarks showing substantial speedups across diverse circuits.

## x86_64 vs Other Architectures

On x86_64, `halo2curves` uses ASM for ~50% performance boost (enabled automatically). On other architectures or WASM, pure Rust fallback is used. If build fails on non-x86, try `--no-default-features`.

## Docs

Benchmark results and optimization notes are stored in `.claude/docs/`. See `bench-one-hot-commit.md` for one-hot batch commitment performance data.
