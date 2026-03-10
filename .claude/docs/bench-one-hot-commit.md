# One-Hot Batch Commitment Benchmark Results

**Branch:** `feat/batch-commit-one-hot`
**Curve:** BN254 (Bn256EngineKZG)
**Block size:** K=16
**Machine:** Apple Silicon (results may vary by platform)

## Commit-Major vs Block-Major

Baseline commit: `6a3f01f` (commit-major `batch_add_one_hot`)
Optimization commit: `fb24a49` (block-major `batch_add_one_hot_block`)

### 2^16 blocks (1M elements)

| Batch | commit-major | block-major | Speedup |
|-------|-------------|-------------|---------|
| x1    | 1.95ms      | 1.93ms      | 1.01×   |
| x16   | 40.1ms      | 25.5ms      | **1.57×** |
| x64   | 157ms       | 94.0ms      | **1.67×** |
| x256  | 626ms       | 388ms       | **1.61×** |

### 2^18 blocks (4M elements)

| Batch | commit-major | block-major | Speedup |
|-------|-------------|-------------|---------|
| x1    | 8.64ms      | 8.62ms      | 1.00×   |
| x16   | 154ms       | 98.6ms      | **1.56×** |
| x64   | 616ms       | 400ms       | **1.54×** |
| x256  | 2.48s       | 1.49s       | **1.66×** |

### 2^20 blocks (16M elements)

| Batch | commit-major | block-major | Speedup |
|-------|-------------|-------------|---------|
| x1    | 37.1ms      | 35.8ms      | 1.04×   |
| x16   | 568ms       | 366ms       | **1.55×** |
| x64   | 2.26s       | 1.38s       | **1.63×** |
| x256  | 9.06s       | 5.49s       | **1.65×** |

## Summary

- Block-major consistently **1.5–1.67× faster** for batch ≥ 16
- x1 performance identical (no regression)
- Speedup slightly increases with larger batch sizes (better cache reuse)
- Block-major works by transposing iteration: generators loaded once per block, distributed to all commits' accumulators
