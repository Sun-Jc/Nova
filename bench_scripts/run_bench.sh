#!/usr/bin/env bash
set -euo pipefail

# ── Defaults ──────────────────────────────────────────────────────
SAMPLES=10
WARMUP=3          # seconds
MEASUREMENT=5     # seconds
SIZES="8,1024,32768,1048576,16777216"
BRANCHES="bench/before-opt,bench/after-opt,bench/before-opt-pasta,bench/after-opt-pasta"
EXTRA_FILTER=""

# ── Parse args ────────────────────────────────────────────────────
usage() {
  cat <<EOF
Usage: $0 [OPTIONS]

Options:
  --fast              Quick mode: 10 samples, 1s warmup, 1s measurement,
                      sizes 8+1024, BN254 branches only
  --samples N         Criterion sample size          (default: 10)
  --warmup  S         Criterion warm-up seconds      (default: 3)
  --measurement S     Criterion measurement seconds  (default: 5)
  --sizes   S1,S2,..  Comma-separated sizes          (default: all five)
  --branches B1,B2,.. Comma-separated branches       (default: all four)
  --filter  REGEX     Extra criterion filter regex    (default: none)
  -h, --help          Show this help
EOF
  exit 0
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --fast)
      SAMPLES=10; WARMUP=1; MEASUREMENT=1; SIZES="8,1024"
      BRANCHES="bench/before-opt,bench/after-opt"
      shift ;;
    --samples)      SAMPLES="$2";      shift 2 ;;
    --warmup)       WARMUP="$2";       shift 2 ;;
    --measurement)  MEASUREMENT="$2";  shift 2 ;;
    --sizes)        SIZES="$2";        shift 2 ;;
    --branches)     BRANCHES="$2";     shift 2 ;;
    --filter)       EXTRA_FILTER="$2"; shift 2 ;;
    -h|--help)      usage ;;
    *) echo "Unknown option: $1"; usage ;;
  esac
done

# ── Derived ───────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
RESULTS_DIR="$PROJECT_DIR/bench_results"
mkdir -p "$RESULTS_DIR"

ORIGINAL_REF=$(git -C "$PROJECT_DIR" rev-parse --abbrev-ref HEAD 2>/dev/null || \
               git -C "$PROJECT_DIR" rev-parse --short HEAD)

# Build size filter regex:  _(8|1024|32768)$
IFS=',' read -ra SIZE_ARR <<< "$SIZES"
SIZE_PATTERN=$(IFS='|'; echo "${SIZE_ARR[*]}")
FILTER="_(${SIZE_PATTERN})$"

# Combine with extra filter if provided
if [[ -n "$EXTRA_FILTER" ]]; then
  FILTER="${EXTRA_FILTER}.*${FILTER}"
fi

echo "============================================"
echo "Benchmark Configuration"
echo "  samples:      $SAMPLES"
echo "  warmup:       ${WARMUP}s"
echo "  measurement:  ${MEASUREMENT}s"
echo "  sizes:        $SIZES"
echo "  branches:     $BRANCHES"
echo "  filter:       $FILTER"
echo "  results:      $RESULTS_DIR"
echo "============================================"
echo ""

# ── Run ───────────────────────────────────────────────────────────
IFS=',' read -ra BRANCH_ARR <<< "$BRANCHES"

for branch in "${BRANCH_ARR[@]}"; do
  safe_name="${branch#bench/}"
  outfile="$RESULTS_DIR/${safe_name}.txt"

  echo "============================================"
  echo "Branch: $branch"
  echo "Output: $outfile"
  echo "============================================"

  git -C "$PROJECT_DIR" checkout "$branch"

  (cd "$PROJECT_DIR" && \
    cargo bench --features test-utils --bench commit -- \
      --sample-size "$SAMPLES" \
      --warm-up-time "$WARMUP" \
      --measurement-time "$MEASUREMENT" \
      "$FILTER" \
  2>&1) | tee "$outfile"

  echo ""
  echo "Done: $branch → $outfile"
  echo ""
done

# ── Restore ───────────────────────────────────────────────────────
git -C "$PROJECT_DIR" checkout "$ORIGINAL_REF" 2>/dev/null || true

echo "============================================"
echo "All benchmarks complete."
echo "Results:"
ls -lh "$RESULTS_DIR"/*.txt
echo "============================================"
