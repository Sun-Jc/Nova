#!/usr/bin/env python3
"""Parse criterion benchmark output and generate throughput / speedup plots."""

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

FILE_MAP = {
    "before-opt": ("BN254", "old"),
    "after-opt": ("BN254", "new"),
    "before-opt-pasta": ("Pallas", "old"),
    "after-opt-pasta": ("Pallas", "new"),
}

BENCH_RE = re.compile(
    r"^(?P<prefix>halo2curves|nova_generic|nova_specialized|nova_batch_add)"
    r"_commit_(?P<bw>[a-z0-9]+)_(?P<size>\d+)(?=\s|$)"
)

PREFIX_TO_METHOD = {
    "halo2curves": "Halo2",
    "nova_generic": "CE-General",
    "nova_specialized": "CE-Small",
    "nova_batch_add": "BA",
}

TIME_RE = re.compile(
    r"time:\s+\[\s*"
    r"[\d.]+ \w+\s+"
    r"([\d.]+) (\w+)\s+"
    r"[\d.]+ \w+\s*\]"
)

UNIT_TO_SEC = {
    "ps": 1e-12,
    "ns": 1e-9,
    "µs": 1e-6,
    "us": 1e-6,
    "ms": 1e-3,
    "s": 1.0,
}

METHOD_STYLES = {
    "Halo2": {"marker": "s", "linestyle": "-"},
    "CE-General": {"marker": "o", "linestyle": "--"},
    "CE-Small": {"marker": "^", "linestyle": ":"},
    "BA": {"marker": "D", "linestyle": "-."},
}

VERSION_COLOR = {"old": "tab:red", "new": "tab:blue"}

# BitWidth display order
BW_ORDER = ["u1", "u10", "u16", "u32", "u64", "random"]

# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------


def parse_criterion_output(filepath: Path) -> dict:
    """Return {bench_name: time_in_seconds}."""
    results = {}
    text = filepath.read_text()
    lines = text.splitlines()

    current_bench = None
    for line in lines:
        stripped = line.strip()
        # Try to match a benchmark name (a line that matches our pattern)
        m = BENCH_RE.match(stripped)
        if m:
            bench_name = m.group(0)
            # Check if time is on the same line (e.g. "halo2curves_commit_u1_8 time:   [...]")
            tm = TIME_RE.search(stripped)
            if tm:
                value = float(tm.group(1))
                unit = tm.group(2)
                if unit in UNIT_TO_SEC:
                    results[bench_name] = value * UNIT_TO_SEC[unit]
                else:
                    print(f"  WARNING: unknown unit '{unit}' for {bench_name}")
                current_bench = None
            else:
                current_bench = bench_name
            continue
        # Try to match time line
        if current_bench:
            tm = TIME_RE.search(line)
            if tm:
                value = float(tm.group(1))
                unit = tm.group(2)
                if unit in UNIT_TO_SEC:
                    results[current_bench] = value * UNIT_TO_SEC[unit]
                else:
                    print(f"  WARNING: unknown unit '{unit}' for {current_bench}")
                current_bench = None
    return results


def load_all_results(results_dir: Path):
    """Return nested dict: data[curve][bitwidth][(method, size, version)] = time_sec."""
    data = defaultdict(lambda: defaultdict(dict))

    for stem, (curve, version) in FILE_MAP.items():
        filepath = results_dir / f"{stem}.txt"
        if not filepath.exists():
            print(f"  SKIP: {filepath} not found")
            continue
        print(f"  Parsing {filepath.name} -> {curve}/{version}")
        results = parse_criterion_output(filepath)
        print(f"    Found {len(results)} benchmarks")

        for bench_name, time_sec in results.items():
            m = BENCH_RE.match(bench_name)
            if not m:
                continue
            method = PREFIX_TO_METHOD[m.group("prefix")]
            bw = m.group("bw")
            size = int(m.group("size"))
            data[curve][bw][(method, size, version)] = time_sec

    return data


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------


def plot_throughput(curve, bw, bw_data, ax):
    """Figure 1: throughput (size / time) vs size."""
    ax.set_title("Throughput")
    ax.set_xlabel("Size (number of scalars)")
    ax.set_ylabel("Throughput (elements / sec)")
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")

    combos = defaultdict(lambda: defaultdict(dict))
    for (method, size, version), t in bw_data.items():
        combos[method][version][size] = size / t

    for method in sorted(combos, key=lambda m: list(METHOD_STYLES).index(m)):
        style = METHOD_STYLES[method]
        for version in ("old", "new"):
            if version not in combos[method]:
                continue
            pts = sorted(combos[method][version].items())
            ax.plot(
                [p[0] for p in pts],
                [p[1] for p in pts],
                linestyle=style["linestyle"],
                color=VERSION_COLOR[version],
                marker=style["marker"],
                markersize=5,
                label=f"{method} ({version})",
            )

    ax.legend(fontsize=7, loc="best")
    ax.grid(True, which="both", alpha=0.3)


def plot_speedup(curve, bw, bw_data, ax):
    """Figure 2: speedup (old_time / new_time) vs size."""
    ax.set_title("Speedup (old / new)")
    ax.set_xlabel("Size (number of scalars)")
    ax.set_ylabel("Speedup")
    ax.set_xscale("log", base=2)
    ax.axhline(y=1.0, color="black", linewidth=0.8, linestyle=":")

    pairs = defaultdict(dict)
    for (method, size, version), t in bw_data.items():
        pairs[(method, size)][version] = t

    combos = defaultdict(dict)
    for (method, size), vt in pairs.items():
        if "old" in vt and "new" in vt:
            combos[method][size] = vt["old"] / vt["new"]

    for method in sorted(combos, key=lambda m: list(METHOD_STYLES).index(m)):
        style = METHOD_STYLES[method]
        pts = sorted(combos[method].items())
        if not pts:
            continue
        ax.plot(
            [p[0] for p in pts],
            [p[1] for p in pts],
            linestyle=style["linestyle"],
            color="tab:purple",
            marker=style["marker"],
            markersize=6,
            label=method,
        )

    ax.legend(fontsize=8, loc="best")
    ax.grid(True, which="both", alpha=0.3)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(description="Plot MSM benchmark results")
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "bench_results",
        help="Directory with criterion output .txt files",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory for plots (default: <input-dir>/plots)",
    )
    parser.add_argument(
        "--format",
        default="png",
        choices=["png", "svg", "pdf"],
        help="Output image format",
    )
    parser.add_argument("--dpi", type=int, default=150)
    parser.add_argument(
        "--curves",
        default=None,
        help="Comma-separated curves to plot (default: all found)",
    )
    parser.add_argument(
        "--bitwidths",
        default=None,
        help="Comma-separated bitwidths to plot (default: all found)",
    )
    args = parser.parse_args()

    output_dir = args.output_dir or (args.input_dir / "plots")
    output_dir.mkdir(parents=True, exist_ok=True)

    print("Loading results...")
    data = load_all_results(args.input_dir)
    if not data:
        print("No data found. Run run_bench.sh first.")
        sys.exit(1)

    curve_filter = set(args.curves.split(",")) if args.curves else None
    bw_filter = set(args.bitwidths.split(",")) if args.bitwidths else None

    for curve in sorted(data.keys()):
        if curve_filter and curve not in curve_filter:
            continue
        # Sort bitwidths by predefined order
        bws = sorted(
            data[curve].keys(),
            key=lambda b: BW_ORDER.index(b) if b in BW_ORDER else 999,
        )
        for bw in bws:
            if bw_filter and bw not in bw_filter:
                continue
            bw_data = data[curve][bw]

            fig, axes = plt.subplots(1, 2, figsize=(14, 5))
            fig.suptitle(f"{curve} — BitWidth: {bw}", fontsize=14)

            plot_throughput(curve, bw, bw_data, axes[0])
            plot_speedup(curve, bw, bw_data, axes[1])

            fig.tight_layout()
            outpath = output_dir / f"{curve}_{bw}.{args.format}"
            fig.savefig(outpath, dpi=args.dpi)
            plt.close(fig)
            print(f"  Saved: {outpath}")

    print(f"\nAll plots saved to {output_dir}")


if __name__ == "__main__":
    main()
