#!/usr/bin/env python3
"""
Run sparse commitment benchmark and generate report with diagrams.

Usage:
    python benches/run_sparse_bench.py [--skip-bench]

Options:
    --skip-bench    Skip running the benchmark, use existing results
"""

import subprocess
import re
import json
import sys
from pathlib import Path
from collections import defaultdict

try:
    import matplotlib.pyplot as plt
    import numpy as np
except ImportError:
    print("Error: matplotlib and numpy are required.")
    print("Install with: pip install matplotlib numpy")
    sys.exit(1)


def run_benchmark():
    """Run the cargo benchmark and stream output in realtime."""
    print("Running benchmark (this may take a while)...")
    print("=" * 60)
    
    output_lines = []
    
    process = subprocess.Popen(
        ["cargo", "bench", "--bench", "sparse-commit", "--", "--noplot"],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
        cwd=Path(__file__).parent.parent
    )
    
    # Stream output in realtime
    for line in process.stdout:
        print(line, end='', flush=True)
        output_lines.append(line)
    
    process.wait()
    
    return ''.join(output_lines)


def parse_criterion_output(output: str) -> dict:
    """
    Parse criterion benchmark output to extract timing data.
    
    Returns dict like:
    {
        (size, density): {
            'make_vec': time_ns,
            'sparse_commit': time_ns,
            'dense_commit': time_ns,
        }
    }
    """
    results = defaultdict(dict)
    
    # Pattern to match benchmark results
    # Example: "make_vec_64_d10        time:   [1.2345 µs 1.2456 µs 1.2567 µs]"
    # Or: "sparse_commit_1048576_d50"
    pattern = r'(\w+)_(\d+)_d(\d+)\s+time:\s+\[[\d.]+ \w+\s+([\d.]+) (µs|ms|ns|s)'
    
    for match in re.finditer(pattern, output):
        bench_type = match.group(1)  # make_vec, sparse_commit, dense_commit
        size = int(match.group(2))
        density = int(match.group(3))
        time_val = float(match.group(4))
        time_unit = match.group(5)
        
        # Convert to nanoseconds
        if time_unit == 'ns':
            time_ns = time_val
        elif time_unit == 'µs':
            time_ns = time_val * 1_000
        elif time_unit == 'ms':
            time_ns = time_val * 1_000_000
        elif time_unit == 's':
            time_ns = time_val * 1_000_000_000
        else:
            time_ns = time_val
        
        results[(size, density)][bench_type] = time_ns
    
    return dict(results)


def load_criterion_json() -> dict:
    """
    Load results from criterion's JSON output files as fallback.
    """
    results = defaultdict(dict)
    criterion_dir = Path(__file__).parent.parent / "target" / "criterion"
    
    if not criterion_dir.exists():
        return {}
    
    for bench_dir in criterion_dir.iterdir():
        if not bench_dir.is_dir():
            continue
        
        # Parse directory name like "make_vec_64_d10"
        match = re.match(r'(\w+)_(\d+)_d(\d+)', bench_dir.name)
        if not match:
            continue
        
        bench_type = match.group(1)
        size = int(match.group(2))
        density = int(match.group(3))
        
        # Look for estimates.json
        estimates_file = bench_dir / "new" / "estimates.json"
        if estimates_file.exists():
            with open(estimates_file) as f:
                data = json.load(f)
                # Get median time in nanoseconds
                time_ns = data.get("median", {}).get("point_estimate", 0)
                results[(size, density)][bench_type] = time_ns
    
    return dict(results)


def generate_report(results: dict, output_dir: Path):
    """Generate diagrams and report from benchmark results."""
    
    if not results:
        print("No benchmark results found!")
        return
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Organize data by density
    densities = sorted(set(d for _, d in results.keys()))
    sizes = sorted(set(s for s, _ in results.keys()))
    
    print(f"\nFound results for densities: {densities}")
    print(f"Found results for sizes: {sizes}")
    
    # Create one figure per density
    for density in densities:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
        fig.suptitle(f'Sparse vs Dense Commitment Benchmark (Density: {density}%)', fontsize=14)
        
        # Filter data for this density
        density_data = {s: results[(s, density)] for s, d in results.keys() if d == density}
        
        if not density_data:
            continue
        
        plot_sizes = sorted(density_data.keys())
        
        # Extract times
        make_vec_times = []
        sparse_times = []
        dense_times = []
        
        for size in plot_sizes:
            data = density_data[size]
            make_vec_times.append(data.get('make_vec', 0))
            sparse_times.append(data.get('sparse_commit', 0))
            dense_times.append(data.get('dense_commit', 0))
        
        make_vec_times = np.array(make_vec_times)
        sparse_times = np.array(sparse_times)
        dense_times = np.array(dense_times)
        
        # Convert to milliseconds for display
        make_vec_ms = make_vec_times / 1_000_000
        sparse_ms = sparse_times / 1_000_000
        dense_ms = dense_times / 1_000_000
        
        # Calculate metrics
        make_vec_percent = np.where(dense_times > 0, (make_vec_times / dense_times) * 100, 0)
        speedup = np.where(sparse_times > 0, dense_times / sparse_times, 0)
        
        # X-axis labels (log scale sizes)
        x_labels = [f'2^{int(np.log2(s))}' if s & (s-1) == 0 else str(s) for s in plot_sizes]
        x_pos = np.arange(len(plot_sizes))
        
        # === Plot 1: Percentage of make_vec in dense commit ===
        ax1.bar(x_pos, make_vec_percent, color='steelblue', alpha=0.8)
        ax1.set_xlabel('Vector Size')
        ax1.set_ylabel('make_vec as % of dense_commit')
        ax1.set_title('Vector Construction Overhead in Dense Commit')
        ax1.set_xticks(x_pos)
        ax1.set_xticklabels(x_labels, rotation=45, ha='right')
        ax1.axhline(y=50, color='red', linestyle='--', alpha=0.5, label='50%')
        ax1.legend()
        ax1.grid(axis='y', alpha=0.3)
        
        # Add percentage labels on bars
        for i, (pct, mv, dc) in enumerate(zip(make_vec_percent, make_vec_ms, dense_ms)):
            if pct > 0:
                ax1.annotate(f'{pct:.1f}%', xy=(i, pct), ha='center', va='bottom', fontsize=8)
        
        # === Plot 2: Speedup of sparse over dense ===
        colors = ['green' if s > 1 else 'red' for s in speedup]
        bars = ax2.bar(x_pos, speedup, color=colors, alpha=0.8)
        ax2.set_xlabel('Vector Size')
        ax2.set_ylabel('Speedup (dense_time / sparse_time)')
        ax2.set_title('Sparse Commit Speedup vs Dense Commit')
        ax2.set_xticks(x_pos)
        ax2.set_xticklabels(x_labels, rotation=45, ha='right')
        ax2.axhline(y=1, color='black', linestyle='-', alpha=0.5, label='Break-even')
        ax2.legend()
        ax2.grid(axis='y', alpha=0.3)
        
        # Add speedup labels on bars
        for i, s in enumerate(speedup):
            if s > 0:
                ax2.annotate(f'{s:.2f}x', xy=(i, s), ha='center', va='bottom', fontsize=8)
        
        plt.tight_layout()
        
        # Save figure
        fig_path = output_dir / f'sparse_bench_d{density}.png'
        plt.savefig(fig_path, dpi=150, bbox_inches='tight')
        print(f"Saved: {fig_path}")
        plt.close()
    
    # === Generate combined timing plot ===
    fig, axes = plt.subplots(1, len(densities), figsize=(6 * len(densities), 6))
    if len(densities) == 1:
        axes = [axes]
    
    fig.suptitle('Absolute Timing Comparison', fontsize=14)
    
    for ax, density in zip(axes, densities):
        density_data = {s: results[(s, density)] for s, d in results.keys() if d == density}
        plot_sizes = sorted(density_data.keys())
        
        make_vec_ms = np.array([density_data[s].get('make_vec', 0) for s in plot_sizes]) / 1_000_000
        sparse_ms = np.array([density_data[s].get('sparse_commit', 0) for s in plot_sizes]) / 1_000_000
        dense_ms = np.array([density_data[s].get('dense_commit', 0) for s in plot_sizes]) / 1_000_000
        
        x_labels = [f'2^{int(np.log2(s))}' if s & (s-1) == 0 else str(s) for s in plot_sizes]
        
        ax.semilogy(x_labels, make_vec_ms, 'o-', label='make_vec', color='blue', alpha=0.7)
        ax.semilogy(x_labels, sparse_ms, 's-', label='sparse_commit', color='green', alpha=0.7)
        ax.semilogy(x_labels, dense_ms, '^-', label='dense_commit', color='red', alpha=0.7)
        
        ax.set_xlabel('Vector Size')
        ax.set_ylabel('Time (ms, log scale)')
        ax.set_title(f'Density: {density}%')
        ax.legend()
        ax.grid(True, alpha=0.3)
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')
    
    plt.tight_layout()
    fig_path = output_dir / 'sparse_bench_timing.png'
    plt.savefig(fig_path, dpi=150, bbox_inches='tight')
    print(f"Saved: {fig_path}")
    plt.close()
    
    # === Generate text report ===
    report_path = output_dir / 'sparse_bench_report.txt'
    with open(report_path, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("SPARSE COMMITMENT BENCHMARK REPORT\n")
        f.write("=" * 80 + "\n\n")
        
        for density in densities:
            f.write(f"\n{'=' * 40}\n")
            f.write(f"DENSITY: {density}%\n")
            f.write(f"{'=' * 40}\n\n")
            
            density_data = {s: results[(s, density)] for s, d in results.keys() if d == density}
            plot_sizes = sorted(density_data.keys())
            
            f.write(f"{'Size':<12} {'make_vec':<12} {'sparse':<12} {'dense':<12} {'vec%':<10} {'speedup':<10}\n")
            f.write("-" * 70 + "\n")
            
            for size in plot_sizes:
                data = density_data[size]
                mv = data.get('make_vec', 0) / 1_000_000  # ms
                sp = data.get('sparse_commit', 0) / 1_000_000
                dn = data.get('dense_commit', 0) / 1_000_000
                
                vec_pct = (data.get('make_vec', 0) / data.get('dense_commit', 1)) * 100 if data.get('dense_commit', 0) > 0 else 0
                spdup = data.get('dense_commit', 0) / data.get('sparse_commit', 1) if data.get('sparse_commit', 0) > 0 else 0
                
                size_str = f"2^{int(np.log2(size))}" if size & (size-1) == 0 else str(size)
                f.write(f"{size_str:<12} {mv:<12.3f} {sp:<12.3f} {dn:<12.3f} {vec_pct:<10.1f} {spdup:<10.2f}x\n")
            
            f.write("\n")
    
    print(f"Saved: {report_path}")
    print("\n" + "=" * 60)
    print("REPORT SUMMARY")
    print("=" * 60)
    
    with open(report_path) as f:
        print(f.read())


def main():
    skip_bench = "--skip-bench" in sys.argv
    
    project_root = Path(__file__).parent.parent
    output_dir = project_root / "target" / "sparse_bench_report"
    
    if skip_bench:
        print("Skipping benchmark, loading existing results...")
        output = ""
    else:
        output = run_benchmark()
    
    # Try to parse from stdout first, then fall back to JSON files
    results = parse_criterion_output(output)
    
    if not results:
        print("Parsing criterion JSON files...")
        results = load_criterion_json()
    
    if not results:
        print("ERROR: No benchmark results found!")
        print("Make sure to run the benchmark first: cargo bench --bench sparse-commit")
        sys.exit(1)
    
    print(f"\nParsed {len(results)} benchmark results")
    generate_report(results, output_dir)
    
    print(f"\nReport generated in: {output_dir}")


if __name__ == "__main__":
    main()
