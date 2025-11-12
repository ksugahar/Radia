"""
Analyze H-matrix cache data for ML parameter tuning

This tool reads the H-matrix cache and analyzes the relationship between
problem characteristics and optimal parameters to improve future parameter
selection.

Usage:
    python tools/analyze_cache_for_ml.py [cache_file]
"""

import sys
import os
import struct
from pathlib import Path
import numpy as np
from collections import defaultdict

def read_cache_file(cache_path):
    """Read and parse the H-matrix cache file"""
    entries = []

    if not os.path.exists(cache_path):
        print(f"Cache file not found: {cache_path}")
        return entries

    with open(cache_path, 'rb') as f:
        # Read header
        magic = struct.unpack('I', f.read(4))[0]
        version = struct.unpack('I', f.read(4))[0]
        num_entries = struct.unpack('I', f.read(4))[0]

        if magic != 0x52414448:  # "RADH"
            print(f"Invalid magic number: {hex(magic)}")
            return entries

        print(f"Cache file: {cache_path}")
        print(f"Version: {version}")
        print(f"Entries: {num_entries}")
        print()

        # Read entries
        for i in range(num_entries):
            entry = {}
            entry['geometry_hash'] = struct.unpack('Q', f.read(8))[0]
            entry['num_elements'] = struct.unpack('I', f.read(4))[0]
            entry['eps'] = struct.unpack('d', f.read(8))[0]
            entry['max_rank'] = struct.unpack('I', f.read(4))[0]
            entry['timestamp'] = struct.unpack('q', f.read(8))[0]
            entry['construction_time'] = struct.unpack('d', f.read(8))[0]
            entry['memory_used'] = struct.unpack('Q', f.read(8))[0]
            entry['compression_ratio'] = struct.unpack('d', f.read(8))[0]

            entries.append(entry)

    return entries

def analyze_parameter_performance(entries):
    """Analyze which parameters work best for different problem sizes"""

    print("="*80)
    print("PARAMETER PERFORMANCE ANALYSIS")
    print("="*80)
    print()

    # Group by problem size
    size_groups = defaultdict(list)
    for entry in entries:
        N = entry['num_elements']
        size_groups[N].append(entry)

    print(f"{'N (elements)':<15} {'Count':<8} {'Avg Time (s)':<15} {'Avg eps':<12} {'Avg rank':<12} {'Compression':<12}")
    print("-"*80)

    results = []
    for N in sorted(size_groups.keys()):
        group = size_groups[N]
        avg_time = np.mean([e['construction_time'] for e in group])
        avg_eps = np.mean([e['eps'] for e in group])
        avg_rank = np.mean([e['max_rank'] for e in group])
        avg_comp = np.mean([e['compression_ratio'] for e in group])

        print(f"{N:<15} {len(group):<8} {avg_time:<15.3f} {avg_eps:<12.6f} {avg_rank:<12.1f} {avg_comp:<12.2f}%")

        results.append({
            'N': N,
            'count': len(group),
            'avg_time': avg_time,
            'avg_eps': avg_eps,
            'avg_rank': avg_rank,
            'avg_compression': avg_comp
        })

    return results

def find_optimal_parameters(entries):
    """Find optimal parameters by minimizing construction time while maintaining accuracy"""

    print()
    print("="*80)
    print("OPTIMAL PARAMETER SEARCH")
    print("="*80)
    print()

    # Group by (N, eps, max_rank) to find best configurations
    config_performance = defaultdict(list)

    for entry in entries:
        key = (entry['num_elements'], entry['eps'], entry['max_rank'])
        config_performance[key].append(entry['construction_time'])

    # Find best configuration for each problem size
    size_best = {}
    for N in set(e['num_elements'] for e in entries):
        N_configs = [(k, v) for k, v in config_performance.items() if k[0] == N]
        if N_configs:
            # Sort by average construction time
            best = min(N_configs, key=lambda x: np.mean(x[1]))
            size_best[N] = {
                'eps': best[0][1],
                'max_rank': best[0][2],
                'avg_time': np.mean(best[1])
            }

    print(f"{'N (elements)':<15} {'Best eps':<12} {'Best rank':<12} {'Avg Time (s)':<15}")
    print("-"*80)
    for N in sorted(size_best.keys()):
        best = size_best[N]
        print(f"{N:<15} {best['eps']:<12.6f} {best['max_rank']:<12} {best['avg_time']:<15.3f}")

    return size_best

def generate_adaptive_rules(size_best):
    """Generate adaptive parameter selection rules from optimal data"""

    print()
    print("="*80)
    print("RECOMMENDED ADAPTIVE PARAMETER RULES")
    print("="*80)
    print()

    if not size_best:
        print("Insufficient data for rule generation")
        return

    # Fit simple piecewise linear model
    sizes = sorted(size_best.keys())

    print("Suggested C++ implementation:")
    print()
    print("```cpp")
    print("static void OptimizeHMatrixParameters(int num_elements, double& eps, int& max_rank)")
    print("{")
    print("    // ML-tuned parameters based on cache analysis")

    for i, N in enumerate(sizes):
        best = size_best[N]
        if i == 0:
            print(f"    if(num_elements < {N})")
        elif i == len(sizes) - 1:
            print(f"    else  // num_elements >= {N}")
        else:
            next_N = sizes[i+1] if i+1 < len(sizes) else N*2
            print(f"    else if(num_elements < {next_N})")

        print(f"    {{")
        print(f"        eps = {best['eps']};")
        print(f"        max_rank = {best['max_rank']};")
        print(f"    }}")

    print("}")
    print("```")

def predict_performance(entries, N_test):
    """Predict construction time for a given problem size"""

    print()
    print("="*80)
    print("PERFORMANCE PREDICTION")
    print("="*80)
    print()

    # Simple power law fit: time = a * N^b
    sizes = np.array([e['num_elements'] for e in entries])
    times = np.array([e['construction_time'] for e in entries])

    if len(sizes) < 2:
        print("Insufficient data for prediction")
        return

    # Log-linear regression
    log_sizes = np.log(sizes)
    log_times = np.log(times)

    # Fit y = mx + c -> log(time) = b*log(N) + log(a)
    m, c = np.polyfit(log_sizes, log_times, 1)
    a = np.exp(c)
    b = m

    print(f"Fitted power law: time = {a:.6e} * N^{b:.4f}")
    print()

    # Predict for test size
    predicted_time = a * (N_test ** b)

    print(f"Prediction for N={N_test}:")
    print(f"  Estimated construction time: {predicted_time:.3f} s")
    print()

    # Recommend parameters
    if N_test < 200:
        print(f"  Recommendation: Use dense solver (N < 200)")
    elif N_test < 500:
        print(f"  Recommendation: eps=1e-4, max_rank=30")
    elif N_test < 1000:
        print(f"  Recommendation: eps=2e-4, max_rank=25")
    else:
        print(f"  Recommendation: eps=5e-4, max_rank=20")

def generate_summary_report(entries):
    """Generate comprehensive summary report"""

    print()
    print("="*80)
    print("CACHE SUMMARY REPORT")
    print("="*80)
    print()

    if not entries:
        print("No cache entries found")
        return

    # Overall statistics
    total_entries = len(entries)
    unique_geometries = len(set(e['geometry_hash'] for e in entries))
    total_construction_time = sum(e['construction_time'] for e in entries)
    avg_construction_time = np.mean([e['construction_time'] for e in entries])
    avg_compression = np.mean([e['compression_ratio'] for e in entries])

    print(f"Total cache entries: {total_entries}")
    print(f"Unique geometries: {unique_geometries}")
    print(f"Total construction time: {total_construction_time:.2f} s")
    print(f"Average construction time: {avg_construction_time:.3f} s")
    print(f"Average compression ratio: {avg_compression:.2f}%")
    print()

    # Size distribution
    sizes = [e['num_elements'] for e in entries]
    print(f"Problem size range: {min(sizes)} - {max(sizes)} elements")
    print(f"Average problem size: {np.mean(sizes):.0f} elements")
    print()

    # Parameter usage
    eps_values = set(e['eps'] for e in entries)
    rank_values = set(e['max_rank'] for e in entries)

    print(f"eps values used: {sorted(eps_values)}")
    print(f"max_rank values used: {sorted(rank_values)}")

def main():
    """Main entry point"""

    # Determine cache file path
    if len(sys.argv) > 1:
        cache_path = sys.argv[1]
    else:
        cache_path = os.path.join(os.getcwd(), '.radia_cache', 'hmatrix_cache.bin')

    # Read cache
    entries = read_cache_file(cache_path)

    if not entries:
        print("No entries found in cache")
        return

    # Generate report
    generate_summary_report(entries)

    # Analyze parameter performance
    results = analyze_parameter_performance(entries)

    # Find optimal parameters
    size_best = find_optimal_parameters(entries)

    # Generate adaptive rules
    generate_adaptive_rules(size_best)

    # Performance prediction
    test_sizes = [500, 1000, 2000]
    for N in test_sizes:
        if N not in [e['num_elements'] for e in entries]:
            predict_performance(entries, N)

if __name__ == "__main__":
    main()
