"""
ML-based parameter optimization for H-matrix

This script uses machine learning to predict optimal H-matrix parameters
(eps, max_rank) based on problem size and characteristics.

Usage:
    python tools/optimize_parameters_ml.py [--train] [--predict N]
"""

import sys
import os
import struct
import argparse
import numpy as np
from pathlib import Path

def read_cache_file(cache_path):
    """Read and parse the H-matrix cache file"""
    entries = []

    if not os.path.exists(cache_path):
        return entries

    with open(cache_path, 'rb') as f:
        # Read header
        magic = struct.unpack('I', f.read(4))[0]
        version = struct.unpack('I', f.read(4))[0]
        num_entries = struct.unpack('I', f.read(4))[0]

        if magic != 0x52414448:  # "RADH"
            return entries

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

def train_ml_model(entries):
    """Train ML model to predict optimal parameters"""

    print("="*80)
    print("ML MODEL TRAINING")
    print("="*80)
    print()

    if len(entries) < 5:
        print("Insufficient training data (need at least 5 entries)")
        return None

    # Extract features and targets
    # For each unique geometry (N), find the configuration with minimum construction time
    from collections import defaultdict

    geom_configs = defaultdict(list)
    for entry in entries:
        N = entry['num_elements']
        geom_configs[N].append(entry)

    # Find best configuration for each problem size
    training_data = []
    for N, configs in geom_configs.items():
        # Find config with minimum construction time
        best = min(configs, key=lambda x: x['construction_time'])
        training_data.append({
            'N': N,
            'best_eps': best['eps'],
            'best_max_rank': best['max_rank'],
            'best_time': best['construction_time']
        })

    # Sort by problem size
    training_data.sort(key=lambda x: x['N'])

    print("Training data:")
    print(f"{'N':<10} {'Best eps':<12} {'Best rank':<12} {'Time (s)':<12}")
    print("-"*50)
    for d in training_data:
        print(f"{d['N']:<10} {d['best_eps']:<12.6f} {d['best_max_rank']:<12} {d['best_time']:<12.3f}")

    print()

    # Fit piecewise linear model for eps and max_rank
    sizes = np.array([d['N'] for d in training_data])
    eps_values = np.array([d['best_eps'] for d in training_data])
    rank_values = np.array([d['best_max_rank'] for d in training_data])

    # Simple regression: eps = a * log(N) + b
    log_sizes = np.log(sizes)
    eps_coeffs = np.polyfit(log_sizes, eps_values, 1)  # Linear fit in log space
    rank_coeffs = np.polyfit(log_sizes, rank_values, 1)

    print("Regression models:")
    print(f"  eps = {eps_coeffs[0]:.6f} * log(N) + {eps_coeffs[1]:.6f}")
    print(f"  max_rank = {rank_coeffs[0]:.2f} * log(N) + {rank_coeffs[1]:.2f}")
    print()

    model = {
        'training_data': training_data,
        'eps_coeffs': eps_coeffs,
        'rank_coeffs': rank_coeffs
    }

    return model

def predict_parameters(model, N):
    """Predict optimal parameters for a given problem size"""

    if model is None:
        # Fallback to rule-based (current implementation)
        if N < 200:
            return 1e-4, 30
        elif N < 500:
            return 1e-4, 30
        elif N < 1000:
            return 2e-4, 25
        else:
            return 5e-4, 20

    # Use regression model
    log_N = np.log(N)
    eps = model['eps_coeffs'][0] * log_N + model['eps_coeffs'][1]
    max_rank = model['rank_coeffs'][0] * log_N + model['rank_coeffs'][1]

    # Clamp to reasonable ranges
    eps = np.clip(eps, 1e-5, 1e-3)
    max_rank = int(np.clip(max_rank, 10, 50))

    return eps, max_rank

def generate_cpp_code(model):
    """Generate C++ code for parameter selection"""

    print("="*80)
    print("C++ CODE GENERATION (ML-BASED)")
    print("="*80)
    print()

    if model is None:
        print("No model available - using rule-based approach")
        return

    training_data = model['training_data']

    print("Suggested C++ implementation:")
    print()
    print("```cpp")
    print("static void OptimizeHMatrixParameters(int num_elements, double& eps, int& max_rank)")
    print("{")
    print("    // ML-optimized parameters (trained from cache data)")
    print("    // Trade-off: construction time vs accuracy")
    print()

    # Generate piecewise rules based on training data
    for i, data in enumerate(training_data):
        N = data['N']
        best_eps = data['best_eps']
        best_rank = data['best_max_rank']

        if i == 0:
            print(f"    if(num_elements < {N})")
        elif i == len(training_data) - 1:
            print(f"    else  // num_elements >= {N}")
        else:
            next_N = training_data[i+1]['N']
            mid_N = (N + next_N) // 2
            print(f"    else if(num_elements < {mid_N})")

        print(f"    {{")
        print(f"        eps = {best_eps};  // Optimized for N~{N}")
        print(f"        max_rank = {best_rank};")
        print(f"    }}")

    print("}")
    print("```")
    print()

    # Also print regression-based version
    print("Alternative: Continuous regression-based prediction:")
    print()
    print("```cpp")
    print("static void OptimizeHMatrixParameters(int num_elements, double& eps, int& max_rank)")
    print("{")
    print("    // Regression model: eps = a*log(N) + b, rank = c*log(N) + d")
    print(f"    double log_N = std::log(num_elements);")
    print(f"    eps = {model['eps_coeffs'][0]:.6f} * log_N + {model['eps_coeffs'][1]:.6f};")
    print(f"    max_rank = static_cast<int>({model['rank_coeffs'][0]:.2f} * log_N + {model['rank_coeffs'][1]:.2f});")
    print()
    print("    // Clamp to reasonable ranges")
    print("    eps = std::max(1e-5, std::min(1e-3, eps));")
    print("    max_rank = std::max(10, std::min(50, max_rank));")
    print("}")
    print("```")

def main():
    """Main entry point"""

    parser = argparse.ArgumentParser(description="ML-based H-matrix parameter optimization")
    parser.add_argument('--train', action='store_true', help='Train ML model from cache')
    parser.add_argument('--predict', type=int, metavar='N', help='Predict parameters for problem size N')
    parser.add_argument('--cache', default='.radia_cache/hmatrix_cache.bin', help='Cache file path')

    args = parser.parse_args()

    # Default: train and generate code
    if not args.train and args.predict is None:
        args.train = True

    # Read cache
    cache_path = os.path.join(os.getcwd(), args.cache)
    entries = read_cache_file(cache_path)

    if not entries:
        print(f"No cache entries found in {cache_path}")
        print("Run generate_training_data.py first to populate the cache")
        return

    print(f"Loaded {len(entries)} cache entries from {cache_path}")
    print()

    # Train model
    if args.train:
        model = train_ml_model(entries)
        if model:
            generate_cpp_code(model)

    # Predict for specific size
    if args.predict is not None:
        model = train_ml_model(entries)
        eps, max_rank = predict_parameters(model, args.predict)

        print("="*80)
        print(f"PREDICTION FOR N={args.predict}")
        print("="*80)
        print()
        print(f"Optimal parameters:")
        print(f"  eps = {eps:.6f}")
        print(f"  max_rank = {max_rank}")
        print()

        # Estimate construction time
        sizes = np.array([e['num_elements'] for e in entries])
        times = np.array([e['construction_time'] for e in entries])

        if len(sizes) >= 2:
            log_sizes = np.log(sizes)
            log_times = np.log(times)
            m, c = np.polyfit(log_sizes, log_times, 1)
            a = np.exp(c)
            b = m

            predicted_time = a * (args.predict ** b)
            print(f"Estimated construction time: {predicted_time:.3f} s")
            print(f"  (based on power law: time = {a:.6e} * N^{b:.4f})")

if __name__ == "__main__":
    main()
