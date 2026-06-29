#!/usr/bin/env python
"""
Run benchmarks based on benchmark_config.json

This script reads the configuration and runs specified benchmarks
in subprocesses for accurate memory measurement.

Usage:
    # Run all nonlinear benchmarks
    python run_benchmarks.py nonlinear

    # Run specific element type
    python run_benchmarks.py nonlinear/hexahedron
    python run_benchmarks.py nonlinear/tetrahedron

    # Run full suite
    python run_benchmarks.py full

    # Run specific solvers only
    python run_benchmarks.py nonlinear/hexahedron --solvers lu hacapk

    # Run specific sizes only
    python run_benchmarks.py nonlinear/hexahedron --n_div 5 10
    python run_benchmarks.py nonlinear/tetrahedron --maxh 0.25 0.15
"""

import os
import sys
import json
import argparse
import subprocess
from pathlib import Path


def load_config(config_path: str) -> dict:
    """Load benchmark configuration from JSON file."""
    with open(config_path, 'r') as f:
        return json.load(f)


def run_benchmark(script_path: str, solver: str, sizes: list, size_type: str,
                  params: dict, verbose: bool = True) -> list:
    """Run a single benchmark script with specified parameters.

    Args:
        script_path: Path to benchmark script
        solver: Solver name ('lu', 'bicgstab', 'hacapk')
        sizes: List of sizes (n_div for hex, maxh for tetra)
        size_type: 'n_div' or 'maxh'
        params: Default parameters dict
        verbose: Print output

    Returns:
        List of result dictionaries
    """
    if not os.path.exists(script_path):
        print('[SKIP] Script not found: %s' % script_path)
        return []

    # Build command
    cmd = [sys.executable, script_path, '--%s' % solver]

    # Add parameters
    cmd.extend(['--hmat_eps', str(params.get('hmat_eps', 1e-4))])
    cmd.extend(['--bicg_tol', str(params.get('bicg_tol', 1e-4))])
    cmd.extend(['--nonl_tol', str(params.get('nonl_tol', 0.001))])
    cmd.extend(['--hmat_leaf_size', str(params.get('hmat_leaf_size', 10))])
    cmd.extend(['--hmat_eta', str(params.get('hmat_eta', 2.0))])

    # Add sizes
    for size in sizes:
        cmd.append(str(size))

    if verbose:
        print('Running: %s' % ' '.join(cmd))

    # Run benchmark
    result = subprocess.run(cmd, capture_output=False)

    if result.returncode != 0:
        print('[ERROR] Benchmark failed: %s' % script_path)
        return []

    return []  # Results are saved to JSON files by the benchmark script


def run_suite(config: dict, suite_name: str, script_dir: str, params: dict,
              solver_filter: list = None, size_filter: list = None,
              verbose: bool = True):
    """Run a benchmark suite.

    Args:
        config: Configuration dictionary
        suite_name: 'nonlinear', 'linear', 'nonlinear/hexahedron', etc.
        script_dir: Base directory for scripts
        params: Default parameters
        solver_filter: If provided, only run these solvers
        size_filter: If provided, only run these sizes
        verbose: Print output
    """
    parts = suite_name.split('/')

    if len(parts) == 1:
        # Run all element types in category (e.g., 'nonlinear')
        category = parts[0]
        if category not in config:
            print('[ERROR] Unknown category: %s' % category)
            return

        for elem_type, elem_config in config[category].items():
            run_element_benchmark(elem_config, script_dir, params,
                                  solver_filter, size_filter, verbose)

    elif len(parts) == 2:
        # Run specific element type (e.g., 'nonlinear/hexahedron')
        category, elem_type = parts
        if category not in config or elem_type not in config[category]:
            print('[ERROR] Unknown benchmark: %s' % suite_name)
            return

        elem_config = config[category][elem_type]
        run_element_benchmark(elem_config, script_dir, params,
                              solver_filter, size_filter, verbose)
    else:
        print('[ERROR] Invalid suite name: %s' % suite_name)


def run_element_benchmark(elem_config: dict, script_dir: str, params: dict,
                          solver_filter: list = None, size_filter: list = None,
                          verbose: bool = True):
    """Run benchmark for a specific element configuration."""
    script_path = os.path.join(script_dir, elem_config['script'])
    solvers = solver_filter or elem_config.get('solvers', ['lu'])

    # Determine size type and values
    if 'n_div' in elem_config:
        size_type = 'n_div'
        sizes = size_filter or elem_config['n_div']
    elif 'maxh' in elem_config:
        size_type = 'maxh'
        sizes = size_filter or elem_config['maxh']
    else:
        print('[ERROR] No sizes specified in config')
        return

    # Run each solver
    for solver in solvers:
        print()
        print('=' * 70)
        print('Solver: %s, Sizes: %s' % (solver, sizes))
        print('=' * 70)
        run_benchmark(script_path, solver, sizes, size_type, params, verbose)


def run_full_suite(config: dict, script_dir: str, params: dict,
                   verbose: bool = True):
    """Run the full benchmark suite."""
    if 'full_suite' not in config:
        print('[ERROR] No full_suite defined in config')
        return

    benchmarks = config['full_suite'].get('benchmarks', [])

    for bench in benchmarks:
        bench_type = bench.get('type', '')
        solvers = bench.get('solvers', ['lu'])
        sizes = bench.get('n_div') or bench.get('maxh', [])

        print()
        print('=' * 70)
        print('Benchmark: %s' % bench_type)
        print('Solvers: %s, Sizes: %s' % (solvers, sizes))
        print('=' * 70)

        # Determine script path
        parts = bench_type.split('/')
        if len(parts) == 2:
            category, elem_type = parts
            if category in config and elem_type in config[category]:
                script_path = os.path.join(script_dir, config[category][elem_type]['script'])
                size_type = 'n_div' if 'n_div' in bench else 'maxh'

                for solver in solvers:
                    run_benchmark(script_path, solver, sizes, size_type, params, verbose)


def main():
    parser = argparse.ArgumentParser(
        description='Run Radia benchmarks from config file',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument('suite', nargs='?', default='nonlinear',
                        help='Benchmark suite: nonlinear, linear, nonlinear/hexahedron, '
                             'nonlinear/tetrahedron, full (default: nonlinear)')
    parser.add_argument('--config', default='benchmark_config.json',
                        help='Path to config file (default: benchmark_config.json)')
    parser.add_argument('--solvers', nargs='+',
                        help='Filter: only run these solvers (lu, bicgstab, hacapk)')
    parser.add_argument('--n_div', nargs='+', type=int,
                        help='Filter: only run these n_div values (for hexahedron)')
    parser.add_argument('--maxh', nargs='+', type=float,
                        help='Filter: only run these maxh values (for tetrahedron)')
    parser.add_argument('--hmat_eps', type=float, default=1e-4,
                        help='ACA tolerance (default: 1e-4)')
    parser.add_argument('--quiet', action='store_true',
                        help='Suppress verbose output')
    args = parser.parse_args()

    # Find script directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    config_path = os.path.join(script_dir, args.config)

    if not os.path.exists(config_path):
        print('[ERROR] Config file not found: %s' % config_path)
        sys.exit(1)

    # Load config
    config = load_config(config_path)

    # Build params
    params = config.get('default_params', {})
    params['hmat_eps'] = args.hmat_eps

    # Determine size filter
    size_filter = None
    if args.n_div:
        size_filter = args.n_div
    elif args.maxh:
        size_filter = args.maxh

    verbose = not args.quiet

    print('=' * 70)
    print('RADIA BENCHMARK RUNNER')
    print('=' * 70)
    print('Config: %s' % config_path)
    print('Suite: %s' % args.suite)
    print('Params: hmat_eps=%.0e' % params['hmat_eps'])
    print()

    if args.suite == 'full':
        run_full_suite(config, script_dir, params, verbose)
    else:
        run_suite(config, args.suite, script_dir, params,
                  args.solvers, size_filter, verbose)

    print()
    print('=' * 70)
    print('BENCHMARK COMPLETE')
    print('=' * 70)


if __name__ == '__main__':
    main()
