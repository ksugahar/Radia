"""
Verify H-matrix parameter accuracy

This script compares H-matrix results with different parameter settings
to verify that ML-optimized parameters maintain acceptable accuracy.

Usage:
    python tools/verify_parameter_accuracy.py
"""

import sys
import os
import time
import numpy as np

# Add paths for Radia module
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'build', 'Release'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'dist'))

import radia as rad

print("="*80)
print("H-MATRIX PARAMETER ACCURACY VERIFICATION")
print("="*80)
print()

def create_test_geometry(n):
    """Create a test geometry for accuracy verification"""
    mat = rad.MatSatIsoFrm([1596.3, 1.1488], [133.11, 0.4268], [18.713, 0.4759])
    cube = rad.ObjRecMag([0, 0, 0], [100, 100, 100], [0, 0, 0])
    rad.MatApl(cube, mat)
    rad.ObjDivMag(cube, [n, n, n])
    return cube

def solve_with_params(cube, eps, max_rank, desc):
    """Solve with specific parameters and return result"""
    rad.SolverHMatrixEnable(1, eps, max_rank)
    t0 = time.time()
    result = rad.Solve(cube, 0.0001, 1000)
    elapsed = time.time() - t0
    return result, elapsed

# Test configurations
test_cases = [
    {
        'name': 'Medium problem (N=343)',
        'subdivision': 7,
        'params': [
            (1e-4, 30, "Current default (conservative)"),
            (2e-4, 25, "Current relaxed (medium)"),
            (5e-4, 20, "ML-optimized (aggressive)"),
        ]
    },
    {
        'name': 'Large problem (N=512)',
        'subdivision': 8,
        'params': [
            (1e-4, 30, "Current default (conservative)"),
            (2e-4, 25, "Current relaxed (medium)"),
            (5e-4, 20, "ML-optimized (aggressive)"),
        ]
    }
]

print("This test compares:")
print("  1. Current default parameters (conservative, high accuracy)")
print("  2. Current relaxed parameters (balanced)")
print("  3. ML-optimized parameters (aggressive, fastest)")
print()
print("Goal: Verify ML-optimized parameters maintain acceptable accuracy")
print()

for test_case in test_cases:
    print("="*80)
    print(test_case['name'])
    print("="*80)
    print()

    n = test_case['subdivision']
    N = n**3

    results = []

    for eps, max_rank, desc in test_case['params']:
        rad.UtiDelAll()
        cube = create_test_geometry(n)

        print(f"Testing: {desc}")
        print(f"  eps={eps}, max_rank={max_rank}")

        try:
            result, elapsed = solve_with_params(cube, eps, max_rank, desc)

            print(f"  Time: {elapsed*1000:.1f} ms")
            print(f"  Result: {result}")

            results.append({
                'desc': desc,
                'eps': eps,
                'max_rank': max_rank,
                'time': elapsed,
                'result': result
            })

        except Exception as e:
            print(f"  ERROR: {e}")

        print()

    # Compare results
    if len(results) >= 2:
        print("Comparison:")
        print(f"{'Config':<40} {'Time (ms)':<12} {'Speedup':<10}")
        print("-"*65)

        baseline_time = results[0]['time']
        for r in results:
            speedup = baseline_time / r['time']
            print(f"{r['desc']:<40} {r['time']*1000:<12.1f} {speedup:<10.2f}x")

        print()

        # Check convergence consistency
        baseline_result = results[0]['result']
        print("Convergence check:")
        for r in results:
            if r['result'] == baseline_result:
                print(f"  {r['desc']}: ✓ Same convergence as baseline")
            else:
                print(f"  {r['desc']}: ✗ Different convergence")

    print()

print("="*80)
print("SUMMARY")
print("="*80)
print()
print("Key findings:")
print("  - All parameter sets should converge to similar results")
print("  - ML-optimized parameters (eps=5e-4, max_rank=20) trade accuracy for speed")
print("  - If convergence differs significantly, more conservative parameters are safer")
print()
print("Recommendation:")
print("  - For production: Use current adaptive parameters (conservative)")
print("  - For exploration: ML-optimized parameters acceptable if accuracy sufficient")
print("  - For benchmarking: Test both and compare field accuracy")
print("="*80)
