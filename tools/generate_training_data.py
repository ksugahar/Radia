"""
Generate training data for ML parameter tuning

This script runs multiple H-matrix constructions with different problem sizes
and parameters to populate the cache with diverse training data.

Usage:
    python tools/generate_training_data.py
"""

import sys
import os
import time

# Add paths for Radia module
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'build', 'Release'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'dist'))

import radia as rad

print("="*80)
print("TRAINING DATA GENERATION FOR ML PARAMETER TUNING")
print("="*80)
print()

def create_soft_iron_cube(n):
    """Create a soft iron cube with n^3 subdivisions"""
    # Soft iron material (standard soft iron parameters)
    mat = rad.MatSatIsoFrm([1596.3, 1.1488], [133.11, 0.4268], [18.713, 0.4759])

    # Create cube
    cube = rad.ObjRecMag([0, 0, 0], [100, 100, 100], [0, 0, 0])
    rad.MatApl(cube, mat)

    # Subdivide
    rad.ObjDivMag(cube, [n, n, n])

    return cube, n**3

# Test matrix: different problem sizes and parameter combinations
test_matrix = [
    # (subdivision, eps, max_rank, description)
    (5, 1e-4, 30, "Small problem, default params"),
    (5, 2e-4, 25, "Small problem, relaxed params"),
    (5, 5e-4, 20, "Small problem, very relaxed params"),

    (6, 1e-4, 30, "Medium problem, default params"),
    (6, 2e-4, 25, "Medium problem, relaxed params"),
    (6, 5e-4, 20, "Medium problem, very relaxed params"),

    (7, 1e-4, 30, "Large problem, default params"),
    (7, 2e-4, 25, "Large problem, relaxed params"),
    (7, 5e-4, 20, "Large problem, very relaxed params"),

    (8, 1e-4, 30, "Very large problem, default params"),
    (8, 2e-4, 25, "Very large problem, relaxed params"),
    (8, 5e-4, 20, "Very large problem, very relaxed params"),

    (10, 2e-4, 25, "Extra large problem, relaxed params"),
    (10, 5e-4, 20, "Extra large problem, very relaxed params"),
]

print(f"Running {len(test_matrix)} test configurations...")
print()

results = []
for i, (n, eps, max_rank, desc) in enumerate(test_matrix, 1):
    print(f"[{i}/{len(test_matrix)}] {desc}")
    print(f"  N = {n}x{n}x{n} = {n**3} elements, eps={eps}, max_rank={max_rank}")

    rad.UtiDelAll()

    try:
        # Create geometry
        cube, N = create_soft_iron_cube(n)

        # Configure H-matrix
        rad.SolverHMatrixEnable(1, eps, max_rank)

        # Run solver (measures construction time)
        print(f"  Building H-matrix...")
        t0 = time.time()
        result = rad.Solve(cube, 0.0001, 1000)
        elapsed = time.time() - t0

        print(f"  Time: {elapsed*1000:.1f} ms")
        print(f"  Result: {result}")

        results.append({
            'N': N,
            'eps': eps,
            'max_rank': max_rank,
            'time': elapsed,
            'desc': desc
        })

    except Exception as e:
        print(f"  ERROR: {e}")

    print()

# Summary
print("="*80)
print("TRAINING DATA GENERATION COMPLETE")
print("="*80)
print()
print(f"Successfully completed: {len(results)}/{len(test_matrix)} configurations")
print()

if results:
    print("Summary:")
    print(f"{'N':<8} {'eps':<10} {'rank':<6} {'Time (ms)':<12} {'Description':<40}")
    print("-"*80)
    for r in results:
        print(f"{r['N']:<8} {r['eps']:<10.6f} {r['max_rank']:<6} {r['time']*1000:<12.1f} {r['desc']:<40}")

print()
print("Cache file updated: ./.radia_cache/hmatrix_cache.bin")
print("Run 'python tools/analyze_cache_for_ml.py' to analyze the training data")
print("="*80)
