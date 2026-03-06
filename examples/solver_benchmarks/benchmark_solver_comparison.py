#!/usr/bin/env python3
"""
Comprehensive Solver Comparison Benchmark
==========================================

Compares three solver methods for magnetostatic problems:
1. Method 0 (LU) - Direct solver, O(N^3) complexity
2. Method 1 (BiCGSTAB) - Iterative solver, O(N^2) per iteration
3. Method 2 (HACApK) - H-matrix accelerated, O(N^2 log N) with O(N log N) memory

This benchmark demonstrates performance trade-offs between methods.

Problem: Nonlinear magnetic material (soft iron) with applied background field
"""

import sys
import time
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../build/Release'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

import radia as rad
import numpy as np

# Set unit system to meters
mm = 1e-3  # 1 mm in meters


def hex_vertices(cx, cy, cz, dx, dy, dz):
    """Generate hexahedron vertices from center and dimensions."""
    hx, hy, hz = dx/2, dy/2, dz/2
    return [
        [cx-hx, cy-hy, cz-hz], [cx+hx, cy-hy, cz-hz],
        [cx+hx, cy+hy, cz-hz], [cx-hx, cy+hy, cz-hz],
        [cx-hx, cy-hy, cz+hz], [cx+hx, cy-hy, cz+hz],
        [cx+hx, cy+hy, cz+hz], [cx-hx, cy+hy, cz+hz]
    ]


def create_geometry(n, size, mat_params, B_bg):
    """Create hexahedral mesh with material and background field."""
    elem_size = size / n
    mat = rad.MatSatIsoFrm(*mat_params)

    elements = []
    for i in range(n):
        for j in range(n):
            for k in range(n):
                x = -size/2 + (i + 0.5) * elem_size
                y = -size/2 + (j + 0.5) * elem_size
                z = -size/2 + (k + 0.5) * elem_size
                vertices = hex_vertices(x, y, z, elem_size, elem_size, elem_size)
                elem = rad.ObjHexahedron(vertices, [0, 0, 0])
                rad.MatApl(elem, mat)
                elements.append(elem)

    bg_field = rad.ObjBckg(lambda p: B_bg)
    return rad.ObjCnt(elements + [bg_field])


print("=" * 80)
print("Comprehensive Solver Comparison Benchmark")
print("Method 0 (LU) vs Method 1 (BiCGSTAB) vs Method 2 (HACApK)")
print("=" * 80)

# Test configurations
test_cases = [
    {"n": 3, "desc": "Small (N=27)"},
    {"n": 5, "desc": "Medium (N=125)"},
    {"n": 7, "desc": "Large (N=343)"},
]

# Solver parameters
precision = 0.0001
max_iter = 1000

# Background field
B_bg = [1.0, 0, 0]

# Material parameters (soft iron)
mat_params = [[1596.3, 1.1488], [133.11, 0.4268], [18.713, 0.4759]]

# Observation point
obs_point = [0, 0, 50 * mm]

print("\nProblem Setup:")
print("  Material: Nonlinear soft iron (MatSatIsoFrm)")
print("  Background field: {} T".format(B_bg))
print("  Observation point: [{}, {}, {}] mm".format(obs_point[0]/mm, obs_point[1]/mm, obs_point[2]/mm))
print("  Solver precision: {}".format(precision))
print("  Max iterations: {}".format(max_iter))

results = []

for test in test_cases:
    n = test["n"]
    desc = test["desc"]
    n_elem = n ** 3
    size = 20.0 * mm

    print("\n" + "=" * 80)
    print("Test Case: {} ({} elements)".format(desc, n_elem))
    print("=" * 80)

    result = {"n": n, "n_elem": n_elem}

    # ========================================
    # Method 0: LU (Direct)
    # ========================================
    print("\n[Method 0] LU Direct Solver")
    print("-" * 40)
    rad.UtiDelAll()
    container = create_geometry(n, size, mat_params, B_bg)

    t_start = time.perf_counter()
    solve_result = rad.Solve(container, precision, max_iter, 0)
    t_lu = time.perf_counter() - t_start

    H_lu = rad.Fld(container, 'h', obs_point)
    iter_lu = int(solve_result[1]) if isinstance(solve_result, (list, tuple)) else 0

    print("  Time: {:.4f} s, Iterations: {}".format(t_lu, iter_lu))
    print("  |H|: {:.6e} A/m".format(np.linalg.norm(H_lu)))
    result['t_lu'] = t_lu
    result['H_lu'] = H_lu

    # ========================================
    # Method 1: BiCGSTAB (Iterative)
    # ========================================
    print("\n[Method 1] BiCGSTAB Iterative Solver")
    print("-" * 40)
    rad.UtiDelAll()
    container = create_geometry(n, size, mat_params, B_bg)

    t_start = time.perf_counter()
    solve_result = rad.Solve(container, precision, max_iter, 1)
    t_bicg = time.perf_counter() - t_start

    H_bicg = rad.Fld(container, 'h', obs_point)
    iter_bicg = int(solve_result[1]) if isinstance(solve_result, (list, tuple)) else 0

    print("  Time: {:.4f} s, Iterations: {}".format(t_bicg, iter_bicg))
    print("  |H|: {:.6e} A/m".format(np.linalg.norm(H_bicg)))
    result['t_bicg'] = t_bicg
    result['H_bicg'] = H_bicg

    # ========================================
    # Method 2: HACApK (H-matrix)
    # ========================================
    print("\n[Method 2] HACApK H-matrix Solver")
    print("-" * 40)
    try:
        rad.UtiDelAll()
        container = create_geometry(n, size, mat_params, B_bg)
        rad.SolverConfig(hacapk_eps=1e-4, hacapk_leaf=10, hacapk_eta=2.0)

        t_start = time.perf_counter()
        solve_result = rad.Solve(container, precision, max_iter, 2)
        t_hacapk = time.perf_counter() - t_start

        H_hacapk = rad.Fld(container, 'h', obs_point)
        iter_hacapk = int(solve_result[1]) if isinstance(solve_result, (list, tuple)) else 0

        print("  Time: {:.4f} s, Iterations: {}".format(t_hacapk, iter_hacapk))
        print("  |H|: {:.6e} A/m".format(np.linalg.norm(H_hacapk)))
        result['t_hacapk'] = t_hacapk
        result['H_hacapk'] = H_hacapk
    except Exception as e:
        print("  HACApK not available: {}".format(e))
        result['t_hacapk'] = float('inf')
        result['H_hacapk'] = None

    # ========================================
    # Comparison
    # ========================================
    print("\n[Comparison]")
    print("-" * 40)

    best_time = min(result['t_lu'], result['t_bicg'], result.get('t_hacapk', float('inf')))

    print("  Speedup vs LU:")
    print("    BiCGSTAB: {:.2f}x".format(result['t_lu'] / result['t_bicg']))
    if result.get('t_hacapk', float('inf')) < float('inf'):
        print("    HACApK:   {:.2f}x".format(result['t_lu'] / result['t_hacapk']))

    # Accuracy comparison (vs LU as reference)
    H_ref = np.array(result['H_lu'])
    H_ref_mag = np.linalg.norm(H_ref)

    err_bicg = np.linalg.norm(np.array(result['H_bicg']) - H_ref) / (H_ref_mag + 1e-15) * 100
    print("\n  Accuracy vs LU:")
    print("    BiCGSTAB: {:.6f}%".format(err_bicg))

    if result.get('H_hacapk') is not None:
        err_hacapk = np.linalg.norm(np.array(result['H_hacapk']) - H_ref) / (H_ref_mag + 1e-15) * 100
        print("    HACApK:   {:.6f}%".format(err_hacapk))

    results.append(result)

# ========================================
# Summary
# ========================================
print("\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)

print("\n{:>6s} {:>8s}  {:>10s} {:>10s} {:>10s}  {:>10s}".format(
    "N", "Elements", "LU (s)", "BiCG (s)", "HACApK (s)", "Best"))
print("-" * 70)

for r in results:
    t_hacapk = r.get('t_hacapk', float('inf'))
    best_time = min(r['t_lu'], r['t_bicg'], t_hacapk)

    if best_time == r['t_lu']:
        best_name = "LU"
    elif best_time == r['t_bicg']:
        best_name = "BiCGSTAB"
    else:
        best_name = "HACApK"

    hacapk_str = "{:.4f}".format(t_hacapk) if t_hacapk < float('inf') else "N/A"

    print("{:>6d} {:>8d}  {:>10.4f} {:>10.4f} {:>10s}  {:>10s}".format(
        r["n"], r["n_elem"],
        r['t_lu'], r['t_bicg'], hacapk_str, best_name))

print("\n" + "=" * 80)
print("RECOMMENDATIONS")
print("=" * 80)

print("""
Solver Selection:

  Method 0 (LU):
    - Best for: Small problems (N < 200)
    - Pros: Direct solution, guaranteed convergence
    - Cons: O(N^3) time, O(N^2) memory

  Method 1 (BiCGSTAB):
    - Best for: Medium problems (200 < N < 2000)
    - Pros: General purpose, good convergence
    - Cons: O(N^2) per iteration

  Method 2 (HACApK):
    - Best for: Large problems (N > 1000)
    - Pros: O(N log N) memory, scales well
    - Cons: Setup overhead for small problems

Usage:
  rad.Solve(obj, precision, max_iterations, method)

Examples:
  rad.Solve(obj, 0.0001, 1000, 0)  # LU for small problem
  rad.Solve(obj, 0.0001, 1000, 1)  # BiCGSTAB for medium
  rad.Solve(obj, 0.0001, 1000, 2)  # HACApK for large
""")

print("=" * 80)
rad.UtiDelAll()
