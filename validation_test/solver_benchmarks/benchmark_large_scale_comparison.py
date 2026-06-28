#!/usr/bin/env python3
"""
Large-Scale Solver Comparison: LU vs BiCGSTAB vs HACApK
========================================================

Comprehensive benchmark comparing three solver methods across a wide range of problem sizes:
1. LU Decomposition (Method 0) - Direct solver, O(N^3)
2. BiCGSTAB (Method 1) - Iterative solver, O(N^2) per iteration
3. HACApK (Method 2) - H-matrix accelerated BiCGSTAB, O(N^2 log N) per iteration

This benchmark demonstrates:
- Performance scaling with problem size
- Memory usage comparison
- Crossover points where each method becomes optimal
- Convergence behavior

Problem: Nonlinear magnetic material (soft iron) with background field
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
    """Create hexahedral mesh with material and background field.

    Args:
        B_bg: Background field in Tesla (ObjBckg callback returns B, not H).
    """
    elem_size = size / n
    mat = rad.MatSatIsoFrm(mat_params)

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
print("Large-Scale Solver Comparison Benchmark")
print("LU (Method 0) vs BiCGSTAB (Method 1) vs HACApK (Method 2)")
print("=" * 80)

# Test configurations - progressively larger problems
test_cases = [
    {"n": 3, "desc": "Tiny (N=27)", "lu_enabled": True},
    {"n": 4, "desc": "Very Small (N=64)", "lu_enabled": True},
    {"n": 5, "desc": "Small (N=125)", "lu_enabled": True},
    {"n": 6, "desc": "Medium-Small (N=216)", "lu_enabled": True},
    {"n": 7, "desc": "Medium (N=343)", "lu_enabled": False},  # LU too slow
    {"n": 8, "desc": "Medium-Large (N=512)", "lu_enabled": False},
    {"n": 10, "desc": "Large (N=1000)", "lu_enabled": False},
]

# Solver parameters
precision = 0.0001
max_iter = 1000

# Background field (applied to soft magnetic material)
B_bg = [1.0, 0, 0]  # 1.0 T background field in X direction

# Observation point
obs_point = [0, 0, 50 * mm]  # 50mm above center

# Nonlinear material (soft iron)
mat_params = [[1596.3, 1.1488], [133.11, 0.4268], [18.713, 0.4759]]

print("\nProblem Setup:")
print("  Material: Nonlinear soft iron (MatSatIsoFrm)")
print("  Background field: [{}, {}, {}] T".format(*B_bg))
print("  Observation point: [{}, {}, {}] mm".format(obs_point[0]/mm, obs_point[1]/mm, obs_point[2]/mm))
print("  Solver precision: {}".format(precision))
print("  Max iterations: {}".format(max_iter))
print("\nTest cases: {} problem sizes".format(len(test_cases)))

results = []

for test in test_cases:
    n = test["n"]
    desc = test["desc"]
    lu_enabled = test["lu_enabled"]
    n_elem = n ** 3
    size = 20.0 * mm  # 20mm (total cube size)

    print("\n" + "=" * 80)
    print("Test Case: {} ({} elements)".format(desc, n_elem))
    print("=" * 80)

    result_entry = {
        "n": n,
        "n_elem": n_elem,
        "desc": desc,
    }

    # ========================================
    # Method 0: LU Decomposition
    # ========================================
    if lu_enabled:
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

        result_entry.update({
            "t_lu": t_lu,
            "iter_lu": iter_lu,
            "H_lu": H_lu,
        })
    else:
        print("\n[Method 0] LU - SKIPPED (too slow for N > 250)")
        result_entry.update({
            "t_lu": None,
            "iter_lu": None,
            "H_lu": None,
        })

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

    result_entry.update({
        "t_bicg": t_bicg,
        "iter_bicg": iter_bicg,
        "H_bicg": H_bicg,
    })

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

        result_entry.update({
            "t_hacapk": t_hacapk,
            "iter_hacapk": iter_hacapk,
            "H_hacapk": H_hacapk,
        })
    except Exception as e:
        print("  HACApK not available: {}".format(e))
        result_entry.update({
            "t_hacapk": float('inf'),
            "iter_hacapk": None,
            "H_hacapk": None,
        })

    # ========================================
    # Comparison
    # ========================================
    print("\n[Comparison]")
    print("-" * 40)

    # Use BiCGSTAB as reference
    H_ref = np.array(result_entry["H_bicg"])
    H_ref_mag = np.linalg.norm(H_ref)

    # Accuracy comparison
    print("  Accuracy vs BiCGSTAB:")
    if lu_enabled and result_entry["H_lu"] is not None:
        err_lu = np.linalg.norm(np.array(result_entry["H_lu"]) - H_ref) / (H_ref_mag + 1e-15) * 100
        print("    LU:     {:.6f}%".format(err_lu))
        result_entry["err_lu"] = err_lu
    else:
        result_entry["err_lu"] = None

    if result_entry.get("H_hacapk") is not None:
        err_hacapk = np.linalg.norm(np.array(result_entry["H_hacapk"]) - H_ref) / (H_ref_mag + 1e-15) * 100
        print("    HACApK: {:.6f}%".format(err_hacapk))
        result_entry["err_hacapk"] = err_hacapk
    else:
        result_entry["err_hacapk"] = None

    # Speedup comparison
    print("\n  Speedup vs BiCGSTAB:")
    if lu_enabled and result_entry["t_lu"] is not None:
        speedup_lu = t_bicg / result_entry["t_lu"]
        print("    LU:     {:.2f}x".format(speedup_lu))
    if result_entry.get("t_hacapk", float('inf')) < float('inf'):
        speedup_hacapk = t_bicg / result_entry["t_hacapk"]
        print("    HACApK: {:.2f}x".format(speedup_hacapk))

    results.append(result_entry)

# ========================================
# Summary
# ========================================
print("\n" + "=" * 80)
print("SUMMARY - Performance vs Problem Size")
print("=" * 80)

print("\n{:>6s} {:>8s}  {:>10s} {:>10s} {:>10s}  {:>10s}".format(
    "N", "Elements", "LU (s)", "BiCG (s)", "HACApK (s)", "Best"))
print("-" * 70)

for r in results:
    lu_str = "{:.4f}".format(r["t_lu"]) if r["t_lu"] is not None else "---"
    bicg_str = "{:.4f}".format(r["t_bicg"])
    hacapk_str = "{:.4f}".format(r["t_hacapk"]) if r.get("t_hacapk", float('inf')) < float('inf') else "N/A"

    # Determine best method
    times = []
    if r["t_lu"] is not None:
        times.append(("LU", r["t_lu"]))
    times.append(("BiCGSTAB", r["t_bicg"]))
    if r.get("t_hacapk", float('inf')) < float('inf'):
        times.append(("HACApK", r["t_hacapk"]))
    best_method = min(times, key=lambda x: x[1])[0]

    print("{:>6d} {:>8d}  {:>10s} {:>10s} {:>10s}  {:>10s}".format(
        r["n"], r["n_elem"],
        lu_str, bicg_str, hacapk_str, best_method))

# ========================================
# Scaling Analysis
# ========================================
print("\n" + "=" * 80)
print("SCALING ANALYSIS")
print("=" * 80)

# Extract data for scaling analysis
n_values = np.array([r["n_elem"] for r in results])
log_n = np.log(n_values)

# BiCGSTAB scaling
bicg_times = np.array([r["t_bicg"] for r in results])
log_bicg = np.log(bicg_times)
A = np.vstack([log_n, np.ones(len(log_n))]).T
alpha_bicg, log_a_bicg = np.linalg.lstsq(A, log_bicg, rcond=None)[0]

# HACApK scaling (filter valid results)
valid_hacapk = [r for r in results if r.get("t_hacapk", float('inf')) < float('inf')]
if len(valid_hacapk) >= 2:
    n_hac = np.array([r["n_elem"] for r in valid_hacapk])
    t_hac = np.array([r["t_hacapk"] for r in valid_hacapk])
    log_n_hac = np.log(n_hac)
    log_t_hac = np.log(t_hac)
    A_hac = np.vstack([log_n_hac, np.ones(len(log_n_hac))]).T
    alpha_hacapk, log_a_hacapk = np.linalg.lstsq(A_hac, log_t_hac, rcond=None)[0]
else:
    alpha_hacapk = None

print("\nComplexity Analysis:")
print("  BiCGSTAB: t = {:.6e} * N^{:.3f}".format(np.exp(log_a_bicg), alpha_bicg))
if alpha_hacapk is not None:
    print("  HACApK:   t = {:.6e} * N^{:.3f}".format(np.exp(log_a_hacapk), alpha_hacapk))

print("\nTheoretical complexity:")
print("  BiCGSTAB: O(N^2) per iteration -> expected alpha ~= 2.0")
print("  HACApK:   O(N^2 log N) per iteration -> expected alpha ~= 2.0-2.5")

print("\nMeasured complexity:")
print("  BiCGSTAB: alpha = {:.3f} -> {}".format(
    alpha_bicg,
    "Matches O(N^2)" if 1.8 <= alpha_bicg <= 2.2 else "Unexpected scaling"))
if alpha_hacapk is not None:
    print("  HACApK:   alpha = {:.3f} -> {}".format(
        alpha_hacapk,
        "Matches O(N^2 log N)" if 1.5 <= alpha_hacapk <= 2.5 else "Unexpected scaling"))

# ========================================
# Recommendations
# ========================================
print("\n" + "=" * 80)
print("RECOMMENDATIONS")
print("=" * 80)

print("""
Solver Selection Guidelines:

  Method 0 (LU):
    - Best for: Small problems (N < 200)
    - Complexity: O(N^3)
    - Pros: Direct solution, guaranteed convergence
    - Cons: Memory O(N^2), slow for large N

  Method 1 (BiCGSTAB):
    - Best for: Medium problems (200 < N < 2000)
    - Complexity: O(N^2) per iteration
    - Pros: Good general-purpose iterative solver
    - Cons: May need tuning for ill-conditioned systems

  Method 2 (HACApK):
    - Best for: Large problems (N > 1000)
    - Complexity: O(N log N) memory, O(N^2 log N) time
    - Pros: Memory efficient, scales well
    - Cons: Setup overhead, best for large N

Usage:
  rad.Solve(obj, precision, max_iterations, method)

  # Example
  rad.Solve(container, 0.0001, 1000, 2)  # HACApK for large problem
""")

print("=" * 80)
print("Benchmark complete!")
print("=" * 80)

rad.UtiDelAll()
