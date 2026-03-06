#!/usr/bin/env python
"""
Solver Methods Benchmark: LU vs BiCGSTAB vs HACApK

Compares three solver methods for magnetic field problems:
1. Method 0 (LU): Direct solver
2. Method 1 (BiCGSTAB): Iterative solver
3. Method 2 (HACApK): H-matrix accelerated solver

Problem: Nonlinear magnetic material with applied background field

Performance metrics:
- Computation time
- Field accuracy (comparing all methods)
"""

import sys
import os
import time
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
        B_bg: Background field in Tesla (NOT H in A/m).
              ObjBckg callback returns B, not H.
    """
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


# MKL/OpenMP warmup - first Solve() includes library initialization overhead
# Running a small warmup problem before benchmarks ensures accurate timing
print("Warming up MKL/OpenMP...")
rad.UtiDelAll()
_warmup_size = 0.01 * mm
_warmup_verts = [
    [-_warmup_size/2, -_warmup_size/2, -_warmup_size/2],
    [_warmup_size/2, -_warmup_size/2, -_warmup_size/2],
    [_warmup_size/2, _warmup_size/2, -_warmup_size/2],
    [-_warmup_size/2, _warmup_size/2, -_warmup_size/2],
    [-_warmup_size/2, -_warmup_size/2, _warmup_size/2],
    [_warmup_size/2, -_warmup_size/2, _warmup_size/2],
    [_warmup_size/2, _warmup_size/2, _warmup_size/2],
    [-_warmup_size/2, _warmup_size/2, _warmup_size/2]
]
_warmup_elem = rad.ObjHexahedron(_warmup_verts, [0, 0, 0])
_warmup_mat = rad.MatLin(100)
rad.MatApl(_warmup_elem, _warmup_mat)
_warmup_bg = rad.ObjBckg(lambda p: [0, 0, 0.01])
_warmup_cnt = rad.ObjCnt([_warmup_elem, _warmup_bg])
rad.Solve(_warmup_cnt, 0.001, 10, 0)  # LU warmup
rad.UtiDelAll()
print("Warmup complete.\n")

print("=" * 80)
print("Solver Methods Benchmark")
print("LU (Method 0) / BiCGSTAB (Method 1) / HACApK (Method 2)")
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

# Background field (uniform field applied to material)
# ObjBckg callback returns B in Tesla, not H in A/m
B_bg = [1.0, 0, 0]  # 1.0 T background field in X direction

# Observation point
obs_point = [0, 0, 50 * mm]  # 50mm above center

# Nonlinear material (soft iron)
mat_params = [[1596.3, 1.1488], [133.11, 0.4268], [18.713, 0.4759]]

print("\nProblem Setup:")
print("  Material: Nonlinear (Steel37)")
print("  Background field: [{}, {}, {}] T".format(*B_bg))
print("  Observation point: [{}, {}, {}] mm".format(obs_point[0]/mm, obs_point[1]/mm, obs_point[2]/mm))
print("  Solver precision: {}".format(precision))
print("  Max iterations: {}".format(max_iter))

results = []

for test in test_cases:
    n = test["n"]
    desc = test["desc"]
    n_elem = n ** 3
    size = 20.0 * mm  # 20mm (total cube size)

    print("\n" + "=" * 80)
    print("Test Case: {} ({} elements)".format(desc, n_elem))
    print("=" * 80)

    result_entry = {"n": n, "n_elem": n_elem}

    # ========================================
    # Method 0: LU (Direct)
    # ========================================
    print("\n[Method 0] LU Direct Solver")
    print("-" * 80)
    rad.UtiDelAll()
    container = create_geometry(n, size, mat_params, B_bg)

    t0 = time.perf_counter()
    result_lu = rad.Solve(container, precision, max_iter, 0)
    t_lu = time.perf_counter() - t0

    H_lu = rad.Fld(container, 'h', obs_point)
    iter_lu = int(result_lu[1]) if isinstance(result_lu, (list, tuple)) else 0

    print("  Time: {:.6f} s, Iterations: {}".format(t_lu, iter_lu))
    print("  H: [{:.6e}, {:.6e}, {:.6e}] A/m".format(*H_lu))
    H_mag_lu = np.linalg.norm(H_lu)
    print("  |H|: {:.6e} A/m".format(H_mag_lu))

    result_entry["t_lu"] = t_lu
    result_entry["H_lu"] = H_lu

    # ========================================
    # Method 1: BiCGSTAB (Iterative)
    # ========================================
    print("\n[Method 1] BiCGSTAB Iterative Solver")
    print("-" * 80)
    rad.UtiDelAll()
    container = create_geometry(n, size, mat_params, B_bg)

    t0 = time.perf_counter()
    result_bicg = rad.Solve(container, precision, max_iter, 1)
    t_bicg = time.perf_counter() - t0

    H_bicg = rad.Fld(container, 'h', obs_point)
    iter_bicg = int(result_bicg[1]) if isinstance(result_bicg, (list, tuple)) else 0

    print("  Time: {:.6f} s, Iterations: {}".format(t_bicg, iter_bicg))
    print("  H: [{:.6e}, {:.6e}, {:.6e}] A/m".format(*H_bicg))
    H_mag_bicg = np.linalg.norm(H_bicg)
    print("  |H|: {:.6e} A/m".format(H_mag_bicg))

    result_entry["t_bicg"] = t_bicg
    result_entry["H_bicg"] = H_bicg

    # ========================================
    # Method 2: HACApK (H-matrix)
    # ========================================
    print("\n[Method 2] HACApK H-matrix Accelerated Solver")
    print("-" * 80)
    try:
        rad.UtiDelAll()
        container = create_geometry(n, size, mat_params, B_bg)
        rad.SolverConfig(hacapk_eps=1e-4, hacapk_leaf=10, hacapk_eta=2.0)

        t0 = time.perf_counter()
        result_hmat = rad.Solve(container, precision, max_iter, 2)
        t_hmat = time.perf_counter() - t0

        H_hmat = rad.Fld(container, 'h', obs_point)
        iter_hmat = int(result_hmat[1]) if isinstance(result_hmat, (list, tuple)) else 0

        print("  Time: {:.6f} s, Iterations: {}".format(t_hmat, iter_hmat))
        print("  H: [{:.6e}, {:.6e}, {:.6e}] A/m".format(*H_hmat))
        H_mag_hmat = np.linalg.norm(H_hmat)
        print("  |H|: {:.6e} A/m".format(H_mag_hmat))

        result_entry["t_hmat"] = t_hmat
        result_entry["H_hmat"] = H_hmat
    except Exception as e:
        print("  HACApK not available: {}".format(e))
        result_entry["t_hmat"] = float('inf')
        result_entry["H_hmat"] = None

    # ========================================
    # Comparison
    # ========================================
    print("\n[Comparison]")
    print("-" * 80)

    # Use LU as reference
    H_ref = np.array(H_lu)
    H_ref_mag = np.linalg.norm(H_ref)

    # BiCGSTAB vs LU error
    diff_bicg = np.array(H_bicg) - H_ref
    err_bicg = np.linalg.norm(diff_bicg) / (H_ref_mag + 1e-15) * 100

    # HACApK vs LU error
    if result_entry.get("H_hmat") is not None:
        diff_hmat = np.array(result_entry["H_hmat"]) - H_ref
        err_hmat = np.linalg.norm(diff_hmat) / (H_ref_mag + 1e-15) * 100
    else:
        err_hmat = float('inf')

    print("  Reference: LU |H| = {:.6e} A/m".format(H_ref_mag))
    print("")
    print("  Method           Time (s)    |H| (A/m)        Error (%)   Speedup")
    print("  " + "-" * 72)
    print("  LU             {:9.6f}  {:14.6e}  {:11.4f}     1.00x".format(
        t_lu, H_mag_lu, 0.0))
    print("  BiCGSTAB       {:9.6f}  {:14.6e}  {:11.4f}  {:7.2f}x".format(
        t_bicg, H_mag_bicg, err_bicg, t_lu/t_bicg if t_bicg > 0 else 0))
    if result_entry.get("t_hmat", float('inf')) < float('inf'):
        print("  HACApK         {:9.6f}  {:14.6e}  {:11.4f}  {:7.2f}x".format(
            result_entry["t_hmat"], H_mag_hmat, err_hmat, t_lu/result_entry["t_hmat"] if result_entry["t_hmat"] > 0 else 0))

    result_entry["err_bicg"] = err_bicg
    result_entry["err_hmat"] = err_hmat
    result_entry["speedup_bicg"] = t_lu/t_bicg if t_bicg > 0 else 0
    result_entry["speedup_hmat"] = t_lu/result_entry.get("t_hmat", float('inf')) if result_entry.get("t_hmat", float('inf')) < float('inf') else 0

    results.append(result_entry)

# ========================================
# Summary
# ========================================
print("\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)
print("\n{:>6s} {:>8s}  {:>12s} {:>12s} {:>12s}  {:>10s} {:>10s}".format(
    "N", "Elements",
    "LU (s)", "BiCGSTAB (s)", "HACApK (s)",
    "BiCG Spdup", "HACApK Spdup"))
print("-" * 90)

for r in results:
    hmat_str = "{:.6f}".format(r.get("t_hmat", float('inf'))) if r.get("t_hmat", float('inf')) < float('inf') else "N/A"
    hmat_spdup = "{:.2f}x".format(r.get("speedup_hmat", 0)) if r.get("speedup_hmat", 0) > 0 else "N/A"

    print("{:>6d} {:>8d}  {:>12.6f} {:>12.6f} {:>12s}  {:>10.2f}x {:>10s}".format(
        r["n"], r["n_elem"],
        r["t_lu"], r["t_bicg"], hmat_str,
        r.get("speedup_bicg", 0), hmat_spdup))

print("=" * 80)
print("\nNotes:")
print("  - LU: Direct solver, best for small problems (N < 200)")
print("  - BiCGSTAB: Iterative solver, good for medium problems")
print("  - HACApK: H-matrix accelerated, best for large problems (N > 1000)")
print("  - Error: Relative difference from LU method")
print("=" * 80)

rad.UtiDelAll()
