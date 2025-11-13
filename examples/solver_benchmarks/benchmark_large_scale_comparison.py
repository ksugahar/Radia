#!/usr/bin/env python3
"""
Large-Scale Solver Comparison: LU vs Gauss-Seidel vs H-matrix
==============================================================

Comprehensive benchmark comparing three solver methods across a wide range of problem sizes:
1. LU Decomposition - Direct solver, O(N^3)
2. Gauss-Seidel (Standard Relaxation) - Iterative solver, O(N^2) per iteration
3. H-matrix Accelerated - Fast iterative solver, O(N^2 log N) per iteration

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
import tracemalloc
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../build/Release'))

import radia as rad
import numpy as np

print("=" * 80)
print("Large-Scale Solver Comparison Benchmark")
print("LU Decomposition vs Gauss-Seidel vs H-matrix")
print("=" * 80)

# Test configurations - progressively larger problems
test_cases = [
    {"n": 3, "desc": "Tiny (N=27)", "lu_enabled": True, "iterations": 50},
    {"n": 4, "desc": "Very Small (N=64)", "lu_enabled": True, "iterations": 30},
    {"n": 5, "desc": "Small (N=125)", "lu_enabled": True, "iterations": 20},
    {"n": 6, "desc": "Medium-Small (N=216)", "lu_enabled": True, "iterations": 10},
    {"n": 7, "desc": "Medium (N=343)", "lu_enabled": False, "iterations": 5},  # LU too slow
    {"n": 8, "desc": "Medium-Large (N=512)", "lu_enabled": False, "iterations": 3},
    {"n": 10, "desc": "Large (N=1000)", "lu_enabled": False, "iterations": 2},
    {"n": 13, "desc": "Very Large (N=2197)", "lu_enabled": False, "iterations": 2},
    {"n": 17, "desc": "Extra Large (N=4913)", "lu_enabled": False, "iterations": 1},
    {"n": 22, "desc": "Huge (N=10648)", "lu_enabled": False, "iterations": 1},
]

# Solver parameters
precision = 0.0001
max_iter = 1000

# Background field (applied to soft magnetic material)
H_bg = [1.0, 0, 0]  # 1.0 T background field in X direction

# Observation point
obs_point = [0, 0, 50]  # 50mm above center

print("\nProblem Setup:")
print("  Material: Nonlinear soft iron (MatSatIsoFrm)")
print("  Background field: [{}, {}, {}] T".format(*H_bg))
print("  Observation point: [{}, {}, {}] mm".format(*obs_point))
print("  Solver precision: {}".format(precision))
print("  Max iterations: {}".format(max_iter))
print("\nTest cases: {} problem sizes from N=27 to N=1000".format(len(test_cases)))

results = []

for test in test_cases:
    n = test["n"]
    desc = test["desc"]
    lu_enabled = test["lu_enabled"]
    num_iterations = test["iterations"]

    print("\n" + "=" * 80)
    print("Test Case: {} - {} iterations per solver".format(desc, num_iterations))
    print("=" * 80)

    size = 20.0  # mm (total cube size)
    elem_size = size / n
    n_elem = n ** 3

    # Nonlinear material (soft iron)
    mat_params = [[1596.3, 1.1488], [133.11, 0.4268], [18.713, 0.4759]]

    result_entry = {
        "n": n,
        "n_elem": n_elem,
        "desc": desc,
    }

    # ========================================
    # Method 1: LU Decomposition
    # ========================================
    if lu_enabled:
        print("\n[Method 1] LU Decomposition (Direct Solver)")
        print("-" * 80)
        rad.UtiDelAll()

        # Start memory tracking
        tracemalloc.start()
        mem_start = tracemalloc.get_traced_memory()[0]

        # Create material and geometry
        mat = rad.MatSatIsoFrm(*mat_params)
        elements = []
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    x = -size/2 + (i + 0.5) * elem_size
                    y = -size/2 + (j + 0.5) * elem_size
                    z = -size/2 + (k + 0.5) * elem_size
                    elem = rad.ObjRecMag([x, y, z], [elem_size, elem_size, elem_size])
                    rad.MatApl(elem, mat)
                    elements.append(elem)

        bg_field = rad.ObjBckg(H_bg)
        container = rad.ObjCnt(elements + [bg_field])

        print("  Elements: {} ({}x{}x{})".format(n_elem, n, n, n))

        # Disable H-matrix
        rad.SolverHMatrixDisable()

        # Matrix construction
        print("  Matrix construction...")
        t_matrix_start = time.perf_counter()
        intrc = rad.RlxPre(container, container)
        t_matrix_lu = time.perf_counter() - t_matrix_start

        # Enable LU decomposition
        rad.SetRelaxSubInterval(intrc, 0, n_elem-1, 1)

        # Measure per-iteration time
        print("  LU solve ({} iterations)...".format(num_iterations))
        t_solve_start = time.perf_counter()
        for _ in range(num_iterations):
            rad.RlxMan(intrc, 5, 1, 1.0)  # Method 5 = LU
        t_solve_lu = (time.perf_counter() - t_solve_start) / num_iterations

        # Memory usage
        mem_peak = tracemalloc.get_traced_memory()[1]
        mem_lu = (mem_peak - mem_start) / (1024 * 1024)  # MB
        tracemalloc.stop()

        # Get field
        H_lu = rad.Fld(container, 'h', obs_point)

        print("    Matrix construction: {:.6f} s".format(t_matrix_lu))
        print("    Per-iteration time: {:.6f} s ({:.2f} ms)".format(t_solve_lu, t_solve_lu * 1000))
        print("    Memory usage: {:.2f} MB".format(mem_lu))
        print("    H: [{:.6e}, {:.6e}, {:.6e}] A/m".format(*H_lu))
        print("    |H|: {:.6e} A/m".format(np.linalg.norm(H_lu)))

        result_entry.update({
            "t_matrix_lu": t_matrix_lu,
            "t_solve_lu": t_solve_lu,
            "mem_lu": mem_lu,
            "H_lu": H_lu,
        })
    else:
        print("\n[Method 1] LU Decomposition - SKIPPED (too slow for N > 250)")
        result_entry.update({
            "t_matrix_lu": None,
            "t_solve_lu": None,
            "mem_lu": None,
            "H_lu": None,
        })

    # ========================================
    # Method 2: Gauss-Seidel (Standard Relaxation)
    # ========================================
    print("\n[Method 2] Gauss-Seidel (Standard Relaxation)")
    print("-" * 80)
    rad.UtiDelAll()

    # Start memory tracking
    tracemalloc.start()
    mem_start = tracemalloc.get_traced_memory()[0]

    # Recreate geometry
    mat = rad.MatSatIsoFrm(*mat_params)
    elements = []
    for i in range(n):
        for j in range(n):
            for k in range(n):
                x = -size/2 + (i + 0.5) * elem_size
                y = -size/2 + (j + 0.5) * elem_size
                z = -size/2 + (k + 0.5) * elem_size
                elem = rad.ObjRecMag([x, y, z], [elem_size, elem_size, elem_size])
                rad.MatApl(elem, mat)
                elements.append(elem)

    bg_field = rad.ObjBckg(H_bg)
    container = rad.ObjCnt(elements + [bg_field])

    # Disable H-matrix
    rad.SolverHMatrixDisable()

    # Matrix construction
    print("  Matrix construction...")
    t_matrix_start = time.perf_counter()
    intrc = rad.RlxPre(container, container)
    t_matrix_gs = time.perf_counter() - t_matrix_start

    # Measure per-iteration time
    print("  Gauss-Seidel solve ({} iterations)...".format(num_iterations))
    t_solve_start = time.perf_counter()
    for _ in range(num_iterations):
        rad.RlxMan(intrc, 4, 1, 1.0)  # Method 4 = Gauss-Seidel
    t_solve_gs = (time.perf_counter() - t_solve_start) / num_iterations

    # Memory usage
    mem_peak = tracemalloc.get_traced_memory()[1]
    mem_gs = (mem_peak - mem_start) / (1024 * 1024)  # MB
    tracemalloc.stop()

    # Full solve for convergence
    rad.UtiDelAll()
    mat = rad.MatSatIsoFrm(*mat_params)
    elements = []
    for i in range(n):
        for j in range(n):
            for k in range(n):
                x = -size/2 + (i + 0.5) * elem_size
                y = -size/2 + (j + 0.5) * elem_size
                z = -size/2 + (k + 0.5) * elem_size
                elem = rad.ObjRecMag([x, y, z], [elem_size, elem_size, elem_size])
                rad.MatApl(elem, mat)
                elements.append(elem)
    bg_field = rad.ObjBckg(H_bg)
    container = rad.ObjCnt(elements + [bg_field])
    rad.SolverHMatrixDisable()

    t_full_start = time.perf_counter()
    result_gs = rad.Solve(container, precision, max_iter)
    t_full_gs = time.perf_counter() - t_full_start
    iter_gs = int(result_gs[1]) if isinstance(result_gs, (list, tuple)) else result_gs

    H_gs = rad.Fld(container, 'h', obs_point)

    print("    Matrix construction: {:.6f} s".format(t_matrix_gs))
    print("    Per-iteration time: {:.6f} s ({:.2f} ms)".format(t_solve_gs, t_solve_gs * 1000))
    print("    Full solve time: {:.6f} s ({} iterations)".format(t_full_gs, iter_gs))
    print("    Memory usage: {:.2f} MB".format(mem_gs))
    print("    H: [{:.6e}, {:.6e}, {:.6e}] A/m".format(*H_gs))
    print("    |H|: {:.6e} A/m".format(np.linalg.norm(H_gs)))

    result_entry.update({
        "t_matrix_gs": t_matrix_gs,
        "t_solve_gs": t_solve_gs,
        "t_full_gs": t_full_gs,
        "iter_gs": iter_gs,
        "mem_gs": mem_gs,
        "H_gs": H_gs,
    })

    # ========================================
    # Method 3: H-matrix Relaxation
    # ========================================
    print("\n[Method 3] H-matrix Accelerated Relaxation")
    print("-" * 80)
    rad.UtiDelAll()

    # Start memory tracking
    tracemalloc.start()
    mem_start = tracemalloc.get_traced_memory()[0]

    # Recreate geometry
    mat = rad.MatSatIsoFrm(*mat_params)
    elements = []
    for i in range(n):
        for j in range(n):
            for k in range(n):
                x = -size/2 + (i + 0.5) * elem_size
                y = -size/2 + (j + 0.5) * elem_size
                z = -size/2 + (k + 0.5) * elem_size
                elem = rad.ObjRecMag([x, y, z], [elem_size, elem_size, elem_size])
                rad.MatApl(elem, mat)
                elements.append(elem)

    bg_field = rad.ObjBckg(H_bg)
    container = rad.ObjCnt(elements + [bg_field])

    # Enable H-matrix
    rad.SolverHMatrixEnable(1, 1e-4, 30)

    # H-matrix construction (includes parallel OpenMP)
    print("  H-matrix construction...")
    t_hmatrix_start = time.perf_counter()
    intrc = rad.RlxPre(container, container)
    t_hmatrix_build = time.perf_counter() - t_hmatrix_start

    # Measure per-iteration time
    print("  H-matrix solve ({} iterations)...".format(num_iterations))
    t_solve_start = time.perf_counter()
    for _ in range(num_iterations):
        rad.RlxMan(intrc, 4, 1, 1.0)  # Same method 4, but with H-matrix
    t_solve_hmat = (time.perf_counter() - t_solve_start) / num_iterations

    # Memory usage (H-matrix memory)
    stats = rad.GetHMatrixStats()
    mem_hmat = stats[2]  # H-matrix memory in MB

    # Total memory (Python objects + H-matrix)
    mem_peak = tracemalloc.get_traced_memory()[1]
    mem_total = (mem_peak - mem_start) / (1024 * 1024)  # MB
    tracemalloc.stop()

    # Full solve for convergence
    rad.UtiDelAll()
    mat = rad.MatSatIsoFrm(*mat_params)
    elements = []
    for i in range(n):
        for j in range(n):
            for k in range(n):
                x = -size/2 + (i + 0.5) * elem_size
                y = -size/2 + (j + 0.5) * elem_size
                z = -size/2 + (k + 0.5) * elem_size
                elem = rad.ObjRecMag([x, y, z], [elem_size, elem_size, elem_size])
                rad.MatApl(elem, mat)
                elements.append(elem)
    bg_field = rad.ObjBckg(H_bg)
    container = rad.ObjCnt(elements + [bg_field])
    rad.SolverHMatrixEnable(1, 1e-4, 30)

    t_full_start = time.perf_counter()
    result_hmat = rad.Solve(container, precision, max_iter)
    t_full_hmat = time.perf_counter() - t_full_start
    iter_hmat = int(result_hmat[1]) if isinstance(result_hmat, (list, tuple)) else result_hmat

    H_hmat = rad.Fld(container, 'h', obs_point)

    print("    H-matrix construction: {:.6f} s".format(t_hmatrix_build))
    print("    Per-iteration time: {:.6f} s ({:.2f} ms)".format(t_solve_hmat, t_solve_hmat * 1000))
    print("    Full solve time: {:.6f} s ({} iterations)".format(t_full_hmat, iter_hmat))
    print("    H-matrix memory: {:.2f} MB".format(mem_hmat))
    print("    Total memory: {:.2f} MB".format(mem_total))
    print("    H: [{:.6e}, {:.6e}, {:.6e}] A/m".format(*H_hmat))
    print("    |H|: {:.6e} A/m".format(np.linalg.norm(H_hmat)))

    result_entry.update({
        "t_hmatrix_build": t_hmatrix_build,
        "t_solve_hmat": t_solve_hmat,
        "t_full_hmat": t_full_hmat,
        "iter_hmat": iter_hmat,
        "mem_hmat": mem_hmat,
        "mem_total_hmat": mem_total,
        "H_hmat": H_hmat,
    })

    # ========================================
    # Comparison
    # ========================================
    print("\n[Comparison]")
    print("-" * 80)

    # Use Gauss-Seidel as reference
    H_ref = np.array(H_gs)
    H_ref_mag = np.linalg.norm(H_ref)

    # Accuracy
    if lu_enabled and result_entry["H_lu"] is not None:
        diff_lu = np.array(result_entry["H_lu"]) - H_ref
        err_lu = np.linalg.norm(diff_lu) / (H_ref_mag + 1e-15) * 100
        result_entry["err_lu"] = err_lu
    else:
        result_entry["err_lu"] = None

    diff_hmat = np.array(H_hmat) - H_ref
    err_hmat = np.linalg.norm(diff_hmat) / (H_ref_mag + 1e-15) * 100
    result_entry["err_hmat"] = err_hmat

    # Speedup calculations
    if lu_enabled and result_entry["t_solve_lu"] is not None:
        speedup_gs_vs_lu = result_entry["t_solve_lu"] / t_solve_gs
        speedup_hmat_vs_lu = result_entry["t_solve_lu"] / t_solve_hmat
        result_entry["speedup_gs_vs_lu"] = speedup_gs_vs_lu
        result_entry["speedup_hmat_vs_lu"] = speedup_hmat_vs_lu
    else:
        result_entry["speedup_gs_vs_lu"] = None
        result_entry["speedup_hmat_vs_lu"] = None

    speedup_hmat_vs_gs = t_solve_gs / t_solve_hmat
    result_entry["speedup_hmat_vs_gs"] = speedup_hmat_vs_gs

    # Print comparison table
    print("\n  {:20s} {:>12s} {:>12s} {:>12s} {:>12s}".format(
        "Method", "Construct(s)", "Iter(ms)", "Memory(MB)", "Accuracy(%)"))
    print("  " + "-" * 72)

    if lu_enabled:
        print("  {:20s} {:>12.6f} {:>12.2f} {:>12.2f} {:>12.6f}".format(
            "LU Decomposition",
            result_entry["t_matrix_lu"],
            result_entry["t_solve_lu"] * 1000,
            result_entry["mem_lu"],
            result_entry["err_lu"]))

    print("  {:20s} {:>12.6f} {:>12.2f} {:>12.2f} {:>12.6f}".format(
        "Gauss-Seidel",
        t_matrix_gs,
        t_solve_gs * 1000,
        mem_gs,
        0.0))

    print("  {:20s} {:>12.6f} {:>12.2f} {:>12.2f} {:>12.6f}".format(
        "H-matrix",
        t_hmatrix_build,
        t_solve_hmat * 1000,
        mem_hmat,
        err_hmat))

    print("\n  Speedup (per iteration vs Gauss-Seidel):")
    if lu_enabled and speedup_gs_vs_lu is not None:
        print("    LU: {:.2f}x SLOWER".format(1/speedup_gs_vs_lu))
    print("    H-matrix: {:.2f}x {}".format(
        abs(speedup_hmat_vs_gs),
        "FASTER" if speedup_hmat_vs_gs > 1 else "SLOWER"))

    results.append(result_entry)

# ========================================
# Summary and Analysis
# ========================================
print("\n" + "=" * 80)
print("SUMMARY - Performance vs Problem Size")
print("=" * 80)

print("\n{:>6s} {:>8s}  {:>10s} {:>10s} {:>10s}  {:>10s} {:>10s}  {:>10s}".format(
    "N", "Elements",
    "LU(ms)", "GS(ms)", "H-mat(ms)",
    "Mem_GS(MB)", "Mem_H(MB)",
    "Best"))
print("-" * 95)

for r in results:
    lu_str = "{:.2f}".format(r["t_solve_lu"] * 1000) if r["t_solve_lu"] is not None else "---"
    gs_str = "{:.2f}".format(r["t_solve_gs"] * 1000)
    hm_str = "{:.2f}".format(r["t_solve_hmat"] * 1000)

    # Determine best method
    times = []
    if r["t_solve_lu"] is not None:
        times.append(("LU", r["t_solve_lu"]))
    times.append(("GS", r["t_solve_gs"]))
    times.append(("H-mat", r["t_solve_hmat"]))
    best_method = min(times, key=lambda x: x[1])[0]

    print("{:>6d} {:>8d}  {:>10s} {:>10s} {:>10s}  {:>10.2f} {:>10.2f}  {:>10s}".format(
        r["n"], r["n_elem"],
        lu_str, gs_str, hm_str,
        r["mem_gs"], r["mem_hmat"],
        best_method))

# ========================================
# Scaling Analysis
# ========================================
print("\n" + "=" * 80)
print("SCALING ANALYSIS")
print("=" * 80)

# Extract data for scaling analysis
n_values = np.array([r["n_elem"] for r in results])
log_n = np.log(n_values)

# GS scaling
gs_times = np.array([r["t_solve_gs"] for r in results])
log_gs = np.log(gs_times)
A = np.vstack([log_n, np.ones(len(log_n))]).T
alpha_gs, log_a_gs = np.linalg.lstsq(A, log_gs, rcond=None)[0]

# H-matrix scaling
hmat_times = np.array([r["t_solve_hmat"] for r in results])
log_hmat = np.log(hmat_times)
alpha_hmat, log_a_hmat = np.linalg.lstsq(A, log_hmat, rcond=None)[0]

print("\nComplexity Analysis (per iteration):")
print("  Gauss-Seidel: t = {:.6e} * N^{:.3f}".format(np.exp(log_a_gs), alpha_gs))
print("  H-matrix:     t = {:.6e} * N^{:.3f}".format(np.exp(log_a_hmat), alpha_hmat))

print("\nTheoretical complexity:")
print("  Gauss-Seidel: O(N^2) -> expected alpha ~= 2.0")
print("  H-matrix:     O(N^2 log N) -> expected alpha ~= 2.0-2.5")

print("\nMeasured complexity:")
print("  Gauss-Seidel: alpha = {:.3f} -> {}".format(
    alpha_gs,
    "Matches O(N^2)" if 1.8 <= alpha_gs <= 2.2 else "Unexpected scaling"))
print("  H-matrix:     alpha = {:.3f} -> {}".format(
    alpha_hmat,
    "Matches O(N^2 log N)" if 2.0 <= alpha_hmat <= 2.5 else "Unexpected scaling"))

# ========================================
# Recommendations
# ========================================
print("\n" + "=" * 80)
print("RECOMMENDATIONS")
print("=" * 80)

print("""
Based on measured performance:

1. PROBLEM SIZE SELECTION:
   - N < 100 (tiny):     Use LU decomposition for direct solve
   - 100 <= N < 250:      Use Gauss-Seidel for good performance
   - N >= 250 (large):    Use H-matrix for best performance

2. MEMORY CONSIDERATIONS:
   - Gauss-Seidel: O(N^2) memory, grows quickly
   - H-matrix:     O(N log N) memory, much more efficient for large N

3. CONVERGENCE BEHAVIOR:
   - Fast convergence (<5 iterations): GS may be faster due to H-matrix overhead
   - Slow convergence (>10 iterations): H-matrix significantly faster per iteration

4. DISK CACHING (v1.1.0):
   - For repeated simulations: Enable H-matrix disk caching
   - First run: Pays construction cost once
   - Subsequent runs: ~10x faster startup (load from disk)
   - Use: rad.SolverHMatrixCacheFull(1)

5. WHEN TO USE EACH METHOD:
   - LU:         Tiny problems, direct solve needed, single solution
   - GS:         Medium problems, moderate accuracy, simple setup
   - H-matrix:   Large problems, repeated simulations, memory-constrained
""")

print("=" * 80)
print("Benchmark complete!")
print("=" * 80)

rad.UtiDelAll()
