#!/usr/bin/env python
"""
Solver Methods Benchmark: Direct vs Relaxation vs Relaxation+H-matrix

Compares three solver methods for magnetic field problems:
1. Direct calculation (no solver, uses initial magnetization)
2. Standard relaxation solver
3. H-matrix accelerated relaxation solver

Problem: Nonlinear magnetic material with applied background field

Performance metrics:
- Computation time
- Field accuracy (comparing all methods)
"""

import sys
import os
import time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../build/Release'))

import radia as rad
import numpy as np

print("=" * 80)
print("Solver Methods Benchmark")
print("Direct / Relaxation / Relaxation+H-matrix Comparison")
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
# Note: In Radia, ObjBckg() uses Tesla, and field values are in Radia's internal units
# where H and B are numerically equal. Use value ~1.0 for moderate field.
H_bg = [1.0, 0, 0]  # 1.0 T background field in X direction

# Observation point
obs_point = [0, 0, 50]  # 50mm above center

print("\nProblem Setup:")
print("  Material: Nonlinear (Steel37)")
print("  Background field: [{}, {}, {}] T".format(*H_bg))
print("  Observation point: [{}, {}, {}] mm".format(*obs_point))
print("  Solver precision: {}".format(precision))
print("  Max iterations: {}".format(max_iter))

results = []

for test in test_cases:
	n = test["n"]
	desc = test["desc"]

	print("\n" + "=" * 80)
	print("Test Case: {}".format(desc))
	print("=" * 80)

	size = 20.0  # mm (total cube size)
	elem_size = size / n

	# ========================================
	# Method 1: Direct Calculation
	# ========================================
	print("\n[Method 1] Direct Calculation (no solver)")
	print("-" * 80)
	rad.UtiDelAll()

	# Create material
	mat = rad.MatSatIsoFrm([1596.3, 1.1488], [133.11, 0.4268], [18.713, 0.4759])

	# Create geometry with initial magnetization
	elements = []
	for i in range(n):
		for j in range(n):
			for k in range(n):
				x = -size/2 + (i + 0.5) * elem_size
				y = -size/2 + (j + 0.5) * elem_size
				z = -size/2 + (k + 0.5) * elem_size
				# No initial magnetization (material will respond to background field)
				elem = rad.ObjRecMag([x, y, z], [elem_size, elem_size, elem_size])
				rad.MatApl(elem, mat)
				elements.append(elem)

	# Add background field source to container
	bg_field = rad.ObjBckg(H_bg)
	container = rad.ObjCnt(elements + [bg_field])

	print("  Elements: {} ({}×{}×{})".format(len(elements), n, n, n))

	# Direct field calculation (using initial magnetization only)
	t0 = time.perf_counter()
	H_direct = rad.Fld(container, 'h', obs_point)
	t_direct = time.perf_counter() - t0

	print("  Time: {:.6f} s".format(t_direct))
	print("  H: [{:.6e}, {:.6e}, {:.6e}] A/m".format(*H_direct))
	H_mag_direct = np.linalg.norm(H_direct)
	print("  |H|: {:.6e} A/m".format(H_mag_direct))

	# ========================================
	# Method 2: Standard Relaxation
	# ========================================
	print("\n[Method 2] Standard Relaxation (no H-matrix)")
	print("-" * 80)
	rad.UtiDelAll()

	# Recreate geometry
	mat = rad.MatSatIsoFrm([1596.3, 1.1488], [133.11, 0.4268], [18.713, 0.4759])

	elements = []
	for i in range(n):
		for j in range(n):
			for k in range(n):
				x = -size/2 + (i + 0.5) * elem_size
				y = -size/2 + (j + 0.5) * elem_size
				z = -size/2 + (k + 0.5) * elem_size
				# No initial magnetization - pure soft magnetic material
				elem = rad.ObjRecMag([x, y, z], [elem_size, elem_size, elem_size])
				rad.MatApl(elem, mat)
				elements.append(elem)

	# Add background field source to container
	bg_field = rad.ObjBckg(H_bg)
	container = rad.ObjCnt(elements + [bg_field])

	# Disable H-matrix for standard relaxation
	rad.SolverHMatrixDisable()

	# Run relaxation solver
	t0 = time.perf_counter()
	result_relax = rad.Solve(container, precision, max_iter)
	t_relax = time.perf_counter() - t0

	H_relax = rad.Fld(container, 'h', obs_point)

	print("  Time: {:.6f} s".format(t_relax))
	print("  H: [{:.6e}, {:.6e}, {:.6e}] A/m".format(*H_relax))
	H_mag_relax = np.linalg.norm(H_relax)
	print("  |H|: {:.6e} A/m".format(H_mag_relax))
	print("  Convergence: {}".format(result_relax))

	# ========================================
	# Method 3: H-matrix Relaxation
	# ========================================
	print("\n[Method 3] H-matrix Accelerated Relaxation")
	print("-" * 80)
	rad.UtiDelAll()

	# Recreate geometry
	mat = rad.MatSatIsoFrm([1596.3, 1.1488], [133.11, 0.4268], [18.713, 0.4759])

	elements = []
	for i in range(n):
		for j in range(n):
			for k in range(n):
				x = -size/2 + (i + 0.5) * elem_size
				y = -size/2 + (j + 0.5) * elem_size
				z = -size/2 + (k + 0.5) * elem_size
				# No initial magnetization - pure soft magnetic material
				elem = rad.ObjRecMag([x, y, z], [elem_size, elem_size, elem_size])
				rad.MatApl(elem, mat)
				elements.append(elem)

	# Add background field source to container
	bg_field = rad.ObjBckg(H_bg)
	container = rad.ObjCnt(elements + [bg_field])

	# Enable H-matrix for relaxation
	rad.SolverHMatrixEnable()

	# Run relaxation solver with H-matrix
	t0 = time.perf_counter()
	result_hmat = rad.Solve(container, precision, max_iter)
	t_hmat = time.perf_counter() - t0

	H_hmat = rad.Fld(container, 'h', obs_point)

	# Get H-matrix stats
	stats = rad.GetHMatrixStats()

	print("  Time: {:.6f} s".format(t_hmat))
	print("  H: [{:.6e}, {:.6e}, {:.6e}] A/m".format(*H_hmat))
	H_mag_hmat = np.linalg.norm(H_hmat)
	print("  |H|: {:.6e} A/m".format(H_mag_hmat))
	print("  Convergence: {}".format(result_hmat))
	print("  H-matrix memory: {:.3f} MB".format(stats[2]))

	# ========================================
	# Comparison
	# ========================================
	print("\n[Comparison]")
	print("-" * 80)

	# Use standard relaxation as reference
	H_ref = np.array(H_relax)
	H_ref_mag = np.linalg.norm(H_ref)

	# Direct vs Relaxation error
	diff_direct = np.array(H_direct) - H_ref
	err_direct = np.linalg.norm(diff_direct) / (H_ref_mag + 1e-15) * 100

	# H-matrix vs Relaxation error
	diff_hmat = np.array(H_hmat) - H_ref
	err_hmat = np.linalg.norm(diff_hmat) / (H_ref_mag + 1e-15) * 100

	print("  Reference: Standard Relaxation |H| = {:.6e} A/m".format(H_ref_mag))
	print("")
	print("  Method           Time (s)    |H| (A/m)        Error (%)   Speedup")
	print("  " + "-" * 72)
	print("  Direct         {:9.6f}  {:14.6e}  {:11.4f}        -".format(
		t_direct, H_mag_direct, err_direct))
	print("  Relaxation     {:9.6f}  {:14.6e}  {:11.4f}     1.00x".format(
		t_relax, H_mag_relax, 0.0))
	print("  H-matrix       {:9.6f}  {:14.6e}  {:11.4f}  {:7.2f}x".format(
		t_hmat, H_mag_hmat, err_hmat, t_relax/t_hmat if t_hmat > 0 else 0))

	results.append({
		"n": n,
		"n_elem": n**3,
		"t_direct": t_direct,
		"t_relax": t_relax,
		"t_hmat": t_hmat,
		"err_direct": err_direct,
		"err_hmat": err_hmat,
		"speedup": t_relax/t_hmat if t_hmat > 0 else 0
	})

# ========================================
# Summary
# ========================================
print("\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)
print("\n{:>6s} {:>8s}  {:>12s} {:>12s} {:>12s}  {:>10s} {:>10s}  {:>10s}".format(
	"N", "Elements",
	"Direct (s)", "Relax (s)", "H-mat (s)",
	"Dir Err(%)", "H-mat Err(%)",
	"Speedup"))
print("-" * 80)

for r in results:
	print("{:>6d} {:>8d}  {:>12.6f} {:>12.6f} {:>12.6f}  {:>10.4f} {:>10.4f}  {:>10.2f}x".format(
		r["n"], r["n_elem"],
		r["t_direct"], r["t_relax"], r["t_hmat"],
		r["err_direct"], r["err_hmat"],
		r["speedup"]))

print("=" * 80)
print("\nNotes:")
print("  - Direct: Uses initial magnetization only (no iterations)")
print("  - Relaxation: Standard solver (reference for accuracy)")
print("  - H-matrix: Accelerated solver (best for large problems N > 100)")
print("  - Error: Relative difference from Relaxation method")
print("=" * 80)

rad.UtiDelAll()
