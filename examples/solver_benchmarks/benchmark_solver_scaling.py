#!/usr/bin/env python3
"""
Benchmark: Solver method comparison
Compare LU (Method 0), BiCGSTAB (Method 1), and HACApK (Method 2)
"""

import sys
import os
import numpy as np
from time import perf_counter
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../build/Release'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

import radia as rad

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


def create_geometry(nx, ny, nz, cube_size, mat):
	"""Create hexahedral mesh geometry."""
	elem_size = cube_size / nx
	elements = []
	for i in range(nx):
		for j in range(ny):
			for k in range(nz):
				x = (i - nx/2 + 0.5) * elem_size
				y = (j - ny/2 + 0.5) * elem_size
				z = (k - nz/2 + 0.5) * elem_size
				vertices = hex_vertices(x, y, z, elem_size, elem_size, elem_size)
				elem = rad.ObjHexahedron(vertices, [0, 0, 0.1])
				rad.MatApl(elem, mat)
				elements.append(elem)
	return rad.ObjCnt(elements)


print("=" * 70)
print("Solver Method Comparison Benchmark")
print("=" * 70)

# Test cases
test_cases = [
	(3, 3, 3),    # 27
	(4, 4, 4),    # 64
	(5, 5, 5),    # 125
	(6, 6, 6),    # 216
	(7, 7, 7),    # 343
	(8, 8, 8),    # 512
]

print("\nConfiguration:")
print("  Material: Nonlinear (soft iron)")
print("  Geometry: 100x100x100 mm cube")
print("  Solver methods:")
print("    - Method 0: LU (direct)")
print("    - Method 1: BiCGSTAB (iterative)")
print("    - Method 2: HACApK (H-matrix accelerated)")

# Nonlinear material (soft iron)
MH_data = [[0, 0], [200, 0.7], [600, 1.2], [1200, 1.4], [2000, 1.5],
           [3500, 1.54], [6000, 1.56], [12000, 1.57]]

results = []

for nx, ny, nz in test_cases:
	n_elem = nx * ny * nz
	cube_size = 100.0 * mm

	print(f"\nN = {nx}x{ny}x{nz} = {n_elem:4d} elements")
	print("-" * 70)

	result = {'n': n_elem}

	# Method 0: LU
	mat = rad.MatSatIsoTab(MH_data)
	grp = create_geometry(nx, ny, nz, cube_size, mat)
	t_start = perf_counter()
	rad.Solve(grp, 0.0001, 100, 0)
	t_lu = perf_counter() - t_start
	result['t_lu'] = t_lu
	print(f"  Method 0 (LU):      {t_lu*1000:8.2f} ms")
	rad.UtiDelAll()

	# Method 1: BiCGSTAB
	mat = rad.MatSatIsoTab(MH_data)
	grp = create_geometry(nx, ny, nz, cube_size, mat)
	t_start = perf_counter()
	rad.Solve(grp, 0.0001, 100, 1)
	t_bicg = perf_counter() - t_start
	result['t_bicg'] = t_bicg
	print(f"  Method 1 (BiCGSTAB): {t_bicg*1000:8.2f} ms")
	rad.UtiDelAll()

	# Method 2: HACApK (if available)
	try:
		mat = rad.MatSatIsoTab(MH_data)
		grp = create_geometry(nx, ny, nz, cube_size, mat)
		rad.SolverConfig(hacapk_eps=1e-4, hacapk_leaf=10, hacapk_eta=2.0)
		t_start = perf_counter()
		rad.Solve(grp, 0.0001, 100, 2)
		t_hacapk = perf_counter() - t_start
		result['t_hacapk'] = t_hacapk
		print(f"  Method 2 (HACApK):  {t_hacapk*1000:8.2f} ms")
		rad.UtiDelAll()
	except Exception as e:
		result['t_hacapk'] = float('inf')
		print(f"  Method 2 (HACApK):  N/A ({e})")
		rad.UtiDelAll()

	# Speedup ratios
	print(f"  Speedup vs LU:  BiCGSTAB {t_lu/t_bicg:.2f}x, HACApK {t_lu/t_hacapk:.2f}x")

	results.append(result)

#============================================================================
# SCALING ANALYSIS
#============================================================================
print("\n" + "=" * 70)
print("Scaling Analysis")
print("=" * 70)

n_values = np.array([r['n'] for r in results])
log_n = np.log(n_values)
A = np.vstack([log_n, np.ones(len(log_n))]).T

# Power law fits
t_lu_values = np.array([r['t_lu'] for r in results])
t_bicg_values = np.array([r['t_bicg'] for r in results])
t_hacapk_values = np.array([r['t_hacapk'] for r in results])

alpha_lu, _ = np.linalg.lstsq(A, np.log(t_lu_values), rcond=None)[0]
alpha_bicg, _ = np.linalg.lstsq(A, np.log(t_bicg_values), rcond=None)[0]
alpha_hacapk, _ = np.linalg.lstsq(A, np.log(t_hacapk_values), rcond=None)[0]

print(f"\nPower law fits: t = a * N^alpha")
print(f"  Method 0 (LU):       O(N^{alpha_lu:.2f})")
print(f"  Method 1 (BiCGSTAB): O(N^{alpha_bicg:.2f})")
print(f"  Method 2 (HACApK):   O(N^{alpha_hacapk:.2f})")

# Summary table
print(f"\n{'N':>6}  {'LU (ms)':>10}  {'BiCG (ms)':>10}  {'HACApK (ms)':>12}  {'Best':>8}")
print("-" * 60)

for r in results:
	n = r['n']
	t_lu = r['t_lu'] * 1000
	t_bicg = r['t_bicg'] * 1000
	t_hacapk = r['t_hacapk'] * 1000

	best = min(t_lu, t_bicg, t_hacapk)
	if best == t_lu:
		best_name = "LU"
	elif best == t_bicg:
		best_name = "BiCGSTAB"
	else:
		best_name = "HACApK"

	print(f"{n:>6}  {t_lu:>10.2f}  {t_bicg:>10.2f}  {t_hacapk:>12.2f}  {best_name:>8}")

print("\n" + "=" * 70)
print("Recommendations")
print("=" * 70)

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

print("=" * 70)
