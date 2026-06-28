#!/usr/bin/env python3
"""
Benchmark: Solver total time scaling
Measures rad.Solve() performance across problem sizes
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

print("=" * 70)
print("Solver Scaling Benchmark")
print("=" * 70)

# Test cases
test_cases = [
	(2, 2, 2),    # 8
	(3, 3, 3),    # 27
	(4, 4, 4),    # 64
	(5, 5, 5),    # 125
	(6, 6, 6),    # 216
	(7, 7, 7),    # 343
	(8, 8, 8),    # 512
	(10, 10, 10), # 1000
]

print("\nConfiguration:")
print("  Material: Nonlinear (soft iron)")
print("  Geometry: 100x100x100 mm cube")
print("  Solver: rad.Solve() with Method 0 (LU)")

# Nonlinear material (soft iron)
MH_data = [[0, 0], [200, 0.7], [600, 1.2], [1200, 1.4], [2000, 1.5],
           [3500, 1.54], [6000, 1.56], [12000, 1.57]]

results = []

for nx, ny, nz in test_cases:
	n_elem = nx * ny * nz

	print(f"\nN = {nx}x{ny}x{nz} = {n_elem:4d} elements ... ", end='', flush=True)

	cube_size = 100.0 * mm
	elem_size = cube_size / nx

	mat = rad.MatSatIsoTab(MH_data)

	# Build geometry
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

	grp = rad.ObjCnt(elements)

	# Measure solve time with LU (Method 0)
	t_start = perf_counter()
	rad.Solve(grp, 0.0001, 100, 0)  # Method 0 = LU
	t_solve = perf_counter() - t_start

	print(f"t = {t_solve*1000:8.2f} ms")

	results.append({
		'n': n_elem,
		't_solve': t_solve,
	})

	rad.UtiDelAll()

#============================================================================
# SCALING ANALYSIS
#============================================================================
print("\n" + "=" * 70)
print("Scaling Analysis")
print("=" * 70)

n_values = np.array([r['n'] for r in results])
t_values = np.array([r['t_solve'] for r in results])

# Power law fit
log_n = np.log(n_values)
log_t = np.log(t_values)
A = np.vstack([log_n, np.ones(len(log_n))]).T
alpha, log_a = np.linalg.lstsq(A, log_t, rcond=None)[0]

print(f"\nPower law fit: t = a * N^alpha")
print(f"  alpha = {alpha:.4f}")
print(f"  a = {np.exp(log_a):.6e}")

# Detailed table
print(f"\n{'N':>6}  {'Time (ms)':>12}  {'t/N^2':>12}  {'t/N^3':>12}")
print("-" * 50)

for r in results:
	n = r['n']
	t_ms = r['t_solve'] * 1000
	t_per_n2 = t_ms / (n * n)
	t_per_n3 = t_ms / (n * n * n)

	print(f"{n:>6}  {t_ms:>12.2f}  {t_per_n2:>12.6f}  {t_per_n3:>12.9f}")

print("\n" + "=" * 70)
print("Interpretation")
print("=" * 70)

print(f"  Measured exponent: alpha = {alpha:.3f}")
if alpha < 2.5:
	print("  -> Scaling better than O(N^3)")
elif alpha < 3.5:
	print("  -> Scaling approximately O(N^3) (LU decomposition)")
else:
	print("  -> Scaling worse than expected")

print("\n" + "=" * 70)
print("Solver Method Comparison")
print("=" * 70)

print("""
Solver methods available in rad.Solve(obj, prec, maxiter, method):

  Method 0 (LU):      Direct solver, O(N^3), best for small problems
  Method 1 (BiCGSTAB): Iterative, O(N^2) per iteration, general purpose
  Method 2 (HACApK):  H-matrix accelerated, O(N log N) memory, large problems

Recommended:
  - N < 500:   Method 0 (LU)
  - N < 2000:  Method 1 (BiCGSTAB)
  - N > 2000:  Method 2 (HACApK)
""")

print("=" * 70)
