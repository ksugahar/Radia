#!/usr/bin/env python3
"""
Benchmark: Matrix construction time only
Verify O(N^2) scaling
"""

import sys
import os
import numpy as np
from time import perf_counter
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../build/Release'))

import radia as rad

print("=" * 70)
print("Matrix Construction Benchmark")
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
	(12, 12, 12), # 1728
	(15, 15, 15), # 3375
]

print("\nConfiguration:")
print("  Material: Nonlinear (soft iron)")
print("  Geometry: 100x100x100 mm cube")
print("  Operation: rad.RlxPre() - Build interaction matrix")
print("  Expected: O(N^2)")

# Nonlinear material (soft iron)
MH_data = [[0, 0], [200, 0.7], [600, 1.2], [1200, 1.4], [2000, 1.5],
           [3500, 1.54], [6000, 1.56], [12000, 1.57]]

results = []

for nx, ny, nz in test_cases:
	n_elem = nx * ny * nz

	print(f"\nN = {nx}x{ny}x{nz} = {n_elem:4d} elements ... ", end='', flush=True)

	cube_size = 100.0
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

				elem = rad.ObjRecMag([x, y, z], [elem_size, elem_size, elem_size], [0, 0, 0.1])
				rad.MatApl(elem, mat)
				elements.append(elem)

	grp = rad.ObjCnt(elements)

	# Disable H-matrix
	rad.SolverHMatrixDisable()

	# Measure matrix construction time
	t_start = perf_counter()
	intrc = rad.RlxPre(grp, grp)
	t_matrix = perf_counter() - t_start

	print(f"t = {t_matrix*1000:8.2f} ms")

	results.append({
		'n': n_elem,
		't_matrix': t_matrix,
	})

	rad.UtiDelAll()

#============================================================================
# SCALING ANALYSIS
#============================================================================
print("\n" + "=" * 70)
print("Scaling Analysis")
print("=" * 70)

n_values = np.array([r['n'] for r in results])
t_values = np.array([r['t_matrix'] for r in results])

# Power law fit
log_n = np.log(n_values)
log_t = np.log(t_values)
A = np.vstack([log_n, np.ones(len(log_n))]).T
alpha, log_a = np.linalg.lstsq(A, log_t, rcond=None)[0]

print(f"\nPower law fit: t = a * N^alpha")
print(f"  alpha = {alpha:.4f}")
print(f"  a = {np.exp(log_a):.6e}")

# Detailed table
print(f"\n{'N':>6}  {'Time (ms)':>12}  {'t/N':>12}  {'t/N^2':>12}")
print("-" * 50)

for r in results:
	n = r['n']
	t_ms = r['t_matrix'] * 1000
	t_per_n = t_ms / n
	t_per_n2 = t_ms / (n * n)

	print(f"{n:>6}  {t_ms:>12.2f}  {t_per_n:>12.6f}  {t_per_n2:>12.6f}")

# Check if t/N^2 is constant
t_per_n2_values = [r['t_matrix'] * 1000 / (r['n']**2) for r in results[3:]]
mean_n2 = np.mean(t_per_n2_values)
std_n2 = np.std(t_per_n2_values)
cv = std_n2 / mean_n2 if mean_n2 > 0 else 0

print(f"\nt/N^2 statistics (N >= {results[3]['n']}):")
print(f"  Mean: {mean_n2:.6f}")
print(f"  Std:  {std_n2:.6f}")
print(f"  CV:   {cv:.3f}")

print("\n" + "=" * 70)
print("Interpretation")
print("=" * 70)

if 1.7 <= alpha <= 2.3:
	print(f"  alpha = {alpha:.3f} -> O(N^2) CONFIRMED")
	if cv < 0.2:
		print(f"  t/N^2 approximately constant (CV={cv:.3f})")
		print("  -> Consistent O(N^2) behavior across range")
	else:
		print(f"  t/N^2 has variation (CV={cv:.3f})")
		print("  -> May indicate transition region or cache effects")
elif alpha < 1.7:
	print(f"  alpha = {alpha:.3f} -> LESS than O(N^2)")
	print("  -> Better than expected, check implementation")
else:
	print(f"  alpha = {alpha:.3f} -> MORE than O(N^2)")
	print("  -> Worse than expected, possible memory/cache issues")

print("\n" + "=" * 70)
print("Conclusion")
print("=" * 70)

if 1.7 <= alpha <= 2.3 and cv < 0.3:
	print(f"Matrix construction correctly shows O(N^2) scaling")
	print(f"  Measured exponent: {alpha:.3f}")
	print(f"  This is the expected cost for building N x N interaction matrix")
else:
	print(f"Unexpected scaling: alpha = {alpha:.3f}, CV = {cv:.3f}")
	print(f"Further investigation may be needed")

print("=" * 70)
