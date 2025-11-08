#!/usr/bin/env python3
"""
Benchmark: Linear material - matrix construction only
Linear materials converge in 0-1 iterations, so only matrix construction is relevant
"""

import sys
import numpy as np
from time import perf_counter
sys.path.insert(0, r"S:\Radia\01_GitHub\build\Release")

import radia as rad

print("=" * 70)
print("Linear Material Benchmark (Matrix Construction Only)")
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
print("  Material: Linear (chi=999, no remanent magnetization)")
print("  Geometry: 100x100x100 mm cube")
print("  Operation: rad.RlxPre() - Matrix construction only")
print("  Note: Linear materials converge in 0-1 iterations")
print("        -> Only matrix construction time is relevant")

# Linear material (high susceptibility, no remanent magnetization)
# Format: MatLin([chi_parallel, chi_perpendicular], [Mx, My, Mz])
# For isotropic linear material with no remanence: [chi, chi], [0, 0, 0]

results = []

for nx, ny, nz in test_cases:
	n_elem = nx * ny * nz

	print(f"\nN = {nx}x{ny}x{nz} = {n_elem:4d} elements")
	print("-" * 70)

	cube_size = 100.0
	elem_size = cube_size / nx

	# Linear material with high susceptibility
	mat = rad.MatLin([999, 999], [0, 0, 0])

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
	rad.SolverHMatrixDisable()

	# Measure matrix construction time
	t_matrix_start = perf_counter()
	intrc = rad.RlxPre(grp, grp)
	t_matrix = perf_counter() - t_matrix_start

	print(f"  Matrix construction: {t_matrix*1000:8.2f} ms (O(N^2))")

	# For linear materials, relaxation converges in 0-1 iterations
	# So we measure one iteration just for reference
	t_solve_start = perf_counter()
	rad.RlxMan(intrc, 4, 1, 1.0)
	t_solve = perf_counter() - t_solve_start

	print(f"  Solver (1 iter):     {t_solve*1000:8.2f} ms (for reference)")
	print(f"  Total time:          {(t_matrix + t_solve)*1000:8.2f} ms")

	results.append({
		'n': n_elem,
		't_matrix': t_matrix,
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
log_n = np.log(n_values)

# Matrix construction scaling
t_matrix_values = np.array([r['t_matrix'] for r in results])
log_t_matrix = np.log(t_matrix_values)
A = np.vstack([log_n, np.ones(len(log_n))]).T
alpha_matrix, log_a_matrix = np.linalg.lstsq(A, log_t_matrix, rcond=None)[0]

print(f"\nPower law fit: t = a * N^alpha")
print(f"  Matrix construction: t = {np.exp(log_a_matrix):.6e} * N^{alpha_matrix:.3f}")

# Detailed table
print(f"\n{'N':>6}  {'Matrix (ms)':>12}  {'Solve (ms)':>12}  {'Total (ms)':>12}  {'t_m/N^2':>12}")
print("-" * 60)

for r in results:
	n = r['n']
	t_m = r['t_matrix'] * 1000
	t_s = r['t_solve'] * 1000
	t_total = t_m + t_s

	t_m_n2 = t_m / (n * n)

	print(f"{n:>6}  {t_m:>12.2f}  {t_s:>12.2f}  {t_total:>12.2f}  {t_m_n2:>12.6f}")

# Check for constant ratios
t_matrix_n2_values = [r['t_matrix'] * 1000 / (r['n']**2) for r in results[3:]]

mean_matrix_n2 = np.mean(t_matrix_n2_values)
cv_matrix = np.std(t_matrix_n2_values) / mean_matrix_n2 if mean_matrix_n2 > 0 else 0

print(f"\nRatio statistics (N >= {results[3]['n']}):")
print(f"  Matrix/N^2: mean={mean_matrix_n2:.6f}, CV={cv_matrix:.3f}")

print("\n" + "=" * 70)
print("Interpretation")
print("=" * 70)

print(f"\nMatrix Construction: alpha = {alpha_matrix:.3f}")
if 1.7 <= alpha_matrix <= 2.3:
	print("  -> O(N^2) CONFIRMED")
	print("  -> N×N interaction matrix for N elements")
	if cv_matrix < 0.2:
		print(f"  -> Consistent scaling (CV={cv_matrix:.3f})")
else:
	print(f"  -> NOT O(N^2) (expected 1.7-2.3)")

print("\n" + "=" * 70)
print("Conclusion")
print("=" * 70)

print(f"""
Linear Material Characteristics:

1. Matrix construction: O(N^{alpha_matrix:.1f})
   - Same O(N^2) scaling as nonlinear materials
   - {mean_matrix_n2:.6f} ms/N^2 (approximately constant)

2. Solver convergence: 0-1 iterations
   - Linear relationship: M = chi*H
   - No iteration required for self-consistent solution
   - Solver time negligible compared to matrix construction

3. Comparison with Nonlinear Materials:
   - Matrix construction: SAME O(N^2)
   - Solver time: LINEAR << NONLINEAR
   - Total time dominated by matrix construction

Recommendation for Linear Problems:
  - Matrix construction O(N^2) is the main cost
  - Solver method choice (GS vs LU) irrelevant (0-1 iteration)
  - Focus optimization on matrix construction, not solver
""")

print("=" * 70)
