#!/usr/bin/env python3
"""
Comprehensive benchmark: LU decomposition vs Gauss-Seidel
Separated timing: Matrix construction vs Solver execution
"""

import sys
import numpy as np
from time import perf_counter
sys.path.insert(0, r"S:\Radia\01_GitHub\build\Release")

import radia as rad

print("=" * 80)
print("Solver Benchmark: LU Decomposition vs Gauss-Seidel")
print("=" * 80)

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
print("  Timing: Separated matrix construction and solver execution")
print("  Repetitions: 1 iteration per measurement")

# Nonlinear material (soft iron)
MH_data = [[0, 0], [200, 0.7], [600, 1.2], [1200, 1.4], [2000, 1.5],
           [3500, 1.54], [6000, 1.56], [12000, 1.57]]

results = []

for nx, ny, nz in test_cases:
	n_elem = nx * ny * nz

	print(f"\n{'='*80}")
	print(f"N = {nx}x{ny}x{nz} = {n_elem:4d} elements")
	print(f"{'='*80}")

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

	#========================================================================
	# MATRIX CONSTRUCTION (same for all methods)
	#========================================================================
	print("\n[Matrix Construction]")
	print("-" * 80)

	rad.SolverHMatrixDisable()

	t_matrix_start = perf_counter()
	intrc = rad.RlxPre(grp, grp)
	t_matrix = perf_counter() - t_matrix_start

	print(f"  Time: {t_matrix*1000:8.2f} ms")
	print(f"  Expected: O(N^2)")

	#========================================================================
	# GAUSS-SEIDEL (Method 4)
	#========================================================================
	print("\n[Gauss-Seidel Solver - Method 4]")
	print("-" * 80)

	# Run 1 iteration
	t_gs_start = perf_counter()
	rad.RlxMan(intrc, 4, 1, 1.0)
	t_gs = perf_counter() - t_gs_start

	print(f"  Time (1 iter): {t_gs*1000:8.2f} ms")
	print(f"  Expected: O(N^2)")

	#========================================================================
	# LU DECOMPOSITION (Method 5)
	#========================================================================
	print("\n[LU Decomposition Solver - Method 5]")
	print("-" * 80)

	# Rebuild geometry for LU test
	rad.UtiDelAll()
	mat2 = rad.MatSatIsoTab(MH_data)
	elements2 = []
	for i in range(nx):
		for j in range(ny):
			for k in range(nz):
				x = (i - nx/2 + 0.5) * elem_size
				y = (j - ny/2 + 0.5) * elem_size
				z = (k - nz/2 + 0.5) * elem_size

				elem = rad.ObjRecMag([x, y, z], [elem_size, elem_size, elem_size], [0, 0, 0.1])
				rad.MatApl(elem, mat2)
				elements2.append(elem)

	grp2 = rad.ObjCnt(elements2)

	# Matrix construction (should be same as above)
	t_matrix2_start = perf_counter()
	intrc2 = rad.RlxPre(grp2, grp2)
	t_matrix2 = perf_counter() - t_matrix2_start

	# Enable LU decomposition
	rad.SetRelaxSubInterval(intrc2, 0, n_elem-1, 1)

	# Run 1 iteration with Method 5 (LU decomposition)
	t_lu_start = perf_counter()
	rad.RlxMan(intrc2, 5, 1, 1.0)  # Method 5 now supported
	t_lu = perf_counter() - t_lu_start

	print(f"  Time (1 iter): {t_lu*1000:8.2f} ms")
	print(f"  Expected: O(N^3)")

	#========================================================================
	# COMPARISON
	#========================================================================
	print("\n[Comparison]")
	print("-" * 80)

	ratio = t_lu / t_gs if t_gs > 0 else 0
	print(f"  Matrix construction: {t_matrix*1000:8.2f} ms (both methods)")
	print(f"  Gauss-Seidel:        {t_gs*1000:8.2f} ms (O(N^2))")
	print(f"  LU decomposition:    {t_lu*1000:8.2f} ms (O(N^3))")
	print(f"  LU / GS ratio:       {ratio:8.2f}x")

	results.append({
		'n': n_elem,
		't_matrix': t_matrix,
		't_gs': t_gs,
		't_lu': t_lu,
	})

	rad.UtiDelAll()

#============================================================================
# SCALING ANALYSIS
#============================================================================
print("\n" + "=" * 80)
print("SCALING ANALYSIS")
print("=" * 80)

n_values = np.array([r['n'] for r in results])
log_n = np.log(n_values)

# Matrix construction scaling
t_matrix_values = np.array([r['t_matrix'] for r in results])
log_t_matrix = np.log(t_matrix_values)
A = np.vstack([log_n, np.ones(len(log_n))]).T
alpha_matrix, log_a_matrix = np.linalg.lstsq(A, log_t_matrix, rcond=None)[0]

# Gauss-Seidel scaling
t_gs_values = np.array([r['t_gs'] for r in results])
log_t_gs = np.log(t_gs_values)
alpha_gs, log_a_gs = np.linalg.lstsq(A, log_t_gs, rcond=None)[0]

# LU decomposition scaling
t_lu_values = np.array([r['t_lu'] for r in results])
log_t_lu = np.log(t_lu_values)
alpha_lu, log_a_lu = np.linalg.lstsq(A, log_t_lu, rcond=None)[0]

print("\nPower law fits: t = a * N^alpha")
print("-" * 80)
print(f"Matrix construction: t = {np.exp(log_a_matrix):.6e} * N^{alpha_matrix:.3f}")
print(f"Gauss-Seidel:        t = {np.exp(log_a_gs):.6e} * N^{alpha_gs:.3f}")
print(f"LU decomposition:    t = {np.exp(log_a_lu):.6e} * N^{alpha_lu:.3f}")

# Detailed table
print(f"\n{'N':>6}  {'Matrix':>10}  {'GS':>10}  {'LU':>10}  {'LU/GS':>8}  "
      f"{'t_m/N^2':>10}  {'t_gs/N^2':>10}  {'t_lu/N^3':>10}")
print("-" * 90)

for r in results:
	n = r['n']
	t_m = r['t_matrix'] * 1000
	t_gs = r['t_gs'] * 1000
	t_lu = r['t_lu'] * 1000
	ratio = t_lu / t_gs if t_gs > 0 else 0

	t_m_n2 = t_m / (n * n)
	t_gs_n2 = t_gs / (n * n)
	t_lu_n3 = t_lu / (n * n * n) * 1000

	print(f"{n:>6}  {t_m:>10.2f}  {t_gs:>10.2f}  {t_lu:>10.2f}  {ratio:>8.2f}  "
	      f"{t_m_n2:>10.6f}  {t_gs_n2:>10.6f}  {t_lu_n3:>10.6f}")

# Check for constant ratios
t_matrix_n2 = [r['t_matrix'] * 1000 / (r['n']**2) for r in results[3:]]
t_gs_n2 = [r['t_gs'] * 1000 / (r['n']**2) for r in results[3:]]
t_lu_n3 = [r['t_lu'] * 1000 / (r['n']**3) * 1000 for r in results[3:]]

print(f"\nRatio statistics (N >= {results[3]['n']}):")
print(f"  Matrix/N^2: mean={np.mean(t_matrix_n2):.6f}, CV={np.std(t_matrix_n2)/np.mean(t_matrix_n2):.3f}")
print(f"  GS/N^2:     mean={np.mean(t_gs_n2):.6f}, CV={np.std(t_gs_n2)/np.mean(t_gs_n2):.3f}")
print(f"  LU/N^3:     mean={np.mean(t_lu_n3):.6f}, CV={np.std(t_lu_n3)/np.mean(t_lu_n3):.3f}")

print("\n" + "=" * 80)
print("INTERPRETATION")
print("=" * 80)

print(f"\nMatrix Construction: alpha = {alpha_matrix:.3f}")
if 1.7 <= alpha_matrix <= 2.3:
	print("  -> O(N^2) CONFIRMED")
else:
	print(f"  -> NOT O(N^2) (expected 1.7-2.3)")

print(f"\nGauss-Seidel: alpha = {alpha_gs:.3f}")
if 1.7 <= alpha_gs <= 2.3:
	print("  -> O(N^2) CONFIRMED")
else:
	print(f"  -> NOT O(N^2) (expected 1.7-2.3)")

print(f"\nLU Decomposition: alpha = {alpha_lu:.3f}")
if 2.7 <= alpha_lu <= 3.5:
	print("  -> O(N^3) CONFIRMED")
	print("  -> Direct matrix inversion as expected")
elif alpha_lu < 2.7:
	print(f"  -> LESS than O(N^3)")
	print("  -> Check implementation or test with larger N")
else:
	print(f"  -> MORE than O(N^3)")
	print("  -> May indicate cache/memory bottleneck")

# Crossover point
if len(results) >= 3:
	print("\n" + "=" * 80)
	print("CROSSOVER ANALYSIS")
	print("=" * 80)

	lu_faster = []
	gs_faster = []

	for r in results:
		if r['t_lu'] < r['t_gs']:
			lu_faster.append(r['n'])
		else:
			gs_faster.append(r['n'])

	if lu_faster:
		print(f"\nLU faster for N: {lu_faster}")
	if gs_faster:
		print(f"GS faster for N: {gs_faster}")

	if gs_faster and lu_faster:
		print(f"\nCrossover: LU becomes slower around N = {max(lu_faster)}-{min(gs_faster)}")
	elif gs_faster:
		print(f"\nGS is faster for all tested N (LU becomes competitive at larger N)")
	else:
		print(f"\nLU is faster for all tested N")

print("\n" + "=" * 80)
print("CONCLUSION")
print("=" * 80)

print(f"""
1. Matrix construction: O(N^{alpha_matrix:.1f}) - Same cost for all methods
2. Gauss-Seidel:        O(N^{alpha_gs:.1f}) - Per-iteration cost
3. LU decomposition:    O(N^{alpha_lu:.1f}) - Direct solve (higher cost, no iteration)

For nonlinear problems requiring many iterations:
  - Use Gauss-Seidel for small-to-medium N
  - Consider LU for large N with few unique configurations

For linear problems (converge in 1 iteration):
  - LU may be competitive despite O(N^3) cost
  - Matrix construction O(N^2) dominates for small N
""")

print("=" * 80)
