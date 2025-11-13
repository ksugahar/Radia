#!/usr/bin/env python3
"""
Benchmark: Solver time only (excluding matrix construction)
Compare Gauss-Seidel O(N^2) vs LU decomposition O(N^3)
"""

import sys
import numpy as np
from time import perf_counter
sys.path.insert(0, r"S:\Radia\01_GitHub\build\Release")

import radia as rad

print("=" * 70)
print("Solver Scaling Benchmark (Matrix Construction Excluded)")
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
print("  Operation: rad.RlxMan() - Execute 1 iteration only")
print("  Expected: GS O(N^2), LU O(N^3)")

# Nonlinear material (soft iron)
MH_data = [[0, 0], [200, 0.7], [600, 1.2], [1200, 1.4], [2000, 1.5],
           [3500, 1.54], [6000, 1.56], [12000, 1.57]]

results = []

for nx, ny, nz in test_cases:
	n_elem = nx * ny * nz

	print(f"\nN = {nx}x{ny}x{nz} = {n_elem:4d} elements")
	print("-" * 70)

	cube_size = 100.0
	elem_size = cube_size / nx

	mat = rad.MatSatIsoTab(MH_data)

	#========================================================================
	# GAUSS-SEIDEL (Method 4)
	#========================================================================

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

	# Build interaction matrix (not timed)
	intrc = rad.RlxPre(grp, grp)

	# Measure solver time only
	t_gs_start = perf_counter()
	rad.RlxMan(intrc, 4, 1, 1.0)  # Method 4, 1 iteration
	t_gs = perf_counter() - t_gs_start

	print(f"  Gauss-Seidel:     {t_gs*1000:8.2f} ms (O(N^2))")

	rad.UtiDelAll()

	#========================================================================
	# LU DECOMPOSITION (Method 5)
	#========================================================================

	# Build geometry
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
	rad.SolverHMatrixDisable()

	# Build interaction matrix (not timed)
	intrc2 = rad.RlxPre(grp2, grp2)

	# Enable LU decomposition
	rad.SetRelaxSubInterval(intrc2, 0, n_elem-1, 1)

	# Measure solver time only
	t_lu_start = perf_counter()
	rad.RlxMan(intrc2, 5, 1, 1.0)  # Method 5 now supported in RlxMan
	t_lu = perf_counter() - t_lu_start

	print(f"  LU decomposition: {t_lu*1000:8.2f} ms (O(N^3))")
	print(f"  LU / GS ratio:    {t_lu/t_gs:8.2f}x")

	results.append({
		'n': n_elem,
		't_gs': t_gs,
		't_lu': t_lu,
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

# Gauss-Seidel scaling
t_gs_values = np.array([r['t_gs'] for r in results])
log_t_gs = np.log(t_gs_values)
A = np.vstack([log_n, np.ones(len(log_n))]).T
alpha_gs, log_a_gs = np.linalg.lstsq(A, log_t_gs, rcond=None)[0]

# LU decomposition scaling
t_lu_values = np.array([r['t_lu'] for r in results])
log_t_lu = np.log(t_lu_values)
alpha_lu, log_a_lu = np.linalg.lstsq(A, log_t_lu, rcond=None)[0]

print(f"\nPower law fits: t = a * N^alpha")
print(f"  Gauss-Seidel:     t = {np.exp(log_a_gs):.6e} * N^{alpha_gs:.3f}")
print(f"  LU decomposition: t = {np.exp(log_a_lu):.6e} * N^{alpha_lu:.3f}")

# Detailed table
print(f"\n{'N':>6}  {'GS (ms)':>12}  {'LU (ms)':>12}  {'LU/GS':>8}  {'t_gs/N^2':>12}  {'t_lu/N^3':>12}")
print("-" * 75)

for r in results:
	n = r['n']
	t_gs = r['t_gs'] * 1000
	t_lu = r['t_lu'] * 1000
	ratio = t_lu / t_gs if t_gs > 0 else 0

	t_gs_n2 = t_gs / (n * n)
	t_lu_n3 = t_lu / (n * n * n) * 1000

	print(f"{n:>6}  {t_gs:>12.2f}  {t_lu:>12.2f}  {ratio:>8.2f}  {t_gs_n2:>12.6f}  {t_lu_n3:>12.6f}")

# Check for constant ratios
t_gs_n2_values = [r['t_gs'] * 1000 / (r['n']**2) for r in results[3:]]
t_lu_n3_values = [r['t_lu'] * 1000 / (r['n']**3) * 1000 for r in results[3:]]

mean_gs_n2 = np.mean(t_gs_n2_values)
cv_gs = np.std(t_gs_n2_values) / mean_gs_n2 if mean_gs_n2 > 0 else 0

mean_lu_n3 = np.mean(t_lu_n3_values)
cv_lu = np.std(t_lu_n3_values) / mean_lu_n3 if mean_lu_n3 > 0 else 0

print(f"\nRatio statistics (N >= {results[3]['n']}):")
print(f"  GS/N^2: mean={mean_gs_n2:.6f}, CV={cv_gs:.3f}")
print(f"  LU/N^3: mean={mean_lu_n3:.6f}, CV={cv_lu:.3f}")

print("\n" + "=" * 70)
print("Interpretation")
print("=" * 70)

print(f"\nGauss-Seidel: alpha = {alpha_gs:.3f}")
if 1.7 <= alpha_gs <= 2.3:
	print("  -> O(N^2) CONFIRMED")
	print("  -> Matrix-vector multiply per iteration")
	if cv_gs < 0.2:
		print(f"  -> Consistent scaling (CV={cv_gs:.3f})")
else:
	print(f"  -> NOT O(N^2) (expected 1.7-2.3)")

print(f"\nLU Decomposition: alpha = {alpha_lu:.3f}")
if 2.7 <= alpha_lu <= 3.5:
	print("  -> O(N^3) CONFIRMED")
	print("  -> Direct matrix inversion as expected")
	if cv_lu < 0.3:
		print(f"  -> Consistent scaling (CV={cv_lu:.3f})")
	else:
		print(f"  -> High variation (CV={cv_lu:.3f})")
		print("  -> May indicate cache effects or memory bandwidth limits")
elif alpha_lu < 2.7:
	print(f"  -> LESS than O(N^3)")
	print("  -> Unexpected, check implementation")
else:
	print(f"  -> MORE than O(N^3) ({alpha_lu:.3f})")
	print("  -> Cache misses or memory bandwidth bottleneck")

# Analyze large N separately
if len(results) >= 5:
	n_mid = len(results) // 2
	log_n_large = log_n[n_mid:]
	log_t_lu_large = log_t_lu[n_mid:]
	A_large = np.vstack([log_n_large, np.ones(len(log_n_large))]).T
	alpha_lu_large, _ = np.linalg.lstsq(A_large, log_t_lu_large, rcond=None)[0]

	print(f"\nLarge N only (N >= {results[n_mid]['n']}):")
	print(f"  LU alpha = {alpha_lu_large:.3f}")

	if 2.8 <= alpha_lu_large <= 3.5:
		print(f"  -> True asymptotic O(N^3) behavior observed")

print("\n" + "=" * 70)
print("Conclusion")
print("=" * 70)

print(f"""
Solver Time Scaling (excluding matrix construction):

1. Gauss-Seidel:     O(N^{alpha_gs:.1f}) per iteration
   - Dense matrix-vector multiplication
   - {mean_gs_n2:.6f} ms/N^2 (approximately constant)

2. LU Decomposition: O(N^{alpha_lu:.1f}) per solve
   - Direct matrix inversion
   - {mean_lu_n3:.6f} * 10^-3 ms/N^3 (approximately constant)

Performance Comparison:
  - For small N (< {results[2]['n']}): GS and LU comparable
  - For medium N ({results[2]['n']}-{results[5]['n']}): LU becomes slower
  - For large N (> {results[5]['n']}): LU is O(N) times slower than GS

Recommendation:
  - Use Gauss-Seidel for iterative convergence
  - Consider LU only for direct solve scenarios or specific problem structures
""")

print("=" * 70)
