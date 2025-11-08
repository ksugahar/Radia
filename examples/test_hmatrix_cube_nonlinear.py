#!/usr/bin/env python3
"""
H-Matrix Test with Nonlinear Material and Background Field

Tests H-matrix with:
- Cube subdivided into N×N×N elements
- Nonlinear material (MatSatIsoFrm) → requires multiple iterations
- Background field (ObjBckg) → all tensor components have rank > 0
"""

import sys
from time import perf_counter
sys.path.insert(0, r"S:\Radia\01_GitHub\build\Release")

import radia as rad

print("=" * 70)
print("H-Matrix Test: Nonlinear Material + Background Field")
print("=" * 70)

# Test cases up to 10×10×10 = 1000
test_cases = [
	(4, 4, 4),    # 64 elements
	(5, 5, 5),    # 125 elements
	(6, 6, 6),    # 216 elements
	(7, 7, 7),    # 343 elements
	(8, 8, 8),    # 512 elements
]

print("\nTest cases (Nx × Ny × Nz):")
for nx, ny, nz in test_cases:
	print(f"  {nx}×{ny}×{nz} = {nx*ny*nz} elements")

print("\nSetup:")
print("  Cube: 100×100×100 mm")
print("  Material: Nonlinear (MatSatIsoFrm), Ms=2000 T, mu_r=1000")
print("  Background field: B = [0, 0, 1.0] T")
print("  → Multiple iterations required\n")

results = []

for nx, ny, nz in test_cases:
	n_elem = nx * ny * nz

	print("=" * 70)
	print(f"Testing {nx}×{ny}×{nz} = {n_elem} elements")
	print("=" * 70)

	cube_size = 100.0
	elem_size_x = cube_size / nx
	elem_size_y = cube_size / ny
	elem_size_z = cube_size / nz

	Bext = [0, 0, 1.0]  # 1 Tesla

	#========================================================================
	# Dense solver
	#========================================================================
	print("\n[1/2] Dense solver...")

	# Nonlinear material: Ms=2000 T, high permeability
	mat_dense = rad.MatSatIsoFrm([2000, 1000], [0.1, 1000])

	elements_dense = []
	for i in range(nx):
		for j in range(ny):
			for k in range(nz):
				x = (i - nx/2 + 0.5) * elem_size_x
				y = (j - ny/2 + 0.5) * elem_size_y
				z = (k - nz/2 + 0.5) * elem_size_z

				elem = rad.ObjRecMag([x, y, z], [elem_size_x, elem_size_y, elem_size_z])
				rad.MatApl(elem, mat_dense)
				elements_dense.append(elem)

	bck_dense = rad.ObjBckg(Bext)
	grp_dense = rad.ObjCnt(elements_dense + [bck_dense])
	rad.SolverHMatrixDisable()

	t_start = perf_counter()
	result_dense = rad.Solve(grp_dense, 1e-4, 1000)
	t_dense = perf_counter() - t_start

	center_idx = n_elem // 2
	m_dense = rad.ObjM(elements_dense[center_idx])
	mag_dense = (m_dense[0][0]**2 + m_dense[0][1]**2 + m_dense[0][2]**2)**0.5

	iter_dense = int(result_dense[3])

	print(f"  Time: {t_dense:.6f} s")
	print(f"  Iterations: {iter_dense}")
	print(f"  |M[center]| = {mag_dense:.3e} A/m")

	#========================================================================
	# H-matrix solver
	#========================================================================
	print("\n[2/2] H-matrix solver...")

	# Nonlinear material: Ms=2000 T, high permeability
	mat_hmat = rad.MatSatIsoFrm([2000, 1000], [0.1, 1000])

	elements_hmat = []
	for i in range(nx):
		for j in range(ny):
			for k in range(nz):
				x = (i - nx/2 + 0.5) * elem_size_x
				y = (j - ny/2 + 0.5) * elem_size_y
				z = (k - nz/2 + 0.5) * elem_size_z

				elem = rad.ObjRecMag([x, y, z], [elem_size_x, elem_size_y, elem_size_z])
				rad.MatApl(elem, mat_hmat)
				elements_hmat.append(elem)

	bck_hmat = rad.ObjBckg(Bext)
	grp_hmat = rad.ObjCnt(elements_hmat + [bck_hmat])
	rad.SolverHMatrixEnable(enable=1, eps=1e-6, max_rank=50)

	t_start = perf_counter()
	result_hmat = rad.Solve(grp_hmat, 1e-4, 1000)
	t_hmat = perf_counter() - t_start

	m_hmat = rad.ObjM(elements_hmat[center_idx])
	mag_hmat = (m_hmat[0][0]**2 + m_hmat[0][1]**2 + m_hmat[0][2]**2)**0.5

	iter_hmat = int(result_hmat[3])

	print(f"  Time: {t_hmat:.6f} s")
	print(f"  Iterations: {iter_hmat}")
	print(f"  |M[center]| = {mag_hmat:.3e} A/m")

	#========================================================================
	# Compare
	#========================================================================
	print("\n[3/3] Comparison:")

	speedup = t_dense / t_hmat if t_hmat > 0 else float('inf')
	error = abs(mag_hmat - mag_dense) / mag_dense if mag_dense > 1e-10 else 0

	print(f"  Dense:    {t_dense:.6f} s ({iter_dense} iter)")
	print(f"  H-matrix: {t_hmat:.6f} s ({iter_hmat} iter)")
	print(f"  Speedup:  {speedup:.3f}x")
	print(f"  Error:    {error*100:.4f}%")

	iter_match = (iter_dense == iter_hmat)
	error_ok = (error < 0.10)

	if speedup > 1.0 and iter_match and error_ok:
		status = f"[SUCCESS] {speedup:.2f}x speedup!"
	elif iter_match and error_ok:
		status = "[OK] Accurate"
	elif not iter_match:
		status = f"[BAD] Iter: {iter_dense} vs {iter_hmat}"
	else:
		status = f"[BAD] Error: {error*100:.2f}%"

	print(f"  Status: {status}")

	results.append({
		'nx': nx, 'ny': ny, 'nz': nz, 'n': n_elem,
		'iter_dense': iter_dense,
		'iter_hmat': iter_hmat,
		't_dense': t_dense,
		't_hmat': t_hmat,
		'speedup': speedup,
		'error': error,
		'status': status
	})

	rad.SolverHMatrixDisable()
	rad.UtiDelAll()

#============================================================================
# Summary
#============================================================================
print("\n" + "=" * 70)
print("Summary")
print("=" * 70)
print(f"\n{'Grid':>10}  {'N':>6}  {'Iter':>8}  {'Dense(s)':>10}  {'H-mat(s)':>10}  {'Speedup':>10}  {'Status':>25}")
print("-" * 70)

for r in results:
	grid_str = f"{r['nx']}×{r['ny']}×{r['nz']}"
	iter_str = f"{r['iter_dense']}/{r['iter_hmat']}"
	print(f"{grid_str:>10}  {r['n']:>6}  {iter_str:>8}  {r['t_dense']:>10.6f}  {r['t_hmat']:>10.6f}  {r['speedup']:>9.3f}x  {r['status']:>25}")

print("=" * 70)

# Check for speedup
speedup_found = any(r['speedup'] > 1.0 for r in results)

if speedup_found:
	best = max(results, key=lambda r: r['speedup'])
	print(f"\n[SUCCESS] H-matrix provides speedup!")
	print(f"Best: {best['speedup']:.2f}x at N={best['n']}")
	print(f"高速化しました！(Accelerated!)")
else:
	print(f"\n[INFO] H-matrix still slower for tested sizes")
	best = max(results, key=lambda r: r['speedup'])
	print(f"Best: {best['speedup']:.3f}x at N={best['n']}")

	# Estimate time per iteration
	if best['iter_hmat'] > 0:
		t_per_iter_dense = best['t_dense'] / best['iter_dense']
		t_per_iter_hmat = best['t_hmat'] / best['iter_hmat']
		print(f"\nTime per iteration:")
		print(f"  Dense:    {t_per_iter_dense*1000:.2f} ms/iter")
		print(f"  H-matrix: {t_per_iter_hmat*1000:.2f} ms/iter")

print("=" * 70)
