#!/usr/bin/env python3
"""
H-Matrix Accuracy Test

Directly compares the interaction matrix computed by H-matrix vs dense solver.
Tests matrix-vector multiplication: H = M * M_in
"""

import sys
import numpy as np
sys.path.insert(0, r"S:\Radia\01_GitHub\build\Release")

import radia as rad

print("=" * 70)
print("H-Matrix Accuracy Test - Direct Matrix Comparison")
print("=" * 70)

# Small test case for detailed comparison
n_elem = 10
print(f"\nTest size: N = {n_elem} elements")
print("Comparing H-matrix vs Dense interaction matrix\n")

# Grid dimensions
nx = int(n_elem**0.5)
ny = (n_elem + nx - 1) // nx

size = 10.0  # mm
gap = 2.0    # mm

# Create quasi-linear material
iron = rad.MatSatIsoFrm([10000, 1000], [0.1, 1000])

# Create permanent magnet as excitation
pm = rad.ObjRecMag([0, 0, -50], [100, 100, 10], [0, 0, 1.0])

print("=" * 70)
print("Test 1: Dense Solver - Get Reference Magnetization")
print("=" * 70)

# Create elements for dense solver
elements_dense = []
for i in range(n_elem):
	x = (i % nx) * (size + gap)
	y = (i // nx) * (size + gap)
	z = 0
	elem = rad.ObjRecMag([x, y, z], [size, size, size])
	rad.MatApl(elem, iron)
	elements_dense.append(elem)

grp_dense = rad.ObjCnt(elements_dense + [pm])
rad.SolverHMatrixDisable()

result_dense = rad.Solve(grp_dense, 1e-4, 1000)
print(f"Dense solver result: {result_dense}")
print(f"Iterations: {result_dense[3]}")

# Get magnetization from dense solver
M_dense = []
for elem in elements_dense:
	m = rad.ObjM(elem)
	M_dense.append(m[0])  # [Mx, My, Mz]

print(f"\nFirst 3 magnetizations (Dense):")
for i in range(min(3, n_elem)):
	print(f"  M[{i}] = [{M_dense[i][0]:.6e}, {M_dense[i][1]:.6e}, {M_dense[i][2]:.6e}] A/m")

rad.UtiDelAll()

print("\n" + "=" * 70)
print("Test 2: H-Matrix Solver - Get H-Matrix Magnetization")
print("=" * 70)

# Create elements for H-matrix solver
iron_hmat = rad.MatSatIsoFrm([10000, 1000], [0.1, 1000])
pm_hmat = rad.ObjRecMag([0, 0, -50], [100, 100, 10], [0, 0, 1.0])

elements_hmat = []
for i in range(n_elem):
	x = (i % nx) * (size + gap)
	y = (i // nx) * (size + gap)
	z = 0
	elem = rad.ObjRecMag([x, y, z], [size, size, size])
	rad.MatApl(elem, iron_hmat)
	elements_hmat.append(elem)

grp_hmat = rad.ObjCnt(elements_hmat + [pm_hmat])
rad.SolverHMatrixEnable(enable=1, eps=1e-6, max_rank=50)

result_hmat = rad.Solve(grp_hmat, 1e-4, 1000)
print(f"H-matrix solver result: {result_hmat}")
print(f"Iterations: {result_hmat[3]}")

# Get magnetization from H-matrix solver
M_hmat = []
for elem in elements_hmat:
	m = rad.ObjM(elem)
	M_hmat.append(m[0])  # [Mx, My, Mz]

print(f"\nFirst 3 magnetizations (H-matrix):")
for i in range(min(3, n_elem)):
	print(f"  M[{i}] = [{M_hmat[i][0]:.6e}, {M_hmat[i][1]:.6e}, {M_hmat[i][2]:.6e}] A/m")

print("\n" + "=" * 70)
print("Test 3: Compare Magnetizations")
print("=" * 70)

print(f"\n{'Index':>5}  {'|M_dense|':>15}  {'|M_hmat|':>15}  {'Rel Error':>12}")
print("-" * 70)

max_error = 0.0
for i in range(n_elem):
	mag_dense = np.linalg.norm(M_dense[i])
	mag_hmat = np.linalg.norm(M_hmat[i])

	if mag_dense > 1e-10:
		error = abs(mag_hmat - mag_dense) / mag_dense
		max_error = max(max_error, error)
	else:
		error = 0.0

	if i < 5 or i >= n_elem - 2:  # Show first 5 and last 2
		print(f"{i:>5}  {mag_dense:>15.6e}  {mag_hmat:>15.6e}  {error*100:>11.2f}%")
	elif i == 5:
		print("  ...")

print("-" * 70)
print(f"Maximum relative error: {max_error*100:.2f}%")

print("\n" + "=" * 70)
print("Summary")
print("=" * 70)

if max_error < 0.01:
	print(f"[SUCCESS] H-matrix accuracy is good (<1% error)")
elif max_error < 0.10:
	print(f"[PARTIAL] H-matrix accuracy is acceptable (<10% error)")
else:
	print(f"[FAIL] H-matrix has large errors (>{max_error*100:.1f}%)")
	print(f"\nPossible causes:")
	print(f"  1. ACA approximation too coarse (eps={1e-6})")
	print(f"  2. Some tensor components approximated as zero (rank=0)")
	print(f"  3. Bug in kernel function or symmetry handling")
	print(f"  4. Bug in H-matrix-vector multiplication")

print("=" * 70)

rad.SolverHMatrixDisable()
