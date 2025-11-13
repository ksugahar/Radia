#!/usr/bin/env python
"""
H-matrix accuracy test with 100+ elements and 100+ observation points

This test demonstrates proper H-matrix usage as specified in CLAUDE.md:
- Use rad.SetHMatrixFieldEval() to enable H-matrix explicitly
- Use at least 100 elements and 100 observation points
- Compare H-matrix vs direct calculation accuracy
"""

import sys
sys.path.insert(0, r"S:\Radia\01_GitHub\build\Release")

import radia as rad
import numpy as np

print("=" * 80)
print("H-Matrix Accuracy Test (100+ elements, 100+ points)")
print("=" * 80)

# Create geometry: 10x10 grid = 100 rectangular magnets
rad.UtiDelAll()

magnets = []
grid_size = 10  # 10x10 = 100 magnets

for i in range(grid_size):
	for j in range(grid_size):
		x = i * 15.0  # 15mm spacing
		y = j * 15.0
		mag = rad.ObjRecMag([x, y, 0], [10, 10, 10], [0, 0, 1])
		magnets.append(mag)

container = rad.ObjCnt(magnets)

# Observation points: 10x10 grid = 100 points above the magnets
obs_points = []
for i in range(grid_size):
	for j in range(grid_size):
		x = i * 15.0 + 7.5  # Center of each cell
		y = j * 15.0 + 7.5
		z = 25.0  # 25mm above
		obs_points.append([x, y, z])

print(f"\n[Setup]")
print(f"  Number of magnets: {len(magnets)}")
print(f"  Number of observation points: {len(obs_points)}")
print(f"  Grid configuration: {grid_size}x{grid_size}")

# Test 1: Direct calculation (use_hmatrix=0)
print("\n[Test 1] Direct calculation")
print("-" * 80)

import time
t0 = time.time()
H_direct = rad.FldBatch(container, 'h', obs_points, 0)  # use_hmatrix=0
t_direct = time.time() - t0

H_direct = np.array(H_direct).reshape(-1, 3)

print(f"  Calculation time: {t_direct:.4f} seconds")
print(f"  First point H: [{H_direct[0,0]:.6e}, {H_direct[0,1]:.6e}, {H_direct[0,2]:.6e}] A/m")
print(f"  Mean |H|: {np.mean(np.sqrt(np.sum(H_direct**2, axis=1))):.6e} A/m")

# Test 2: H-matrix (use_hmatrix=1)
print("\n[Test 2] H-matrix calculation")
print("-" * 80)

# Enable H-matrix explicitly as per CLAUDE.md guidelines
rad.SetHMatrixFieldEval(1, 1e-6)

t0 = time.time()
H_hmat = rad.FldBatch(container, 'h', obs_points, 1)  # use_hmatrix=1
t_hmat = time.time() - t0

H_hmat = np.array(H_hmat).reshape(-1, 3)

print(f"  Calculation time: {t_hmat:.4f} seconds")
print(f"  First point H: [{H_hmat[0,0]:.6e}, {H_hmat[0,1]:.6e}, {H_hmat[0,2]:.6e}] A/m")
print(f"  Mean |H|: {np.mean(np.sqrt(np.sum(H_hmat**2, axis=1))):.6e} A/m")

# Compute error
diff = H_hmat - H_direct
abs_errors = np.sqrt(np.sum(diff**2, axis=1))
mag_H = np.sqrt(np.sum(H_direct**2, axis=1))
rel_errors = abs_errors / (mag_H + 1e-15)  # Avoid division by zero

max_rel_error = np.max(rel_errors)
mean_rel_error = np.mean(rel_errors)

print(f"\n[Accuracy]")
print(f"  Max relative error:  {max_rel_error:.6e} ({max_rel_error*100:.4f}%)")
print(f"  Mean relative error: {mean_rel_error:.6e} ({mean_rel_error*100:.4f}%)")
print(f"  Max absolute error:  {np.max(abs_errors):.6e} A/m")

print(f"\n[Performance]")
print(f"  Direct time:  {t_direct:.4f} seconds")
print(f"  H-matrix time: {t_hmat:.4f} seconds")
print(f"  Speedup: {t_direct/t_hmat:.2f}x")

# Get H-matrix stats
stats = rad.GetHMatrixStats()
print(f"\n[H-matrix Stats]")
print(f"  Enabled: {stats[0]}")
print(f"  Cached:  {stats[1]}")
print(f"  Memory:  {stats[2]:.3f} MB")

# Pass/Fail criteria
print("\n" + "=" * 80)
if max_rel_error < 0.01:
	print(f"PASS: H-matrix accuracy < 1% (actual: {max_rel_error*100:.4f}%)")
	exit_code = 0
else:
	print(f"FAIL: H-matrix accuracy > 1% (actual: {max_rel_error*100:.4f}%)")
	exit_code = 1

print("=" * 80)

sys.exit(exit_code)
