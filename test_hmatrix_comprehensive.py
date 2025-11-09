#!/usr/bin/env python
"""
H-Matrix Comprehensive Test

Tests H-matrix field evaluation with proper subdivisions.

Author: Radia development team
Date: 2025-11-09
"""

import sys
sys.path.insert(0, r"S:\Radia\01_GitHub\build\Release")

import radia as rad
import numpy as np
import time

print("=" * 80)
print("H-Matrix Comprehensive Test")
print("=" * 80)

# Test 1: Single element WITHOUT subdivision (baseline)
print("\n[Test 1] Single element WITHOUT subdivision")
print("-" * 80)

m1 = rad.ObjRecMag([0, 0, 0], [10, 10, 10], [0, 0, 1.0])
obs_pts_1 = [[20, 0, 0], [0, 20, 0], [0, 0, 20]]

# Direct
H_direct1 = rad.FldBatch(m1, 'h', obs_pts_1, 0)
print(f"Direct field at [20,0,0]: [{H_direct1[0][0]:.6e}, {H_direct1[0][1]:.6e}, {H_direct1[0][2]:.6e}]")

# H-matrix
rad.SetHMatrixFieldEval(1, 1e-8)
H_hmat1 = rad.FldBatch(m1, 'h', obs_pts_1, 1)
print(f"H-matrix field at [20,0,0]: [{H_hmat1[0][0]:.6e}, {H_hmat1[0][1]:.6e}, {H_hmat1[0][2]:.6e}]")

# Error
H_direct1 = np.array(H_direct1)
H_hmat1 = np.array(H_hmat1)
diff1 = H_hmat1 - H_direct1
rms1 = np.sqrt(np.mean(diff1**2))
max_H1 = np.max(np.abs(H_direct1))
rel_rms1 = rms1 / max_H1 if max_H1 > 0 else 0
print(f"Relative RMS error: {rel_rms1:.2e}")
print(f"Result: {'PASS' if rel_rms1 < 0.01 else 'FAIL'} (expected: 0.00e+00 for single element)")

rad.UtiDelAll()

# Test 2: Single element WITH subdivision
print("\n[Test 2] Single element WITH subdivision (3x3x3)")
print("-" * 80)

m2 = rad.ObjRecMag([0, 0, 0], [10, 10, 10], [0, 0, 1.0])
m2_div = rad.ObjDivMag(m2, [3, 3, 3])
print("Element subdivided into 3x3x3 = 27 sub-elements")

# Direct
H_direct2 = rad.FldBatch(m2_div, 'h', obs_pts_1, 0)
print(f"Direct field at [20,0,0]: [{H_direct2[0][0]:.6e}, {H_direct2[0][1]:.6e}, {H_direct2[0][2]:.6e}]")

# H-matrix
H_hmat2 = rad.FldBatch(m2_div, 'h', obs_pts_1, 1)
print(f"H-matrix field at [20,0,0]: [{H_hmat2[0][0]:.6e}, {H_hmat2[0][1]:.6e}, {H_hmat2[0][2]:.6e}]")

# Error
H_direct2 = np.array(H_direct2)
H_hmat2 = np.array(H_hmat2)
diff2 = H_hmat2 - H_direct2
rms2 = np.sqrt(np.mean(diff2**2))
max_H2 = np.max(np.abs(H_direct2))
rel_rms2 = rms2 / max_H2 if max_H2 > 0 else 0
print(f"Relative RMS error: {rel_rms2:.2e}")
print(f"Result: {'PASS' if rel_rms2 < 0.01 else 'FAIL'} (expected: < 1% with subdivision)")

rad.UtiDelAll()

# Test 3: Multiple elements WITH subdivisions
print("\n[Test 3] Multiple elements WITH subdivisions (5 magnets x 2^3 = 40 sources)")
print("-" * 80)

magnets = []
for i in range(5):
	x = i * 20.0
	y = 0.0
	z = 0.0

	m = rad.ObjRecMag([x, y, z], [8, 8, 8], [0.8, 0.5, 1.0])
	m_div = rad.ObjDivMag(m, [2, 2, 2])  # 2x2x2 = 8 sub-elements each
	magnets.append(m_div)

geometry = rad.ObjCnt(magnets)
print(f"Created 5 magnets, each subdivided 2x2x2 = 8 sub-elements")
print(f"Total sub-elements: 5 x 8 = 40")

# Test points
test_pts = [[50, 20, 30], [30, 15, 25], [70, 10, 20]]

# Direct
t0 = time.time()
H_direct3 = rad.FldBatch(geometry, 'h', test_pts, 0)
t_direct = time.time() - t0
print(f"\nDirect calculation:")
print(f"  Time: {t_direct*1000:.2f} ms")
print(f"  Field at [50,20,30]: [{H_direct3[0][0]:.6e}, {H_direct3[0][1]:.6e}, {H_direct3[0][2]:.6e}]")

# H-matrix
t0 = time.time()
H_hmat3 = rad.FldBatch(geometry, 'h', test_pts, 1)
t_hmat = time.time() - t0
print(f"\nH-matrix calculation:")
print(f"  Time: {t_hmat*1000:.2f} ms")
print(f"  Field at [50,20,30]: [{H_hmat3[0][0]:.6e}, {H_hmat3[0][1]:.6e}, {H_hmat3[0][2]:.6e}]")

# Error
H_direct3 = np.array(H_direct3)
H_hmat3 = np.array(H_hmat3)
diff3 = H_hmat3 - H_direct3
rms3 = np.sqrt(np.mean(diff3**2))
max_H3 = np.max(np.abs(H_direct3))
rel_rms3 = rms3 / max_H3 if max_H3 > 0 else 0

print(f"\nAccuracy:")
print(f"  Relative RMS error: {rel_rms3:.2e}")
print(f"  Result: {'PASS' if rel_rms3 < 0.01 else 'FAIL'} (expected: < 1%)")

rad.UtiDelAll()

# Test 4: Large-scale test with subdivisions
print("\n[Test 4] Large-scale test (20 magnets x 3^3 = 540 sources, 100 obs. points)")
print("-" * 80)

magnets = []
for i in range(20):
	x = (i % 5) * 25.0
	y = (i // 5) * 25.0
	z = 0.0

	angle = i * np.pi / 10.0
	Mx = 0.8 * np.cos(angle)
	My = 0.8 * np.sin(angle)
	Mz = 0.6

	m = rad.ObjRecMag([x, y, z], [10, 10, 10], [Mx, My, Mz])
	m_div = rad.ObjDivMag(m, [3, 3, 3])  # 3x3x3 = 27 sub-elements
	magnets.append(m_div)

geometry = rad.ObjCnt(magnets)
print(f"Created 20 magnets, each subdivided 3x3x3 = 27 sub-elements")
print(f"Total sub-elements: 20 x 27 = 540")

# Observation points
obs_grid = []
for i in range(10):
	for j in range(10):
		x = 50 + i * 5.0
		y = 50 + j * 5.0
		z = 30.0
		obs_grid.append([x, y, z])

print(f"Observation grid: 10 x 10 = 100 points")

# Direct calculation
print("\nDirect calculation...")
t0 = time.time()
H_direct4 = rad.FldBatch(geometry, 'h', obs_grid, 0)
t_direct = time.time() - t0
print(f"  Time: {t_direct*1000:.2f} ms")

# H-matrix calculation
print("\nH-matrix calculation...")
t0 = time.time()
H_hmat4 = rad.FldBatch(geometry, 'h', obs_grid, 1)
t_hmat = time.time() - t0
print(f"  Time: {t_hmat*1000:.2f} ms")

# Analysis
H_direct4 = np.array(H_direct4)
H_hmat4 = np.array(H_hmat4)
diff4 = H_hmat4 - H_direct4
rms4 = np.sqrt(np.mean(diff4**2))
max_H4 = np.max(np.abs(H_direct4))
rel_rms4 = rms4 / max_H4 if max_H4 > 0 else 0

speedup = t_direct / t_hmat if t_hmat > 0 else 0

print(f"\nResults:")
print(f"  Speedup: {speedup:.2f}x")
print(f"  Relative RMS error: {rel_rms4:.2e}")
print(f"  Accuracy: {'PASS' if rel_rms4 < 0.05 else 'FAIL'} (expected: < 5%)")

# Summary
print("\n" + "=" * 80)
print("Summary")
print("=" * 80)

tests = [
	("Test 1: Single element (no subdivision)", rel_rms1, 0.01, "0%"),
	("Test 2: Single element (3x3x3 subdivision)", rel_rms2, 0.01, "< 1%"),
	("Test 3: 5 magnets (2x2x2 subdivision)", rel_rms3, 0.01, "< 1%"),
	("Test 4: 20 magnets (3x3x3 subdivision)", rel_rms4, 0.05, "< 5%"),
]

all_passed = True
for name, error, threshold, expected in tests:
	status = "PASS" if error < threshold else "FAIL"
	if status == "FAIL":
		all_passed = False
	print(f"{name:50s} Error: {error:8.2e}  Expected: {expected:6s}  [{status}]")

print("\n" + "=" * 80)
if all_passed:
	print("OVERALL RESULT: ALL TESTS PASSED")
	print("\nH-matrix implementation is working correctly!")
	print("Remember to subdivide elements using rad.ObjDivMag() for accurate results.")
else:
	print("OVERALL RESULT: SOME TESTS FAILED")
	print("\nPlease check the implementation.")

print("=" * 80)
