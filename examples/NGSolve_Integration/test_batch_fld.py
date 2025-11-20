#!/usr/bin/env python3
"""
Test: Does rad.Fld() support batch evaluation?

Tests whether Radia can evaluate field at multiple points in one call.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "build", "Release"))

import radia as rad
import numpy as np
import time

print("=" * 80)
print("Test: Batch Field Evaluation with rad.Fld()")
print("=" * 80)

# Create simple magnet
rad.UtiDelAll()
rad.FldUnits('m')  # Set units to meters (required for NGSolve integration)
magnet = rad.ObjRecMag([0, 0, 0], [0.020, 0.020, 0.020], [0, 0, 1.2])  # meters

print("\n[Test 1] Single point evaluation")
print("-" * 80)

point = [0.010, 0.0, 0.0]  # meters
B_single = rad.Fld(magnet, 'b', point)
print(f"  Point: {point} m")
print(f"  B = {B_single} T")

print("\n[Test 2] Try passing multiple points as list")
print("-" * 80)

points_list = [
	[0.010, 0.0, 0.0],   # meters
	[0.0, 0.010, 0.0],
	[0.0, 0.0, 0.010],
]

try:
	result = rad.Fld(magnet, 'b', points_list)
	print(f"  Result type: {type(result)}")
	print(f"  Result: {result}")

	if isinstance(result, list) and len(result) == 3:
		# Check if each element is itself a list (field vector)
		if isinstance(result[0], list) and len(result[0]) == 3:
			print("\n  [OK] SUCCESS: rad.Fld() supports batch evaluation!")
			print(f"  Number of points evaluated: {len(result)}")
			for i, (point, field) in enumerate(zip(points_list, result)):
				print(f"    Point {i}: {point} m -> B = {field} T")
		else:
			print("\n  [FAIL] Result is single field vector, not batch")
			print(f"    This means rad.Fld() interpreted the list as a single point")
	else:
		print(f"\n  [FAIL] Unexpected result format")

except Exception as e:
	print(f"  [FAIL] ERROR: {e}")
	print(f"  rad.Fld() does not support batch evaluation")

print("\n[Test 3] Try passing flattened array")
print("-" * 80)

# Flatten: [x1,y1,z1, x2,y2,z2, x3,y3,z3]
points_flat = [coord for point in points_list for coord in point]
print(f"  Flattened array length: {len(points_flat)}")

try:
	result = rad.Fld(magnet, 'b', points_flat)
	print(f"  Result type: {type(result)}")
	print(f"  Result: {result}")
	print(f"  [FAIL] Interpreted as single point with wrong dimension")
except Exception as e:
	print(f"  [FAIL] ERROR: {e}")

print("\n[Test 4] Performance comparison")
print("-" * 80)

N_points = 1000
print(f"  Evaluating {N_points} points...")

# Generate random points in meters
np.random.seed(42)
points = [[np.random.uniform(-0.010, 0.010),
           np.random.uniform(-0.010, 0.010),
           np.random.uniform(-0.010, 0.010)] for _ in range(N_points)]

# Method 1: Loop (current approach)
print("\n  Method 1: Loop over single point calls")
t_start = time.perf_counter()
results_loop = []
for point in points:
	B = rad.Fld(magnet, 'b', point)
	results_loop.append(B)
t_loop = time.perf_counter() - t_start
print(f"    Time: {t_loop*1000:.2f} ms ({t_loop*1e6/N_points:.2f} us/point)")

# Method 2: Try batch (if supported)
print("\n  Method 2: Batch evaluation (if supported)")
try:
	t_start = time.perf_counter()
	results_batch = rad.Fld(magnet, 'b', points)
	t_batch = time.perf_counter() - t_start

	# Verify result format
	if isinstance(results_batch, list) and len(results_batch) == N_points:
		print(f"    Time: {t_batch*1000:.2f} ms ({t_batch*1e6/N_points:.2f} us/point)")
		print(f"    Speedup: {t_loop/t_batch:.2f}x")

		# Verify correctness
		max_error = 0.0
		for i in range(min(10, N_points)):
			error = np.linalg.norm(np.array(results_batch[i]) - np.array(results_loop[i]))
			max_error = max(max_error, error)
		print(f"    Max error (first 10 points): {max_error:.6e} T")
	else:
		print(f"    [FAIL] Unexpected result format (not batch)")
except Exception as e:
	print(f"    [FAIL] Not supported: {e}")

print("\n" + "=" * 80)
print("CONCLUSION")
print("=" * 80)

print("""
If batch evaluation is supported:
  - Can implement SIMD evaluation in rad_ngsolve.cpp
  - Expected speedup: 2-10x for GridFunction.Set()

If not supported:
  - Need to modify Radia Python bindings first
  - Or implement batching at C++ level in rad_ngsolve.cpp
""")

print("=" * 80)
