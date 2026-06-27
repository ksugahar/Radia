#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Analytical Solution Comparison for Quadrupole Field with ObjBckg

Tests magnetizable cube (approximating sphere) in quadrupole background field.
Compares Radia numerical solution with analytical quadrupole field
at points outside the cube.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../build/Release'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

import numpy as np
import radia as rd

mm = 1e-3  # 1 mm in meters

print("=" * 80)
print("Quadrupole Field - Analytical Solution Comparison")
print("=" * 80)

# Parameters
gradient = 10.0  # T/m
R_sphere = 5.0 * mm   # half-size of cube (sphere approximation)

print("\nParameters:")
print(f"  Cube half-size: {R_sphere/mm:.1f} mm")
print(f"  Quadrupole gradient: {gradient} T/m")

# ============================================================================
# Create Geometry
# ============================================================================

print("\n[Step 1] Creating Geometry")
print("-" * 80)

# Simple cubic approximation of sphere: 10mm cube centered at origin
size = 2 * R_sphere  # 10mm cube
half = size / 2
# Hexahedron vertices for cube centered at [0, 0, 0] with dimensions [10, 10, 10] mm
vertices = [
	[-half, -half, -half], [half, -half, -half], [half, half, -half], [-half, half, -half],
	[-half, -half, half], [half, -half, half], [half, half, half], [-half, half, half]
]
cube = rd.ObjHexahedron(vertices, [0, 0, 0])
mat = rd.MatSatIsoFrm([[1596.3, 1.1488], [133.11, 0.4268], [18.713, 0.4759]])
rd.MatApl(cube, mat)
print(f"  Created {size/mm:.0f}x{size/mm:.0f}x{size/mm:.0f} mm cube with MatSatIsoFrm (nonlinear)")

# ============================================================================
# Create Quadrupole Background Field
# ============================================================================

print("\n[Step 2] Creating Quadrupole Background Field")
print("-" * 80)

def quadrupole_field(pos):
	"""Quadrupole field: Bx = g*y, By = g*x, Bz = 0"""
	x, y, z = pos  # Position in meters (Radia always uses meters)
	Bx = gradient * y  # [T]
	By = gradient * x  # [T]
	Bz = 0.0
	return [Bx, By, Bz]

bckg_cf = rd.ObjBckg(quadrupole_field)
print("  Quadrupole field created: Bx = g*y, By = g*x")

# Container with cube and background field
container = rd.ObjCnt([cube, bckg_cf])
print("  Container created")

# ============================================================================
# Solve
# ============================================================================

print("\n[Step 3] Solving Magnetostatic Problem")
print("-" * 80)

print("  Solving...")
solve_result = rd.Solve(container, 1e-5, 5000)
max_abs_M = solve_result[0]  # convergence residual (max |dM|)
n_iter = int(solve_result[3])  # iteration count
print(f"  Solve result: residual={max_abs_M:.2e}, iterations={n_iter}")
if max_abs_M < 1e-5:
	print("  [OK] Solution converged")
else:
	print(f"  [WARNING] Solution may not have converged (max|dM|={max_abs_M:.2e})")

# ============================================================================
# Compare with Analytical Solution
# ============================================================================

print("\n[Step 4] Compare with Analytical Quadrupole Field")
print("-" * 80)

# Test points outside the cube (field ~= background + stray field from cube)
# Note: The magnetized cube produces stray fields that decay with distance.
# Points farther from the cube give closer agreement with the pure background.
test_points = [
	# Along X-axis (y=0, z=0)
	[20*mm, 0, 0],    # r = 20mm
	[30*mm, 0, 0],    # r = 30mm
	[50*mm, 0, 0],    # r = 50mm
	[100*mm, 0, 0],   # r = 100mm
	# Along Y-axis (x=0, z=0)
	[0, 20*mm, 0],    # r = 20mm
	[0, 30*mm, 0],    # r = 30mm
	[0, 50*mm, 0],    # r = 50mm
	[0, 100*mm, 0],   # r = 100mm
	# Diagonal points
	[20*mm, 20*mm, 0],   # r = 28.28mm
	[50*mm, 50*mm, 0],   # r = 70.71mm
	[100*mm, 100*mm, 0], # r = 141.42mm
]

print("\nAnalytical quadrupole field: Bx = g*y, By = g*x, Bz = 0")
print(f"where g = {gradient} T/m")
print(f"\nComparison at points outside cube (r > {R_sphere/mm:.1f} mm):")
print()
print(f"{'Point (mm)':<25} {'r (mm)':>8} | {'B_Radia (T)':^35} | {'B_Analytical (T)':^35} | {'|ΔB| (T)':>12} {'Error (%)':>10}")
print("-" * 140)

errors = []
for pt in test_points:
	# Calculate distance from center
	r = np.sqrt(pt[0]**2 + pt[1]**2 + pt[2]**2)

	# Radia solution
	B_radia = rd.Fld(container, 'b', pt)

	# Analytical quadrupole field (positions already in meters)
	B_analytical = np.array([gradient * pt[1], gradient * pt[0], 0.0])

	# Calculate error
	B_radia_arr = np.array(B_radia)
	delta_B = B_radia_arr - B_analytical
	error_mag = np.linalg.norm(delta_B)
	B_analytical_mag = np.linalg.norm(B_analytical)

	if B_analytical_mag > 1e-10:
		error_pct = error_mag / B_analytical_mag * 100
	else:
		error_pct = 0.0

	errors.append(error_pct)

	# Format output
	pt_mm = [p/mm for p in pt]
	B_radia_str = f"[{B_radia[0]:8.5f}, {B_radia[1]:8.5f}, {B_radia[2]:8.5f}]"
	B_analytical_str = f"[{B_analytical[0]:8.5f}, {B_analytical[1]:8.5f}, {B_analytical[2]:8.5f}]"

	print(f"{str([f'{p:.0f}' for p in pt_mm]):<25} {r/mm:8.2f} | {B_radia_str:^35} | {B_analytical_str:^35} | {error_mag:12.6e} {error_pct:9.4f}%")

# ============================================================================
# Statistics
# ============================================================================

print("\n[Step 5] Error Statistics")
print("-" * 80)

errors_arr = np.array(errors)
print("\nError statistics:")
print(f"  Mean error:    {errors_arr.mean():.4f}%")
print(f"  Median error:  {np.median(errors_arr):.4f}%")
print(f"  Max error:     {errors_arr.max():.4f}%")
print(f"  Min error:     {errors_arr.min():.4f}%")
print(f"  Std deviation: {errors_arr.std():.4f}%")

# Group by distance
print("\nError vs. distance from center:")
distances_mm = [20, 30, 50, 100]
for d_mm in distances_mm:
	d = d_mm * mm
	# Find errors for points at this distance (+-1mm tolerance)
	d_errors = []
	for i, pt in enumerate(test_points):
		r = np.sqrt(pt[0]**2 + pt[1]**2 + pt[2]**2)
		if abs(r - d) < 1.5*mm:  # Tolerance for diagonal points
			d_errors.append(errors[i])

	if d_errors:
		avg_error = np.mean(d_errors)
		print(f"  r ~ {d_mm:2d} mm: {avg_error:6.4f}% average error ({len(d_errors)} points)")

# ============================================================================
# Physical Interpretation
# ============================================================================

print("\n[Step 6] Physical Interpretation")
print("-" * 80)

print(f"\nExpected behavior:")
print(f"  1. Far from cube (r >> {R_sphere/mm:.1f} mm): B_Radia ~ B_Analytical (pure quadrupole)")
print(f"  2. Near cube (r ~ {R_sphere/mm:.1f} mm): Small distortion due to magnetizable material")
print("  3. Error should decrease as 1/r^2 (dipole perturbation)")

# Compute far-field errors (r >= 50mm)
far_field_errors = []
for i, pt in enumerate(test_points):
	r = np.sqrt(pt[0]**2 + pt[1]**2 + pt[2]**2)
	if r >= 50*mm:
		far_field_errors.append(errors[i])
far_field_mean = np.mean(far_field_errors) if far_field_errors else float('inf')

if far_field_mean < 1.0:
	print(f"\n  [OK] Far-field accuracy (r>=50mm): {far_field_mean:.4f}% < 1%")
else:
	print(f"\n  [WARNING] Far-field error higher than expected: {far_field_mean:.4f}%")

if errors_arr[-1] < errors_arr[0]:
	print("  [OK] Error decreases with distance (as expected)")
else:
	print("  [WARNING] Error does not decrease with distance")

# ============================================================================
# Summary
# ============================================================================

print("\n" + "=" * 80)
print("Summary")
print("=" * 80)

print("\n1. ObjBckg successfully implements quadrupole background field")
print(f"2. Radia numerical solution compared with analytical quadrupole at {len(test_points)} points")
print(f"3. Average error: {errors_arr.mean():.4f}%")
print(f"4. Far-field agreement (r>=50mm): {far_field_mean:.4f}%")

if errors_arr.mean() < 5.0:
	print("\n[OK] Good agreement with analytical solution (avg error < 5%)")
elif errors_arr.mean() < 15.0:
	print("\n[OK] Reasonable agreement with analytical solution (avg error < 15%)")
else:
	print("\n[WARNING] Significant deviation from analytical solution")

print("\n" + "=" * 80)
print("Test Complete")
print("=" * 80)

rd.UtiDelAll()
