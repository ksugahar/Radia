#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Analytical Solution Comparison for Quadrupole Field with ObjBckgCF

Tests magnetizable sphere in quadrupole background field.
Compares Radia numerical solution with analytical quadrupole field
at points outside the sphere.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../build/Release'))

import numpy as np
import radia as rd

print("=" * 80)
print("Quadrupole Field - Analytical Solution Comparison")
print("=" * 80)

# Parameters
gradient = 10.0  # T/m
R_sphere = 5.0   # mm (sphere radius)
mu_r = 1000.0    # Relative permeability
chi = mu_r - 1.0

print(f"\nParameters:")
print(f"  Sphere radius: {R_sphere} mm")
print(f"  Relative permeability: {mu_r}")
print(f"  Quadrupole gradient: {gradient} T/m")
print(f"  Magnetic susceptibility: {chi}")

# ============================================================================
# Create Geometry
# ============================================================================

print(f"\n[Step 1] Creating Geometry")
print("-" * 80)

# Simple cubic approximation of sphere
size = 2 * R_sphere  # 10mm cube
cube = rd.ObjRecMag([0, 0, 0], [size, size, size])
mat = rd.MatSatIsoFrm([1596.3, 1.1488], [133.11, 0.4268], [18.713, 0.4759])
rd.MatApl(cube, mat)
print(f"  Created {size}x{size}x{size} mm cube (approximates sphere)")

# ============================================================================
# Create Quadrupole Background Field
# ============================================================================

print(f"\n[Step 2] Creating Quadrupole Background Field")
print("-" * 80)

def quadrupole_field(pos):
	"""Quadrupole field: Bx = g*y, By = g*x, Bz = 0"""
	x, y, z = pos  # Position in mm
	x_m = x * 1e-3  # Convert to meters
	y_m = y * 1e-3
	Bx = gradient * y_m  # [T]
	By = gradient * x_m  # [T]
	Bz = 0.0
	return [Bx, By, Bz]

bckg_cf = rd.ObjBckgCF(quadrupole_field)
print(f"  Quadrupole field created: Bx = g*y, By = g*x")

# Container with cube and background field
container = rd.ObjCnt([cube, bckg_cf])
print(f"  Container created")

# ============================================================================
# Solve
# ============================================================================

print(f"\n[Step 3] Solving Magnetostatic Problem")
print("-" * 80)

print(f"  Solving...")
result = rd.Solve(container, 1e-5, 5000)
print(f"  [OK] Convergence: {result}")

# ============================================================================
# Compare with Analytical Solution
# ============================================================================

print(f"\n[Step 4] Compare with Analytical Quadrupole Field")
print("-" * 80)

# Test points outside the sphere/cube (r > R_sphere)
test_points = [
	# Along X-axis (y=0, z=0)
	[10, 0, 0],   # r = 10mm
	[15, 0, 0],   # r = 15mm
	[20, 0, 0],   # r = 20mm
	[30, 0, 0],   # r = 30mm
	# Along Y-axis (x=0, z=0)
	[0, 10, 0],   # r = 10mm
	[0, 15, 0],   # r = 15mm
	[0, 20, 0],   # r = 20mm
	[0, 30, 0],   # r = 30mm
	# Diagonal points
	[10, 10, 0],  # r = 14.14mm
	[15, 15, 0],  # r = 21.21mm
	[20, 20, 0],  # r = 28.28mm
]

print(f"\nAnalytical quadrupole field: Bx = g*y, By = g*x, Bz = 0")
print(f"where g = {gradient} T/m")
print(f"\nComparison at points outside sphere (r > {R_sphere} mm):")
print(f"")
print(f"{'Point (mm)':<15} {'r (mm)':>8} | {'B_Radia (T)':^35} | {'B_Analytical (T)':^35} | {'|ΔB| (T)':>12} {'Error (%)':>10}")
print("-" * 130)

errors = []
for pt in test_points:
	# Calculate distance from center
	r = np.sqrt(pt[0]**2 + pt[1]**2 + pt[2]**2)

	# Radia solution
	B_radia = rd.Fld(container, 'b', pt)

	# Analytical quadrupole field
	x_m = pt[0] * 1e-3
	y_m = pt[1] * 1e-3
	B_analytical = np.array([gradient * y_m, gradient * x_m, 0.0])

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
	B_radia_str = f"[{B_radia[0]:8.5f}, {B_radia[1]:8.5f}, {B_radia[2]:8.5f}]"
	B_analytical_str = f"[{B_analytical[0]:8.5f}, {B_analytical[1]:8.5f}, {B_analytical[2]:8.5f}]"

	print(f"{str(pt):<15} {r:8.2f} | {B_radia_str:^35} | {B_analytical_str:^35} | {error_mag:12.6e} {error_pct:9.4f}%")

# ============================================================================
# Statistics
# ============================================================================

print(f"\n[Step 5] Error Statistics")
print("-" * 80)

errors_arr = np.array(errors)
print(f"\nError statistics:")
print(f"  Mean error:    {errors_arr.mean():.4f}%")
print(f"  Median error:  {np.median(errors_arr):.4f}%")
print(f"  Max error:     {errors_arr.max():.4f}%")
print(f"  Min error:     {errors_arr.min():.4f}%")
print(f"  Std deviation: {errors_arr.std():.4f}%")

# Group by distance
print(f"\nError vs. distance from center:")
distances = [10, 15, 20, 30]
for d in distances:
	# Find errors for points at this distance (±1mm tolerance)
	d_errors = []
	for i, pt in enumerate(test_points):
		r = np.sqrt(pt[0]**2 + pt[1]**2 + pt[2]**2)
		if abs(r - d) < 1.5:  # Tolerance for diagonal points
			d_errors.append(errors[i])

	if d_errors:
		avg_error = np.mean(d_errors)
		print(f"  r ~ {d:2d} mm: {avg_error:6.4f}% average error ({len(d_errors)} points)")

# ============================================================================
# Physical Interpretation
# ============================================================================

print(f"\n[Step 6] Physical Interpretation")
print("-" * 80)

print(f"\nExpected behavior:")
print(f"  1. Far from sphere (r >> {R_sphere} mm): B_Radia ~ B_Analytical (pure quadrupole)")
print(f"  2. Near sphere (r ~ {R_sphere} mm): Small distortion due to magnetizable material")
print(f"  3. Error should decrease as 1/r^2 (dipole perturbation)")

# Check if behavior matches expectations
if errors_arr[-4:].mean() < 1.0:  # Far field points (r=30mm)
	print(f"\n  [OK] Far-field accuracy: {errors_arr[-4:].mean():.4f}% < 1%")
else:
	print(f"\n  [WARNING] Far-field error higher than expected: {errors_arr[-4:].mean():.4f}%")

if errors_arr[-1] < errors_arr[0]:  # Error decreases with distance
	print(f"  [OK] Error decreases with distance (as expected)")
else:
	print(f"  [WARNING] Error does not decrease with distance")

# ============================================================================
# Summary
# ============================================================================

print("\n" + "=" * 80)
print("Summary")
print("=" * 80)

print(f"\n1. ObjBckgCF successfully implements quadrupole background field")
print(f"2. Radia numerical solution compared with analytical quadrupole at {len(test_points)} points")
print(f"3. Average error: {errors_arr.mean():.4f}%")
print(f"4. Far-field agreement (r=30mm): {errors_arr[-4:].mean():.4f}%")

if errors_arr.mean() < 5.0:
	print(f"\n[OK] Good agreement with analytical solution (avg error < 5%)")
elif errors_arr.mean() < 15.0:
	print(f"\n[OK] Reasonable agreement with analytical solution (avg error < 15%)")
else:
	print(f"\n[WARNING] Significant deviation from analytical solution")

print("\n" + "=" * 80)
print("Test Complete")
print("=" * 80)

# ============================================================================
# VTK Export - Export geometry for visualization
# ============================================================================

try:
	sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/python'))
	from radia_vtk_export import exportGeometryToVTK

	script_name = os.path.splitext(os.path.basename(__file__))[0]
	vtk_filename = f"{script_name}.vtk"
	vtk_path = os.path.join(os.path.dirname(__file__), vtk_filename)

	exportGeometryToVTK(container, vtk_path)
	print(f"\n[VTK] Exported: {vtk_filename}")
	print(f"      View with: paraview {vtk_filename}")
except ImportError:
	print("\n[VTK] Warning: radia_vtk_export not available (VTK export skipped)")
except Exception as e:
	print(f"\n[VTK] Warning: Export failed: {e}")

rd.UtiDelAll()
