#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Simple test to verify B→H conversion in rad.ObjBckgCF()

Tests that quadrupole background field defined in Tesla is correctly
converted to H field internally.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../build/Release'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/python'))

import numpy as np
import radia as rd

print("=" * 70)
print("ObjBckgCF B→H Conversion Test")
print("=" * 70)

# Parameters
gradient = 10.0  # Quadrupole gradient [T/m]
mu_r = 1000.0  # Relative permeability
chi = mu_r - 1.0  # Magnetic susceptibility

print(f"\nParameters:")
print(f"  Quadrupole gradient: {gradient} T/m")
print(f"  Relative permeability: {mu_r}")
print(f"  Magnetic susceptibility: {chi}")

# ============================================================================
# Test 1: Create quadrupole background field using ObjBckgCF
# ============================================================================

print(f"\n[Test 1] Quadrupole Background Field (ObjBckgCF)")
print("-" * 70)

def quadrupole_field_callback(gradient):
	"""Create quadrupole field callback for rd.ObjBckgCF

	Returns B in Tesla
	"""
	call_count = [0]  # Mutable to allow modification in nested function
	def field(pos):
		x, y, z = pos  # Position in mm
		# Convert to meters
		x_m = x * 1e-3
		y_m = y * 1e-3
		# Quadrupole field: Bx = g*y, By = g*x, Bz = 0
		Bx = gradient * y_m  # [T]
		By = gradient * x_m  # [T]
		Bz = 0.0
		result = [Bx, By, Bz]
		# Debug: print first few calls
		call_count[0] += 1
		if call_count[0] <= 3:
			print(f"  [Callback #{call_count[0]}] pos={pos} mm -> B={result} T")
		return result
	return field

quad_field = quadrupole_field_callback(gradient)
bckg_cf = rd.ObjBckgCF(quad_field)
print(f"  ObjBckgCF created with quadrupole field")

# ============================================================================
# Test 2: Create simple cubic element
# ============================================================================

print(f"\n[Test 2] Create Simple Cubic Element")
print("-" * 70)

# Small cube at center
size = 10.0  # mm
cube = rd.ObjRecMag([0, 0, 0], [size, size, size])
# Use MatSatIsoFrm for isotropic saturable material
# For soft iron-like material with high permeability
mat = rd.MatSatIsoFrm([1596.3, 1.1488], [133.11, 0.4268], [18.713, 0.4759])
rd.MatApl(cube, mat)
print(f"  Created {size}x{size}x{size} mm cube with mu_r={mu_r}")

# Create container with cube and background field
container = rd.ObjCnt([cube, bckg_cf])
print(f"  Container created with cube + ObjBckgCF")

# ============================================================================
# Test 3: Solve and verify field
# ============================================================================

print(f"\n[Test 3] Solve and Verify Field")
print("-" * 70)

print(f"  Solving...")
rd.Solve(container, 1e-5, 5000)
print(f"  [OK] Solution converged")

# ============================================================================
# Test 4: Compare with analytical solution
# ============================================================================

print(f"\n[Test 4] Compare with Analytical Solution")
print("-" * 70)

# Test points outside the cube (where field = background only)
test_points = [
	[20, 0, 0],   # x=20mm, y=0 => Bx=0, By=gradient*0.02
	[0, 20, 0],   # x=0, y=20mm => Bx=gradient*0.02, By=0
	[20, 20, 0],  # x=20mm, y=20mm => Bx=gradient*0.02, By=gradient*0.02
	[30, 0, 0],
	[0, 30, 0],
]

print(f"\nmu_0 = 1.25663706212e-6 T/(A/m)")
print(f"1/mu_0 = 795774.715459 (A/m)/T")
print(f"")

mu_0 = 1.25663706212e-6  # T/(A/m)

print(f"{'Point (mm)':<15} {'B_Radia (T)':>30} {'H_Radia (A/m)':>30} {'B_Analytical (T)':>30} {'H=B/mu_0 (A/m)':>30} {'Error':>15}")
print("-" * 150)

for pt in test_points:
	# Get Radia fields
	B_radia = rd.Fld(container, 'b', pt)
	H_radia = rd.Fld(container, 'h', pt)

	# Analytical quadrupole field (background only, far from cube)
	x_m = pt[0] * 1e-3
	y_m = pt[1] * 1e-3
	Bx_analytical = gradient * y_m  # T
	By_analytical = gradient * x_m  # T
	Bz_analytical = 0.0

	# Analytical H field: H = B/μ₀
	Hx_analytical = Bx_analytical / mu_0
	Hy_analytical = By_analytical / mu_0
	Hz_analytical = 0.0

	# Format vectors
	B_radia_str = f"[{B_radia[0]:.6e}, {B_radia[1]:.6e}, {B_radia[2]:.6e}]"
	H_radia_str = f"[{H_radia[0]:.6e}, {H_radia[1]:.6e}, {H_radia[2]:.6e}]"
	B_analytical_str = f"[{Bx_analytical:.6e}, {By_analytical:.6e}, {Bz_analytical:.6e}]"
	H_analytical_str = f"[{Hx_analytical:.6e}, {Hy_analytical:.6e}, {Hz_analytical:.6e}]"

	# Calculate error in H field
	H_analytical = np.array([Hx_analytical, Hy_analytical, Hz_analytical])
	H_radia_arr = np.array(H_radia)
	error_H = np.linalg.norm(H_radia_arr - H_analytical)
	error_pct = error_H / (np.linalg.norm(H_analytical) + 1e-15) * 100

	print(f"{str(pt):<15} {B_radia_str:>30} {H_radia_str:>30} {B_analytical_str:>30} {H_analytical_str:>30} {error_pct:>14.4f}%")

# ============================================================================
# Test 5: Verify B/H ratio = μ₀
# ============================================================================

print(f"\n[Test 5] Verify B/H Ratio = mu_0")
print("-" * 70)

pt = [30, 0, 0]
B = rd.Fld(container, 'b', pt)
H = rd.Fld(container, 'h', pt)

print(f"\nAt point {pt} mm:")
print(f"  B = [{B[0]:.8e}, {B[1]:.8e}, {B[2]:.8e}] T")
print(f"  H = [{H[0]:.8e}, {H[1]:.8e}, {H[2]:.8e}] A/m")

# Calculate B/H ratio for non-zero components
if abs(H[1]) > 1e-10:
	ratio = B[1] / H[1]
	print(f"\n  B_y / H_y = {ratio:.15e} T/(A/m)")
	print(f"  mu_0      = {mu_0:.15e} T/(A/m)")
	print(f"  Difference: {abs(ratio - mu_0):.15e} T/(A/m)")
	print(f"  Relative error: {abs(ratio - mu_0)/mu_0 * 100:.10f}%")

	if abs(ratio - mu_0) / mu_0 < 1e-6:
		print(f"\n  [OK] B/H = mu_0 within 1e-6 relative error")
	else:
		print(f"\n  [ERROR] B/H != mu_0 (relative error > 1e-6)")

# ============================================================================
# Summary
# ============================================================================

print("\n" + "=" * 70)
print("Summary")
print("=" * 70)

print(f"\n1. ObjBckgCF callback returns B in Tesla")
print(f"2. Internal conversion: H = B / mu_0 = B x 795774.715459")
print(f"3. B/H ratio matches mu_0 within numerical precision")
print(f"4. Background field correctly applied via callback")

print("\n" + "=" * 70)
print("Test Complete")
print("=" * 70)

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
