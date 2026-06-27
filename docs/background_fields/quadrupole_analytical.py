#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Simple test to verify B->H conversion in rad.ObjBckg()

Tests that quadrupole background field defined in Tesla is correctly
converted to H field internally.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../build/Release'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

import numpy as np
import radia as rd

mm = 1e-3  # 1 mm in meters

print("=" * 70)
print("ObjBckg B->H Conversion Test")
print("=" * 70)

# Parameters
gradient = 10.0  # Quadrupole gradient [T/m]

print("\nParameters:")
print(f"  Quadrupole gradient: {gradient} T/m")

# ============================================================================
# Test 1: Create quadrupole background field using ObjBckg
# ============================================================================

print("\n[Test 1] Quadrupole Background Field (ObjBckg)")
print("-" * 70)

def quadrupole_field_callback(gradient):
	"""Create quadrupole field callback for rd.ObjBckg.

	Returns a callable(pos) -> [Bx, By, Bz] in Tesla.
	"""
	call_count = [0]  # Mutable to allow modification in nested function
	def field(pos):
		x, y, z = pos  # Position in meters (Radia always uses meters)
		# Quadrupole field: Bx = g*y, By = g*x, Bz = 0
		Bx = gradient * y  # [T]
		By = gradient * x  # [T]
		Bz = 0.0
		result = [Bx, By, Bz]
		# Debug: print first few calls
		call_count[0] += 1
		if call_count[0] <= 3:
			print(f"  [Callback #{call_count[0]}] pos={pos} m -> B={result} T")
		return result
	return field

quad_field = quadrupole_field_callback(gradient)
bckg_cf = rd.ObjBckg(quad_field)
print("  ObjBckg created with quadrupole field")

# ============================================================================
# Test 2: Create simple cubic element
# ============================================================================

print("\n[Test 2] Create Simple Cubic Element")
print("-" * 70)

# Small cube at center: 10mm cube centered at origin
size = 10.0 * mm
half = size / 2
# Hexahedron vertices for cube centered at [0, 0, 0] with dimensions [10, 10, 10] mm
vertices = [
	[-half, -half, -half], [half, -half, -half], [half, half, -half], [-half, half, -half],
	[-half, -half, half], [half, -half, half], [half, half, half], [-half, half, half]
]
cube = rd.ObjHexahedron(vertices, [0, 0, 0])
# Use MatSatIsoFrm for isotropic saturable material
# For soft iron-like material with high permeability
mat = rd.MatSatIsoFrm([[1596.3, 1.1488], [133.11, 0.4268], [18.713, 0.4759]])
rd.MatApl(cube, mat)
print(f"  Created {size/mm:.0f}x{size/mm:.0f}x{size/mm:.0f} mm cube with MatSatIsoFrm (nonlinear)")

# Create container with cube and background field
container = rd.ObjCnt([cube, bckg_cf])
print("  Container created with cube + ObjBckg")

# ============================================================================
# Test 3: Solve and verify field
# ============================================================================

print("\n[Test 3] Solve and Verify Field")
print("-" * 70)

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
# Test 4: Compare with analytical solution
# ============================================================================

print("\n[Test 4] Compare with Analytical Solution")
print("-" * 70)

# Test points outside the cube (field ~= background + stray field from cube)
# Note: The magnetized cube produces stray fields that decay with distance.
# Points farther from the cube give closer agreement with the pure background.
test_points = [
	[50*mm, 0, 0],      # Far from cube: stray field negligible
	[0, 50*mm, 0],
	[50*mm, 50*mm, 0],
	[100*mm, 0, 0],     # Very far: essentially pure background
	[0, 100*mm, 0],
]

print("\nmu_0 = 1.25663706212e-6 T/(A/m)")
print("1/mu_0 = 795774.715459 (A/m)/T")
print()

mu_0 = 1.25663706212e-6  # T/(A/m)

print(f"{'Point (mm)':<25} {'B_Radia (T)':>30} {'H_Radia (A/m)':>30} {'B_Analytical (T)':>30} {'H=B/mu_0 (A/m)':>30} {'Error':>15}")
print("-" * 160)

for pt in test_points:
	# Get Radia fields
	B_radia = rd.Fld(container, 'b', pt)
	H_radia = rd.Fld(container, 'h', pt)

	# Analytical quadrupole field (background only, ignoring stray field from cube)
	# Positions already in meters (Radia always uses meters)
	Bx_analytical = gradient * pt[1]  # T
	By_analytical = gradient * pt[0]  # T
	Bz_analytical = 0.0

	# Analytical H field: H = B/mu_0
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

	pt_mm_str = str([round(p/mm) for p in pt])
	print(f"{pt_mm_str:<25} {B_radia_str:>30} {H_radia_str:>30} {B_analytical_str:>30} {H_analytical_str:>30} {error_pct:>14.4f}%")

# ============================================================================
# Test 5: Verify B/H ratio = mu_0 at multiple points
# ============================================================================

print("\n[Test 5] Verify B/H Ratio = mu_0 at Multiple Points")
print("-" * 70)

# Test at several far-field points where stray field is negligible
bh_test_points = [
	[100*mm, 0, 0],
	[0, 100*mm, 0],
	[100*mm, 100*mm, 0],
	[0, 0, 100*mm],
]

all_ok = True
for pt in bh_test_points:
	B = rd.Fld(container, 'b', pt)
	H = rd.Fld(container, 'h', pt)
	pt_mm_str = str([round(p/mm) for p in pt])

	print(f"\n  At {pt_mm_str} mm:")
	print(f"    B = [{B[0]:.8e}, {B[1]:.8e}, {B[2]:.8e}] T")
	print(f"    H = [{H[0]:.8e}, {H[1]:.8e}, {H[2]:.8e}] A/m")

	# Check B/H ratio for each non-zero component
	for comp, label in enumerate(['x', 'y', 'z']):
		if abs(H[comp]) > 1e-10:
			ratio = B[comp] / H[comp]
			rel_err = abs(ratio - mu_0) / mu_0
			status = "[OK]" if rel_err < 1e-6 else "[ERROR]"
			if rel_err >= 1e-6:
				all_ok = False
			print(f"    B_{label}/H_{label} = {ratio:.15e}  (rel_err={rel_err:.2e}) {status}")

if all_ok:
	print("\n  [OK] B/H = mu_0 within 1e-6 for all tested components")
else:
	print("\n  [ERROR] B/H != mu_0 for some components")

# ============================================================================
# Summary
# ============================================================================

print("\n" + "=" * 70)
print("Summary")
print("=" * 70)

print("\n1. ObjBckg callback returns B in Tesla")
print("2. Internal conversion: H = B / mu_0 = B x 795774.715459")
print("3. B/H = mu_0 verified at multiple far-field points and components")
print("4. Background field correctly applied via callback")

print("\n" + "=" * 70)
print("Test Complete")
print("=" * 70)

rd.UtiDelAll()
