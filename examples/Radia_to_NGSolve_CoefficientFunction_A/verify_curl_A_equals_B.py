#!/usr/bin/env python
"""
Verify that curl(A) = B for Radia background field with vector potential

This example demonstrates:
1. Creating a Radia background field that provides both B and A
2. Passing the vector potential A to NGSolve as a CoefficientFunction
3. Computing curl(A) in NGSolve
4. Verifying that curl(A) = B
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "build", "Release"))

import radia as rad
try:
	from ngsolve import *
	from netgen.occ import *
	import rad_ngsolve
	NGSOLVE_AVAILABLE = True
except ImportError:
	print("NGSolve not available. This example requires NGSolve.")
	NGSOLVE_AVAILABLE = False
	sys.exit(1)

import numpy as np

print("=" * 80)
print("Verification: curl(A) = B for Radia Background Field")
print("=" * 80)

# ============================================================================
# Step 1: Create Radia magnet with known field
# ============================================================================
print("\n[Step 1] Creating Radia rectangular magnet")
print("-" * 80)

rad.UtiDelAll()

# Create a simple rectangular magnet
# Center at origin, dimensions 20x20x30 mm, magnetization in Z direction
magnet = rad.ObjRecMag(
	[0, 0, 0],           # Center (mm)
	[20, 20, 30],        # Dimensions (mm)
	[0, 0, 1000]         # Magnetization (A/m)
)

print(f"  Magnet ID: {magnet}")
print(f"  Center: [0, 0, 0] mm")
print(f"  Dimensions: [20, 20, 30] mm")
print(f"  Magnetization: [0, 0, 1000] A/m")

# ============================================================================
# Step 2: Create Radia background field providing both B and A
# ============================================================================
print("\n[Step 2] Creating background field wrapper for B and A")
print("-" * 80)

def radia_field_with_A(coords):
	"""
	Callback that returns both B field and vector potential A from Radia

	Args:
		coords: [x, y, z] in mm

	Returns:
		dict: {'B': [Bx, By, Bz], 'A': [Ax, Ay, Az]}
	"""
	x, y, z = coords

	# Get B field from Radia (returns Tesla)
	B = rad.Fld(magnet, 'b', [x, y, z])

	# Get A field from Radia (returns T*m)
	A = rad.Fld(magnet, 'a', [x, y, z])

	return {'B': list(B), 'A': list(A)}

# Create background field object
bg_field = rad.ObjBckgCF(radia_field_with_A)
print(f"  Background field ID: {bg_field}")

# Test at a point
test_point = [0.025, 0.015, 0.040]  # meters
test_point_mm = [p * 1000 for p in test_point]  # convert to mm

result = rad.Fld(bg_field, 'ba', test_point_mm)
B_radia = np.array(result[0:3])
A_radia = np.array(result[3:6])

print(f"  Test point: {test_point} m")
print(f"  B from Radia: {B_radia} T")
print(f"  A from Radia: {A_radia} T*m")

# ============================================================================
# Step 3: Create NGSolve mesh and get A as CoefficientFunction
# ============================================================================
print("\n[Step 3] Creating NGSolve mesh and CoefficientFunction")
print("-" * 80)

# Create a simple mesh
box = Box((0.01, 0.01, 0.02), (0.04, 0.03, 0.06))  # meters
geo = OCCGeometry(box)
mesh = Mesh(geo.GenerateMesh(maxh=0.01))

print(f"  Mesh created: {mesh.nv} vertices, {mesh.ne} elements")

# Get vector potential A from Radia as CoefficientFunction
A_cf = rad_ngsolve.RadiaField(bg_field, 'a')

print(f"  A as CoefficientFunction: {type(A_cf)}")

# ============================================================================
# Step 4: Compute curl(A) in NGSolve and compare with B
# ============================================================================
print("\n[Step 4] Computing curl(A) and comparing with B")
print("-" * 80)

# Compute curl(A) using NGSolve
curl_A_cf = curl(A_cf)

# Get B directly from Radia
B_cf = rad_ngsolve.RadiaField(bg_field, 'b')

# Evaluate at several test points
test_points_meters = [
	(0.025, 0.015, 0.030),
	(0.025, 0.015, 0.040),
	(0.020, 0.020, 0.035),
	(0.030, 0.025, 0.045),
]

print("\n  Point (m)                curl(A) (T)              B (T)                 Error")
print("  " + "-" * 76)

max_error = 0.0
for point in test_points_meters:
	mip = mesh(*point)

	curl_A_val = curl_A_cf(mip)
	B_val = B_cf(mip)

	error = np.linalg.norm(np.array(curl_A_val) - np.array(B_val))
	max_error = max(max_error, error)

	print(f"  {point}  [{curl_A_val[0]:8.5f}, {curl_A_val[1]:8.5f}, {curl_A_val[2]:8.5f}]  "
	      f"[{B_val[0]:8.5f}, {B_val[1]:8.5f}, {B_val[2]:8.5f}]  {error:.2e}")

# ============================================================================
# Step 5: Verification result
# ============================================================================
print("\n[Step 5] Verification Result")
print("-" * 80)

tolerance = 1e-4  # Tesla
if max_error < tolerance:
	print(f"  [SUCCESS] curl(A) = B verified!")
	print(f"  Maximum error: {max_error:.6e} T (tolerance: {tolerance:.6e} T)")
else:
	print(f"  [FAILED] curl(A) != B")
	print(f"  Maximum error: {max_error:.6e} T (tolerance: {tolerance:.6e} T)")

# ============================================================================
# Step 6: Visualization (optional - create VTK output)
# ============================================================================
print("\n[Step 6] Creating visualization")
print("-" * 80)

# Create finite element space
fes = HCurl(mesh, order=1)
print(f"  FE space created: {fes.ndof} DOFs")

# Project A onto FE space
A_gf = GridFunction(fes)
A_gf.Set(A_cf)

# Compute curl in FE space
curl_A_gf = curl(A_gf)

# Project B for comparison
B_gf = GridFunction(fes)
B_gf.Set(B_cf)

print("  GridFunctions created for A, curl(A), and B")

# Export to VTK
try:
	vtk = VTKOutput(mesh, coefs=[A_gf, curl_A_gf, B_gf],
	                names=["A_vector_potential", "curl_A", "B_field"],
	                filename="curl_A_verification",
	                subdivision=2)
	vtk.Do()
	print("  [OK] VTK output saved: curl_A_verification.vtu")
	print("       Open in ParaView to visualize A, curl(A), and B")
except Exception as e:
	print(f"  [INFO] VTK export failed (optional): {e}")

print("\n" + "=" * 80)
print("Verification complete!")
print("=" * 80)

# Cleanup
rad.UtiDelAll()
