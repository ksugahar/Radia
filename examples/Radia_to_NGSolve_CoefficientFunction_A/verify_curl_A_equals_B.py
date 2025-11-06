#!/usr/bin/env python
"""
Rigorous verification that curl(A) = B for Radia background field with vector potential

This example performs comprehensive verification:
1. Creates a Radia background field providing both B and A
2. Extracts A as NGSolve CoefficientFunction
3. Computes curl(A) numerically in NGSolve
4. Compares curl(A) with B at multiple points
5. Generates detailed VTK outputs for visualization:
   - Vector potential A field
   - curl(A) computed from A
   - B field from Radia
   - Error field: |curl(A) - B|
6. Statistical error analysis
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
	print("ERROR: NGSolve not available. This example requires NGSolve.")
	NGSOLVE_AVAILABLE = False
	sys.exit(1)

import numpy as np

print("=" * 80)
print("RIGOROUS VERIFICATION: curl(A) = B for Radia Background Field")
print("=" * 80)

# ============================================================================
# Step 1: Create Radia magnet with known field
# ============================================================================
print("\n[Step 1] Creating Radia rectangular magnet")
print("-" * 80)

rad.UtiDelAll()

# Create a rectangular magnet with significant field
# Using larger dimensions and magnetization for better SNR in verification
magnet = rad.ObjRecMag(
	[0, 0, 0],           # Center (mm)
	[40, 40, 60],        # Dimensions (mm) - larger for stronger field
	[0, 0, 1.2]          # Magnetization (T) - NdFeB typical value
)

print(f"  Magnet ID: {magnet}")
print(f"  Center: [0, 0, 0] mm")
print(f"  Dimensions: [40, 40, 60] mm")
print(f"  Magnetization: [0, 0, 1.2] T")

# Calculate field at a reference point
ref_point = [0, 0, 50]  # mm
B_ref = rad.Fld(magnet, 'b', ref_point)
print(f"  Reference B at [0, 0, 50] mm: [{B_ref[0]:.6f}, {B_ref[1]:.6f}, {B_ref[2]:.6f}] T")
print(f"  |B| = {np.linalg.norm(B_ref):.6f} T")

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

# Test at a point to verify callback works
test_point = [30, 20, 40]  # mm
result = rad.Fld(bg_field, 'ba', test_point)
B_test = np.array(result[0:3])
A_test = np.array(result[3:6])

print(f"  Callback test at [{test_point[0]}, {test_point[1]}, {test_point[2]}] mm:")
print(f"    B = [{B_test[0]:.6e}, {B_test[1]:.6e}, {B_test[2]:.6e}] T")
print(f"    A = [{A_test[0]:.6e}, {A_test[1]:.6e}, {A_test[2]:.6e}] T*m")

# ============================================================================
# Step 3: Create NGSolve mesh
# ============================================================================
print("\n[Step 3] Creating NGSolve mesh")
print("-" * 80)

# Create mesh covering region with significant field
# Region: 10mm to 60mm in x,y,z (0.01m to 0.06m)
box = Box((0.01, 0.01, 0.02), (0.06, 0.06, 0.08))  # meters
geo = OCCGeometry(box)
mesh = Mesh(geo.GenerateMesh(maxh=0.008))  # 8mm max element size

print(f"  Mesh created:")
print(f"    Vertices: {mesh.nv}")
print(f"    Elements: {mesh.ne}")
print(f"    Bounding box: [{mesh.GetBoundaries()[0]} to {mesh.GetBoundaries()[1]}]")

# ============================================================================
# Step 4: Extract A and B as CoefficientFunctions
# ============================================================================
print("\n[Step 4] Creating CoefficientFunctions")
print("-" * 80)

# Get vector potential A from Radia as CoefficientFunction
A_cf = rad_ngsolve.RadiaField(bg_field, 'a')
print(f"  A_cf created: {type(A_cf)}")

# Get B directly from Radia
B_cf = rad_ngsolve.RadiaField(bg_field, 'b')
print(f"  B_cf created: {type(B_cf)}")

# ============================================================================
# Step 5: Create GridFunctions for visualization
# ============================================================================
print("\n[Step 5] Creating GridFunctions for VTK export")
print("-" * 80)

# Create H(curl) finite element space
fes_hcurl = HCurl(mesh, order=2)
print(f"  H(curl) FE space: {fes_hcurl.ndof} DOFs")

# Create H1 space for scalar error field
fes_h1 = H1(mesh, order=2)
print(f"  H1 FE space: {fes_h1.ndof} DOFs")

# Project A onto H(curl) space
print("  Projecting A to GridFunction...")
A_gf = GridFunction(fes_hcurl)
A_gf.Set(A_cf)

# Compute curl(A) in FE space
print("  Computing curl(A) in FE space...")
curl_A_gf = curl(A_gf)

# Project B for comparison
print("  Projecting B to GridFunction...")
B_gf = GridFunction(fes_hcurl)
B_gf.Set(B_cf)

# Compute error field: |curl(A) - B|
print("  Computing error field...")
error_cf = sqrt((curl_A_gf[0]-B_gf[0])**2 + (curl_A_gf[1]-B_gf[1])**2 + (curl_A_gf[2]-B_gf[2])**2)
error_gf = GridFunction(fes_h1)
error_gf.Set(error_cf)

# Compute field magnitudes
print("  Computing field magnitudes...")
A_magnitude = GridFunction(fes_h1)
A_magnitude.Set(sqrt(A_cf[0]**2 + A_cf[1]**2 + A_cf[2]**2))

B_magnitude = GridFunction(fes_h1)
B_magnitude.Set(sqrt(B_cf[0]**2 + B_cf[1]**2 + B_cf[2]**2))

curl_A_magnitude = GridFunction(fes_h1)
curl_A_magnitude.Set(sqrt(curl_A_gf[0]**2 + curl_A_gf[1]**2 + curl_A_gf[2]**2))

print("  All GridFunctions created")

# ============================================================================
# Step 6: Point-wise verification at test points
# ============================================================================
print("\n[Step 6] Point-wise verification: curl(A) vs B")
print("-" * 80)

# Define comprehensive test points (in meters)
test_points_meters = [
	(0.030, 0.020, 0.030),  # Near magnet surface
	(0.030, 0.020, 0.040),  # Near center
	(0.030, 0.020, 0.050),  # Center region
	(0.020, 0.020, 0.040),  # Corner region
	(0.040, 0.040, 0.050),  # Off-center
	(0.050, 0.030, 0.060),  # Edge region
	(0.035, 0.035, 0.045),  # Diagonal
	(0.025, 0.045, 0.055),  # Asymmetric
]

print("\n  Point (m)                curl(A) (T)                    B (T)                       Error (T)")
print("  " + "-" * 100)

errors = []
for point in test_points_meters:
	try:
		mip = mesh(*point)

		curl_A_val = np.array(curl_A_gf(mip))
		B_val = np.array(B_cf(mip))

		error_vec = curl_A_val - B_val
		error_norm = np.linalg.norm(error_vec)
		errors.append(error_norm)

		# Relative error
		B_norm = np.linalg.norm(B_val)
		rel_error = error_norm / B_norm if B_norm > 1e-10 else 0.0

		print(f"  {point}  [{curl_A_val[0]:9.6f}, {curl_A_val[1]:9.6f}, {curl_A_val[2]:9.6f}]  "
		      f"[{B_val[0]:9.6f}, {B_val[1]:9.6f}, {B_val[2]:9.6f}]  "
		      f"{error_norm:.3e} ({rel_error*100:.4f}%)")

	except Exception as e:
		print(f"  {point}  [Point outside mesh: {e}]")

# ============================================================================
# Step 7: Statistical error analysis
# ============================================================================
print("\n[Step 7] Statistical Error Analysis")
print("-" * 80)

if len(errors) > 0:
	errors = np.array(errors)
	mean_error = np.mean(errors)
	std_error = np.std(errors)
	max_error = np.max(errors)
	min_error = np.min(errors)

	print(f"  Number of test points: {len(errors)}")
	print(f"  Mean error:   {mean_error:.6e} T")
	print(f"  Std deviation: {std_error:.6e} T")
	print(f"  Min error:    {min_error:.6e} T")
	print(f"  Max error:    {max_error:.6e} T")

	# Verification criterion
	tolerance = 1e-4  # Tesla (0.1 mT)
	if max_error < tolerance:
		print(f"\n  [SUCCESS] curl(A) = B verified!")
		print(f"    Maximum error {max_error:.6e} T < tolerance {tolerance:.6e} T")
	else:
		print(f"\n  [WARNING] Large errors detected!")
		print(f"    Maximum error {max_error:.6e} T >= tolerance {tolerance:.6e} T")
else:
	print("  No valid test points evaluated")

# ============================================================================
# Step 8: Export VTK files
# ============================================================================
print("\n[Step 8] Exporting VTK visualization files")
print("-" * 80)

try:
	# Export vector fields (A, curl(A), B)
	vtk_vectors = VTKOutput(
		mesh,
		coefs=[A_gf, curl_A_gf, B_gf],
		names=["A_vector_potential", "curl_A", "B_field"],
		filename="vector_fields",
		subdivision=2
	)
	vtk_vectors.Do()
	print("  [OK] vector_fields.vtu exported (A, curl(A), B)")

	# Export scalar fields (magnitudes and error)
	vtk_scalars = VTKOutput(
		mesh,
		coefs=[A_magnitude, B_magnitude, curl_A_magnitude, error_gf],
		names=["A_magnitude", "B_magnitude", "curl_A_magnitude", "error_magnitude"],
		filename="scalar_fields",
		subdivision=2
	)
	vtk_scalars.Do()
	print("  [OK] scalar_fields.vtu exported (magnitudes and error)")

	# Export error field separately for detailed analysis
	vtk_error = VTKOutput(
		mesh,
		coefs=[error_gf],
		names=["curl_A_minus_B_error"],
		filename="error_field",
		subdivision=3  # Higher subdivision for error visualization
	)
	vtk_error.Do()
	print("  [OK] error_field.vtu exported (|curl(A) - B|)")

	print("\n  Visualization files created:")
	print("    - vector_fields.vtu: A, curl(A), B vector fields")
	print("    - scalar_fields.vtu: Field magnitudes and error")
	print("    - error_field.vtu: Detailed error distribution")
	print("\n  Open in ParaView:")
	print("    1. Load vector_fields.vtu")
	print("    2. Apply 'Glyph' filter to visualize vector fields")
	print("    3. Load error_field.vtu to see error distribution")
	print("    4. Use 'Threshold' filter to highlight regions with high error")

except Exception as e:
	print(f"  [ERROR] VTK export failed: {e}")
	import traceback
	traceback.print_exc()

# ============================================================================
# Step 9: Summary and recommendations
# ============================================================================
print("\n[Step 9] Verification Summary")
print("=" * 80)

if len(errors) > 0 and max_error < tolerance:
	print("STATUS: [VERIFICATION PASSED]")
	print(f"\ncurl(A) = B relationship verified to within {tolerance:.6e} T")
	print(f"Maximum error: {max_error:.6e} T")
	print(f"Mean error: {mean_error:.6e} T")
	print("\nThe vector potential A correctly represents the magnetic field B.")
	print("The relationship B = curl(A) is satisfied numerically.")
else:
	print("STATUS: [VERIFICATION INCOMPLETE]")
	if len(errors) == 0:
		print("\nNo valid test points could be evaluated.")
	else:
		print(f"\nMaximum error {max_error:.6e} T exceeds tolerance {tolerance:.6e} T")
	print("\nPossible causes:")
	print("  - Numerical differentiation error in curl computation")
	print("  - Mesh resolution insufficient")
	print("  - Field evaluation outside valid region")

print("\nNext steps:")
print("  1. Open VTK files in ParaView")
print("  2. Visualize curl(A) and B side-by-side")
print("  3. Check error_field.vtu for spatial distribution of errors")
print("  4. Verify that errors are concentrated in low-field regions")

print("\n" + "=" * 80)
print("Verification complete!")
print("=" * 80)

# Cleanup
rad.UtiDelAll()
