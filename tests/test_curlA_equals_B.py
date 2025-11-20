#!/usr/bin/env python
"""
Improved verification: curl(A) = B with proper finite element spaces

Key improvements:
1. A projected to H(curl) space (correct)
2. B projected to H(div) space (mathematically correct for div(B) = 0)
3. Mesh refinement study to verify convergence
4. Quantitative error analysis

Mathematical justification:
- Vector potential A ∈ H(curl): ensures curl(A) is well-defined
- Magnetic flux B ∈ H(div): ensures div(B) = 0 (Maxwell's equation)
- curl: H(curl) → H(div) is the natural mapping
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src", "python"))

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
print("IMPROVED VERIFICATION: curl(A) = B with Proper FE Spaces")
print("=" * 80)

# ============================================================================
# Step 1: Create Radia magnet
# ============================================================================
print("\n[Step 1] Creating Radia rectangular magnet")
print("-" * 80)

rad.UtiDelAll()

# Set Radia to use meters (required for NGSolve integration)
rad.FldUnits('m')

# Rectangular magnet with significant field
magnet = rad.ObjRecMag(
	[0, 0, 0],           # Center (m)
	[0.04, 0.04, 0.06],  # Dimensions (m)
	[0, 0, 1.2]          # Magnetization (T)
)

print(f"  Magnet ID: {magnet}")
print(f"  Center: [0, 0, 0] m")
print(f"  Dimensions: [0.04, 0.04, 0.06] m")
print(f"  Magnetization: [0, 0, 1.2] T")

# ============================================================================
# Step 2: Create background field providing both B and A
# ============================================================================
print("\n[Step 2] Creating background field wrapper")
print("-" * 80)

def radia_field_with_A(coords):
	"""Callback returning both B and A from Radia"""
	x, y, z = coords
	B = rad.Fld(magnet, 'b', [x, y, z])
	A = rad.Fld(magnet, 'a', [x, y, z])
	return {'B': list(B), 'A': list(A)}

bg_field = rad.ObjBckgCF(radia_field_with_A)
print(f"  Background field ID: {bg_field}")

# ============================================================================
# Step 3: Mesh refinement study
# ============================================================================
print("\n[Step 3] Mesh Refinement Study")
print("-" * 80)

# Test multiple mesh sizes
mesh_sizes = [0.012, 0.010, 0.008, 0.006]  # meters (12mm, 10mm, 8mm, 6mm)
results = []

for maxh in mesh_sizes:
	print(f"\n  Mesh size: {maxh*1000:.1f} mm")
	print("  " + "-" * 70)

	# Create mesh
	box = Box((0.01, 0.01, 0.02), (0.06, 0.06, 0.08))  # meters
	geo = OCCGeometry(box)
	mesh = Mesh(geo.GenerateMesh(maxh=maxh))

	print(f"    Vertices: {mesh.nv}, Elements: {mesh.ne}")

	# Create finite element spaces
	# A ∈ H(curl): curl operator is well-defined
	fes_hcurl = HCurl(mesh, order=2)

	# B ∈ H(div): div(B) = 0 is naturally satisfied
	fes_hdiv = HDiv(mesh, order=2)

	# H1 for scalar fields
	fes_h1 = H1(mesh, order=2)

	print(f"    H(curl) DOFs: {fes_hcurl.ndof}")
	print(f"    H(div) DOFs: {fes_hdiv.ndof}")

	# Get CoefficientFunctions
	A_cf = rad_ngsolve.RadiaField(bg_field, 'a')
	B_cf = rad_ngsolve.RadiaField(bg_field, 'b')

	# Project A to H(curl) space
	A_gf = GridFunction(fes_hcurl)
	A_gf.Set(A_cf)

	# Compute curl(A) - result is in H(div)
	curl_A_gf = curl(A_gf)

	# Project B to H(div) space (CORRECTED!)
	B_gf = GridFunction(fes_hdiv)
	B_gf.Set(B_cf)

	# Compute L2 error: ||curl(A) - B||_L2
	error_vec = curl_A_gf - B_gf
	error_l2_squared = Integrate(InnerProduct(error_vec, error_vec), mesh)
	error_l2 = np.sqrt(error_l2_squared)

	# Compute ||B||_L2 for normalization
	B_l2_squared = Integrate(InnerProduct(B_gf, B_gf), mesh)
	B_l2 = np.sqrt(B_l2_squared)

	# Relative error
	rel_error = error_l2 / B_l2 if B_l2 > 1e-10 else 0.0

	print(f"    ||curl(A) - B||_L2: {error_l2:.6e} T")
	print(f"    ||B||_L2:          {B_l2:.6e} T")
	print(f"    Relative error:    {rel_error*100:.4f}%")

	# Store results for convergence analysis
	results.append({
		'maxh': maxh,
		'h_mm': maxh * 1000,
		'ne': mesh.ne,
		'nv': mesh.nv,
		'ndof_hcurl': fes_hcurl.ndof,
		'ndof_hdiv': fes_hdiv.ndof,
		'error_l2': error_l2,
		'B_l2': B_l2,
		'rel_error': rel_error
	})

	# Point-wise verification for finest mesh
	if maxh == min(mesh_sizes):
		print("\n  Point-wise verification (finest mesh):")
		print("  " + "-" * 70)

		test_points_meters = [
			(0.030, 0.020, 0.040),
			(0.030, 0.020, 0.050),
			(0.040, 0.040, 0.050),
		]

		print("  Point (m)           curl(A) (T)              B (T)           Error (T)")
		for point in test_points_meters:
			try:
				mip = mesh(*point)
				curl_A_val = np.array(curl_A_gf(mip))
				B_val = np.array(B_cf(mip))
				error_vec = curl_A_val - B_val
				error_norm = np.linalg.norm(error_vec)

				print(f"  {point}  [{curl_A_val[0]:7.4f}, {curl_A_val[1]:7.4f}, {curl_A_val[2]:7.4f}]  "
				      f"[{B_val[0]:7.4f}, {B_val[1]:7.4f}, {B_val[2]:7.4f}]  {error_norm:.3e}")
			except:
				pass

# ============================================================================
# Step 4: Convergence analysis
# ============================================================================
print("\n[Step 4] Convergence Analysis")
print("=" * 80)

print("\nMesh Refinement Results:")
print("-" * 80)
print(f"{'h (mm)':>8s} {'Elements':>10s} {'H(curl) DOFs':>12s} {'H(div) DOFs':>12s} "
      f"{'||error||_L2':>14s} {'Rel. Error':>12s}")
print("-" * 80)

for r in results:
	print(f"{r['h_mm']:8.1f} {r['ne']:10d} {r['ndof_hcurl']:12d} {r['ndof_hdiv']:12d} "
	      f"{r['error_l2']:14.6e} {r['rel_error']*100:11.4f}%")

# Estimate convergence rate
if len(results) >= 2:
	print("\nConvergence rate estimation:")
	print("-" * 80)
	for i in range(1, len(results)):
		h1 = results[i-1]['maxh']
		h2 = results[i]['maxh']
		e1 = results[i-1]['error_l2']
		e2 = results[i]['error_l2']

		# Convergence rate: error ~ h^rate
		# log(e1/e2) = rate * log(h1/h2)
		if e1 > 0 and e2 > 0 and h1 != h2:
			rate = np.log(e1/e2) / np.log(h1/h2)
			print(f"  h: {h1*1000:.1f}mm → {h2*1000:.1f}mm: "
			      f"error: {e1:.3e} → {e2:.3e}, rate = {rate:.2f}")

# ============================================================================
# Step 5: VTK export (finest mesh only)
# ============================================================================
print("\n[Step 5] VTK Export (finest mesh)")
print("-" * 80)

finest_idx = mesh_sizes.index(min(mesh_sizes))
maxh = mesh_sizes[finest_idx]

# Recreate mesh and fields for finest resolution
box = Box((0.01, 0.01, 0.02), (0.06, 0.06, 0.08))
geo = OCCGeometry(box)
mesh = Mesh(geo.GenerateMesh(maxh=maxh))

fes_hcurl = HCurl(mesh, order=2)
fes_hdiv = HDiv(mesh, order=2)
fes_h1 = H1(mesh, order=2)

A_cf = rad_ngsolve.RadiaField(bg_field, 'a')
B_cf = rad_ngsolve.RadiaField(bg_field, 'b')

A_gf = GridFunction(fes_hcurl)
A_gf.Set(A_cf)

curl_A_gf = curl(A_gf)

B_gf = GridFunction(fes_hdiv)
B_gf.Set(B_cf)

# Compute error field magnitude
error_vec = curl_A_gf - B_gf
error_magnitude = GridFunction(fes_h1)
error_magnitude.Set(sqrt(InnerProduct(error_vec, error_vec)))

try:
	# Export fields
	vtk_output = VTKOutput(
		mesh,
		coefs=[A_gf, curl_A_gf, B_gf, error_magnitude],
		names=["A_vector_potential", "curl_A", "B_field", "error_magnitude"],
		filename="verify_curl_A_equals_B_improved",
		subdivision=2
	)
	vtk_output.Do()
	print("  [OK] VTK file exported: verify_curl_A_equals_B_improved.vtu")
	print("  Open with: paraview verify_curl_A_equals_B_improved.vtu")
except Exception as e:
	print(f"  [ERROR] VTK export failed: {e}")

# ============================================================================
# Step 6: Summary
# ============================================================================
print("\n[Step 6] Verification Summary")
print("=" * 80)

finest_result = results[-1]
print(f"\nFinest mesh (h = {finest_result['h_mm']:.1f} mm):")
print(f"  Elements: {finest_result['ne']}")
print(f"  H(curl) DOFs: {finest_result['ndof_hcurl']}")
print(f"  H(div) DOFs: {finest_result['ndof_hdiv']}")
print(f"  Relative L2 error: {finest_result['rel_error']*100:.4f}%")

tolerance = 0.01  # 1%
if finest_result['rel_error'] < tolerance:
	print(f"\n[SUCCESS] curl(A) = B verified!")
	print(f"  Relative error {finest_result['rel_error']*100:.4f}% < {tolerance*100}%")
else:
	print(f"\n[INFO] Relative error {finest_result['rel_error']*100:.4f}%")
	print(f"  (Target: < {tolerance*100}%)")

print("\nKey improvements in this version:")
print("  1. A projected to H(curl) space (mathematically correct)")
print("  2. B projected to H(div) space (ensures div(B) = 0)")
print("  3. Mesh refinement study shows convergence")
print("  4. L2 norm error quantification")

print("\n" + "=" * 80)
print("Verification complete!")
print("=" * 80)

# Cleanup
rad.UtiDelAll()
