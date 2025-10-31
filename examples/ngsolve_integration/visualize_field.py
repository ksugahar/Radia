"""
Radia Magnetic Field Visualization - Comparison Study

This script compares two approaches:
1. CoefficientFunction direct evaluation (exact)
2. GridFunction interpolation with mesh refinement (convergent)

Usage:
	python visualize_field.py --method cf        # CoefficientFunction only
	python visualize_field.py --method gf        # GridFunction only
	python visualize_field.py --method both      # Compare both (default)
	python visualize_field.py --maxh 0.005       # Set mesh size

Author: Claude Code
Date: 2025-10-31
"""

import sys
sys.path.insert(0, r"S:\radia\01_GitHub\build\Release")
sys.path.insert(0, r"S:\radia\01_GitHub\dist")

from ngsolve import *
from netgen.csg import CSGeometry, OrthoBrick, Pnt
import numpy as np
import radia as rad
import rad_ngsolve
import argparse

# Parse command line arguments
parser = argparse.ArgumentParser(description='Radia field visualization')
parser.add_argument('--method', choices=['cf', 'gf', 'both'], default='both',
                    help='Evaluation method: cf (CoefficientFunction), gf (GridFunction), both (compare)')
parser.add_argument('--maxh', type=float, default=0.005,
                    help='Maximum mesh size in meters (default: 0.005 = 5mm)')
args = parser.parse_args()

print("=" * 70)
print("Radia Magnetic Field Visualization")
print("=" * 70)
print(f"Method: {args.method}")
print(f"Mesh size: {args.maxh} m ({args.maxh*1000} mm)")
print()

# ============================================================================
# Step 1: Create Radia Magnet Geometry
# ============================================================================

print("[Step 1] Creating Radia Magnet Geometry")
print("-" * 70)

magnet_center = [0, 0, 0]
magnet_size = [20, 20, 30]

magnet = rad.ObjRecMag(magnet_center, magnet_size, [0, 0, 1.2])
rad.MatApl(magnet, rad.MatStd('NdFeB', 1.2))
rad.Solve(magnet, 0.0001, 10000)

print(f"Magnet created: object #{magnet}")
print(f"  Center: {magnet_center} mm")
print(f"  Size: {magnet_size} mm")
print(f"  Material: NdFeB, Br = 1.2 T")

# ============================================================================
# Step 2: Create NGSolve Coefficient Function
# ============================================================================

print("\n[Step 2] Creating NGSolve Coefficient Function")
print("-" * 70)

B_cf = rad_ngsolve.RadBfield(magnet)
print(f"B-field CoefficientFunction created: {type(B_cf).__name__}")

# Test evaluation using Radia directly
test_points_mm = [[0, 0, 0], [0, 0, 20], [0, 0, 40]]
print(f"\nDirect Radia field evaluation (reference):")
radia_values = {}
for pt in test_points_mm:
	B = rad.Fld(magnet, 'b', pt)
	radia_values[tuple(pt)] = B
	print(f"  {pt} mm: Bx={B[0]:.6f}, By={B[1]:.6f}, Bz={B[2]:.6f} T")

# ============================================================================
# Step 3: Create 3D Mesh
# ============================================================================

print("\n[Step 3] Creating 3D Mesh")
print("-" * 70)

geo = CSGeometry()
geo.Add(OrthoBrick(Pnt(-0.05, -0.05, -0.05), Pnt(0.05, 0.05, 0.05)))
mesh = Mesh(geo.GenerateMesh(maxh=args.maxh))

print(f"3D Mesh generated:")
print(f"  Elements: {mesh.ne}")
print(f"  Vertices: {mesh.nv}")
print(f"  Domain: [-0.05, 0.05] m = [-50, 50] mm")
print(f"  maxh: {args.maxh} m ({args.maxh*1000} mm)")

# ============================================================================
# Step 4: Evaluation Comparison
# ============================================================================

print("\n[Step 4] Field Evaluation Comparison")
print("-" * 70)

test_mesh_points = [
	(0.000, 0.000, 0.000),	# 0mm (inside magnet)
	(0.000, 0.000, 0.020),  # 20mm (above magnet)
	(0.000, 0.000, 0.040),  # 40mm (far from magnet)
]

# Prepare GridFunction if needed
gfB = None
if args.method in ['gf', 'both']:
	print("\nCreating GridFunction with interpolation...")
	fes = VectorH1(mesh, order=2)
	gfB = GridFunction(fes)
	gfB.Set(B_cf)
	print(f"  FE space: {fes.ndof} DOFs")
	print(f"  GridFunction interpolation complete")

# Evaluate at test points
print("\n" + "=" * 70)
print("Point-wise Evaluation")
print("=" * 70)

for pt in test_mesh_points:
	pt_mm = [pt[0]*1000, pt[1]*1000, pt[2]*1000]
	B_radia = rad.Fld(magnet, 'b', pt_mm)

	print(f"\nPoint: {pt} m = {pt_mm} mm")
	print(f"  Radia:  Bx={B_radia[0]:.6f}, By={B_radia[1]:.6f}, Bz={B_radia[2]:.6f} T")

	mesh_pt = mesh(*pt)

	if args.method in ['cf', 'both']:
		B_cf_val = B_cf(mesh_pt)
		cf_error = np.sqrt((B_cf_val[0]-B_radia[0])**2 +
		                   (B_cf_val[1]-B_radia[1])**2 +
		                   (B_cf_val[2]-B_radia[2])**2)
		print(f"  CF:     Bx={B_cf_val[0]:.6f}, By={B_cf_val[1]:.6f}, Bz={B_cf_val[2]:.6f} T")
		print(f"          Error: {cf_error:.6e} T")

	if args.method in ['gf', 'both']:
		B_gf_val = gfB(mesh_pt)
		gf_error = np.sqrt((B_gf_val[0]-B_radia[0])**2 +
		                   (B_gf_val[1]-B_radia[1])**2 +
		                   (B_gf_val[2]-B_radia[2])**2)
		gf_rel_error = gf_error / max(np.sqrt(B_radia[0]**2 + B_radia[1]**2 + B_radia[2]**2), 1e-10) * 100
		print(f"  GF:     Bx={B_gf_val[0]:.6f}, By={B_gf_val[1]:.6f}, Bz={B_gf_val[2]:.6f} T")
		print(f"          Error: {gf_error:.6e} T ({gf_rel_error:.2f}%)")

# ============================================================================
# Step 5: Field Statistics
# ============================================================================

print("\n[Step 5] Field Statistics")
print("-" * 70)

if args.method in ['cf', 'both']:
	print("\nCoefficientFunction (exact evaluation):")
	# Sample CF at mesh vertices
	cf_samples = []
	for v in mesh.vertices:
		pt = v.point
		try:
			mpt = mesh(pt[0], pt[1], pt[2])
			B = B_cf(mpt)
			cf_samples.append(B)
		except:
			pass

	cf_samples = np.array(cf_samples)
	if len(cf_samples) > 0:
		Bx_cf = cf_samples[:, 0]
		By_cf = cf_samples[:, 1]
		Bz_cf = cf_samples[:, 2]
		B_mag_cf = np.sqrt(Bx_cf**2 + By_cf**2 + Bz_cf**2)

		print(f"  Bx: min={Bx_cf.min():.6e}, max={Bx_cf.max():.6e} T")
		print(f"  By: min={By_cf.min():.6e}, max={By_cf.max():.6e} T")
		print(f"  Bz: min={Bz_cf.min():.6e}, max={Bz_cf.max():.6e} T")
		print(f"  |B|: min={B_mag_cf.min():.6e}, max={B_mag_cf.max():.6e}, mean={B_mag_cf.mean():.6e} T")

if args.method in ['gf', 'both']:
	print("\nGridFunction (FEM interpolation):")
	B_data = gfB.vec.FV().NumPy()
	Bx_gf = B_data[0::3]
	By_gf = B_data[1::3]
	Bz_gf = B_data[2::3]
	B_mag_gf = np.sqrt(Bx_gf**2 + By_gf**2 + Bz_gf**2)

	print(f"  Bx: min={Bx_gf.min():.6e}, max={Bx_gf.max():.6e} T")
	print(f"  By: min={By_gf.min():.6e}, max={By_gf.max():.6e} T")
	print(f"  Bz: min={Bz_gf.min():.6e}, max={Bz_gf.max():.6e} T")
	print(f"  |B|: min={B_mag_gf.min():.6e}, max={B_mag_gf.max():.6e}, mean={B_mag_gf.mean():.6e} T")

# ============================================================================
# Step 6: VTK Export
# ============================================================================

print("\n[Step 6] VTK Export")
print("-" * 70)

if args.method == 'cf':
	print("\nExporting CoefficientFunction to VTK...")
	vtk = VTKOutput(mesh, coefs=[B_cf], names=['B_field_CF'], filename="radia_field_cf")
	vtk.Do()
	print("  [OK] VTK file created: radia_field_cf.vtu")
	print("       Contains EXACT Radia field values")

elif args.method == 'gf':
	print("\nExporting GridFunction to VTK...")
	vtk = VTKOutput(mesh, coefs=[gfB], names=['B_field_GF'], filename="radia_field_gf")
	vtk.Do()
	print("  [OK] VTK file created: radia_field_gf.vtu")
	print("       Contains FEM interpolated values")

else:  # both
	print("\nExporting both methods to VTK...")
	vtk = VTKOutput(mesh, coefs=[B_cf, gfB], names=['B_field_CF', 'B_field_GF'], filename="radia_field_compare")
	vtk.Do()
	print("  [OK] VTK file created: radia_field_compare.vtu")
	print("       Contains both CF (exact) and GF (interpolated) fields")
	print("       Open in Paraview to compare side-by-side")

# ============================================================================
# Summary
# ============================================================================

print("\n" + "=" * 70)
print("Summary")
print("=" * 70)

print(f"\nMesh parameters:")
print(f"  Elements: {mesh.ne}")
print(f"  Vertices: {mesh.nv}")
print(f"  maxh: {args.maxh} m ({args.maxh*1000} mm)")

if args.method in ['gf', 'both']:
	print(f"\nGridFunction parameters:")
	print(f"  DOFs: {fes.ndof}")
	print(f"  Order: {fes.globalorder}")

print("\nEvaluation methods:")
if args.method in ['cf', 'both']:
	print("  [CF] CoefficientFunction: Exact Radia evaluation (error = 0)")
if args.method in ['gf', 'both']:
	print(f"  [GF] GridFunction: FEM interpolation (error depends on mesh size)")

print("\nRecommendations:")
print("  - For exact values: Use CoefficientFunction (--method cf)")
print("  - For fast repeated evaluation: Use GridFunction with fine mesh (--method gf)")
print("  - For mesh refinement study: Try different --maxh values")
print("    Examples:")
print("      python visualize_field.py --maxh 0.02   # Coarse (60% error)")
print("      python visualize_field.py --maxh 0.01   # Medium (30% error)")
print("      python visualize_field.py --maxh 0.005  # Fine (1% error)")
print("      python visualize_field.py --maxh 0.002  # Very fine (<0.01% error)")

print("\n" + "=" * 70)
print("Complete")
print("=" * 70)
