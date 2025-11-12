"""
Radia Field Types Demo - NGSolve Integration

This script demonstrates how to use the new unified RadiaField interface
to access different physical field quantities:
- 'b': Magnetic flux density (Tesla)
- 'h': Magnetic field (A/m)
- 'a': Vector potential (T*m)
- 'm': Magnetization (A/m)

Usage:
	python demo_field_types.py


Date: 2025-11-01
"""

import sys
sys.path.insert(0, r"S:\radia\01_GitHub\build\Release")
sys.path.insert(0, r"S:\radia\01_GitHub\dist")

from ngsolve import *
from netgen.csg import CSGeometry, OrthoBrick, Pnt
import numpy as np
import radia as rad
import rad_ngsolve

print("=" * 70)
print("Radia Field Types Demo")
print("=" * 70)

# ============================================================================
# Step 1: Create Radia Magnet Geometry
# ============================================================================

print("\n[Step 1] Creating Radia Magnet")
print("-" * 70)

magnet_center = [0, 0, 0]
magnet_size = [20, 20, 30]

magnet = rad.ObjRecMag(magnet_center, magnet_size, [0, 0, 1.2])
rad.MatApl(magnet, rad.MatPM(1.2, 900000, [0, 0, 1]))  # NdFeB
rad.Solve(magnet, 0.0001, 10000)

print(f"Magnet created: object #{magnet}")
print(f"  Center: {magnet_center} mm")
print(f"  Size: {magnet_size} mm")
print(f"  Material: NdFeB, Br = 1.2 T")

# ============================================================================
# Step 2: Create CoefficientFunctions for Different Field Types
# ============================================================================

print("\n[Step 2] Creating CoefficientFunctions for All Field Types")
print("-" * 70)

# Method 1: New unified interface (recommended)
B_cf = rad_ngsolve.RadiaField(magnet, 'b')  # Magnetic flux density
H_cf = rad_ngsolve.RadiaField(magnet, 'h')  # Magnetic field
A_cf = rad_ngsolve.RadiaField(magnet, 'a')  # Vector potential
M_cf = rad_ngsolve.RadiaField(magnet, 'm')  # Magnetization

print(f"Created CoefficientFunctions:")
print(f"  B (flux density):  {type(B_cf).__name__}, field_type='{B_cf.field_type}'")
print(f"  H (magnetic field): {type(H_cf).__name__}, field_type='{H_cf.field_type}'")
print(f"  A (vector potential): {type(A_cf).__name__}, field_type='{A_cf.field_type}'")
print(f"  M (magnetization):  {type(M_cf).__name__}, field_type='{M_cf.field_type}'")

# ============================================================================
# Step 3: Create 3D Mesh
# ============================================================================

print("\n[Step 3] Creating 3D Mesh")
print("-" * 70)

geo = CSGeometry()
geo.Add(OrthoBrick(Pnt(-0.05, -0.05, -0.05), Pnt(0.05, 0.05, 0.05)))
mesh = Mesh(geo.GenerateMesh(maxh=0.01))

print(f"3D Mesh generated:")
print(f"  Elements: {mesh.ne}")
print(f"  Vertices: {mesh.nv}")
print(f"  Domain: [-0.05, 0.05] m = [-50, 50] mm")

# ============================================================================
# Step 4: Evaluate Fields at Test Points
# ============================================================================

print("\n[Step 4] Field Evaluation at Test Points")
print("-" * 70)

test_points = [
	(0.000, 0.000, 0.000),	# 0mm (center)
	(0.000, 0.000, 0.020),  # 20mm (above magnet)
	(0.000, 0.000, 0.040),  # 40mm (far from magnet)
]

print("\n" + "=" * 70)
print("Field Values at Different Points")
print("=" * 70)

for pt in test_points:
	pt_mm = [pt[0]*1000, pt[1]*1000, pt[2]*1000]
	mesh_pt = mesh(*pt)

	print(f"\nPoint: {pt} m = {pt_mm} mm")

	# B field (Tesla)
	B_val = B_cf(mesh_pt)
	print(f"  B (flux density):    Bx={B_val[0]:8.6f}, By={B_val[1]:8.6f}, Bz={B_val[2]:8.6f} T")

	# H field (A/m)
	H_val = H_cf(mesh_pt)
	print(f"  H (magnetic field):  Hx={H_val[0]:8.2f}, Hy={H_val[1]:8.2f}, Hz={H_val[2]:8.2f} A/m")

	# A field (T*m)
	A_val = A_cf(mesh_pt)
	print(f"  A (vector potential): Ax={A_val[0]:8.6f}, Ay={A_val[1]:8.6f}, Az={A_val[2]:8.6f} T*m")

	# M field (A/m)
	M_val = M_cf(mesh_pt)
	print(f"  M (magnetization):   Mx={M_val[0]:8.2f}, My={M_val[1]:8.2f}, Mz={M_val[2]:8.2f} A/m")

# ============================================================================
# Step 5: Compare with Radia Direct Evaluation
# ============================================================================

print("\n" + "=" * 70)
print("Comparison with Radia Direct Evaluation")
print("=" * 70)

pt_mm = [0, 0, 0]  # Center point in mm
mesh_pt = mesh(0, 0, 0)

print(f"\nAt center ({pt_mm} mm):")

# Radia direct
B_radia = rad.Fld(magnet, 'b', pt_mm)
H_radia = rad.Fld(magnet, 'h', pt_mm)
A_radia = rad.Fld(magnet, 'a', pt_mm)
M_radia = rad.Fld(magnet, 'm', pt_mm)

# NGSolve CoefficientFunction
B_ngsolve = B_cf(mesh_pt)
H_ngsolve = H_cf(mesh_pt)
A_ngsolve = A_cf(mesh_pt)
M_ngsolve = M_cf(mesh_pt)

# Compare
print(f"\nB field:")
print(f"  Radia:   Bz = {B_radia[2]:.6f} T")
print(f"  NGSolve: Bz = {B_ngsolve[2]:.6f} T")
print(f"  Error:   {abs(B_radia[2] - B_ngsolve[2]):.2e} T")

print(f"\nH field:")
print(f"  Radia:   Hz = {H_radia[2]:.2f} A/m")
print(f"  NGSolve: Hz = {H_ngsolve[2]:.2f} A/m")
print(f"  Error:   {abs(H_radia[2] - H_ngsolve[2]):.2e} A/m")

print(f"\nA field:")
print(f"  Radia:   |A| = {np.linalg.norm(A_radia):.6e} T*m")
print(f"  NGSolve: |A| = {np.linalg.norm(A_ngsolve):.6e} T*m")

print(f"\nM field (Magnetization):")
print(f"  Radia:   Mz = {M_radia[2]:.2f} A/m")
print(f"  NGSolve: Mz = {M_ngsolve[2]:.2f} A/m")
print(f"  Error:   {abs(M_radia[2] - M_ngsolve[2]):.2e} A/m")

# ============================================================================
# Summary
# ============================================================================

print("\n" + "=" * 70)
print("Summary")
print("=" * 70)

print("\nNew Unified Interface:")
print("  rad_ngsolve.RadiaField(radia_obj, field_type)")
print("\nField types:")
print("  'b' - Magnetic flux density B (Tesla)")
print("  'h' - Magnetic field H (A/m)")
print("  'a' - Vector potential A (T*m)")
print("  'm' - Magnetization M (A/m)")

print("\nRecommendation:")
print("  Use the new RadiaField interface for all new code.")
print("  Legacy interfaces (RadBfield, RadHfield, etc.) remain")
print("  supported for backward compatibility.")

print("\n" + "=" * 70)
print("Complete")
print("=" * 70)
