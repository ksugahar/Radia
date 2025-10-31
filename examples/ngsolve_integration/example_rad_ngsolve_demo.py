"""
Demonstration: Using rad_ngsolve for FEM coupling with Radia

This example shows how to:
1. Create a Radia magnetic structure
2. Create NGSolve CoefficientFunction from Radia fields
3. Use the field in NGSolve FEM calculations

Author: Generated with Claude Code
Date: 2025-10-30
"""
import sys
sys.path.insert(0, r'S:\radia\01_GitHub\build\Release')
sys.path.insert(0, r'S:\radia\01_GitHub\dist')

# CRITICAL: Import ngsolve BEFORE rad_ngsolve to load libngsolve.dll
import ngsolve
from ngsolve import Mesh, H1, GridFunction, CoefficientFunction, Integrate, dx
from netgen.geom2d import SplineGeometry
import radia as rad
import rad_ngsolve
import numpy as np

print("=" * 70)
print("rad_ngsolve Demonstration")
print("C++ CoefficientFunction for Radia + NGSolve coupling")
print("=" * 70)

# ========================================================================
# Step 1: Create Radia magnetic structure
# ========================================================================
print("\n[Step 1] Creating Radia magnetic structure...")

# Create permanent magnet (20mm cube, magnetized in +z direction)
magnet_size = 20  # mm
magnetization = 1.2  # Tesla (NdFeB)

magnet = rad.ObjRecMag([0, 0, 0], [magnet_size, magnet_size, magnet_size], [0, 0, magnetization] )

# Apply material properties
rad.MatApl(magnet, rad.MatStd('NdFeB', magnetization))

# Solve magnetostatics
rad.Solve(magnet, 0.0001, 10000)

print(f"  [OK] Magnet created: {magnet_size}mm cube, Br={magnetization}T")
print(f"  [OK] Radia object ID: {magnet}")

# Verify field at key points
B_center = rad.Fld(magnet, 'b', [0, 0, 0])
B_above = rad.Fld(magnet, 'b', [0, 0, 30])
print(f"  [OK] Field at center: Bz = {B_center[2]:.4f} T")
print(f"  [OK] Field at 30mm above: Bz = {B_above[2]:.4f} T")

# ========================================================================
# Step 2: Create NGSolve mesh (region above magnet)
# ========================================================================
print("\n[Step 2] Creating NGSolve mesh...")

# Create 2D mesh in x-y plane (will evaluate at different z heights)
geo = SplineGeometry()
geo.AddRectangle(
    p1=(-40, -40),  # mm
    p2=(40, 40),    # mm
    bc="outer"
)
mesh = Mesh(geo.GenerateMesh(maxh=10))

print(f"  [OK] Mesh created: {mesh.ne} elements")
print(f"  [OK] Domain: [-40, 40] x [-40, 40] mm^2")

# ========================================================================
# Step 3: Create CoefficientFunction from Radia
# ========================================================================
print("\n[Step 3] Creating CoefficientFunction from Radia...")

# Create B-field CoefficientFunction
B_cf = rad_ngsolve.RadBfield(magnet)
H_cf = rad_ngsolve.RadHfield(magnet)
A_cf = rad_ngsolve.RadAfield(magnet)

print(f"  [OK] RadBfield created: {type(B_cf)}")
print(f"  [OK] Is CoefficientFunction: {isinstance(B_cf, CoefficientFunction)}")
print(f"  [OK] RadHfield created: {type(H_cf)}")
print(f"  [OK] RadAfield created: {type(A_cf)}")

# ========================================================================
# Step 4: Evaluate field on NGSolve mesh
# ========================================================================
print("\n[Step 4] Evaluating Radia field on NGSolve mesh...")

# Create finite element space (vector-valued for 3D field)
fes = H1(mesh, order=2, dim=3)
print(f"  [OK] FE space created: {fes.ndof} DOFs")

# Evaluate B field at z=0 plane (through magnet center)
gf_B = GridFunction(fes)
gf_B.Set(B_cf)

# Get field data
B_data = gf_B.vec.FV().NumPy()
print(f"  [OK] Field evaluated on mesh")
print(f"  [OK] Field range: [{B_data.min():.6f}, {B_data.max():.6f}]")

# ========================================================================
# Step 5: Compute field integrals using NGSolve
# ========================================================================
print("\n[Step 5] Computing field integrals...")

# Compute magnetic energy density integral: Int B*H dV
# Note: In vacuum, B = mu_0*H, so energy density = B^2/(2mu_0)
mu_0 = 4 * np.pi * 1e-7  # H/m = T*m/A
# Convert to mm units: H/mm = 1e-3 H/m
mu_0_mm = mu_0 * 1e-3

# Integrate |B|^2 over the domain
B_squared = B_cf[0]**2 + B_cf[1]**2 + B_cf[2]**2
integral_B2 = Integrate(B_squared, mesh)
print(f"  [OK] Int |B|^2 dA = {integral_B2:.6e} T^2*mm^2")

# Integrate Bz component
integral_Bz = Integrate(B_cf[2], mesh)
print(f"  [OK] Int Bz dA = {integral_Bz:.6e} T*mm^2")

# ========================================================================
# Step 6: Sample field along a line
# ========================================================================
print("\n[Step 6] Sampling field along z-axis...")

z_points = np.linspace(-30, 50, 17)  # mm
Bz_along_z = []

print("  z [mm]    Bz [T]")
print("  " + "-" * 22)
for z in z_points:
    B = rad.Fld(magnet, 'b', [0, 0, z])
    Bz_along_z.append(B[2])
    if abs(z) < 31 or z > 29:  # Print selected points
        print(f"  {z:6.1f}    {B[2]:8.5f}")

print(f"\n  [OK] Field sampled at {len(z_points)} points")

# ========================================================================
# Summary
# ========================================================================
print("\n" + "=" * 70)
print("SUCCESS: Demonstration completed")
print("=" * 70)
print("\nKey achievements:")
print("  [OK] Radia magnetic structure created and solved")
print("  [OK] NGSolve mesh created for FEM domain")
print("  [OK] C++ CoefficientFunction bridges Radia -> NGSolve")
print("  [OK] Field evaluated on mesh (no Python callback overhead)")
print("  [OK] Field integrals computed using NGSolve")
print("  [OK] Ready for coupled FEM analysis")
print("\nThis enables:")
print("  - Magnetostatic-thermal coupling")
print("  - Magnetostatic-structural coupling")
print("  - Eddy current analysis with Radia background fields")
print("  - Optimization with gradient-based methods")
print("=" * 70)
