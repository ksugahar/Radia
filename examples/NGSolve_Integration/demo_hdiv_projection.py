"""
HDiv Function Space Projection Demo - NGSolve Integration

This example demonstrates the RECOMMENDED way to project Radia magnetic field
onto NGSolve GridFunction using the HDiv function space.

WHY HDiv?
- HDiv (H(div) space) is mathematically appropriate for B field because:
  - div(B) = 0 (Maxwell's equation - no magnetic monopoles)
  - HDiv naturally preserves normal continuity across element boundaries
  - Best accuracy for magnetic field representation

CLAUDE.md Best Practices:
1. ALWAYS use rad.FldUnits('m') for NGSolve integration
2. Use HDiv(mesh, order=2) for best accuracy
3. Evaluate GridFunction at distances > 1 mesh cell from magnet surface
4. Use CoefficientFunction directly for maximum accuracy near boundaries

Output:
- VTK file for visualization in ParaView
- Comparison between HDiv, HCurl, and VectorH1 spaces

Usage:
    python demo_hdiv_projection.py

Date: 2025-12-05
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../build/Release'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

import numpy as np


# Check NGSolve availability
try:
    from ngsolve import *
    from netgen.csg import CSGeometry, OrthoBrick, Pnt
    import radia as rad
    import radia_ngsolve
    NGSOLVE_AVAILABLE = True
except ImportError as e:
    print(f"ERROR: Required module not available: {e}")
    print("This example requires NGSolve and radia_ngsolve.")
    print("Install NGSolve: pip install ngsolve")
    print("Build radia_ngsolve: cmake --build build --target radia_ngsolve")
    sys.exit(1)

print("=" * 70)
print("HDiv Function Space Projection Demo")
print("NGSolve Integration - Best Practice Example")
print("=" * 70)

# ============================================================================
# Step 1: Setup Radia with meters (REQUIRED for NGSolve)
# ============================================================================
print("\n[Step 1] Setting up Radia with meters")
print("-" * 70)

rad.UtiDelAll()
rad.FldUnits('m')  # CRITICAL: Required for NGSolve integration

units_str = rad.FldUnits()
print(f"  Unit system: {units_str.split()[0]}")
assert 'Length:  m' in units_str, "ERROR: Units must be meters for NGSolve!"

# ============================================================================
# Step 2: Create Radia Permanent Magnet
# ============================================================================
print("\n[Step 2] Creating Radia permanent magnet")
print("-" * 70)

magnet_center = [0, 0, 0]  # meters
magnet_size = [0.020, 0.020, 0.030]  # 20mm x 20mm x 30mm in meters
magnetization = [0, 0, 1.2]  # 1.2 T in z-direction (NdFeB typical)

# Create hexahedron vertices centered at [0, 0, 0] with dimensions [0.020, 0.020, 0.030] m
cx, cy, cz = magnet_center
dx, dy, dz = magnet_size[0] / 2, magnet_size[1] / 2, magnet_size[2] / 2
vertices = [
    [cx - dx, cy - dy, cz - dz],  # vertex 1
    [cx + dx, cy - dy, cz - dz],  # vertex 2
    [cx + dx, cy + dy, cz - dz],  # vertex 3
    [cx - dx, cy + dy, cz - dz],  # vertex 4
    [cx - dx, cy - dy, cz + dz],  # vertex 5
    [cx + dx, cy - dy, cz + dz],  # vertex 6
    [cx + dx, cy + dy, cz + dz],  # vertex 7
    [cx - dx, cy + dy, cz + dz],  # vertex 8
]

magnet = rad.ObjHexahedron(vertices, magnetization)
rad.MatApl(magnet, rad.MatPM(1.2, 900000, [0, 0, 1]))  # NdFeB material
rad.Solve(magnet, 0.0001, 10000)

print(f"  Magnet ID: {magnet}")
print(f"  Center: {magnet_center} m")
print(f"  Size: {magnet_size} m ({[s*1000 for s in magnet_size]} mm)")
print(f"  Magnetization: {magnetization} T")

# Reference field calculation
ref_point = [0, 0, 0.025]  # 25mm from center (5mm from top surface)
B_ref = rad.Fld(magnet, 'b', ref_point)
print(f"  B at {ref_point} m: [{B_ref[0]:.6f}, {B_ref[1]:.6f}, {B_ref[2]:.6f}] T")

# ============================================================================
# Step 3: Create NGSolve Mesh (Air Region)
# ============================================================================
print("\n[Step 3] Creating NGSolve mesh (air region)")
print("-" * 70)

# Create mesh OUTSIDE the magnet region
# Magnet extends from -0.01 to 0.01 in x and y, -0.015 to 0.015 in z
# We mesh the region x=[0.03, 0.08], y=[-0.03, 0.03], z=[-0.03, 0.03]
# This is far enough from the magnet surface (> 1 mesh cell)

geo = CSGeometry()
mesh_region = OrthoBrick(
    Pnt(0.025, -0.025, -0.025),  # Min corner (25mm from origin)
    Pnt(0.070, 0.025, 0.025)     # Max corner
)
geo.Add(mesh_region)

mesh_size = 0.008  # 8mm mesh size
ngmesh = geo.GenerateMesh(maxh=mesh_size)
mesh = Mesh(ngmesh)

print(f"  Mesh region: x=[0.025, 0.070] m, y=[-0.025, 0.025] m, z=[-0.025, 0.025] m")
print(f"  Mesh size: {mesh_size*1000} mm")
print(f"  Elements: {mesh.ne}")
print(f"  Vertices: {mesh.nv}")

# ============================================================================
# Step 4: Create HDiv Function Space (RECOMMENDED)
# ============================================================================
print("\n[Step 4] Creating HDiv function space (order=2)")
print("-" * 70)

# HDiv is RECOMMENDED for B field because:
# - div(B) = 0 is naturally satisfied
# - Normal continuity is preserved across elements
# - Best accuracy for magnetic flux density

fes_hdiv = HDiv(mesh, order=2)  # CLAUDE.md recommends order=2
B_gf_hdiv = GridFunction(fes_hdiv)

print(f"  HDiv space DOFs: {fes_hdiv.ndof}")
print(f"  Polynomial order: 2")

# Create RadiaField CoefficientFunction
B_cf = radia_ngsolve.RadiaField(magnet, 'b')
print(f"  RadiaField created: field_type='{B_cf.field_type}'")

# Project to HDiv GridFunction
print("  Projecting B field to HDiv GridFunction...")
B_gf_hdiv.Set(B_cf)
print("  [OK] HDiv projection complete")

# ============================================================================
# Step 5: Compare with Other Function Spaces
# ============================================================================
print("\n[Step 5] Comparing function spaces")
print("-" * 70)

# HCurl space (for vector potential A, but can also be used for B)
fes_hcurl = HCurl(mesh, order=2)
B_gf_hcurl = GridFunction(fes_hcurl)
B_gf_hcurl.Set(B_cf)
print(f"  HCurl space DOFs: {fes_hcurl.ndof}")

# VectorH1 space (continuous vector field)
fes_vh1 = VectorH1(mesh, order=2)
B_gf_vh1 = GridFunction(fes_vh1)
B_gf_vh1.Set(B_cf)
print(f"  VectorH1 space DOFs: {fes_vh1.ndof}")

# ============================================================================
# Step 6: Accuracy Comparison at Test Points
# ============================================================================
print("\n[Step 6] Accuracy comparison at test points")
print("-" * 70)

# Test points (far from magnet surface, > 1 mesh cell)
test_points = [
    (0.035, 0.0, 0.0),   # 25mm from magnet center
    (0.045, 0.0, 0.0),   # 35mm from magnet center
    (0.055, 0.0, 0.0),   # 45mm from magnet center
    (0.035, 0.015, 0.0), # Off-axis point
    (0.035, 0.0, 0.015), # Off-axis point
]

print(f"{'Point (m)':<25s} {'Radia |B|':>12s} {'HDiv |B|':>12s} {'HCurl |B|':>12s} {'VH1 |B|':>12s}")
print("-" * 75)

errors_hdiv = []
errors_hcurl = []
errors_vh1 = []

for pt in test_points:
    # Direct Radia evaluation
    B_radia = rad.Fld(magnet, 'b', list(pt))
    B_radia_mag = np.linalg.norm(B_radia)

    # NGSolve GridFunction evaluations
    mip = mesh(*pt)
    B_hdiv = B_gf_hdiv(mip)
    B_hcurl = B_gf_hcurl(mip)
    B_vh1 = B_gf_vh1(mip)

    B_hdiv_mag = np.linalg.norm(B_hdiv)
    B_hcurl_mag = np.linalg.norm(B_hcurl)
    B_vh1_mag = np.linalg.norm(B_vh1)

    # Calculate relative errors
    if B_radia_mag > 1e-6:
        err_hdiv = abs(B_hdiv_mag - B_radia_mag) / B_radia_mag * 100
        err_hcurl = abs(B_hcurl_mag - B_radia_mag) / B_radia_mag * 100
        err_vh1 = abs(B_vh1_mag - B_radia_mag) / B_radia_mag * 100
    else:
        err_hdiv = err_hcurl = err_vh1 = 0.0

    errors_hdiv.append(err_hdiv)
    errors_hcurl.append(err_hcurl)
    errors_vh1.append(err_vh1)

    pt_str = f"({pt[0]:.3f},{pt[1]:.3f},{pt[2]:.3f})"
    print(f"{pt_str:<25s} {B_radia_mag:>12.6f} {B_hdiv_mag:>12.6f} {B_hcurl_mag:>12.6f} {B_vh1_mag:>12.6f}")

print("-" * 75)
print(f"{'Max error (%)':<25s} {'':<12s} {max(errors_hdiv):>12.3f} {max(errors_hcurl):>12.3f} {max(errors_vh1):>12.3f}")
print(f"{'Avg error (%)':<25s} {'':<12s} {np.mean(errors_hdiv):>12.3f} {np.mean(errors_hcurl):>12.3f} {np.mean(errors_vh1):>12.3f}")

# ============================================================================
# Step 7: Export to VTK for Visualization
# ============================================================================
print("\n[Step 7] Exporting to VTK")
print("-" * 70)

output_dir = os.path.dirname(os.path.abspath(__file__))
vtk_filename = os.path.join(output_dir, "demo_hdiv_projection")

# Export HDiv GridFunction
vtkopts = VTKOutput(
    mesh,
    names=["B_hdiv"],
    coefs=[B_gf_hdiv],
    filename=vtk_filename
)
vtkopts.Do()

print(f"  VTK file exported: {vtk_filename}.vtu")
print("  Open in ParaView to visualize the magnetic field")

# ============================================================================
# Summary
# ============================================================================
print("\n" + "=" * 70)
print("SUMMARY: HDiv Function Space Projection")
print("=" * 70)
print("""
Best Practices for NGSolve Integration:

1. ALWAYS call rad.FldUnits('m') before creating Radia objects
   - NGSolve uses SI units (meters)
   - Radia default is millimeters

2. Use HDiv(mesh, order=2) for magnetic flux density B
   - Preserves div(B) = 0 (no magnetic monopoles)
   - Normal continuity across element boundaries
   - Best accuracy for magnetic field

3. Evaluate GridFunction at distances > 1 mesh cell from magnet surface
   - Projection interpolation can introduce errors near boundaries
   - Use CoefficientFunction directly for points near magnet

4. For vector potential A, use HCurl space
   - Preserves curl structure
   - Natural for A where B = curl(A)
""")

print("=" * 70)
print("[DONE] HDiv projection demo complete")
print("=" * 70)

# Cleanup
rad.UtiDelAll()
