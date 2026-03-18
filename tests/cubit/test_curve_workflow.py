"""
Test: Complete Cubit-to-NGSolve High-Order Curving Workflow

This test verifies the complete workflow:
1. Cubit: Create geometry -> Export STEP
2. Cubit: Reimport STEP -> Generate mesh
3. Netgen: Load STEP geometry + Import Cubit mesh
4. NGSolve: mesh.Curve(order) for high-order elements

Compares different approaches:
- Linear mesh (no curving)
- SetDeformation (workaround)
- mesh.Curve() with geometry (proper approach)

Run: python test_curve_workflow.py
"""

import sys
import os
import math

cubit_path = os.environ.get("CUBIT_PATH")
if cubit_path:
	sys.path.append(cubit_path)

work_dir = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.dirname(os.path.dirname(work_dir))
sys.path.insert(0, repo_root)

# Use locally built NGSolve (ksugahar fork with SetGeomInfo API)
sys.path.insert(0, "s:/NGSolve/01_GitHub/install_ksugahar/Lib/site-packages")

from netgen.occ import OCCGeometry
from ngsolve import Mesh, Integrate, CF, BND
import cubit
import cubit_mesh_export

print("=" * 60)
print("Test: Complete Cubit-to-NGSolve High-Order Curving Workflow")
print("=" * 60)
print()

# ============================================================
# Parameters
# ============================================================
R = 0.5  # Radius
H = 2.0  # Height
MESH_SIZE = 0.2  # Coarser mesh for clearer error comparison

expected_area = 2*math.pi*R*H + 2*math.pi*R*R
expected_vol = math.pi*R*R*H

print(f"Cylinder: R={R}, H={H}")
print(f"Expected: Area={expected_area:.6f}, Vol={expected_vol:.6f}")
print()

# ============================================================
# Step 1: Create geometry in Cubit and export STEP
# ============================================================
print("Step 1: Create geometry in Cubit")

cubit.init(['cubit', '-nojournal', '-batch'])
cubit.cmd("reset")
cubit.cmd(f"create cylinder height {H} radius {R}")

step_file = os.path.join(work_dir, "test_cylinder.step")
cubit.cmd(f'export step "{step_file}" overwrite')
print(f"  STEP exported: {step_file}")

# ============================================================
# Step 2: Reimport STEP and mesh in Cubit
# ============================================================
print("\nStep 2: Reimport STEP and mesh in Cubit")

cubit.cmd("reset")
cubit.cmd(f'import step "{step_file}" heal')
print(f"  Surfaces after reimport: {cubit.get_surface_count()}")

cubit.cmd("volume all scheme tetmesh")
cubit.cmd(f"volume all size {MESH_SIZE}")
cubit.cmd("mesh volume all")
cubit.cmd("block 1 add tet all")
cubit.cmd('block 1 name "domain"')
cubit.cmd("block 2 add tri all")
cubit.cmd('block 2 name "boundary"')

print(f"  Tets: {cubit.get_tet_count()}")

# ============================================================
# Step 3: Export to Netgen with OCC geometry
# ============================================================
print("\nStep 3: Export to Netgen with OCC geometry")

geo = OCCGeometry(step_file)
ngmesh = cubit_mesh_export.export_netgen(cubit, geometry=geo)
mesh = Mesh(ngmesh)

print(f"  Elements: {mesh.ne}")
print(f"  Boundaries: {mesh.GetBoundaries()}")

# ============================================================
# Test A: Linear mesh (baseline)
# ============================================================
print("\n" + "=" * 60)
print("Test A: Linear mesh (no curving)")
print("=" * 60)

mesh_linear = Mesh(ngmesh)
mesh_linear.Curve(1)

area_linear = Integrate(CF(1), mesh_linear, VOL_or_BND=BND)
vol_linear = Integrate(CF(1), mesh_linear)

print(f"  Area: {area_linear:.6f} (error: {abs(area_linear-expected_area)/expected_area*100:.4f}%)")
print(f"  Vol:  {vol_linear:.6f} (error: {abs(vol_linear-expected_vol)/expected_vol*100:.4f}%)")

# ============================================================
# Test B: mesh.Curve(2) without SetDeformation
# ============================================================
print("\n" + "=" * 60)
print("Test B: mesh.Curve(2) directly (no SetDeformation)")
print("=" * 60)

mesh_curve2 = Mesh(ngmesh)
try:
    mesh_curve2.Curve(2)
    area_curve2 = Integrate(CF(1), mesh_curve2, VOL_or_BND=BND)
    vol_curve2 = Integrate(CF(1), mesh_curve2)
    print(f"  Area: {area_curve2:.6f} (error: {abs(area_curve2-expected_area)/expected_area*100:.4f}%)")
    print(f"  Vol:  {vol_curve2:.6f} (error: {abs(vol_curve2-expected_vol)/expected_vol*100:.4f}%)")
except Exception as e:
    print(f"  Failed: {e}")
    area_curve2 = area_linear
    vol_curve2 = vol_linear

# ============================================================
# Test C: mesh.Curve(3) without SetDeformation
# ============================================================
print("\n" + "=" * 60)
print("Test C: mesh.Curve(3) directly (no SetDeformation)")
print("=" * 60)

mesh_curve3 = Mesh(ngmesh)
try:
    mesh_curve3.Curve(3)
    area_curve3 = Integrate(CF(1), mesh_curve3, VOL_or_BND=BND)
    vol_curve3 = Integrate(CF(1), mesh_curve3)
    print(f"  Area: {area_curve3:.6f} (error: {abs(area_curve3-expected_area)/expected_area*100:.4f}%)")
    print(f"  Vol:  {vol_curve3:.6f} (error: {abs(vol_curve3-expected_vol)/expected_vol*100:.4f}%)")
except Exception as e:
    print(f"  Failed: {e}")
    area_curve3 = area_linear
    vol_curve3 = vol_linear

# ============================================================
# Test D: SetDeformation (workaround approach)
# ============================================================
print("\n" + "=" * 60)
print("Test D: SetDeformation (workaround approach)")
print("=" * 60)

if hasattr(cubit_mesh_export, 'detect_cylinder_boundaries'):
    mesh_deform = Mesh(ngmesh)
    mesh_deform.Curve(1)

    cyl_boundaries = cubit_mesh_export.detect_cylinder_boundaries(mesh_deform, R, axis='z')
    print(f"  Cylinder boundaries: {cyl_boundaries}")

    cubit_mesh_export.apply_cylinder_deformation(
        mesh_deform, radius=R, boundary_names=cyl_boundaries, order=2, axis='z'
    )

    area_deform = Integrate(CF(1), mesh_deform, VOL_or_BND=BND)
    vol_deform = Integrate(CF(1), mesh_deform)

    print(f"  Area: {area_deform:.6f} (error: {abs(area_deform-expected_area)/expected_area*100:.4f}%)")
    print(f"  Vol:  {vol_deform:.6f} (error: {abs(vol_deform-expected_vol)/expected_vol*100:.4f}%)")
else:
    print("  SKIPPED: detect_cylinder_boundaries not yet implemented")
    vol_deform = vol_linear

# ============================================================
# Test E: SetDeformation order=3
# ============================================================
print("\n" + "=" * 60)
print("Test E: SetDeformation order=3")
print("=" * 60)

if hasattr(cubit_mesh_export, 'apply_cylinder_deformation'):
    mesh_deform3 = Mesh(ngmesh)
    mesh_deform3.Curve(1)

    cubit_mesh_export.apply_cylinder_deformation(
        mesh_deform3, radius=R, boundary_names=cyl_boundaries, order=3, axis='z'
    )

    area_deform3 = Integrate(CF(1), mesh_deform3, VOL_or_BND=BND)
    vol_deform3 = Integrate(CF(1), mesh_deform3)

    print(f"  Area: {area_deform3:.6f} (error: {abs(area_deform3-expected_area)/expected_area*100:.4f}%)")
    print(f"  Vol:  {vol_deform3:.6f} (error: {abs(vol_deform3-expected_vol)/expected_vol*100:.4f}%)")
else:
    print("  SKIPPED: apply_cylinder_deformation not yet implemented")
    vol_deform3 = vol_linear

# ============================================================
# Summary
# ============================================================
print("\n" + "=" * 60)
print("Summary: Volume Error Comparison")
print("=" * 60)
print()
print(f"{'Method':<35} {'Vol Error':>12}")
print("-" * 50)
print(f"{'Linear (baseline)':<35} {abs(vol_linear-expected_vol)/expected_vol*100:>11.4f}%")
print(f"{'mesh.Curve(2) direct':<35} {abs(vol_curve2-expected_vol)/expected_vol*100:>11.4f}%")
print(f"{'mesh.Curve(3) direct':<35} {abs(vol_curve3-expected_vol)/expected_vol*100:>11.4f}%")
print(f"{'SetDeformation order=2':<35} {abs(vol_deform-expected_vol)/expected_vol*100:>11.4f}%")
print(f"{'SetDeformation order=3':<35} {abs(vol_deform3-expected_vol)/expected_vol*100:>11.4f}%")
print()

# Cleanup
os.remove(step_file)

print("=" * 60)
print("Test Complete")
print("=" * 60)
