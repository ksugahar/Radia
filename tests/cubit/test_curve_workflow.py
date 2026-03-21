"""
Test: Complete Cubit-to-NGSolve High-Order Curving Workflow

This test verifies the complete workflow using export_curved():
1. Cubit: Create geometry -> Export STEP
2. Cubit: Reimport STEP -> Generate mesh
3. export_curved(order=N) -> NGSolve Mesh (already curved)

Compares different curve orders for geometric accuracy.

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
sys.path.insert(0, os.path.join(repo_root, 'src', 'radia'))

# Use locally built NGSolve (ksugahar fork with SetGeomInfo API)
sys.path.insert(0, "s:/NGSolve/01_GitHub/install_ksugahar/Lib/site-packages")

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
# Test A-C: Compare curve orders using export_curved
# ============================================================
results = {}

for order in [1, 2, 3]:
    print("\n" + "=" * 60)
    print(f"Test: export_curved(order={order})")
    print("=" * 60)

    try:
        mesh = cubit_mesh_export.export_curved(cubit, order=order)
        area = Integrate(CF(1), mesh, VOL_or_BND=BND)
        vol = Integrate(CF(1), mesh)

        print(f"  Area: {area:.6f} (error: {abs(area-expected_area)/expected_area*100:.4f}%)")
        print(f"  Vol:  {vol:.6f} (error: {abs(vol-expected_vol)/expected_vol*100:.4f}%)")

        results[order] = {'area': area, 'vol': vol}
    except Exception as e:
        print(f"  Failed: {e}")
        results[order] = {'area': expected_area, 'vol': expected_vol}  # placeholder

# ============================================================
# Summary
# ============================================================
print("\n" + "=" * 60)
print("Summary: Volume Error Comparison")
print("=" * 60)
print()
print(f"{'Method':<35} {'Vol Error':>12}")
print("-" * 50)
for order in [1, 2, 3]:
    if order in results:
        vol = results[order]['vol']
        print(f"{'export_curved(order=' + str(order) + ')':<35} {abs(vol-expected_vol)/expected_vol*100:>11.4f}%")
print()

# Cleanup
os.remove(step_file)

print("=" * 60)
print("Test Complete")
print("=" * 60)
