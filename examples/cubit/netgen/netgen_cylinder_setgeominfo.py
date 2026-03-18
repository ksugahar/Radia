"""
Netgen Export Example: Cylinder with SetGeomInfo API

This example demonstrates the proper high-order curving workflow
using the SetGeomInfo API instead of SetDeformation.

Workflow:
1. Cubit: Create geometry -> Export STEP
2. Cubit: Reimport STEP -> Generate mesh
3. Netgen: Import mesh with geometry
4. Set UV parameters using set_cylinder_geominfo()
5. NGSolve: mesh.Curve(order) for high-order elements

Run: python netgen_cylinder_setgeominfo.py
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

# ============================================================
# Parameters
# ============================================================
R = 0.5  # Radius
H = 2.0  # Height
ORDER = 2  # Polynomial order for curving

print("=== Netgen Export: Cylinder with SetGeomInfo API ===")
print(f"Parameters: R={R}, H={H}, order={ORDER}")
print()

# ============================================================
# Step 1: Create geometry in Cubit
# ============================================================
print("Step 1: Create geometry in Cubit")

cubit.init(['cubit', '-nojournal', '-batch'])
cubit.cmd("reset")
cubit.cmd(f"create cylinder height {H} radius {R}")
print(f"  Surfaces: {cubit.get_surface_count()}")

# ============================================================
# Step 2: Export to STEP
# ============================================================
print("\nStep 2: Export to STEP")

step_file = os.path.join(work_dir, "cylinder_gi.step")
cubit.cmd(f'export step "{step_file}" overwrite')

# ============================================================
# Step 3: Reimport STEP (seam-aware)
# ============================================================
print("\nStep 3: Reimport STEP")

cubit.cmd("reset")
cubit.cmd(f'import step "{step_file}" heal')
print(f"  Surfaces after reimport: {cubit.get_surface_count()}")

# ============================================================
# Step 4: Mesh
# ============================================================
print("\nStep 4: Mesh")

cubit.cmd("volume all scheme tetmesh")
cubit.cmd("volume all size 0.15")
cubit.cmd("mesh volume all")
cubit.cmd("block 1 add tet all")
cubit.cmd('block 1 name "domain"')
cubit.cmd("block 2 add tri all")
cubit.cmd('block 2 name "boundary"')

print(f"  Tets: {cubit.get_tet_count()}")

# ============================================================
# Step 5: Export to Netgen with geometry
# ============================================================
print("\nStep 5: Export to Netgen with geometry")

geo = OCCGeometry(step_file)
ngmesh = cubit_mesh_export.export_netgen(cubit, geometry=geo)

print(f"  Elements: {ngmesh.ne}")

# ============================================================
# Step 6: Set UV parameters using SetGeomInfo API
# ============================================================
print("\nStep 6: Set UV parameters (SetGeomInfo API)")

modified = cubit_mesh_export.set_cylinder_geominfo(
    ngmesh, radius=R, height=H, center=(0, 0, 0), axis='z'
)
print(f"  Modified {modified} vertex geominfo entries")

# ============================================================
# Step 7: Create NGSolve mesh and curve
# ============================================================
print(f"\nStep 7: mesh.Curve({ORDER})")

mesh = Mesh(ngmesh)
mesh.Curve(ORDER)

print(f"  Boundaries: {mesh.GetBoundaries()}")

# ============================================================
# Step 8: Accuracy comparison
# ============================================================
print("\nStep 8: Accuracy")

expected_area = 2*math.pi*R*H + 2*math.pi*R*R
expected_vol = math.pi*R*R*H

area = Integrate(CF(1), mesh, VOL_or_BND=BND)
vol = Integrate(CF(1), mesh)

print(f"  Expected:  Area={expected_area:.6f}, Vol={expected_vol:.6f}")
print(f"  Computed:  Area={area:.6f}, Vol={vol:.6f}")
print(f"  Error:     Area={abs(area-expected_area)/expected_area*100:.4f}%, "
      f"Vol={abs(vol-expected_vol)/expected_vol*100:.4f}%")

# ============================================================
# Compare with higher order
# ============================================================
print(f"\nCompare with order=3:")
mesh3 = Mesh(cubit_mesh_export.export_netgen(cubit, geometry=geo))
cubit_mesh_export.set_cylinder_geominfo(
    mesh3.ngmesh, radius=R, height=H, center=(0, 0, 0), axis='z'
)
mesh3.Curve(3)

area3 = Integrate(CF(1), mesh3, VOL_or_BND=BND)
vol3 = Integrate(CF(1), mesh3)
print(f"  Error:     Area={abs(area3-expected_area)/expected_area*100:.4f}%, "
      f"Vol={abs(vol3-expected_vol)/expected_vol*100:.4f}%")

# Cleanup
os.remove(step_file)

print("\n=== DONE ===")
