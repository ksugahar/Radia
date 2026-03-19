"""
Netgen Export Example: Hex Cylinder with SetGeomInfo API

Key use case: Cubit's hex meshing + NGSolve high-order elements
This demonstrates that SetGeomInfo works with hexahedral meshes too.

Run: python netgen_hex_cylinder_setgeominfo.py
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
R = 0.5    # Radius
H = 2.0    # Height
ORDER = 2  # Curving order

print("=== Netgen Export: Hex Cylinder with SetGeomInfo API ===")
print(f"Parameters: R={R}, H={H}, order={ORDER}")
print("(Cubit hex meshing + NGSolve high-order)")
print()

# ============================================================
# Step 1: Create geometry in Cubit
# ============================================================
print("Step 1: Create geometry in Cubit")

cubit.init(['cubit', '-nojournal', '-batch'])
cubit.cmd("reset")
cubit.cmd(f"create cylinder height {H} radius {R}")

# ============================================================
# Step 2: Export STEP and reimport
# ============================================================
print("\nStep 2: Export STEP and reimport")

step_file = os.path.join(work_dir, "hex_cylinder_gi.step")
cubit.cmd(f'export step "{step_file}" overwrite')
cubit.cmd("reset")
cubit.cmd(f'import step "{step_file}" heal')
print(f"  Surfaces: {cubit.get_surface_count()}")

# ============================================================
# Step 3: Hex mesh using sweep
# ============================================================
print("\nStep 3: Hex mesh (sweep)")

cubit.cmd("volume 1 scheme sweep")
cubit.cmd("surface all size 0.15")
cubit.cmd("mesh volume 1")

cubit.cmd("block 1 add hex all")
cubit.cmd('block 1 name "domain"')
cubit.cmd("block 2 add face all")
cubit.cmd('block 2 name "boundary"')

print(f"  Hexes: {cubit.get_hex_count()}")
print(f"  Faces: {len(cubit.get_block_faces(2))}")

# ============================================================
# Step 4: Export to Netgen with geometry
# ============================================================
print("\nStep 4: Export to Netgen with geometry")

geo = OCCGeometry(step_file)
ngmesh = cubit_mesh_export.export_netgen(cubit, geometry=geo)
print(f"  Elements: {ngmesh.ne}")

# ============================================================
# Step 5: Set UV parameters (SetGeomInfo API)
# ============================================================
print("\nStep 5: Set UV parameters (SetGeomInfo API)")

modified = cubit_mesh_export.set_cylinder_geominfo(
    ngmesh, radius=R, height=H, center=(0, 0, 0), axis='z'
)
print(f"  Modified {modified} vertex geominfo entries")

# ============================================================
# Step 6: mesh.Curve()
# ============================================================
print(f"\nStep 6: mesh.Curve({ORDER})")

mesh = Mesh(ngmesh)
mesh.Curve(ORDER)
print(f"  Boundaries: {mesh.GetBoundaries()}")

# ============================================================
# Step 7: Accuracy
# ============================================================
print("\nStep 7: Accuracy")

expected_area = 2 * math.pi * R * H + 2 * math.pi * R * R
expected_vol = math.pi * R * R * H

area = Integrate(CF(1), mesh, VOL_or_BND=BND)
vol = Integrate(CF(1), mesh)

print(f"  Expected:  Area={expected_area:.6f}, Vol={expected_vol:.6f}")
print(f"  Computed:  Area={area:.6f}, Vol={vol:.6f}")
print(f"  Error:     Area={abs(area-expected_area)/expected_area*100:.4f}%, "
      f"Vol={abs(vol-expected_vol)/expected_vol*100:.4f}%")

# Compare with order=3
print("\nCompare with order=3:")
mesh2 = Mesh(ngmesh)
mesh2.Curve(3)
area3 = Integrate(CF(1), mesh2, VOL_or_BND=BND)
vol3 = Integrate(CF(1), mesh2)
print(f"  Error:     Area={abs(area3-expected_area)/expected_area*100:.4f}%, "
      f"Vol={abs(vol3-expected_vol)/expected_vol*100:.4f}%")

# Cleanup
os.remove(step_file)

print("\n=== DONE ===")
