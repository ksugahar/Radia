"""
Netgen Export Example: Torus with SetGeomInfo API

This example demonstrates the proper workflow using SetGeomInfo API
for high-order curving of externally imported meshes.

Run: python netgen_torus_setgeominfo.py
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
R_major = 1.0  # Major radius (center of tube to center of torus)
R_minor = 0.3  # Minor radius (tube radius)
ORDER = 2      # Curving order

print("=== Netgen Export: Torus with SetGeomInfo API ===")
print(f"Parameters: R_major={R_major}, R_minor={R_minor}, order={ORDER}")
print()

# ============================================================
# Step 1: Create geometry in Cubit
# ============================================================
print("Step 1: Create geometry in Cubit")

cubit.init(['cubit', '-nojournal', '-batch'])
cubit.cmd("reset")
cubit.cmd(f"create torus major radius {R_major} minor radius {R_minor}")
print(f"  Surfaces: {cubit.get_surface_count()}")

# ============================================================
# Step 2: Export to STEP
# ============================================================
print("\nStep 2: Export to STEP")

step_file = os.path.join(work_dir, "torus_gi.step")
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
cubit.cmd("volume all size 0.08")
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
# Step 6: Set UV parameters (SetGeomInfo API)
# ============================================================
print("\nStep 6: Set UV parameters (SetGeomInfo API)")

modified = cubit_mesh_export.set_torus_geominfo(
    ngmesh, major_radius=R_major, minor_radius=R_minor,
    center=(0, 0, 0), axis='z'
)
print(f"  Modified {modified} vertex geominfo entries")

# ============================================================
# Step 7: mesh.Curve()
# ============================================================
print(f"\nStep 7: mesh.Curve({ORDER})")

mesh = Mesh(ngmesh)
mesh.Curve(ORDER)
print(f"  Boundaries: {mesh.GetBoundaries()}")

# ============================================================
# Step 8: Accuracy
# ============================================================
print("\nStep 8: Accuracy")

# Torus surface area = 4 * pi^2 * R_major * R_minor
# Torus volume = 2 * pi^2 * R_major * R_minor^2
expected_area = 4 * math.pi**2 * R_major * R_minor
expected_vol = 2 * math.pi**2 * R_major * R_minor**2

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
