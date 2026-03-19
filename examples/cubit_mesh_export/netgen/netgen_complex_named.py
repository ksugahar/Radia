"""
Netgen Export Example: Complex Geometry with Name-based Face Mapping

This example demonstrates the NEW workflow using name-based face mapping
for complex geometries (Boolean operations). This achieves Netgen-native
accuracy for any geometry.

Workflow:
1. OCC: Create geometry
2. OCC: Name faces with name_occ_faces()
3. OCC: Export STEP
4. Cubit: Import STEP (names preserved)
5. Cubit: Mesh
6. export_netgen_with_names(): Export with correct face mapping
7. SetGeomInfo: Set UV for curved surfaces
8. mesh.Curve(): High-order curving

Run: python netgen_complex_named.py
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

from netgen.occ import OCCGeometry, Box, Cylinder, gp_Ax2, gp_Dir, gp_Pnt
from ngsolve import Mesh, Integrate, CF, BND
import cubit
import cubit_mesh_export

# ============================================================
# Parameters
# ============================================================
BRICK_SIZE = 2.0   # Brick dimension
R_HOLE = 0.3       # Cylindrical hole radius
ORDER = 2          # Curving order

print("=== Netgen Export: Complex Geometry (Name-based) ===")
print(f"(Brick {BRICK_SIZE}x{BRICK_SIZE}x{BRICK_SIZE} with cylindrical hole R={R_HOLE})")
print()

# ============================================================
# Step 1: Create geometry in OCC (NOT Cubit)
# ============================================================
print("Step 1: Create geometry in OCC")

brick = Box(gp_Pnt(-BRICK_SIZE/2, -BRICK_SIZE/2, -BRICK_SIZE/2),
            gp_Pnt(BRICK_SIZE/2, BRICK_SIZE/2, BRICK_SIZE/2))
cyl = Cylinder(gp_Ax2(gp_Pnt(0, 0, -BRICK_SIZE), gp_Dir(0, 0, 1)),
               R_HOLE, 2 * BRICK_SIZE)
shape = brick - cyl

print(f"  OCC faces: {len(shape.faces)}")

# ============================================================
# Step 2: Name OCC faces
# ============================================================
print("\nStep 2: Name OCC faces")

num_named = cubit_mesh_export.name_occ_faces(shape)
print(f"  Named {num_named} faces")

# ============================================================
# Step 3: Export STEP from OCC
# ============================================================
print("\nStep 3: Export STEP from OCC")

step_file = os.path.join(work_dir, "complex_named.step")
shape.WriteStep(step_file)
print(f"  Exported: {step_file}")

# ============================================================
# Step 4: Load OCCGeometry (for mesh reference)
# ============================================================
print("\nStep 4: Load OCCGeometry")

geo = OCCGeometry(step_file)
print(f"  Faces: {len(geo.shape.faces)}")

# ============================================================
# Step 5: Import STEP into Cubit and mesh
# ============================================================
print("\nStep 5: Import and mesh in Cubit")

cubit.init(['cubit', '-nojournal', '-batch'])
cubit.cmd("reset")
cubit.cmd(f'import step "{step_file}" noheal')

# Verify names are preserved
surface_ids = cubit.get_entities('surface')
print(f"  Cubit surfaces: {len(surface_ids)}")
for sid in surface_ids:
    name = cubit.get_entity_name('surface', sid)
    print(f"    Surface {sid}: {name}")

# Mesh
cubit.cmd("volume all scheme tetmesh")
cubit.cmd("volume all size 0.15")
cubit.cmd("mesh volume all")
cubit.cmd("block 1 add tet all")
cubit.cmd('block 1 name "domain"')
cubit.cmd("block 2 add tri all")
cubit.cmd('block 2 name "boundary"')

print(f"  Tets: {cubit.get_tet_count()}")

# ============================================================
# Step 6: Export to Netgen with name-based mapping
# ============================================================
print("\nStep 6: Export to Netgen (name-based)")

ngmesh = cubit_mesh_export.export_netgen_with_names(cubit, geo)
print(f"  Elements: {ngmesh.ne}")

face_indices = sorted(set(el.index for el in ngmesh.Elements2D()))
print(f"  Face indices: {face_indices}")

# ============================================================
# Step 7: Set UV parameters for cylindrical hole
# ============================================================
print("\nStep 7: Set UV parameters (SetGeomInfo)")

modified = cubit_mesh_export.set_cylinder_geominfo(
    ngmesh, radius=R_HOLE, height=BRICK_SIZE, center=(0, 0, 0), axis='z'
)
print(f"  Modified {modified} entries for cylindrical surface")

# ============================================================
# Step 8: mesh.Curve() and accuracy test
# ============================================================
print(f"\nStep 8: mesh.Curve({ORDER}) and accuracy")

mesh = Mesh(ngmesh)
mesh.Curve(ORDER)

expected_vol = BRICK_SIZE**3 - math.pi * R_HOLE**2 * BRICK_SIZE
vol = Integrate(CF(1), mesh)

print(f"  Expected volume: {expected_vol:.6f}")
print(f"  Computed volume: {vol:.6f}")
print(f"  Error: {abs(vol-expected_vol)/expected_vol*100:.4f}%")

# Compare with order=3
print("\nCompare with order=3:")
mesh3 = Mesh(ngmesh)
mesh3.Curve(3)
vol3 = Integrate(CF(1), mesh3)
print(f"  Error: {abs(vol3-expected_vol)/expected_vol*100:.4f}%")

# Cleanup
os.remove(step_file)

print("\n=== DONE ===")
