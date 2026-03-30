"""
Netgen Export Example: Complex Geometry with High-Order Curving

Uses export_NGSolveCurvedMesh() for automatic geometry association and curving.

Workflow:
1. OCC: Create geometry (brick with cylindrical hole)
2. OCC: Export STEP -> Cubit imports STEP (geometry transfer)
3. Cubit: Mesh
4. export_NGSolveCurvedMesh(order=N): Export with curving (uses Cubit's ACIS kernel)
5. Accuracy test via volume integration

Run: python netgen_complex_named.py
"""

import sys
import os
import math

# Auto-detect Cubit installation
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'src', 'radia'))
from install_panels import find_cubit_bin as _fcb
_cubit_path = _fcb()
if _cubit_path and _cubit_path not in sys.path:
    sys.path.append(_cubit_path)

work_dir = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.dirname(os.path.dirname(work_dir))
sys.path.insert(0, repo_root)
sys.path.insert(0, os.path.join(repo_root, 'src', 'radia'))

# Use locally built NGSolve (ksugahar fork with SetGeomInfo API)

from netgen.occ import Box, Cylinder, gp_Ax2, gp_Dir, gp_Pnt
from ngsolve import Mesh, Integrate, CF, BND
import cubit
import cubit_mesh_export

# ============================================================
# Parameters
# ============================================================
BRICK_SIZE = 2.0   # Brick dimension
R_HOLE = 0.3       # Cylindrical hole radius
ORDER = 2          # Curving order

print("=== Netgen Export: Complex Geometry (export_NGSolveCurvedMesh) ===")
print(f"(Brick {BRICK_SIZE}x{BRICK_SIZE}x{BRICK_SIZE} with cylindrical hole R={R_HOLE})")
print()

# ============================================================
# Step 1: Create geometry in OCC
# ============================================================
print("Step 1: Create geometry in OCC")

brick = Box(gp_Pnt(-BRICK_SIZE/2, -BRICK_SIZE/2, -BRICK_SIZE/2),
            gp_Pnt(BRICK_SIZE/2, BRICK_SIZE/2, BRICK_SIZE/2))
cyl = Cylinder(gp_Ax2(gp_Pnt(0, 0, -BRICK_SIZE), gp_Dir(0, 0, 1)),
               R_HOLE, 2 * BRICK_SIZE)
shape = brick - cyl

print(f"  OCC faces: {len(shape.faces)}")

# ============================================================
# Step 2: Export STEP from OCC
# ============================================================
print("\nStep 2: Export STEP from OCC")

step_file = os.path.join(work_dir, "complex_named.step")
shape.WriteStep(step_file)
print(f"  Exported: {step_file}")

# ============================================================
# Step 3: Import STEP into Cubit and mesh
# ============================================================
print("\nStep 3: Import and mesh in Cubit")

cubit.init(['cubit', '-nojournal', '-batch'])
cubit.cmd("reset")
cubit.cmd(f'import step "{step_file}" noheal')

# Verify surfaces
surface_ids = cubit.get_entities('surface')
print(f"  Cubit surfaces: {len(surface_ids)}")

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
# Step 4: Export with curving and accuracy test
# ============================================================
print(f"\nStep 4: export_NGSolveCurvedMesh(order={ORDER}) and accuracy")

mesh = cubit_mesh_export.export_NGSolveCurvedMesh(cubit, order=ORDER)

expected_vol = BRICK_SIZE**3 - math.pi * R_HOLE**2 * BRICK_SIZE
vol = Integrate(CF(1), mesh)

print(f"  Expected volume: {expected_vol:.6f}")
print(f"  Computed volume: {vol:.6f}")
print(f"  Error: {abs(vol-expected_vol)/expected_vol*100:.4f}%")

# Compare with order=3
print("\nCompare with order=3:")
mesh3 = cubit_mesh_export.export_NGSolveCurvedMesh(cubit, order=3)
vol3 = Integrate(CF(1), mesh3)
print(f"  Error: {abs(vol3-expected_vol)/expected_vol*100:.4f}%")

# Cleanup
os.remove(step_file)

print("\n=== DONE ===")
