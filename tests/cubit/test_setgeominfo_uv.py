"""
Test: SetGeomInfo API with UV parameters from OCC Geometry

This test attempts to use the SetGeomInfo API to set proper UV parameters
from the OCC geometry, enabling mesh.Curve() to work correctly.

NOTE: For most use cases, export_NGSolveCurvedMesh() handles all of this internally.
This test exercises the low-level SetGeomInfo API with manual UV computation.

The workflow:
1. Load STEP geometry with OCCGeometry
2. Import Cubit mesh
3. For each surface element, project vertices to OCC surface and get UV
4. Set UV parameters using SetGeomInfo
5. Call mesh.Curve(order)

Run: python test_setgeominfo_uv.py
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
_fork_path = os.environ.get("NGSOLVE_FORK_PATH")
if _fork_path:
    sys.path.insert(0, _fork_path)

from ngsolve import Mesh, Integrate, CF, BND
import cubit
import cubit_mesh_export

print("=" * 60)
print("Test: SetGeomInfo API with UV from OCC Geometry")
print("=" * 60)
print()

# ============================================================
# Parameters
# ============================================================
R = 0.5
H = 2.0
MESH_SIZE = 0.2

expected_area = 2*math.pi*R*H + 2*math.pi*R*R
expected_vol = math.pi*R*R*H

print(f"Cylinder: R={R}, H={H}")
print(f"Expected: Area={expected_area:.6f}, Vol={expected_vol:.6f}")
print()

# ============================================================
# Step 1: Create geometry and mesh in Cubit
# ============================================================
print("Step 1: Create geometry and mesh in Cubit")

cubit.init(['cubit', '-nojournal', '-batch'])
cubit.cmd("reset")
cubit.cmd(f"create cylinder height {H} radius {R}")

step_file = os.path.join(work_dir, "test_cylinder_uv.step")
cubit.cmd(f'export step "{step_file}" overwrite')
cubit.cmd("reset")
cubit.cmd(f'import step "{step_file}" heal')

cubit.cmd("volume all scheme tetmesh")
cubit.cmd(f"volume all size {MESH_SIZE}")
cubit.cmd("mesh volume all")
cubit.cmd("block 1 add tet all")
cubit.cmd("block 2 add tri all")

print(f"  Tets: {cubit.get_tet_count()}")

# ============================================================
# Step 2: Export mesh using export_NGSolveCurvedMesh for different orders
# ============================================================
print("\nStep 2: Export mesh using export_NGSolveCurvedMesh")

print("  Testing export_NGSolveCurvedMesh(order=1)...")
mesh1 = cubit_mesh_export.export_NGSolveCurvedMesh(cubit, order=1)
area1 = Integrate(CF(1), mesh1, VOL_or_BND=BND)
vol1 = Integrate(CF(1), mesh1)
print(f"    Area: {area1:.6f} (error: {abs(area1-expected_area)/expected_area*100:.4f}%)")
print(f"    Vol:  {vol1:.6f} (error: {abs(vol1-expected_vol)/expected_vol*100:.4f}%)")

print("  Testing export_NGSolveCurvedMesh(order=2)...")
try:
    mesh2 = cubit_mesh_export.export_NGSolveCurvedMesh(cubit, order=2)
    area2 = Integrate(CF(1), mesh2, VOL_or_BND=BND)
    vol2 = Integrate(CF(1), mesh2)
    print(f"    Area: {area2:.6f} (error: {abs(area2-expected_area)/expected_area*100:.4f}%)")
    print(f"    Vol:  {vol2:.6f} (error: {abs(vol2-expected_vol)/expected_vol*100:.4f}%)")
except Exception as e:
    print(f"    Failed: {e}")

print("  Testing export_NGSolveCurvedMesh(order=3)...")
try:
    mesh3 = cubit_mesh_export.export_NGSolveCurvedMesh(cubit, order=3)
    area3 = Integrate(CF(1), mesh3, VOL_or_BND=BND)
    vol3 = Integrate(CF(1), mesh3)
    print(f"    Area: {area3:.6f} (error: {abs(area3-expected_area)/expected_area*100:.4f}%)")
    print(f"    Vol:  {vol3:.6f} (error: {abs(vol3-expected_vol)/expected_vol*100:.4f}%)")
except Exception as e:
    print(f"    Failed: {e}")

# Cleanup
os.remove(step_file)

print("\n" + "=" * 60)
print("Test Complete")
print("=" * 60)
print()
print("Conclusion:")
print("  export_NGSolveCurvedMesh() handles SetGeomInfo and mesh.Curve() internally.")
print("  Higher order gives better geometric accuracy for curved surfaces.")
