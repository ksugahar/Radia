"""
Test: SetGeomInfo API with UV parameters from OCC Geometry

This test attempts to use the SetGeomInfo API to set proper UV parameters
from the OCC geometry, enabling mesh.Curve() to work correctly.

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

# Use locally built NGSolve (ksugahar fork with SetGeomInfo API)
sys.path.insert(0, "s:/NGSolve/01_GitHub/install_ksugahar/Lib/site-packages")

from netgen.occ import OCCGeometry
from netgen.meshing import Mesh as NetgenMesh, MeshPoint, Element2D, FaceDescriptor
from netgen.csg import Pnt
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
# Step 2: Load OCC geometry and examine faces
# ============================================================
print("\nStep 2: Load OCC geometry")

geo = OCCGeometry(step_file)
print(f"  Geometry loaded from: {step_file}")

# Get the shape and faces
shape = geo.shape
print(f"  Shape type: {type(shape)}")

# Try to get faces from the shape
try:
    faces = shape.faces
    print(f"  Number of faces: {len(faces)}")
    for i, face in enumerate(faces):
        print(f"    Face {i}: {face.name if hasattr(face, 'name') else 'unnamed'}")
except Exception as e:
    print(f"  Could not get faces: {e}")
    faces = []

# ============================================================
# Step 3: Export mesh to Netgen
# ============================================================
print("\nStep 3: Export mesh to Netgen")

ngmesh = cubit_mesh_export.export_netgen(cubit, geometry=geo)
print(f"  Elements: {ngmesh.ne}")

# ============================================================
# Step 4: Examine surface elements and their geominfo
# ============================================================
print("\nStep 4: Examine surface elements")

el2d_list = list(ngmesh.Elements2D())
print(f"  Surface elements: {len(el2d_list)}")

# Check current geominfo state
sample_el = el2d_list[0] if el2d_list else None
if sample_el:
    print(f"  Sample element vertices: {sample_el.vertices}")
    print(f"  Sample element index (face): {sample_el.index}")
    gi = sample_el.geominfo
    print(f"  Sample geominfo: {gi}")

    # Check if geominfo looks valid
    is_valid = all(abs(g[1]) < 1e10 and abs(g[2]) < 1e10 for g in gi)
    print(f"  Geominfo appears valid: {is_valid}")

# ============================================================
# Step 5: Try to set UV from geometry (experimental)
# ============================================================
print("\nStep 5: Attempt to set UV from geometry")

# For a cylinder, we can compute UV analytically:
# - For cylindrical surface: u = angle (theta), v = height (z)
# - For circular caps: u, v are based on x, y coordinates

def compute_cylinder_uv(x, y, z, R, H, face_type):
    """
    Compute UV parameters for a cylinder surface.
    face_type: 'side', 'top', 'bottom'
    """
    if face_type == 'side':
        # Cylindrical surface: u = angle, v = normalized height
        theta = math.atan2(y, x)
        if theta < 0:
            theta += 2 * math.pi
        u = theta / (2 * math.pi)  # Normalize to [0, 1]
        v = (z + H/2) / H  # Normalize height to [0, 1]
    elif face_type == 'top':
        # Top cap: use x, y as UV
        u = (x + R) / (2 * R)
        v = (y + R) / (2 * R)
    else:  # bottom
        u = (x + R) / (2 * R)
        v = (y + R) / (2 * R)
    return u, v

# Classify elements by face
def classify_face(vertices, points, R, H):
    """Determine if a face is on the side, top, or bottom of cylinder."""
    zs = [points[v.nr-1][2] for v in vertices]  # Use .nr to get integer
    avg_z = sum(zs) / len(zs)

    # Check if all vertices are at approximately z = H/2 or z = -H/2
    if all(abs(z - H/2) < 0.01 for z in zs):
        return 'top'
    elif all(abs(z + H/2) < 0.01 for z in zs):
        return 'bottom'
    else:
        return 'side'

# Get mesh points
points = [(p.p[0], p.p[1], p.p[2]) for p in ngmesh.Points()]

# Set geominfo for surface elements
modified_count = 0
for el in ngmesh.Elements2D():
    verts = list(el.vertices)
    face_type = classify_face(verts, points, R, H)

    for i, v in enumerate(verts):
        x, y, z = points[v.nr - 1]  # Use .nr to get integer
        u, v_param = compute_cylinder_uv(x, y, z, R, H, face_type)
        try:
            el.SetGeomInfo(i, u, v_param)
            modified_count += 1
        except Exception as e:
            print(f"  Error setting geominfo: {e}")
            break

print(f"  Modified {modified_count} vertex geominfo entries")

# Check geominfo after modification
if sample_el:
    gi_after = el2d_list[0].geominfo
    print(f"  Sample geominfo after: {gi_after}")

# ============================================================
# Step 6: Test mesh.Curve() with modified geominfo
# ============================================================
print("\nStep 6: Test mesh.Curve() with modified geominfo")

mesh_with_gi = Mesh(ngmesh)

print("  Testing mesh.Curve(1)...")
mesh_with_gi.Curve(1)
area1 = Integrate(CF(1), mesh_with_gi, VOL_or_BND=BND)
vol1 = Integrate(CF(1), mesh_with_gi)
print(f"    Area: {area1:.6f} (error: {abs(area1-expected_area)/expected_area*100:.4f}%)")
print(f"    Vol:  {vol1:.6f} (error: {abs(vol1-expected_vol)/expected_vol*100:.4f}%)")

print("  Testing mesh.Curve(2)...")
mesh_curve2 = Mesh(ngmesh)
try:
    mesh_curve2.Curve(2)
    area2 = Integrate(CF(1), mesh_curve2, VOL_or_BND=BND)
    vol2 = Integrate(CF(1), mesh_curve2)
    print(f"    Area: {area2:.6f} (error: {abs(area2-expected_area)/expected_area*100:.4f}%)")
    print(f"    Vol:  {vol2:.6f} (error: {abs(vol2-expected_vol)/expected_vol*100:.4f}%)")
except Exception as e:
    print(f"    Failed: {e}")

# ============================================================
# Step 7: Compare with SetDeformation
# ============================================================
print("\nStep 7: Compare with SetDeformation (baseline)")

mesh_deform = Mesh(cubit_mesh_export.export_netgen(cubit, geometry=geo))
mesh_deform.Curve(1)

if hasattr(cubit_mesh_export, 'detect_cylinder_boundaries'):
    cyl_boundaries = cubit_mesh_export.detect_cylinder_boundaries(mesh_deform, R, axis='z')
    cubit_mesh_export.apply_cylinder_deformation(
        mesh_deform, radius=R, boundary_names=cyl_boundaries, order=2, axis='z'
    )

    area_deform = Integrate(CF(1), mesh_deform, VOL_or_BND=BND)
    vol_deform = Integrate(CF(1), mesh_deform)
    print(f"  Area: {area_deform:.6f} (error: {abs(area_deform-expected_area)/expected_area*100:.4f}%)")
    print(f"  Vol:  {vol_deform:.6f} (error: {abs(vol_deform-expected_vol)/expected_vol*100:.4f}%)")
else:
    print("  SKIPPED: detect_cylinder_boundaries not yet implemented")

# Cleanup
os.remove(step_file)

print("\n" + "=" * 60)
print("Test Complete")
print("=" * 60)
print()
print("Conclusion:")
print("  SetGeomInfo API allows setting UV parameters, but proper UV")
print("  computation from OCC geometry requires more work.")
print("  SetDeformation remains the practical solution for now.")
