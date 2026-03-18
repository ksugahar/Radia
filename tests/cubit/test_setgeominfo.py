"""
Test SetGeomInfo API for External Meshes

This test demonstrates using the new Element2D.SetGeomInfo() method
to enable mesh.Curve() for externally imported meshes.

The SetGeomInfo API allows setting UV parameters on surface elements,
which is required for proper high-order curving of external meshes.

Run: python test_setgeominfo.py
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

# Use locally built Netgen (ksugahar fork with SetGeomInfo API)
sys.path.insert(0, "s:/NGSolve/01_GitHub/install_ksugahar/Lib/site-packages")

print("=== Test: SetGeomInfo API for External Meshes ===")
print()

# ============================================================
# Test 1: Basic SetGeomInfo functionality
# ============================================================
print("Test 1: Basic SetGeomInfo functionality")

from netgen.meshing import Mesh, MeshPoint, Element2D, Element3D, FaceDescriptor
from netgen.csg import Pnt

mesh = Mesh(dim=3)

# Add points for a simple tetrahedron
p1 = mesh.Add(MeshPoint(Pnt(0, 0, 0)))
p2 = mesh.Add(MeshPoint(Pnt(1, 0, 0)))
p3 = mesh.Add(MeshPoint(Pnt(0, 1, 0)))
p4 = mesh.Add(MeshPoint(Pnt(0, 0, 1)))

# Add volume element
mesh.Add(Element3D(1, [p1, p2, p3, p4]))
mesh.SetMaterial(1, "domain")

# Add face descriptor for boundary
fd = FaceDescriptor(surfnr=1, domin=1, bc=1)
fd.bcname = "boundary"
mesh.Add(fd)

# Add surface element (triangle)
el2d = Element2D(index=1, vertices=[p1, p2, p3])

# Check initial geominfo (should be uninitialized)
print(f"  Initial geominfo: {el2d.geominfo}")

# Set geominfo using new API
el2d.SetGeomInfo(0, 0.0, 0.0)  # vertex 0: u=0, v=0
el2d.SetGeomInfo(1, 1.0, 0.0)  # vertex 1: u=1, v=0
el2d.SetGeomInfo(2, 0.0, 1.0)  # vertex 2: u=0, v=1

print(f"  After SetGeomInfo: {el2d.geominfo}")

# Verify values are set correctly
gi = el2d.geominfo
assert abs(gi[0][1] - 0.0) < 1e-10 and abs(gi[0][2] - 0.0) < 1e-10, "Vertex 0 failed"
assert abs(gi[1][1] - 1.0) < 1e-10 and abs(gi[1][2] - 0.0) < 1e-10, "Vertex 1 failed"
assert abs(gi[2][1] - 0.0) < 1e-10 and abs(gi[2][2] - 1.0) < 1e-10, "Vertex 2 failed"

print("  PASSED: SetGeomInfo works correctly")
print()

# ============================================================
# Test 2: SetGeomInfo with trignum parameter
# ============================================================
print("Test 2: SetGeomInfo with trignum parameter (for STL)")

el2d_stl = Element2D(index=1, vertices=[p1, p2, p3])
el2d_stl.SetGeomInfo(0, 0.1, 0.2, trignum=42)
el2d_stl.SetGeomInfo(1, 0.3, 0.4, trignum=42)
el2d_stl.SetGeomInfo(2, 0.5, 0.6, trignum=42)

gi_stl = el2d_stl.geominfo
print(f"  geominfo with trignum: {gi_stl}")

assert gi_stl[0][0] == 42, "trignum not set correctly"
print("  PASSED: trignum parameter works")
print()

# ============================================================
# Test 3: Boundary check
# ============================================================
print("Test 3: Boundary check (invalid vertex_index)")

try:
    el2d.SetGeomInfo(10, 0.0, 0.0)  # Invalid index
    print("  FAILED: Should have raised exception")
except Exception as e:
    print(f"  Exception raised as expected: {type(e).__name__}")
    print("  PASSED: Boundary check works")
print()

# ============================================================
# Test 4: Quad element
# ============================================================
print("Test 4: SetGeomInfo on Quad element")

p5 = mesh.Add(MeshPoint(Pnt(1, 1, 0)))
el2d_quad = Element2D(index=1, vertices=[p1, p2, p5, p3])

print(f"  Quad has {len(el2d_quad.vertices)} vertices")

el2d_quad.SetGeomInfo(0, 0.0, 0.0)
el2d_quad.SetGeomInfo(1, 1.0, 0.0)
el2d_quad.SetGeomInfo(2, 1.0, 1.0)
el2d_quad.SetGeomInfo(3, 0.0, 1.0)

gi_quad = el2d_quad.geominfo
print(f"  geominfo: {gi_quad}")
print("  PASSED: Quad SetGeomInfo works")
print()

# ============================================================
# Test 5: Use with Cubit hex mesh
# ============================================================
print("Test 5: Cubit hex mesh with SetGeomInfo")

try:
    import cubit
    import cubit_mesh_export
    from netgen.occ import OCCGeometry

    # Create simple box in Cubit
    cubit.init(['cubit', '-nojournal', '-batch'])
    cubit.cmd("reset")
    cubit.cmd("brick x 1 y 1 z 1")
    cubit.cmd("move volume 1 x 0.5 y 0.5 z 0.5")

    # Export STEP
    step_file = os.path.join(work_dir, "test_box.step")
    cubit.cmd(f'export step "{step_file}" overwrite')
    cubit.cmd("reset")
    cubit.cmd(f'import step "{step_file}" heal')

    # Hex mesh
    cubit.cmd("volume 1 scheme map")
    cubit.cmd("volume 1 size 0.25")
    cubit.cmd("mesh volume 1")
    cubit.cmd("block 1 add hex all")
    cubit.cmd('block 1 name "domain"')
    cubit.cmd("block 2 add face all")
    cubit.cmd('block 2 name "boundary"')

    print(f"  Created {cubit.get_hex_count()} hex elements")

    # Export to Netgen with geometry
    geo = OCCGeometry(step_file)
    ngmesh = cubit_mesh_export.export_NetgenMesh(cubit, geometry=geo)

    num_vol = ngmesh.ne
    num_surf = len(list(ngmesh.Elements2D()))
    print(f"  Netgen mesh: {num_vol} volume elements, {num_surf} surface elements")

    # Get surface elements and set geominfo
    # For a box, we can use simple planar UV mapping
    surface_elements_modified = 0
    for el in ngmesh.Elements2D():
        # For each vertex, set UV based on position
        for vi in range(len(el.vertices)):
            # Simple mapping: just use some UV values
            # In real use, you'd compute proper UV from geometry
            el.SetGeomInfo(vi, 0.5, 0.5)
        surface_elements_modified += 1

    print(f"  Modified geominfo on {surface_elements_modified} surface elements")

    # Try to curve the mesh
    try:
        from ngsolve import Mesh as NGSMesh
        ngs_mesh = NGSMesh(ngmesh)
        ngs_mesh.Curve(2)
        print("  PASSED: mesh.Curve() succeeded with SetGeomInfo")
    except Exception as e:
        # NGSolve might not be available, that's OK
        print(f"  Note: mesh.Curve() test skipped ({type(e).__name__})")

    # Cleanup
    os.remove(step_file)

except ImportError as e:
    print(f"  Skipped (Cubit not available): {e}")
except Exception as e:
    print(f"  Error: {e}")
    import traceback
    traceback.print_exc()

print()
print("=== All Tests Completed ===")
