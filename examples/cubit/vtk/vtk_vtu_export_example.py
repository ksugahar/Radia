"""
Example: VTK and VTU mesh export

This example demonstrates:
1. Creating 1st and 2nd order meshes in Cubit
2. Exporting to VTK Legacy format (.vtk)
3. Exporting to VTK XML format (.vtu)
4. Auto-detection of element order
5. BlockID and NodeID data arrays

Output files can be visualized in ParaView.

Run this script:
    python vtk_vtu_export_example.py
"""

import sys
import os
cubit_path = os.environ.get("CUBIT_PATH")
if cubit_path:
	sys.path.append(cubit_path)

# Add parent directory to path for cubit_mesh_export
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))), 'src', 'radia'))

import cubit
import cubit_mesh_export

# Initialize Cubit
cubit.init(['cubit', '-nojournal', '-batch'])

# ============================================================
# Example 1: 1st Order Tet Mesh
# ============================================================
print("=" * 60)
print("Example 1: 1st Order Tet Mesh")
print("=" * 60)

cubit.cmd("reset")
cubit.cmd("create sphere radius 1")
cubit.cmd("volume 1 scheme tetmesh")
cubit.cmd("volume 1 size 0.3")
cubit.cmd("mesh volume 1")

cubit.cmd("block 1 add tet all in volume 1")
cubit.cmd("block 1 name 'sphere_volume'")
cubit.cmd("block 2 add tri all in surface all")
cubit.cmd("block 2 name 'sphere_surface'")

num_tets = len(cubit.get_block_tets(1))
print(f"  Mesh: {num_tets} tetrahedra (1st order, 4 nodes each)")

# Export to VTK Legacy
cubit_mesh_export.export_vtk(cubit, "sphere_1st_order.vtk")
print("  Exported: sphere_1st_order.vtk")

# Export to VTK XML (VTU)
cubit_mesh_export.export_vtu(cubit, "sphere_1st_order.vtu")
print("  Exported: sphere_1st_order.vtu")

# ============================================================
# Example 2: 2nd Order Tet Mesh
# ============================================================
print("\n" + "=" * 60)
print("Example 2: 2nd Order Tet Mesh (Auto-detected)")
print("=" * 60)

cubit.cmd("reset")
cubit.cmd("create sphere radius 1")
cubit.cmd("volume 1 scheme tetmesh")
cubit.cmd("volume 1 size 0.4")
cubit.cmd("mesh volume 1")

cubit.cmd("block 1 add tet all in volume 1")
cubit.cmd("block 1 name 'sphere_volume'")
cubit.cmd("block 1 element type tetra10")  # Convert to 2nd order

cubit.cmd("block 2 add tri all in surface all")
cubit.cmd("block 2 name 'sphere_surface'")
cubit.cmd("block 2 element type tri6")  # Convert to 2nd order

num_tets = len(cubit.get_block_tets(1))
print(f"  Mesh: {num_tets} tetrahedra (2nd order, 10 nodes each)")

# Export - order is auto-detected
cubit_mesh_export.export_vtk(cubit, "sphere_2nd_order.vtk")
print("  Exported: sphere_2nd_order.vtk (auto-detected 2nd order)")

cubit_mesh_export.export_vtu(cubit, "sphere_2nd_order.vtu")
print("  Exported: sphere_2nd_order.vtu (auto-detected 2nd order)")

# ============================================================
# Example 3: Mixed Element Types
# ============================================================
print("\n" + "=" * 60)
print("Example 3: Mixed Element Types")
print("=" * 60)

cubit.cmd("reset")

# Hex region
cubit.cmd("brick x 1 y 1 z 1")
cubit.cmd("volume 1 move 0 0 0.5")
cubit.cmd("volume 1 scheme map")
cubit.cmd("volume 1 size 0.3")
cubit.cmd("mesh volume 1")

# Tet region
cubit.cmd("brick x 1 y 1 z 1")
cubit.cmd("volume 2 move 0 0 1.5")
cubit.cmd("volume 2 scheme tetmesh")
cubit.cmd("volume 2 size 0.3")
cubit.cmd("mesh volume 2")

cubit.cmd("block 1 add hex all")
cubit.cmd("block 1 name 'hex_region'")
cubit.cmd("block 2 add tet all")
cubit.cmd("block 2 name 'tet_region'")
cubit.cmd("block 3 add quad all in surface all")
cubit.cmd("block 3 name 'quad_boundary'")
cubit.cmd("block 4 add tri all in surface all")
cubit.cmd("block 4 name 'tri_boundary'")

num_hex = len(cubit.get_block_hexes(1))
num_tet = len(cubit.get_block_tets(2))
print(f"  Mesh: {num_hex} hexahedra + {num_tet} tetrahedra")

cubit_mesh_export.export_vtk(cubit, "mixed_elements.vtk")
print("  Exported: mixed_elements.vtk")

cubit_mesh_export.export_vtu(cubit, "mixed_elements.vtu")
print("  Exported: mixed_elements.vtu")

# ============================================================
# Example 4: 2D Surface Mesh
# ============================================================
print("\n" + "=" * 60)
print("Example 4: 2D Surface Mesh")
print("=" * 60)

cubit.cmd("reset")
cubit.cmd("create surface circle radius 1 zplane")
cubit.cmd("surface 1 scheme trimesh")
cubit.cmd("surface 1 size 0.2")
cubit.cmd("mesh surface 1")

cubit.cmd("block 1 add tri all")
cubit.cmd("block 1 name 'circle'")
cubit.cmd("block 2 add edge all in curve all")
cubit.cmd("block 2 name 'boundary'")

num_tri = len(cubit.get_block_tris(1))
print(f"  Mesh: {num_tri} triangles")

cubit_mesh_export.export_vtk(cubit, "circle_2d.vtk")
print("  Exported: circle_2d.vtk")

cubit_mesh_export.export_vtu(cubit, "circle_2d.vtu")
print("  Exported: circle_2d.vtu")

# ============================================================
# Summary
# ============================================================
print("\n" + "=" * 60)
print("Summary")
print("=" * 60)
print("Generated files:")
print("  - sphere_1st_order.vtk / .vtu  (1st order tet)")
print("  - sphere_2nd_order.vtk / .vtu  (2nd order tet)")
print("  - mixed_elements.vtk / .vtu    (hex + tet)")
print("  - circle_2d.vtk / .vtu         (2D triangles)")
print()
print("VTK Legacy (.vtk):")
print("  - ASCII format, human-readable")
print("  - Compatible with older VTK versions")
print()
print("VTK XML (.vtu):")
print("  - Modern XML format")
print("  - Includes BlockID (cell data) and NodeID (point data)")
print("  - Better ParaView integration")
print()
print("Open in ParaView to visualize!")
print("=" * 60)
