"""
MEG Export Example

Demonstrates exporting Cubit mesh to ELF/MAGIC MEG format.

Output files:
    - cube_tet.meg : 3D tetrahedral mesh (DIM='T')
    - sphere_hex.meg : 3D hexahedral mesh (DIM='T')
    - plate_2d.meg : 2D planar mesh (DIM='K')
    - axisym.meg : Axisymmetric mesh (DIM='R')
    - with_spatial_nodes.meg : 3D mesh with MGR2 spatial nodes

ELF element names (magnetic material elements):
    3D (DIM='T'):
        MMB8T - 8-node hexahedron
        MMB6T - 6-node wedge/prism
        MMB4T - 4-node tetrahedron
    2D planar (DIM='K'):
        MMB4K - 4-node quadrilateral
        MMB3K - 3-node triangle
    Axisymmetric (DIM='R'):
        MMB4R - 4-node quadrilateral
        MMB3R - 3-node triangle
"""

import sys
import os

cubit_path = os.environ.get("CUBIT_PATH")
if cubit_path:
	sys.path.append(cubit_path)

import cubit
cubit.init(['cubit', '-nojournal', '-batch'])

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))), 'src', 'radia'))
import cubit_mesh_export

# ============================================================
# 3D Tetrahedral Mesh Export
# ============================================================
print("Creating 3D tetrahedral mesh...")
cubit.cmd("reset")
cubit.cmd("create brick x 2 y 2 z 2")
cubit.cmd("volume 1 scheme tetmesh")
cubit.cmd("volume 1 size 0.5")
cubit.cmd("mesh volume 1")

cubit.cmd("block 1 add tet all")
cubit.cmd("block 1 name 'MMB4T'")  # ELF magnetic material element (4-node tet)

print("\nExporting to MEG (3D tet, DIM='T')...")
cubit_mesh_export.export_meg(cubit, "cube_tet.meg", DIM='T')
print("  Created: cube_tet.meg")

# ============================================================
# 3D Hexahedral Mesh Export (Sphere)
# ============================================================
print("\nCreating 3D hexahedral mesh (sphere)...")
cubit.cmd("reset")
cubit.cmd("create sphere radius 1")
cubit.cmd("volume 1 scheme polyhedron")
cubit.cmd("volume 1 size 0.3")
cubit.cmd("mesh volume 1")

cubit.cmd("block 1 add hex all")
cubit.cmd("block 1 name 'MMB8T'")  # ELF magnetic material element (8-node hex)

print("\nExporting to MEG (3D hex, DIM='T')...")
cubit_mesh_export.export_meg(cubit, "sphere_hex.meg", DIM='T')
print("  Created: sphere_hex.meg")

# ============================================================
# 2D Planar Mesh Export (XY plane)
# ============================================================
print("\nCreating 2D planar mesh (XY plane)...")
cubit.cmd("reset")
cubit.cmd("create surface rectangle width 2 height 2 zplane")
cubit.cmd("surface 1 scheme trimesh")
cubit.cmd("surface 1 size 0.3")
cubit.cmd("mesh surface 1")

cubit.cmd("block 1 add tri all")
cubit.cmd("block 1 name 'MMB3K'")  # ELF magnetic material element (3-node tri, 2D)

print("\nExporting to MEG (2D planar, DIM='K')...")
cubit_mesh_export.export_meg(cubit, "plate_2d.meg", DIM='K')
print("  Created: plate_2d.meg")

# ============================================================
# Axisymmetric Mesh Export (XZ plane, Y=0)
# ============================================================
print("\nCreating axisymmetric mesh (XZ plane)...")
cubit.cmd("reset")
cubit.cmd("create surface rectangle width 2 height 2 yplane")  # XZ plane (Y=0)
cubit.cmd("surface 1 move 1 0 1")  # Move to positive X quadrant (r > 0)
cubit.cmd("surface 1 scheme trimesh")
cubit.cmd("surface 1 size 0.3")
cubit.cmd("mesh surface 1")

cubit.cmd("block 1 add tri all")
cubit.cmd("block 1 name 'MMB3R'")  # ELF magnetic material element (3-node tri, axisym)

print("\nExporting to MEG (axisymmetric, DIM='R')...")
cubit_mesh_export.export_meg(cubit, "axisym.meg", DIM='R')
print("  Created: axisym.meg")

# ============================================================
# With Spatial Nodes (MGR2)
# ============================================================
print("\nCreating 3D hex mesh with spatial nodes...")
cubit.cmd("reset")
cubit.cmd("create sphere radius 1")
cubit.cmd("volume 1 scheme polyhedron")
cubit.cmd("volume 1 size 0.3")
cubit.cmd("mesh volume 1")

cubit.cmd("block 1 add hex all")
cubit.cmd("block 1 name 'MMB8T'")  # ELF magnetic material element (8-node hex)

spatial_nodes = [
    [0.0, 0.0, 2.0],    # Spatial node (near sphere)
    [0.0, 0.0, -2.0],   # Spatial node (near sphere)
]
cubit_mesh_export.export_meg(cubit, "with_spatial_nodes.meg", DIM='T', MGR2=spatial_nodes)
print("  Created: with_spatial_nodes.meg (with MGR2 spatial nodes)")

print("\nDone!")
