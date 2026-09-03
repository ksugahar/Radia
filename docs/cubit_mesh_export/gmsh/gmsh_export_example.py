"""
Gmsh Export Example

Demonstrates exporting Cubit mesh to Gmsh v4.1 format
using the export gmsh command.

Output files:
    - cube.msh            : Gmsh v4.1 format (1st order)
    - cube_2nd_order.msh  : Gmsh v4.1 format (2nd order)
"""

import sys

# Auto-detect Cubit installation
from radia.install_panels import find_cubit_bin as _fcb
_cubit_path = _fcb()
if _cubit_path and _cubit_path not in sys.path:
    sys.path.append(_cubit_path)

import cubit
cubit.init(['cubit', '-nojournal', '-batch'])

# Create geometry and mesh
print("Creating geometry and mesh...")
cubit.cmd("reset")
cubit.cmd("create brick x 2 y 2 z 2")
cubit.cmd("volume 1 scheme tetmesh")
cubit.cmd("volume 1 size 0.5")
cubit.cmd("mesh volume 1")

# Define blocks
cubit.cmd("block 1 add tet all")
cubit.cmd("block 1 name 'solid'")
cubit.cmd("block 2 add tri all")
cubit.cmd("block 2 name 'boundary'")

# Export to Gmsh v4.1 (default)
print("\nExporting to Gmsh v4.1...")
cubit.cmd('export gmsh "cube.msh" overwrite')
print("  Created: cube.msh")

# 2nd order example
print("\nCreating 2nd order mesh...")
cubit.cmd("block 1 element type tetra10")
cubit.cmd("block 2 element type tri6")
cubit.cmd('export gmsh "cube_2nd_order.msh" overwrite')
print("  Created: cube_2nd_order.msh")

print("\nDone!")
