"""
Gmsh Export Example

Demonstrates exporting Cubit mesh to Gmsh format (v2.2 and v4.1)
using the radia_export gmsh command.

Output files:
    - cube_v2.msh : Gmsh v2.2 format (simple flat structure)
    - cube_v4.msh : Gmsh v4.1 format (with $Entities section)
"""

import sys
import os

# Auto-detect Cubit installation
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'src', 'radia'))
from install_panels import find_cubit_bin as _fcb
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

# Export to Gmsh v2.2 (default)
print("\nExporting to Gmsh v2.2...")
cubit.cmd('radia_export gmsh "cube_v2.msh" version 2 overwrite')
print("  Created: cube_v2.msh")

# Export to Gmsh v4.1
print("\nExporting to Gmsh v4.1...")
cubit.cmd('radia_export gmsh "cube_v4.msh" version 4 overwrite')
print("  Created: cube_v4.msh")

# 2nd order example
print("\nCreating 2nd order mesh...")
cubit.cmd("block 1 element type tetra10")
cubit.cmd("block 2 element type tri6")
cubit.cmd('radia_export gmsh "cube_2nd_order.msh" version 2 overwrite')
print("  Created: cube_2nd_order.msh")

print("\nDone!")
