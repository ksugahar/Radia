"""
Nastran Export Example

Demonstrates exporting Cubit mesh to Nastran BDF format.

Output files:
    - cube_3d.bdf : 3D solid mesh
    - plate_2d.bdf : 2D shell mesh
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

# ============================================================
# 3D Mesh Export
# ============================================================
print("Creating 3D mesh...")
cubit.cmd("reset")
cubit.cmd("create brick x 2 y 2 z 2")
cubit.cmd("volume 1 scheme tetmesh")
cubit.cmd("volume 1 size 0.5")
cubit.cmd("mesh volume 1")

cubit.cmd("block 1 add tet all")
cubit.cmd("block 1 name 'solid'")

print("\nExporting to Nastran (3D)...")
cubit.cmd('export jmag_nastran "cube_3d.bdf" dimension 3 overwrite')
print("  Created: cube_3d.bdf")

# ============================================================
# 2D Mesh Export
# ============================================================
print("\nCreating 2D mesh...")
cubit.cmd("reset")
cubit.cmd("create surface rectangle width 2 height 2 zplane")
cubit.cmd("surface 1 scheme trimesh")
cubit.cmd("surface 1 size 0.3")
cubit.cmd("mesh surface 1")

cubit.cmd("block 1 add tri all")
cubit.cmd("block 1 name 'plate'")

print("\nExporting to Nastran (2D)...")
cubit.cmd('export jmag_nastran "plate_2d.bdf" dimension 2 overwrite')
print("  Created: plate_2d.bdf")

# ============================================================
# Mixed Hex/Tet with Pyramid handling
# ============================================================
print("\nCreating mixed hex/tet mesh...")
cubit.cmd("reset")
cubit.cmd("create brick x 1 y 1 z 1")
cubit.cmd("volume 1 move 0.5 0 0")
cubit.cmd("create brick x 1 y 1 z 1")
cubit.cmd("volume 2 move -0.5 0 0")
cubit.cmd("volume 1 scheme map")
cubit.cmd("volume 2 scheme tetmesh")
cubit.cmd("volume all size 0.25")
cubit.cmd("mesh volume all")

cubit.cmd("block 1 add hex all")
cubit.cmd("block 1 add tet all")
cubit.cmd("block 1 add pyramid all")
cubit.cmd("block 1 name 'mixed'")

# Export with pyramid as CPYRAM
print("\nExporting with pyramids...")
cubit.cmd('export jmag_nastran "mixed_with_pyramid.bdf" dimension 3 overwrite')
print("  Created: mixed_with_pyramid.bdf")

# Export with pyramid as degenerate hex (for JMAG compatibility)
print("\nExporting without pyramids (for JMAG)...")
cubit.cmd('export jmag_nastran "mixed_for_jmag.bdf" dimension 3 nopyramid overwrite')
print("  Created: mixed_for_jmag.bdf")

print("\nDone!")
