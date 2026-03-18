"""
Exodus Export Example

Demonstrates exporting Cubit mesh to Exodus II format (.exo).
Exodus II is Cubit's native format and supports all features including:
- All element types and orders
- Nodesets and Sidesets
- Block/material definitions
- Large model support (64-bit)

Output files:
    - cube_tet.exo : Simple tetrahedral mesh
    - cube_hex.exo : Simple hexahedral mesh
    - sphere_2nd_order.exo : 2nd order mesh with curved geometry
    - mixed_with_bc.exo : Mixed elements with boundary conditions
"""

import sys
import os

cubit_path = os.environ.get("CUBIT_PATH")
if cubit_path:
	sys.path.append(cubit_path)

import cubit
cubit.init(['cubit', '-nojournal', '-batch'])

# Add parent directory for cubit_mesh_export
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))), 'src', 'radia'))
import cubit_mesh_export

# ============================================================
# Example 1: Simple Tetrahedral Mesh
# ============================================================
print("=" * 60)
print("Example 1: Simple Tetrahedral Mesh")
print("=" * 60)

cubit.cmd("reset")
cubit.cmd("create brick x 2 y 2 z 2")
cubit.cmd("volume 1 scheme tetmesh")
cubit.cmd("volume 1 size 0.5")
cubit.cmd("mesh volume 1")

# Define block with mesh elements
cubit.cmd("block 1 add tet all")
cubit.cmd("block 1 name 'solid'")

cubit_mesh_export.export_exodus(cubit, "cube_tet.exo")
print("  Created: cube_tet.exo\n")

# ============================================================
# Example 2: Hexahedral Mesh
# ============================================================
print("=" * 60)
print("Example 2: Hexahedral Mesh")
print("=" * 60)

cubit.cmd("reset")
cubit.cmd("create brick x 2 y 2 z 2")
cubit.cmd("volume 1 scheme map")
cubit.cmd("volume 1 size 0.5")
cubit.cmd("mesh volume 1")

cubit.cmd("block 1 add hex all")
cubit.cmd("block 1 name 'solid'")

cubit_mesh_export.export_exodus(cubit, "cube_hex.exo")
print("  Created: cube_hex.exo\n")

# ============================================================
# Example 3: 2nd Order Mesh (Sphere)
# ============================================================
print("=" * 60)
print("Example 3: 2nd Order Mesh (Sphere)")
print("=" * 60)

cubit.cmd("reset")
cubit.cmd("create sphere radius 1")
cubit.cmd("volume 1 scheme tetmesh")
cubit.cmd("volume 1 size 0.3")
cubit.cmd("mesh volume 1")

# Add mesh elements to block
cubit.cmd("block 1 add tet all")
cubit.cmd("block 1 name 'sphere'")

# Convert to 2nd order (TET10)
cubit.cmd("block 1 element type tetra10")

# Add boundary surface
cubit.cmd("block 2 add tri all")
cubit.cmd("block 2 name 'surface'")
cubit.cmd("block 2 element type tri6")

cubit_mesh_export.export_exodus(cubit, "sphere_2nd_order.exo")
print("  Created: sphere_2nd_order.exo\n")

# ============================================================
# Example 4: Mixed Elements with Boundary Conditions
# ============================================================
print("=" * 60)
print("Example 4: Mixed Elements with Nodesets/Sidesets")
print("=" * 60)

cubit.cmd("reset")

# Create two adjacent volumes
cubit.cmd("create brick x 1 y 1 z 1")
cubit.cmd("volume 1 move 0.5 0 0")
cubit.cmd("create brick x 1 y 1 z 1")
cubit.cmd("volume 2 move -0.5 0 0")

# Different mesh schemes
cubit.cmd("volume 1 scheme map")
cubit.cmd("volume 2 scheme tetmesh")
cubit.cmd("volume all size 0.25")
cubit.cmd("mesh volume all")

# Define material blocks
cubit.cmd("block 1 add hex all")
cubit.cmd("block 1 name 'material_1'")
cubit.cmd("block 2 add tet all")
cubit.cmd("block 2 name 'material_2'")
cubit.cmd("block 3 add pyramid all")
cubit.cmd("block 3 name 'interface'")

# Define nodesets for boundary conditions
cubit.cmd("nodeset 1 add node in surface 1")  # Left face
cubit.cmd("nodeset 1 name 'fixed'")
cubit.cmd("nodeset 2 add node in surface 7")  # Right face
cubit.cmd("nodeset 2 name 'load'")

# Define sidesets for surface loads
cubit.cmd("sideset 1 add surface 1")
cubit.cmd("sideset 1 name 'inlet'")
cubit.cmd("sideset 2 add surface 7")
cubit.cmd("sideset 2 name 'outlet'")

cubit_mesh_export.export_exodus(cubit, "mixed_with_bc.exo")
print("  Created: mixed_with_bc.exo\n")

# ============================================================
# Example 5: Large Model Option
# ============================================================
print("=" * 60)
print("Example 5: Large Model (64-bit integers)")
print("=" * 60)

cubit.cmd("reset")
cubit.cmd("create brick x 5 y 5 z 5")
cubit.cmd("volume 1 scheme tetmesh")
cubit.cmd("volume 1 size 0.3")
cubit.cmd("mesh volume 1")

cubit.cmd("block 1 add tet all")
cubit.cmd("block 1 name 'large_model'")

# Use large_model=True for meshes with >2^31 elements/nodes
cubit_mesh_export.export_exodus(cubit, "large_model.exo", large_model=True)
print("  Created: large_model.exo (with 64-bit integers)\n")

print("=" * 60)
print("All Exodus exports completed!")
print("=" * 60)
print("\nGenerated files:")
print("  - cube_tet.exo         : Simple tet mesh")
print("  - cube_hex.exo         : Simple hex mesh")
print("  - sphere_2nd_order.exo : 2nd order curved mesh")
print("  - mixed_with_bc.exo    : Mixed mesh with nodesets/sidesets")
print("  - large_model.exo      : Large model with 64-bit integers")
print("\nThese files can be opened with:")
print("  - ParaView (with Exodus plugin)")
print("  - Cubit (import mesh)")
print("  - SIERRA tools")
print("  - VisIt")
