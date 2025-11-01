import os, sys
import numpy as np
sys.path.append("C:/Program Files/Coreform Cubit 2025.3/bin")

import cubit
cubit.init(['cubit','-nojournal','-batch'])

# Strategy: Decompose sphere into convex sectors for quadrupole symmetry
# Quadrupole field has 4-fold rotational symmetry
# We'll divide into 16 sectors: 8 azimuthal × 2 polar (upper/lower)

cubit.cmd('reset')

# Create sphere
cubit.cmd('create sphere radius 10')

# Method: Create 8 octants by orthogonal cuts
# This provides sufficient resolution for quadrupole field analysis
# Each octant is a convex sector

# Step 1: Cut by X plane (creates 2 halves)
cubit.cmd('webcut volume all with plane xplane offset 0')

# Step 2: Cut by Y plane (creates 4 quadrants)
cubit.cmd('webcut volume all with plane yplane offset 0')

# Step 3: Cut by Z plane (creates 8 octants)
cubit.cmd('webcut volume all with plane zplane offset 0')

# Check how many volumes we have
num_volumes = cubit.get_volume_count()
print(f"\nTotal number of octants created: {num_volumes}")
print(f"Note: 8 octants provide sufficient symmetry for quadrupole field")
# Set mesh size (element edge length ~1 mm for better accuracy)
cubit.cmd('volume all size 1')

# Use tetrahedral mesh (tetrahedra are always convex, required for Radia)
cubit.cmd('volume all scheme tetmesh')
cubit.cmd('mesh volume all')

# Create separate blocks for each octant (each octant = one material group)
# This allows each octant to be treated as a separate convex polyhedron
num_volumes = cubit.get_volume_count()
print(f"Number of volumes (octants): {num_volumes}")

for i in range(1, num_volumes + 1):
	cubit.cmd(f'block {i} add tri in volume {i}')
	cubit.cmd(f'block {i} name "octant_{i}"')

FileName = 'sphere'
import cubit_mesh_export

#cubit_mesh_export.export_Gmsh_ver2(cubit, FileName + '.msh')
#
cubit_mesh_export.export_Nastran(cubit, FileName + '.bdf', DIM='3D', PYRAM=False)
#
#cubit_mesh_export.export_meg(cubit, FileName + '.meg', DIM='T', MGR2=[])
#
#cubit_mesh_export.export_vtk(cubit, FileName + '.vtk', ORDER="2nd")
