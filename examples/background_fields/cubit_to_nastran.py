#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Generate tetrahedral mesh of a sphere using Coreform Cubit.

Creates a sphere (radius 10 mm), decomposes it into 8 convex octants,
meshes with tetrahedra (~1 mm element size), and exports to Nastran .bdf format.

Requires:
    - Coreform Cubit (adjust CUBIT_PATH below for your installation)
    - cubit_mesh_export module
"""

import os
import sys

CUBIT_PATH = "C:/Program Files/Coreform Cubit 2025.3/bin"
sys.path.append(CUBIT_PATH)

import cubit
cubit.init(['cubit', '-nojournal', '-batch'])

cubit.cmd('reset')

# Create sphere
cubit.cmd('create sphere radius 10')

# Decompose into 8 convex octants by orthogonal cuts.
# Each octant is convex, which is required for Radia polyhedra.
cubit.cmd('webcut volume all with plane xplane offset 0')
cubit.cmd('webcut volume all with plane yplane offset 0')
cubit.cmd('webcut volume all with plane zplane offset 0')

num_volumes = cubit.get_volume_count()
print(f"\nTotal number of octants created: {num_volumes}")

# Set mesh size (element edge length ~1 mm for better accuracy)
cubit.cmd('volume all size 1')

# Use tetrahedral mesh (tetrahedra are always convex, required for Radia)
cubit.cmd('volume all scheme tetmesh')
cubit.cmd('mesh volume all')

# Create separate blocks for each octant (each octant = one material group)
for i in range(1, num_volumes + 1):
	cubit.cmd(f'block {i} add tri in volume {i}')
	cubit.cmd(f'block {i} name "octant_{i}"')

# Export to Nastran format
output_dir = os.path.dirname(os.path.abspath(__file__))
filename = os.path.join(output_dir, 'sphere')

import cubit_mesh_export
cubit_mesh_export.export_Nastran(cubit, filename + '.bdf', DIM='3D', PYRAM=False)
print(f"Exported: {filename}.bdf")
