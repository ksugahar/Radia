#!/usr/bin/env python
"""
Debug single tetrahedron with centroid charge implementation
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../build/Release'))

import numpy as np
import radia as rad

rad.FldUnits('m')

# Create single tetrahedron manually
# Vertices of a simple tetrahedron (0-indexed)
v0 = [-0.05, -0.05, -0.05]
v1 = [0.05, -0.05, -0.05]
v2 = [0, 0.05, -0.05]
v3 = [0, 0, 0.05]

vertices = [v0, v1, v2, v3]

# Face definitions (1-indexed for Radia, from netgen_mesh_import.py)
# TETRA_FACES: correct winding for outward normals
faces = [
    [1, 3, 2],  # Face 0: v0-v2-v1
    [1, 2, 4],  # Face 1: v0-v1-v3
    [2, 3, 4],  # Face 2: v1-v2-v3
    [3, 1, 4]   # Face 3: v2-v0-v3
]

# Magnetization
M = [0, 0, 1.2]

# Create tetrahedron using polyhedron
print("Creating single tetrahedron...")
try:
    tetra = rad.ObjPolyhdr(vertices, faces, M)
    print(f"  Created: ID = {tetra}")

    # Test point
    test_pt = [0.08, 0, 0]

    print(f"\nComputing field at {test_pt}...")
    H = rad.Fld(tetra, 'h', test_pt)
    print(f"  H = {H}")
    print(f"  |H| = {np.linalg.norm(H)}")

except Exception as e:
    print(f"Error: {e}")
    import traceback
    traceback.print_exc()

print("\nTest complete!")
