#!/usr/bin/env python
"""
Case 3: Polyhedron (Cube) Creation
Converted from Mathematica/Wolfram Language to Python
"""

import sys
import os
import math
import numpy as np

# Add parent directory to path to import radia
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'dist'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'build', 'lib', 'Release'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "build", "Release"))

import radia as rad

# Clear all objects
rad.UtiDelAll()

# Arc parameters (defined but not used in this case)
rmin = 100
rmax = 150
phimin = 0
phimax = 2 * math.pi
h = 20
nseg = 20
j = 10

# Define vertices of a cube
p1 = [-1, 1, -1]
p2 = [-1, -1, -1]
p3 = [1, -1, -1]
p4 = [1, 1, -1]
p5 = [-1, 1, 1]
p6 = [-1, -1, 1]
p7 = [1, -1, 1]
p8 = [1, 1, 1]

# Define vertices list
vertices = [p1, p2, p3, p4, p5, p6, p7, p8]

# Define faces (1-indexed in Mathematica, 0-indexed in Python)
# Each face is defined by vertex indices
faces = [
	[1, 2, 3, 4],  # Bottom face
	[5, 6, 7, 8],  # Top face
	[1, 2, 6, 5],  # Front face
	[2, 3, 7, 6],  # Right face
	[3, 4, 8, 7],  # Back face
	[4, 1, 5, 8]   # Left face
]

# Create polyhedron with magnetization [0, 0, 0]
g1 = rad.ObjPolyhdr(vertices, faces, [0, 0, 0])

# Print object ID
print(f"Polyhedron object ID: {g1}")

# Note: 3D visualization requires additional libraries
# For now, we skip the Graphics3D export

# Note: g2 is not defined in the original script, so the field calculation would fail
# Commenting out the field calculation as it references undefined g2
# field = rad.Fld(g2, 'b', [0, 0, 0])
# print(f"Magnetic field at origin: Bx={field[0]:.6e}, By={field[1]:.6e}, Bz={field[2]:.6e} T")

print("Cube polyhedron created successfully")
print(f"Vertices: {len(vertices)}")
print(f"Faces: {len(faces)}")
print(f"Note: Field calculation requires defining g2 object")
