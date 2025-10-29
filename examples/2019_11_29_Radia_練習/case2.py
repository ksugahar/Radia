#!/usr/bin/env python
"""
Case 2: Multiple Extrusion with Chamfer
Converted from Mathematica/Wolfram Language to Python
"""

import sys
import os
import math
import numpy as np

# Add parent directory to path to import radia
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'dist'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'build', 'lib', 'Release'))

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

# Geometry parameters
gap = 10
thick = 50
width = 40
chamfer = 8
current = -2000

# Dimensions
lx1 = thick / 2
ly1 = width
lz1 = 20
l1 = [lx1, ly1, lz1]

# Define cross-sections for multiple extrusion
# k1: first chamfered section
k1 = [
	[thick/4 - chamfer/2, 0, gap/2],
	[thick/2 - chamfer, ly1 - 2*chamfer]
]

# k2: second section at chamfer height
k2 = [
	[thick/4, 0, gap/2 + chamfer],
	[thick/2, ly1]
]

# k3: final section
k3 = [
	[thick/4, 0, gap/2 + lz1],
	[thick/2, ly1]
]

# Create multiple extrusion object
g1 = rad.ObjMltExtRtg([k1, k2, k3])

# Subdivide the magnet
n1 = [2, 3, 2]
rad.ObjDivMag(g1, n1)

# Print object ID
print(f"Object ID: {g1}")

# Note: 3D visualization requires additional libraries
# For now, we skip the Graphics3D export

# Note: g2 is not defined in the original script, so the field calculation would fail
# Commenting out the field calculation as it references undefined g2
# field = rad.Fld(g2, 'b', [0, 0, 0])
# print(f"Magnetic field at origin: Bx={field[0]:.6e}, By={field[1]:.6e}, Bz={field[2]:.6e} T")

print("Geometry created and subdivided successfully")
print(f"Note: Field calculation requires defining g2 object")
