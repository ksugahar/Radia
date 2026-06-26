#!/usr/bin/env python
"""
Case 1: Arc Current with Two Rectangular Magnets
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
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'src', 'radia'))

import radia as rad

# Clear all objects
rad.UtiDelAll()

# Parameters (in meters)
rmin = 0.100       # 100 mm
rmax = 0.150       # 150 mm
phimin = 0
phimax = 2 * math.pi
h = 0.020          # 20 mm
nseg = 20
j = 10e6           # 10 A/mm^2 = 10e6 A/m^2

# Create arc with current
g1 = rad.ObjArcCur([0, 0, 0], [rmin, rmax], [phimin, phimax], h, nseg, 'man', 'z', j)

# Create two hexahedral magnets with magnetization
# Note: Radia magnetization unit is Tesla (T), not A/m
# For permanent magnets, set magnetization directly (no material needed)
# 300x300x5 mm centered at [0, 0, -0.050 m]
vertices1 = [[-0.150, -0.150, -0.0525], [0.150, -0.150, -0.0525], [0.150, 0.150, -0.0525], [-0.150, 0.150, -0.0525],
             [-0.150, -0.150, -0.0475], [0.150, -0.150, -0.0475], [0.150, 0.150, -0.0475], [-0.150, 0.150, -0.0475]]
g2 = rad.ObjHexahedron(vertices1, [0, 0, 1.0])
# 200x200x5 mm centered at [0, 0, 0.050 m]
vertices2 = [[-0.100, -0.100, 0.0475], [0.100, -0.100, 0.0475], [0.100, 0.100, 0.0475], [-0.100, 0.100, 0.0475],
             [-0.100, -0.100, 0.0525], [0.100, -0.100, 0.0525], [0.100, 0.100, 0.0525], [-0.100, 0.100, 0.0525]]
g3 = rad.ObjHexahedron(vertices2, [0, 0, 0.8])

# Combine magnets into a container
g2 = rad.ObjCnt([g2, g3])

# Note: Material properties (MatLin, MatSatIso) are for soft magnetic materials
# like iron yokes, NOT for permanent magnets with fixed magnetization

# Create final container with arc and magnets
g = rad.ObjCnt([g1, g2])

# Print object ID
print(f"Container object ID: {g}")

# Note: 3D visualization requires additional libraries
# For now, we skip the Graphics3D export

# Calculate magnetic field at origin
field = rad.Fld(g2, 'b', [0, 0, 0])
print(f"Magnetic field at origin: Bx={field[0]:.6e}, By={field[1]:.6e}, Bz={field[2]:.6e} T")

print("Calculation complete.")
