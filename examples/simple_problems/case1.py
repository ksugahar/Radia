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

import radia as rad

# Clear all objects
rad.UtiDelAll()

# Parameters
rmin = 100
rmax = 150
phimin = 0
phimax = 2 * math.pi
h = 20
nseg = 20
j = 10

# Create arc with current
g1 = rad.ObjArcCur([0, 0, 0], [rmin, rmax], [phimin, phimax], h, nseg, j)

# Define linear material
m1 = rad.MatLin([0.001, 0.001], [0, 0, 0])

# Create two rectangular magnets with magnetization
# Note: Radia magnetization unit is Tesla (T), not A/m
g2 = rad.ObjRecMag([0, 0, -50], [300, 300, 5], [0, 0, 1.0])
g3 = rad.ObjRecMag([0, 0, 50], [200, 200, 5], [0, 0, 0.8])

# Combine magnets into a container
g2 = rad.ObjCnt([g2, g3])

# Apply material to magnet container
rad.MatApl(g2, m1)

# Set drawing attributes
rad.ObjDrwAtr(g1, [1, 0, 0], 0.001)  # Red for arc current
rad.ObjDrwAtr(g2, [0, 0, 1], 0.001)  # Blue for magnets

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
