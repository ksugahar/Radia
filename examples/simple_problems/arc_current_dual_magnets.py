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

# Create two hexahedral magnets with magnetization
# Note: Radia magnetization unit is Tesla (T), not A/m
# For permanent magnets, set magnetization directly (no material needed)
# 300x300x5 mm centered at [0, 0, -50]
vertices1 = [[-150, -150, -52.5], [150, -150, -52.5], [150, 150, -52.5], [-150, 150, -52.5],
             [-150, -150, -47.5], [150, -150, -47.5], [150, 150, -47.5], [-150, 150, -47.5]]
g2 = rad.ObjHexahedron(vertices1, [0, 0, 1.0])
# 200x200x5 mm centered at [0, 0, 50]
vertices2 = [[-100, -100, 47.5], [100, -100, 47.5], [100, 100, 47.5], [-100, 100, 47.5],
             [-100, -100, 52.5], [100, -100, 52.5], [100, 100, 52.5], [-100, 100, 52.5]]
g3 = rad.ObjHexahedron(vertices2, [0, 0, 0.8])

# Combine magnets into a container
g2 = rad.ObjCnt([g2, g3])

# Note: Material properties (MatLin, MatSatIso) are for soft magnetic materials
# like iron yokes, NOT for permanent magnets with fixed magnetization

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

# VTS Export - Export field distribution with same filename as script
try:
	# Get script basename without extension
	script_name = os.path.splitext(os.path.basename(__file__))[0]
	vts_filename = f"{script_name}.vts"
	vts_path = os.path.join(os.path.dirname(__file__), vts_filename)

	# Geometry: arc current rmin=100, rmax=150, two magnets at z=-50 and z=50
	# Extend range to cover geometry
	x_range = [-200, 200]
	y_range = [-200, 200]
	z_range = [-100, 100]

	rad.FldVTS(g, vts_path, x_range, y_range, z_range, 21, 21, 21, 1, 0, 1.0)
	print(f"\n[VTS] Exported: {vts_filename}")
	print(f"      View with: paraview {vts_filename}")
except Exception as e:
	print(f"\n[VTS] Warning: Export failed: {e}")
