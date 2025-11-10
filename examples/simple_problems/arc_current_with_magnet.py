#!/usr/bin/env python
"""
Case 0: Arc Current with Rectangular Magnet
Converted from Mathematica/Wolfram Language to Python
"""

import sys
import os
import math
import numpy as np

# Add parent directory to path to import radia
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'dist'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'build', 'lib', 'Release'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'build', 'Release'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'src', 'python'))

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

# Create rectangular magnet with magnetization [0,0,1.0] T
# Note: Radia magnetization unit is Tesla (T), not A/m
# For permanent magnets, set magnetization directly (no material needed)
g2 = rad.ObjRecMag([0, 0, -50], [300, 300, 5], [0, 0, 1.0])

# Note: Material properties (MatLin, MatSatIso) are for soft magnetic materials
# like iron yokes, NOT for permanent magnets with fixed magnetization

# Set drawing attributes
rad.ObjDrwAtr(g1, [1, 0, 0], 0.001)  # Red for g1
rad.ObjDrwAtr(g2, [0, 0, 1], 0.001)  # Blue for g2

# Create container with both objects
g = rad.ObjCnt([g1, g2])

# Print object ID
print(f"Container object ID: {g}")

# Note: 3D visualization requires additional libraries (matplotlib with mplot3d)
# For now, we skip the Graphics3D export

# Calculate magnetic field at origin
field = rad.Fld(g2, 'b', [0, 0, 0])
print(f"Magnetic field at origin: Bx={field[0]:.6e}, By={field[1]:.6e}, Bz={field[2]:.6e} T")

print("Calculation complete.")

# VTK Export - Export geometry with same filename as script
try:
	from radia_vtk_export import exportGeometryToVTK

	# Get script basename without extension
	script_name = os.path.splitext(os.path.basename(__file__))[0]
	vtk_filename = f"{script_name}.vtk"
	vtk_path = os.path.join(os.path.dirname(__file__), vtk_filename)

	# Export geometry
	exportGeometryToVTK(g, vtk_path)
	print(f"\n[VTK] Exported: {vtk_filename}")
	print(f"      View with: paraview {vtk_filename}")
except ImportError:
	print("\n[VTK] Warning: radia_vtk_export not available (VTK export skipped)")
except Exception as e:
	print(f"\n[VTK] Warning: Export failed: {e}")
