#!/usr/bin/env python
"""
Case 3: Hexahedron (Cube) Magnet with Field Calculation
Converted from Mathematica/Wolfram Language to Python

This example demonstrates:
- Creating a cubic magnet using ObjHexahedron
- Applying magnetization to a hexahedron
- Calculating magnetic field at various points
"""

import sys
import os
import math
import numpy as np

# Add parent directory to path to import radia
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'src', 'radia'))

import radia as rad

# Clear all objects
rad.UtiDelAll()

print("=" * 70)
print("Case 3: Cubic Magnet using Polyhedron")
print("=" * 70)

# Define vertices of a cube (20mm x 20mm x 20mm centered at origin)
# Coordinates in mm
size = 10  # Half-size: 10mm -> 20mm cube
p1 = [-size, -size, -size]  # Bottom front-left
p2 = [size, -size, -size]   # Bottom front-right
p3 = [size, size, -size]    # Bottom back-right
p4 = [-size, size, -size]   # Bottom back-left
p5 = [-size, -size, size]   # Top front-left
p6 = [size, -size, size]    # Top front-right
p7 = [size, size, size]     # Top back-right
p8 = [-size, size, size]    # Top back-left

# Define vertices list (8 vertices for hexahedron)
vertices = [p1, p2, p3, p4, p5, p6, p7, p8]

# Create hexahedron with magnetization [0, 0, 1.2] T (NdFeB typical value)
# ObjHexahedron automatically generates the correct face topology
# Magnetization in Z direction
magnetization = [0, 0, 1.2]  # Tesla
g1 = rad.ObjHexahedron(vertices, magnetization)

print(f"\nCube magnet created:")
print(f"  Object ID: {g1}")
print(f"  Size: {2*size} x {2*size} x {2*size} mm")
print(f"  Magnetization: {magnetization} T")
print(f"  Vertices: {len(vertices)}")

# Set drawing attributes (blue color)
rad.ObjDrwAtr(g1, [0, 0, 1], 0.001)

# Calculate magnetic field at various points
print("\n" + "=" * 70)
print("Magnetic Field Calculation")
print("=" * 70)

test_points = [
	[0, 0, 0],      # Center of cube
	[0, 0, 20],     # 20mm above cube
	[0, 0, -20],    # 20mm below cube
	[20, 0, 0],     # 20mm to the right
	[0, 20, 0],     # 20mm to the back
]

print(f"\n{'Point (mm)':<20} {'Bx (mT)':<12} {'By (mT)':<12} {'Bz (mT)':<12} {'|B| (mT)':<12}")
print("-" * 70)

for point in test_points:
	field = rad.Fld(g1, 'b', point)
	Bx_mT = field[0] * 1000
	By_mT = field[1] * 1000
	Bz_mT = field[2] * 1000
	B_mag = math.sqrt(Bx_mT**2 + By_mT**2 + Bz_mT**2)

	point_str = f"({point[0]:5.1f}, {point[1]:5.1f}, {point[2]:5.1f})"
	print(f"{point_str:<20} {Bx_mT:<12.3f} {By_mT:<12.3f} {Bz_mT:<12.3f} {B_mag:<12.3f}")

# Additional test: Verify symmetry
print("\n" + "=" * 70)
print("Symmetry Verification (Bz component)")
print("=" * 70)

symmetric_points = [
	([0, 0, 15], "Above center"),
	([0, 0, -15], "Below center"),
	([10, 0, 15], "Above right"),
	([-10, 0, 15], "Above left"),
]

print(f"\n{'Location':<20} {'Point (mm)':<20} {'Bz (mT)':<12}")
print("-" * 55)

for point, desc in symmetric_points:
	field = rad.Fld(g1, 'b', point)
	Bz_mT = field[2] * 1000
	point_str = f"({point[0]:5.1f}, {point[1]:5.1f}, {point[2]:5.1f})"
	print(f"{desc:<20} {point_str:<20} {Bz_mT:<12.3f}")

print("\n" + "=" * 70)
print("Calculation complete.")
print("=" * 70)

# VTS Export - Export field distribution with same filename as script
try:
	# Get script basename without extension
	script_name = os.path.splitext(os.path.basename(__file__))[0]
	vts_filename = f"{script_name}.vts"
	vts_path = os.path.join(os.path.dirname(__file__), vts_filename)

	# Cube is 20mm centered at origin, extend range to 40mm for far-field
	x_range = [-40, 40]
	y_range = [-40, 40]
	z_range = [-40, 40]

	rad.FldVTS(g1, vts_path, x_range, y_range, z_range, 21, 21, 21, 1, 0, 1.0)
	print(f"\n[VTS] Exported: {vts_filename}")
	print(f"      View with: paraview {vts_filename}")
except Exception as e:
	print(f"\n[VTS] Warning: Export failed: {e}")
