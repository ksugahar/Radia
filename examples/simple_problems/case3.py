#!/usr/bin/env python
"""
Case 3: Polyhedron (Cube) Magnet with Field Calculation
Converted from Mathematica/Wolfram Language to Python

This example demonstrates:
- Creating a cubic magnet using ObjPolyhdr
- Applying magnetization to a polyhedron
- Calculating magnetic field at various points
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

# Define vertices list
vertices = [p1, p2, p3, p4, p5, p6, p7, p8]

# Define faces using 1-based indexing (Radia convention)
# Each face is defined by vertex indices
faces = [
	[1, 4, 3, 2],  # Bottom face (-Z)
	[5, 6, 7, 8],  # Top face (+Z)
	[1, 2, 6, 5],  # Front face (-Y)
	[2, 3, 7, 6],  # Right face (+X)
	[3, 4, 8, 7],  # Back face (+Y)
	[4, 1, 5, 8]   # Left face (-X)
]

# Create polyhedron with magnetization [0, 0, 1.2] T (NdFeB typical value)
# Magnetization in Z direction
magnetization = [0, 0, 1.2]  # Tesla
g1 = rad.ObjPolyhdr(vertices, faces, magnetization)

print(f"\nCube magnet created:")
print(f"  Object ID: {g1}")
print(f"  Size: {2*size} x {2*size} x {2*size} mm")
print(f"  Magnetization: {magnetization} T")
print(f"  Vertices: {len(vertices)}")
print(f"  Faces: {len(faces)}")

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
