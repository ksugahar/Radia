#!/usr/bin/env python
"""
Case 2: Multiple Extrusion with Chamfer - Pole Piece Magnet
Converted from Mathematica/Wolfram Language to Python

This example demonstrates:
- Creating complex 3D geometry using ObjMltExtRtg (multiple extrusion)
- Applying chamfer to edges
- Subdividing geometry for accurate field calculation
- Calculating magnetic field from extruded magnet
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
print("Case 2: Chamfered Pole Piece using Multiple Extrusion")
print("=" * 70)

# Geometry parameters (all dimensions in mm)
gap = 10           # Gap between pole pieces
thick = 50         # Total thickness
width = 40         # Width
chamfer = 8        # Chamfer size
lz1 = 20           # Height of pole piece

# Magnetization (T)
magnetization = [0, 0, 1.0]  # 1.0 T in Z direction

print(f"\nGeometry parameters:")
print(f"  Gap: {gap} mm")
print(f"  Thickness: {thick} mm")
print(f"  Width: {width} mm")
print(f"  Chamfer: {chamfer} mm")
print(f"  Height: {lz1} mm")
print(f"  Magnetization: {magnetization} T")

# Define cross-sections for multiple extrusion
# This creates a pole piece with chamfered edges
# Each section is defined as [[x, y, z], [dx, dy]]

# k1: first chamfered section (narrower due to chamfer)
k1 = [
	[thick/4 - chamfer/2, 0, gap/2],
	[thick/2 - chamfer, width - 2*chamfer]
]

# k2: second section at chamfer height (full width)
k2 = [
	[thick/4, 0, gap/2 + chamfer],
	[thick/2, width]
]

# k3: final section (same width as k2)
k3 = [
	[thick/4, 0, gap/2 + lz1],
	[thick/2, width]
]

print(f"\nCross-section levels:")
print(f"  Level 1 (z={gap/2} mm): size=[{thick/2 - chamfer}, {width - 2*chamfer}] mm (chamfered)")
print(f"  Level 2 (z={gap/2 + chamfer} mm): size=[{thick/2}, {width}] mm (full)")
print(f"  Level 3 (z={gap/2 + lz1} mm): size=[{thick/2}, {width}] mm (full)")

# Create multiple extrusion object with magnetization
g1 = rad.ObjMltExtRtg([k1, k2, k3], magnetization)

print(f"\nMagnet object created: ID = {g1}")

# Subdivide the magnet for accurate field calculation
# [nx, ny, nz] - number of subdivisions in each direction
n1 = [2, 3, 2]
rad.ObjDivMag(g1, n1)

print(f"Subdivisions applied: {n1[0]} × {n1[1]} × {n1[2]} = {n1[0]*n1[1]*n1[2]} segments")

# Set drawing attributes (blue color)
rad.ObjDrwAtr(g1, [0, 0, 1], 0.001)

# Calculate magnetic field at various points
print("\n" + "=" * 70)
print("Magnetic Field Calculation")
print("=" * 70)

test_points = [
	[0, 0, gap/2 + lz1/2],  # Center of pole piece
	[0, 0, gap/2 + lz1 + 10],  # 10mm above pole piece
	[thick/2 + 5, 0, gap/2 + lz1/2],  # 5mm to the side
	[0, width/2 + 5, gap/2 + lz1/2],  # 5mm to front
]

print(f"\n{'Point (mm)':<30} {'Bx (mT)':<12} {'By (mT)':<12} {'Bz (mT)':<12} {'|B| (mT)':<12}")
print("-" * 80)

for point in test_points:
	field = rad.Fld(g1, 'b', point)
	Bx_mT = field[0] * 1000
	By_mT = field[1] * 1000
	Bz_mT = field[2] * 1000
	B_mag = math.sqrt(Bx_mT**2 + By_mT**2 + Bz_mT**2)

	point_str = f"({point[0]:6.1f}, {point[1]:6.1f}, {point[2]:6.1f})"
	print(f"{point_str:<30} {Bx_mT:<12.3f} {By_mT:<12.3f} {Bz_mT:<12.3f} {B_mag:<12.3f}")

print("\n" + "=" * 70)
print("Note: This example demonstrates multiple extrusion (ObjMltExtRtg)")
print("      The chamfer creates a smooth transition from narrow to full width")
print("=" * 70)
print("Calculation complete.")
print("=" * 70)
