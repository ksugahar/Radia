#!/usr/bin/env python3
"""
Debug: Check correct field values

Compare direct Radia calculation with rad_ngsolve evaluation.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "build", "Release"))

import radia as rad
try:
	from ngsolve import *
	from netgen.occ import *
	import rad_ngsolve
except ImportError:
	print("ERROR: NGSolve not available")
	sys.exit(1)

import numpy as np

print("=" * 80)
print("Debug: Field Value Verification")
print("=" * 80)

# Create simple magnet
rad.UtiDelAll()
n = 5
cube_size = 100.0  # mm
elem_size = cube_size / n
mag_value = 1.2  # T

print(f"\n[Setup] Creating magnet: {n}x{n}x{n} = {n**3} elements")
elements = []
for i in range(n):
	for j in range(n):
		for k in range(n):
			x = (i - n/2 + 0.5) * elem_size
			y = (j - n/2 + 0.5) * elem_size
			z = (k - n/2 + 0.5) * elem_size
			elem = rad.ObjRecMag([x, y, z], [elem_size, elem_size, elem_size],
			                      [0, 0, mag_value])
			elements.append(elem)

magnet = rad.ObjCnt(elements)
print(f"  Magnet center: [0, 0, 0] mm")
print(f"  Magnet size: {cube_size}x{cube_size}x{cube_size} mm")
print(f"  Magnetization: [0, 0, {mag_value}] T")

# Test point
test_point_m = [0.01, 0.0, 0.0]  # meters
test_point_mm = [p * 1000 for p in test_point_m]  # millimeters

print(f"\n[Test Point]")
print(f"  Position (NGSolve): {test_point_m} m")
print(f"  Position (Radia):   {test_point_mm} mm")

# Method 1: Direct Radia calculation
print(f"\n[Method 1] Direct Radia calculation")
B_direct = rad.Fld(magnet, 'b', test_point_mm)
print(f"  B = [{B_direct[0]:.6e}, {B_direct[1]:.6e}, {B_direct[2]:.6e}] T")

# Method 2: rad_ngsolve.RadiaField evaluation at single point
print(f"\n[Method 2] rad_ngsolve.RadiaField at single point")

# Create mesh with just one element containing our test point
box = Box(Pnt(0.005, -0.01, -0.01), Pnt(0.015, 0.01, 0.01))
geo = OCCGeometry(box)
mesh = Mesh(geo.GenerateMesh(maxh=0.1))  # Very coarse

B_cf = rad_ngsolve.RadiaField(magnet, 'b')

# Evaluate at test point
try:
	mip = mesh(test_point_m[0], test_point_m[1], test_point_m[2])
	B_ngsolve = B_cf(mip)
	print(f"  B = [{B_ngsolve[0]:.6e}, {B_ngsolve[1]:.6e}, {B_ngsolve[2]:.6e}] T")
except Exception as e:
	print(f"  ERROR: {e}")
	B_ngsolve = [0, 0, 0]

# Compare
print(f"\n[Comparison]")
error = np.linalg.norm(np.array(B_direct) - np.array(B_ngsolve))
print(f"  Error: {error:.6e} T")

if error < 1e-6:
	print(f"  [OK] Values match!")
else:
	print(f"  [FAIL] Values do NOT match!")
	print(f"  Difference: [{B_direct[0]-B_ngsolve[0]:.6e}, {B_direct[1]-B_ngsolve[1]:.6e}, {B_direct[2]-B_ngsolve[2]:.6e}]")

print("\n" + "=" * 80)
