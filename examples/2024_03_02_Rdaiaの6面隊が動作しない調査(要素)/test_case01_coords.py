#!/usr/bin/env python
"""
Test case_01 hexahedron with actual coordinates from angle_1.bdf
This tests the problematic hexahedron from webcut at 1 degree angle
"""

import sys
import os

# Add parent directory to path to import radia
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'dist'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'build', 'lib', 'Release'))

import radia as rad

print("="*60)
print("Test: Webcut hexahedron from angle_1.bdf (1 degree cut)")
print("="*60)

# Actual coordinates from angle_1.bdf
# CHEXA   1       1       7       5       6       8       1       3       4       2
# Node order: 7, 5, 6, 8, 1, 3, 4, 2
coords_cubit = [
	[-0.5, -0.5, -0.5],           # Node 7
	[-0.5, -0.5, 0.5],            # Node 5
	[0.5, -0.5, 0.5],             # Node 6
	[0.5, -0.5, -0.5],            # Node 8
	[-0.5, 0.40878846368167, -0.5],   # Node 1
	[-0.5, 0.39133339875345, 0.5],    # Node 3
	[0.5, 0.39133339875345, 0.5],     # Node 4
	[0.5, 0.40878846368167, -0.5]     # Node 2
]

# Standard hexahedron face connectivity
face = [[1,2,3,4],[1,5,6,2],[2,6,7,3],[3,7,8,4],[1,4,8,5],[5,8,7,6]]

# Magnetization
M = [0, 0, 0.5]

print(f"Coordinates (8 nodes):")
for i, coord in enumerate(coords_cubit):
	print(f"  Node {i+1}: [{coord[0]:8.5f}, {coord[1]:8.5f}, {coord[2]:8.5f}]")

print()
print("Attempting to create hexahedron with rad.ObjPolyhdr()...")

try:
	hex1 = rad.ObjPolyhdr(coords_cubit, face, M)
	print(f"SUCCESS: Created hexahedron with ID {hex1}")
	print()

	# Calculate field at various points
	test_points = [
		[0, 0, 0],
		[0, -0.25, 0],
		[0, 0.2, 0]
	]

	print("Field calculations:")
	for pt in test_points:
		field = rad.Fld(hex1, 'b', pt)
		print(f"  Point [{pt[0]:5.2f}, {pt[1]:5.2f}, {pt[2]:5.2f}]: Bx={field[0]:.6e}, By={field[1]:.6e}, Bz={field[2]:.6e} T")

	print()
	print("="*60)
	print("RESULT: No issues detected - hexahedron works correctly!")
	print("="*60)

except Exception as e:
	print(f"ERROR: Failed to create hexahedron")
	print(f"Exception: {e}")
	print()
	print("="*60)
	print("RESULT: Issue detected with this hexahedron geometry")
	print("="*60)
