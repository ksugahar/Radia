#!/usr/bin/env python
"""
Test degenerate hexahedra one by one to identify which causes segfault
"""

import sys
import os

# Add parent directory to path to import radia
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'dist'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'build', 'lib', 'Release'))

import radia as rad

face = [[1,2,3,4],[1,5,6,2],[2,6,7,3],[3,7,8,4],[1,4,8,5],[5,8,7,6]]
M = [0, 0, 0.5]

# Select which test to run
test_num = int(sys.argv[1]) if len(sys.argv) > 1 else 1

if test_num == 1:
	print("Test 1: Very flat hexahedron (thickness = 0.001)")
	coords = [
		[-1, -1, -0.0005], [1, -1, -0.0005], [1, 1, -0.0005], [-1, 1, -0.0005],
		[-1, -1, 0.0005], [1, -1, 0.0005], [1, 1, 0.0005], [-1, 1, 0.0005]
	]

elif test_num == 2:
	print("Test 2: Extremely flat hexahedron (thickness = 0.00001)")
	coords = [
		[-1, -1, -0.000005], [1, -1, -0.000005], [1, 1, -0.000005], [-1, 1, -0.000005],
		[-1, -1, 0.000005], [1, -1, 0.000005], [1, 1, 0.000005], [-1, 1, 0.000005]
	]

elif test_num == 3:
	print("Test 3: Zero thickness hexahedron (all Z = 0)")
	coords = [
		[-1, -1, 0], [1, -1, 0], [1, 1, 0], [-1, 1, 0],
		[-1, -1, 0], [1, -1, 0], [1, 1, 0], [-1, 1, 0]
	]

elif test_num == 4:
	print("Test 4: Inverted hexahedron (inside-out faces)")
	coords = [
		[-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1],
		[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1]
	]

elif test_num == 5:
	print("Test 5: Non-convex (bowtie) hexahedron")
	coords = [
		[-1, -1, -1], [1, 1, -1], [1, -1, -1], [-1, 1, -1],
		[-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]
	]

elif test_num == 6:
	print("Test 6: All nodes at same point (zero volume)")
	coords = [
		[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0],
		[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]
	]

print(f"Attempting to create hexahedron...")
try:
	hex_obj = rad.ObjPolyhdr(coords, face, M)
	print(f"SUCCESS: Created hex ID {hex_obj}")

	field = rad.Fld(hex_obj, 'b', [0, 0, 0])
	print(f"Field calculation SUCCESS: Bz={field[2]:.6e} T")
except Exception as e:
	print(f"FAILED with exception: {e}")
