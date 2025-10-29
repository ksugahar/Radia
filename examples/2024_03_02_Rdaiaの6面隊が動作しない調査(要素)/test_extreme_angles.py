#!/usr/bin/env python
"""
Test hexahedra with extreme angles and distortions
Simulating webcut at various angles from 1 to 30 degrees
"""

import sys
import os
import math

# Add parent directory to path to import radia
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'dist'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'build', 'lib', 'Release'))

import radia as rad

print("="*60)
print("Test: Hexahedra with various cut angles (1-30 degrees)")
print("="*60)

# Standard hexahedron face connectivity
face = [[1,2,3,4],[1,5,6,2],[2,6,7,3],[3,7,8,4],[1,4,8,5],[5,8,7,6]]
M = [0, 0, 0.5]

# Test various angles
test_angles = [1, 5, 10, 15, 20, 25, 30]
successful = []
failed = []

for angle_deg in test_angles:
	angle_rad = math.radians(angle_deg)

	# Simulate a webcut at y=0.4 with rotation about x-axis
	# Bottom face at y=-0.5 (unchanged)
	# Top face at y=0.4 rotated by angle about x-axis through center (0, 0, 0)

	# For rotation about x-axis:
	# y' = y*cos(angle) - z*sin(angle)
	# z' = y*sin(angle) + z*cos(angle)

	y_cut = 0.4
	y_offset = y_cut - 0  # offset from rotation center

	coords = [
		[-0.5, -0.5, -0.5],  # Node 1: bottom face
		[-0.5, -0.5, 0.5],   # Node 2: bottom face
		[0.5, -0.5, 0.5],    # Node 3: bottom face
		[0.5, -0.5, -0.5],   # Node 4: bottom face
		# Top face - rotated
		[-0.5, y_offset * math.cos(angle_rad) - (-0.5) * math.sin(angle_rad),
		       y_offset * math.sin(angle_rad) + (-0.5) * math.cos(angle_rad)],  # Node 5
		[-0.5, y_offset * math.cos(angle_rad) - 0.5 * math.sin(angle_rad),
		       y_offset * math.sin(angle_rad) + 0.5 * math.cos(angle_rad)],     # Node 6
		[0.5, y_offset * math.cos(angle_rad) - 0.5 * math.sin(angle_rad),
		      y_offset * math.sin(angle_rad) + 0.5 * math.cos(angle_rad)],      # Node 7
		[0.5, y_offset * math.cos(angle_rad) - (-0.5) * math.sin(angle_rad),
		      y_offset * math.sin(angle_rad) + (-0.5) * math.cos(angle_rad)]    # Node 8
	]

	try:
		hex_obj = rad.ObjPolyhdr(coords, face, M)
		field = rad.Fld(hex_obj, 'b', [0, 0, 0])
		successful.append(angle_deg)
		print(f"Angle {angle_deg:2d} deg: SUCCESS - Hex ID {hex_obj}, Bz={field[2]:.6e} T")
	except Exception as e:
		failed.append(angle_deg)
		print(f"Angle {angle_deg:2d} deg: FAILED - {str(e)[:50]}")

print()
print("="*60)
print("Summary:")
print(f"  Successful: {len(successful)}/{len(test_angles)} angles")
print(f"  Failed: {len(failed)}/{len(test_angles)} angles")

if len(failed) > 0:
	print(f"  Failed angles: {failed}")
	print()
	print("RESULT: Some hexahedra failed with extreme distortion")
else:
	print()
	print("RESULT: All hexahedra created successfully!")
print("="*60)
