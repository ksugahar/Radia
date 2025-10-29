#!/usr/bin/env python
"""
Test hexahedron creation with rad.ObjPolyhdr()
Testing various cases to identify potential issues
"""

import sys
import os

# Add parent directory to path to import radia
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'dist'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'build', 'lib', 'Release'))

import radia as rad

print("="*60)
print("Test 1: Simple regular cube")
print("="*60)

coords1 = [
	[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],  # bottom
	[-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]       # top
]
face = [[1,2,3,4],[1,5,6,2],[2,6,7,3],[3,7,8,4],[1,4,8,5],[5,8,7,6]]
M = [0, 0, 0.5]

try:
	hex1 = rad.ObjPolyhdr(coords1, face, M)
	print(f'Success: Created hexahedron with ID {hex1}')

	# Calculate field
	field = rad.Fld(hex1, 'b', [0, 0, 0])
	print(f'  Field at origin: Bx={field[0]:.6e}, By={field[1]:.6e}, Bz={field[2]:.6e} T')
except Exception as e:
	print(f'ERROR: {e}')

print()
print("="*60)
print("Test 2: Skewed hexahedron (angled)")
print("="*60)

# Skewed hexahedron similar to webcut result
coords2 = [
	[-5.0, -5.0, -5.0],
	[5.0, -5.0, -5.0],
	[5.0, 5.0, -5.0],
	[-5.0, 5.0, -5.0],
	[-5.0, -4.6, 5.0],  # slightly angled
	[5.0, -4.6, 5.0],
	[5.0, 5.4, 5.0],
	[-5.0, 5.4, 5.0]
]

try:
	hex2 = rad.ObjPolyhdr(coords2, face, M)
	print(f'Success: Created skewed hexahedron with ID {hex2}')

	# Calculate field
	field = rad.Fld(hex2, 'b', [0, 0, 0])
	print(f'  Field at origin: Bx={field[0]:.6e}, By={field[1]:.6e}, Bz={field[2]:.6e} T')
except Exception as e:
	print(f'ERROR: {e}')

print()
print("="*60)
print("Test 3: Highly distorted hexahedron")
print("="*60)

# Highly distorted hexahedron
coords3 = [
	[0.0, 0.0, 0.0],
	[10.0, 0.0, 0.0],
	[10.0, 10.0, 0.0],
	[0.0, 10.0, 0.0],
	[1.0, 1.0, 10.0],   # highly distorted top
	[9.0, 1.0, 10.0],
	[9.0, 9.0, 10.0],
	[1.0, 9.0, 10.0]
]

try:
	hex3 = rad.ObjPolyhdr(coords3, face, M)
	print(f'Success: Created distorted hexahedron with ID {hex3}')

	# Calculate field
	field = rad.Fld(hex3, 'b', [5, 5, 5])
	print(f'  Field at [5,5,5]: Bx={field[0]:.6e}, By={field[1]:.6e}, Bz={field[2]:.6e} T')
except Exception as e:
	print(f'ERROR: {e}')

print()
print("="*60)
print("Test 4: Multiple hexahedra in container")
print("="*60)

try:
	hex_list = []
	# Create 3 cubes
	for i in range(3):
		offset = i * 3
		coords = [
			[-1+offset, -1, -1], [1+offset, -1, -1], [1+offset, 1, -1], [-1+offset, 1, -1],
			[-1+offset, -1, 1], [1+offset, -1, 1], [1+offset, 1, 1], [-1+offset, 1, 1]
		]
		hex_list.append(rad.ObjPolyhdr(coords, face, M))

	g = rad.ObjCnt(hex_list)
	print(f'Success: Created container with {len(hex_list)} hexahedra, ID {g}')

	# Calculate field
	field = rad.Fld(g, 'b', [3, 0, 0])
	print(f'  Field at [3,0,0]: Bx={field[0]:.6e}, By={field[1]:.6e}, Bz={field[2]:.6e} T')
except Exception as e:
	print(f'ERROR: {e}')

print()
print("="*60)
print("Summary: All hexahedron tests completed")
print("="*60)
