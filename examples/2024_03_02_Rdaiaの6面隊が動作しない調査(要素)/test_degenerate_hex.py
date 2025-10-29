#!/usr/bin/env python
"""
Test degenerate and problematic hexahedra
This tests edge cases that might fail:
1. Very flat hexahedra (near-zero volume)
2. Inverted hexahedra (inside-out)
3. Non-convex hexahedra
4. Zero-volume hexahedra
"""

import sys
import os

# Add parent directory to path to import radia
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'dist'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'build', 'lib', 'Release'))

import radia as rad

face = [[1,2,3,4],[1,5,6,2],[2,6,7,3],[3,7,8,4],[1,4,8,5],[5,8,7,6]]
M = [0, 0, 0.5]

print("="*60)
print("Test 1: Very flat hexahedron (thickness = 0.001)")
print("="*60)

coords1 = [
	[-1, -1, -0.0005], [1, -1, -0.0005], [1, 1, -0.0005], [-1, 1, -0.0005],
	[-1, -1, 0.0005], [1, -1, 0.0005], [1, 1, 0.0005], [-1, 1, 0.0005]
]

try:
	hex1 = rad.ObjPolyhdr(coords1, face, M)
	field = rad.Fld(hex1, 'b', [0, 0, 0])
	print(f"SUCCESS: Hex ID {hex1}, Bz={field[2]:.6e} T")
except Exception as e:
	print(f"FAILED: {e}")

print()
print("="*60)
print("Test 2: Extremely flat hexahedron (thickness = 0.00001)")
print("="*60)

coords2 = [
	[-1, -1, -0.000005], [1, -1, -0.000005], [1, 1, -0.000005], [-1, 1, -0.000005],
	[-1, -1, 0.000005], [1, -1, 0.000005], [1, 1, 0.000005], [-1, 1, 0.000005]
]

try:
	hex2 = rad.ObjPolyhdr(coords2, face, M)
	field = rad.Fld(hex2, 'b', [0, 0, 0])
	print(f"SUCCESS: Hex ID {hex2}, Bz={field[2]:.6e} T")
except Exception as e:
	print(f"FAILED: {e}")

print()
print("="*60)
print("Test 3: Zero thickness hexahedron (all Z = 0)")
print("="*60)

coords3 = [
	[-1, -1, 0], [1, -1, 0], [1, 1, 0], [-1, 1, 0],
	[-1, -1, 0], [1, -1, 0], [1, 1, 0], [-1, 1, 0]
]

try:
	hex3 = rad.ObjPolyhdr(coords3, face, M)
	field = rad.Fld(hex3, 'b', [0, 0, 1])
	print(f"SUCCESS: Hex ID {hex3}, Bz={field[2]:.6e} T")
except Exception as e:
	print(f"FAILED: {e}")

print()
print("="*60)
print("Test 4: Inverted hexahedron (inside-out faces)")
print("="*60)

coords4 = [
	[-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1],   # top and bottom swapped
	[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1]
]

try:
	hex4 = rad.ObjPolyhdr(coords4, face, M)
	field = rad.Fld(hex4, 'b', [0, 0, 0])
	print(f"SUCCESS: Hex ID {hex4}, Bz={field[2]:.6e} T")
except Exception as e:
	print(f"FAILED: {e}")

print()
print("="*60)
print("Test 5: Non-convex (bowtie) hexahedron")
print("="*60)

coords5 = [
	[-1, -1, -1], [1, 1, -1], [1, -1, -1], [-1, 1, -1],  # crossed bottom
	[-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]
]

try:
	hex5 = rad.ObjPolyhdr(coords5, face, M)
	field = rad.Fld(hex5, 'b', [0, 0, 0])
	print(f"SUCCESS: Hex ID {hex5}, Bz={field[2]:.6e} T")
except Exception as e:
	print(f"FAILED: {e}")

print()
print("="*60)
print("Test 6: All nodes at same point (zero volume)")
print("="*60)

coords6 = [
	[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0],
	[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]
]

try:
	hex6 = rad.ObjPolyhdr(coords6, face, M)
	field = rad.Fld(hex6, 'b', [1, 1, 1])
	print(f"SUCCESS: Hex ID {hex6}, Bz={field[2]:.6e} T")
except Exception as e:
	print(f"FAILED: {e}")

print()
print("="*60)
print("Summary: Degenerate hexahedra tests completed")
print("="*60)
