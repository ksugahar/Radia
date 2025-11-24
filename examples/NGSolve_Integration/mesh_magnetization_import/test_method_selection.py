#!/usr/bin/env python
"""Verify method selection is working"""
import sys, os

# Test 1: CENTROID method
print("="*70)
print("Test 1: CENTROID Method")
print("="*70)
os.environ['RADIA_TETRA_METHOD'] = 'CENTROID'
print(f"Set environment: RADIA_TETRA_METHOD={os.environ['RADIA_TETRA_METHOD']}")

# Import radia AFTER setting environment variable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../build/Release'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/python'))
import radia as rad
import numpy as np

rad.FldUnits('m')

# Create single tetrahedron
vertices = [
    [0, 0, 0],
    [0.1, 0, 0],
    [0.05, 0.1, 0],
    [0.05, 0.05, 0.1]
]
faces = [[1, 2, 3], [1, 2, 4], [1, 3, 4], [2, 3, 4]]
M = [0, 0, 1.2]

tetra = rad.ObjPolyhdr(vertices, faces, M)
test_pt = [0.08, 0, 0]
H = rad.Fld(tetra, 'h', test_pt)
H_mag = np.linalg.norm(H)

print(f"Single tetrahedron field at {test_pt}:")
print(f"  H = {H}")
print(f"  |H| = {H_mag:.6f}")

print("\n" + "="*70)
print("IMPORTANT: Environment variable is read at C++ level")
print("If method selection happens in C++, we can't see it from Python")
print("Both methods might give similar results if:")
print("  1. Environment variable not being read correctly")
print("  2. Centroid cancellation formula needs adjustment")
print("  3. Single tetrahedron doesn't show the difference")
print("="*70)
