#!/usr/bin/env python
"""
Debug tetrahedral permanent magnet implementation

Simplest possible test: Single tetrahedron with magnetization
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../build/Release'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/python'))

import numpy as np
import radia as rad

rad.FldUnits('m')

print("="*70)
print("Debug: Single Tetrahedron with Permanent Magnetization")
print("="*70)

# Define a simple tetrahedron manually - standard test tetrahedron
# Unit tetrahedron with vertices at standard positions
vertices = [
    [0, 0, 0],          # vertex 0
    [1, 0, 0],          # vertex 1
    [0, 1, 0],          # vertex 2
    [0, 0, 1]           # vertex 3
]

# Faces (4 triangular faces) - CCW ordering when viewed from outside
TETRA_FACES = [[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]]

# Magnetization
M = [0, 0, 1.2]  # T

print(f"\nCreating tetrahedron:")
print(f"  Vertices: {vertices}")
print(f"  Faces: {TETRA_FACES}")
print(f"  Magnetization: M = {M} T")

try:
    tetra = rad.ObjPolyhdr(vertices, TETRA_FACES, M)
    print(f"  [OK] Created polyhedron ID: {tetra}")
except Exception as e:
    print(f"  [FAIL] Error creating polyhedron: {e}")
    sys.exit(1)

# Test field evaluation at a nearby point
test_pt = [0.2, 0, 0]  # Outside tetrahedron

print(f"\nEvaluating field at test point {test_pt}:")

try:
    H = rad.Fld(tetra, 'h', test_pt)
    print(f"  H = {H}")

    B = rad.Fld(tetra, 'b', test_pt)
    print(f"  B = {B}")

    M_eval = rad.Fld(tetra, 'm', test_pt)
    print(f"  M = {M_eval}")

    # Check if results are nan
    if np.any(np.isnan(H)):
        print(f"\n  [ERROR] H field contains NaN!")
    else:
        H_mag = np.linalg.norm(H)
        print(f"  |H| = {H_mag:.6f} A/m")

except Exception as e:
    print(f"  [FAIL] Error evaluating field: {e}")
    import traceback
    traceback.print_exc()

print("\n" + "="*70)
print("Debug complete")
print("="*70)
