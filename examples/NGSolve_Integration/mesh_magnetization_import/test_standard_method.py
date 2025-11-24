#!/usr/bin/env python
"""Test standard polygon method"""
import sys, os

# Set environment variable to STANDARD (or don't set it, as STANDARD is default)
os.environ['RADIA_TETRA_METHOD'] = 'STANDARD'

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../build/Release'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/python'))

import numpy as np
import radia as rad
from netgen.occ import Box, OCCGeometry
from ngsolve import Mesh
from netgen_mesh_import import netgen_mesh_to_radia

rad.FldUnits('m')
cube_size = 0.1
M = np.array([0, 0, 1.2])

print("="*70)
print("Standard Polygon Method Test")
print("="*70)
print(f"Environment: RADIA_TETRA_METHOD={os.environ.get('RADIA_TETRA_METHOD', 'not set')}")

# Reference: Built-in ObjDivMag
print("\n[Reference] Built-in Radia ObjDivMag")
cube_builtin = rad.ObjRecMag([-cube_size/2, -cube_size/2, -cube_size/2],
                              [cube_size, cube_size, cube_size],
                              [M[0], M[1], M[2]])
rad.ObjDivMag(cube_builtin, [5, 5, 5])

test_pt = [0.08, 0, 0]
H_builtin = rad.Fld(cube_builtin, 'h', test_pt)
H_builtin_mag = np.linalg.norm(H_builtin)
print(f"  |H| at {test_pt}: {H_builtin_mag:.6f}")

# Generate tetrahedral mesh
geo = OCCGeometry(Box((-cube_size/2, -cube_size/2, -cube_size/2),
                      (cube_size/2, cube_size/2, cube_size/2)), dim=3)
ngmesh_native = geo.GenerateMesh(maxh=cube_size/5)
ngmesh = Mesh(ngmesh_native)
print(f"\n[Mesh] {ngmesh_native.ne} tetrahedral elements")

# Test standard method
print("\n[Method] Standard Polygon Method")
cube_standard = netgen_mesh_to_radia(ngmesh,
                                      material={'magnetization': [M[0], M[1], M[2]]},
                                      verbose=False)
H_standard = rad.Fld(cube_standard, 'h', test_pt)
H_standard_mag = np.linalg.norm(H_standard)
ratio_standard = H_standard_mag / H_builtin_mag

print(f"  |H| at {test_pt}: {H_standard_mag:.6f}")
print(f"  Ratio vs reference: {ratio_standard:.3f}")
print(f"  Error: {100*(ratio_standard - 1.0):.1f}%")

print("\n" + "="*70)
print("Summary")
print("="*70)
print(f"Reference (ObjDivMag):     {H_builtin_mag:.6f}")
print(f"Standard method:           {H_standard_mag:.6f} ({ratio_standard:.3f}x)")
print("="*70)
