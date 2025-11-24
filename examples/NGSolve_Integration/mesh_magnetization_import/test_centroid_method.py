#!/usr/bin/env python
"""Test centroid charge cancellation method"""
import sys, os

# Set environment variable BEFORE importing radia
os.environ['RADIA_TETRA_METHOD'] = 'CENTROID'

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
print("Centroid Charge Cancellation Method Test")
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

# Test centroid method
print("\n[Method] Centroid Charge Cancellation")
cube_centroid = netgen_mesh_to_radia(ngmesh,
                                      material={'magnetization': [M[0], M[1], M[2]]},
                                      verbose=False)
H_centroid = rad.Fld(cube_centroid, 'h', test_pt)
H_centroid_mag = np.linalg.norm(H_centroid)
ratio_centroid = H_centroid_mag / H_builtin_mag

print(f"  |H| at {test_pt}: {H_centroid_mag:.6f}")
print(f"  Ratio vs reference: {ratio_centroid:.3f}")
print(f"  Error: {100*(ratio_centroid - 1.0):.1f}%")

print("\n" + "="*70)
print("Summary")
print("="*70)
print(f"Reference (ObjDivMag):     {H_builtin_mag:.6f}")
print(f"Centroid method:           {H_centroid_mag:.6f} ({ratio_centroid:.3f}x)")

if 0.8 < ratio_centroid < 1.2:
    print("\n[SUCCESS] Centroid method gives reasonable results (within 20%)")
    print("Interior face cancellation is working!")
elif ratio_centroid > 10:
    print("\n[ISSUE] Ratio too high - centroid cancellation may not be working")
else:
    print(f"\n[INFO] Ratio: {ratio_centroid:.3f}x reference")

print("="*70)
