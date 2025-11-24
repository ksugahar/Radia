#!/usr/bin/env python
"""Test field evaluation at various distances from surface"""
import sys, os
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
print("Far-Field Distance Test")
print("="*70)

# Generate mesh
geo = OCCGeometry(Box((-cube_size/2, -cube_size/2, -cube_size/2),
                      (cube_size/2, cube_size/2, cube_size/2)), dim=3)
ngmesh_native = geo.GenerateMesh(maxh=cube_size/5)
ngmesh = Mesh(ngmesh_native)
mesh_h = cube_size/5

print(f"\nMesh parameters:")
print(f"  Cube size: {cube_size} m")
print(f"  Mesh h: {mesh_h} m")
print(f"  Surface (right): x = {cube_size/2} m")
print(f"  Recommended distance: {2*mesh_h} m from surface")

# Reference
cube_builtin = rad.ObjRecMag([-cube_size/2, -cube_size/2, -cube_size/2],
                              [cube_size, cube_size, cube_size],
                              [M[0], M[1], M[2]])
rad.ObjDivMag(cube_builtin, [5, 5, 5])

# Test mesh
cube_mesh = netgen_mesh_to_radia(ngmesh,
                                  material={'magnetization': [M[0], M[1], M[2]]},
                                  verbose=False)

# Test at multiple distances
distances = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.10]  # from surface
print(f"\n{'Distance':<12} {'x position':<12} {'|H| ref':<12} {'|H| mesh':<12} {'Ratio':<8} {'Status'}")
print("-"*70)

for dist in distances:
    x_pos = cube_size/2 + dist
    test_pt = [x_pos, 0, 0]

    H_ref = rad.Fld(cube_builtin, 'h', test_pt)
    H_mesh = rad.Fld(cube_mesh, 'h', test_pt)

    H_ref_mag = np.linalg.norm(H_ref)
    H_mesh_mag = np.linalg.norm(H_mesh)

    ratio = H_mesh_mag / H_ref_mag if H_ref_mag > 0 else float('inf')

    # Status indicator
    if dist >= 2*mesh_h:
        status = "OK (>= 2h)"
    elif dist >= mesh_h:
        status = "Marginal"
    else:
        status = "Too close"

    print(f"{dist:<12.4f} {x_pos:<12.4f} {H_ref_mag:<12.6f} {H_mesh_mag:<12.6f} {ratio:<8.3f} {status}")

print("\n" + "="*70)
print("Summary:")
print(f"  Mesh h = {mesh_h} m")
print(f"  Recommended min distance: {2*mesh_h} m (>= 2h)")
print(f"  Current test point (0.08): distance = {0.08 - 0.05} m = {(0.08-0.05)/mesh_h:.1f}h")
print("="*70)
