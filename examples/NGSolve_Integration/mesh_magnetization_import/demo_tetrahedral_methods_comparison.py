#!/usr/bin/env python
"""
Demonstrate ANALYTICAL tetrahedral method with NGSolve mesh import

This example shows how to:
1. Generate a tetrahedral mesh with NGSolve/Netgen
2. Import it into Radia using the ANALYTICAL method
3. Compare with Radia's built-in hexahedral method

The ANALYTICAL method uses basis vectors extracted from face vertices
and the proven RadAnalyticalFieldFromPolygonCharge formula.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../build/Release'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/python'))

import numpy as np
import radia as rad
from netgen.occ import Box, OCCGeometry
from ngsolve import Mesh
from netgen_mesh_import import netgen_mesh_to_radia

# Enable ANALYTICAL method for tetrahedral elements
os.environ['RADIA_TETRA_METHOD'] = 'ANALYTICAL'

rad.FldUnits('m')
cube_size = 0.1
M = np.array([0, 0, 1.2])  # Tesla (permanent magnet)

print("="*70)
print("Tetrahedral Mesh Import Demo (ANALYTICAL Method)")
print("="*70)

# Reference: Built-in ObjDivMag (hexahedral)
print("\n[Reference] Built-in Radia ObjDivMag (hexahedral)")
cube_builtin = rad.ObjRecMag([-cube_size/2, -cube_size/2, -cube_size/2],
                              [cube_size, cube_size, cube_size],
                              [M[0], M[1], M[2]])
rad.ObjDivMag(cube_builtin, [5, 5, 5])

test_pt = [0.08, 0, 0]
H_builtin = rad.Fld(cube_builtin, 'h', test_pt)
H_builtin_mag = np.linalg.norm(H_builtin)
print(f"  |H| at {test_pt}: {H_builtin_mag:.6f} A/m")

# Generate tetrahedral mesh with NGSolve/Netgen
print("\n[Mesh Generation] NGSolve/Netgen tetrahedral mesh")
geo = OCCGeometry(Box((-cube_size/2, -cube_size/2, -cube_size/2),
                      (cube_size/2, cube_size/2, cube_size/2)), dim=3)
ngmesh_native = geo.GenerateMesh(maxh=cube_size/5)
ngmesh = Mesh(ngmesh_native)
print(f"  Generated {ngmesh_native.ne} tetrahedral elements")

# Import into Radia with ANALYTICAL method
print("\n[Test] ANALYTICAL Tetrahedral Method")
print("  Method: Extract basis vectors from face vertices")
print("  Formula: RadAnalyticalFieldFromPolygonCharge")

cube_analytical = netgen_mesh_to_radia(ngmesh,
                                        material={'magnetization': [M[0], M[1], M[2]]},
                                        verbose=False)
H_analytical = rad.Fld(cube_analytical, 'h', test_pt)
H_analytical_mag = np.linalg.norm(H_analytical)
error = 100 * abs(H_analytical_mag - H_builtin_mag) / H_builtin_mag

print(f"  |H| at {test_pt}: {H_analytical_mag:.6f} A/m")
print(f"  Error vs reference: {error:.1f}%")

# Summary
print("\n" + "="*70)
print("Summary")
print("="*70)
print(f"Reference (hexahedral):    {H_builtin_mag:.6f} A/m")
print(f"ANALYTICAL (tetrahedral):  {H_analytical_mag:.6f} A/m")
print(f"Error:                     {error:.1f}%")

if error < 20:
    print("\n[PASS] Error < 20% - ANALYTICAL method working correctly")
else:
    print("\n[WARNING] Error > 20% - may need investigation")

print("\n" + "="*70)
print("Implementation Details")
print("="*70)
print("ANALYTICAL method features:")
print("- Extracts orthonormal basis (AA, BB, CC) from face vertices")
print("- Projects 3D vertices to 2D local coordinates")
print("- Transforms magnetization to local coordinate system")
print("- Uses RadAnalyticalFieldFromPolygonCharge (proven formula)")
print("\nSee: ANALYTICAL_METHOD_SUCCESS.md for full implementation report")
print("="*70)
