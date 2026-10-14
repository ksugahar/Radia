"""Basic NGSolve integration example with Radia.

This example demonstrates how to use rad.RadiaField
to create an NGSolve CoefficientFunction from a Radia magnet.

Requirements:
    - NGSolve installed
    - radia v2.5.0+ (RadiaField integrated into main module)
"""

import sys
import os

# Add paths for local development
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../build/Release'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

import numpy as np
import radia as rad


# Radia always uses meters (SI), compatible with NGSolve. No unit setup needed.
rad.UtiDelAll()

print("="*70)
print("NGSolve Integration Demo: Basic Field Evaluation")
print("="*70)

# Create a simple rectangular magnet using ObjHexahedron
# Size: 40mm x 40mm x 60mm (0.04 x 0.04 x 0.06 m), magnetization 1.2 T in z-direction
cx, cy, cz = 0, 0, 0
dx, dy, dz = 0.02, 0.02, 0.03  # Half-dimensions
vertices = [
    [cx - dx, cy - dy, cz - dz],  # vertex 1
    [cx + dx, cy - dy, cz - dz],  # vertex 2
    [cx + dx, cy + dy, cz - dz],  # vertex 3
    [cx - dx, cy + dy, cz - dz],  # vertex 4
    [cx - dx, cy - dy, cz + dz],  # vertex 5
    [cx + dx, cy - dy, cz + dz],  # vertex 6
    [cx + dx, cy + dy, cz + dz],  # vertex 7
    [cx - dx, cy + dy, cz + dz],  # vertex 8
]

# Radia magnetization is in A/m (NOT Tesla): M = Br / mu_0.
# Br = 1.2 T  ->  M = 1.2 / (4*pi*1e-7) = 954930 A/m
MU_0 = 4 * np.pi * 1e-7
Mr = 1.2 / MU_0  # 954930 A/m
magnet = rad.ObjHexahedron(vertices, [0, 0, Mr])

# A permanent magnet has FIXED magnetization -- no rad.Solve() is needed.
# (rad.Solve builds the demagnetization interaction matrix for soft-magnetic
#  materials; a bare PM has no such material, so Solve would raise
#  "Failed to create Interaction Matrix".)

print("\nRadia magnet created:")
print(f"  Size: 40mm x 40mm x 60mm")
print(f"  Magnetization: 1.2 T in z-direction")

# Test field at a point using rad.Fld directly
test_point = [0.05, 0, 0]  # 50mm away from center in x
B_direct = rad.Fld(magnet, 'b', test_point)
print(f"\nField at {test_point} (rad.Fld):")
print(f"  B = [{B_direct[0]:.6f}, {B_direct[1]:.6f}, {B_direct[2]:.6f}] T")

# Use rad.RadiaField (integrated since v2.5.0)
print("\nCreating RadiaField CoefficientFunction...")

# Create CoefficientFunction
B_cf = rad.RadiaField(magnet, 'b', units='m')
print(f"  RadiaField created with units='m'")

# Import NGSolve for mesh creation
try:
    from ngsolve import *
    from ngsolve import TaskManager
    from netgen.csg import unit_cube

    # Create a simple mesh for testing
    with TaskManager():
        mesh = Mesh(unit_cube.GenerateMesh(maxh=0.3))

        # Evaluate field on mesh
        print("\nEvaluating field on NGSolve mesh...")

        # Integrate field magnitude over mesh
        B_mag = sqrt(B_cf[0]**2 + B_cf[1]**2 + B_cf[2]**2)
        integral = Integrate(B_mag, mesh)
        print(f"  Integral of |B| over unit cube: {integral:.6f} T*m^3")

        # Evaluate at mesh center
        mip = mesh(0.5, 0.5, 0.5)
        B_at_center = B_cf(mip)
        print(f"  B at mesh center (0.5, 0.5, 0.5): {B_at_center}")

        print("\nNGSolve integration test PASSED")

except ImportError as e:
    print(f"\nNGSolve not available: {e}")
    print("Skipping mesh-based tests")

print("\n" + "="*70)
print("Demo complete")
print("="*70)

rad.UtiDelAll()
