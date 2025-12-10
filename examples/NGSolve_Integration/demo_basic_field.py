"""Basic NGSolve integration example with Radia.

This example demonstrates how to use radia_ngsolve.RadiaField
to create an NGSolve CoefficientFunction from a Radia magnet.

Requirements:
    - NGSolve installed
    - radia_ngsolve module built and available
"""

import sys
import os

# Add paths for local development
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../build/Release'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

import numpy as np
import radia as rad

# Set units to meters for NGSolve compatibility
rad.FldUnits('m')

print("="*70)
print("NGSolve Integration Demo: Basic Field Evaluation")
print("="*70)

# Create a simple rectangular magnet
# Size: 40mm x 40mm x 60mm, magnetization 1.2 T in z-direction
magnet = rad.ObjRecMag([0, 0, 0], [0.04, 0.04, 0.06], [0, 0, 1.2e6])

# Solve the magnetization problem (for soft magnetic materials)
# For permanent magnets this is optional but doesn't hurt
rad.Solve(magnet, 0.0001, 1000)

print("\nRadia magnet created:")
print(f"  Size: 40mm x 40mm x 60mm")
print(f"  Magnetization: 1.2 T in z-direction")

# Test field at a point using rad.Fld directly
test_point = [0.05, 0, 0]  # 50mm away from center in x
B_direct = rad.Fld(magnet, 'b', test_point)
print(f"\nField at {test_point} (rad.Fld):")
print(f"  B = [{B_direct[0]:.6f}, {B_direct[1]:.6f}, {B_direct[2]:.6f}] T")

# Try to import radia_ngsolve
try:
    import radia_ngsolve
    print("\nradia_ngsolve module loaded successfully")

    # Create CoefficientFunction
    B_cf = radia_ngsolve.RadiaField(magnet, 'b', units='m')
    print(f"  RadiaField created with units='m'")

    # Import NGSolve for mesh creation
    try:
        from ngsolve import *
        from netgen.csg import unit_cube

        # Create a simple mesh for testing
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

except ImportError as e:
    print(f"\nradia_ngsolve not available: {e}")
    print("Build radia_ngsolve module first using Build_NGSolve.ps1")

print("\n" + "="*70)
print("Demo complete")
print("="*70)
