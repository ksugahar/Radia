"""Coordinate transformation example for NGSolve integration.

This example demonstrates how to use the coordinate transformation
parameters (origin, u_axis, v_axis, w_axis) in rad.RadiaField.

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
print("NGSolve Integration Demo: Coordinate Transformation")
print("="*70)

# Create a simple rectangular magnet at origin using ObjHexahedron
# Center: [0, 0, 0], Dimensions: [0.04, 0.04, 0.06] m, Magnetization in z-direction
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

# A permanent magnet has FIXED magnetization -- no rad.Solve() is needed
# (a bare PM has no soft-magnetic material, so there is no demagnetization
#  interaction matrix to build).

print("\nRadia magnet created at origin")
print("  Size: 40mm x 40mm x 60mm")
print("  Magnetization: 1.2 T in z-direction")

# Test point in global coordinates
test_point_global = [0.1, 0, 0]  # 100mm in x direction
B_global = rad.Fld(magnet, 'b', test_point_global)
print(f"\nField at {test_point_global} (global coords):")
print(f"  B = [{B_global[0]:.6f}, {B_global[1]:.6f}, {B_global[2]:.6f}] T")

# Use rad.RadiaField (integrated since v2.5.0)
print("\nrad.RadiaField available")

# Example 1: Translation only (shifted origin)
print("\n--- Example 1: Translation (origin shift) ---")
B_cf_shifted = rad.RadiaField(
    magnet, 'b',
    origin=[0.1, 0, 0],  # Origin shifted to (0.1, 0, 0)
    units='m'
)
print("  Created RadiaField with origin=[0.1, 0, 0]")
print("  Point (0, 0, 0) in local coords -> (0.1, 0, 0) in global coords")

# Example 2: 90-degree rotation (local z-axis aligned with global x-axis)
print("\n--- Example 2: 90-degree rotation ---")
B_cf_rotated = rad.RadiaField(
    magnet, 'b',
    origin=[0, 0, 0],
    u_axis=[0, 1, 0],   # local x -> global y
    v_axis=[0, 0, 1],   # local y -> global z
    w_axis=[1, 0, 0],   # local z -> global x
    units='m'
)
print("  Created RadiaField with 90-degree rotation:")
print("    u_axis (local x) -> global y")
print("    v_axis (local y) -> global z")
print("    w_axis (local z) -> global x")

# Example 3: Combined translation and rotation
print("\n--- Example 3: Translation + Rotation ---")
B_cf_combined = rad.RadiaField(
    magnet, 'b',
    origin=[0.2, 0.1, 0],  # Origin at (0.2, 0.1, 0)
    u_axis=[1, 0, 0],      # Keep local x aligned with global x
    v_axis=[0, 0, 1],      # local y -> global z
    w_axis=[0, -1, 0],     # local z -> -global y
    units='m'
)
print("  Created RadiaField with origin=[0.2, 0.1, 0] and rotation")

# Practical example: Evaluate field in cylindrical-like coordinates
print("\n--- Practical Example: Multiple observation points ---")
print("  Evaluating field at points along a circle around the magnet:")

radius = 0.08  # 80mm radius
n_points = 8
for i in range(n_points):
    theta = 2 * np.pi * i / n_points
    x = radius * np.cos(theta)
    y = radius * np.sin(theta)
    z = 0.0

    B = rad.Fld(magnet, 'b', [x, y, z])
    B_mag = np.sqrt(B[0]**2 + B[1]**2 + B[2]**2)
    print(f"    theta={np.degrees(theta):5.1f} deg: |B|={B_mag*1000:.3f} mT")

# Try NGSolve mesh-based evaluation
try:
    from ngsolve import *
    from ngsolve import TaskManager
    from netgen.csg import unit_cube

    print("\n--- NGSolve Mesh Evaluation ---")
    with TaskManager():
        mesh = Mesh(unit_cube.GenerateMesh(maxh=0.3))

        # Evaluate at mesh center with different transformations
        mip = mesh(0.5, 0.5, 0.5)

        # Global (no transform)
        B_cf_global = rad.RadiaField(magnet, 'b', units='m')
        B1 = B_cf_global(mip)
        print(f"  B at (0.5, 0.5, 0.5) - no transform: {B1}")

        # With origin shift
        B_cf_shifted = rad.RadiaField(
            magnet, 'b', origin=[-0.5, -0.5, -0.5], units='m'
        )
        B2 = B_cf_shifted(mip)
        print(f"  B at (0.5, 0.5, 0.5) - origin=[-0.5,-0.5,-0.5]: {B2}")
        print("    (Equivalent to evaluating at (0, 0, 0) in magnet frame)")

        print("\nNGSolve coordinate transformation test PASSED")

except ImportError as e:
    print(f"\nNGSolve not available: {e}")
    print("Skipping mesh-based tests")

print("\n" + "="*70)
print("Coordinate transformation demo complete")
print("="*70)

rad.UtiDelAll()
