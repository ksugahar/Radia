#!/usr/bin/env python
"""
Demo: Vector Potential A from a Permanent Magnet using RadiaField

This example demonstrates how to:
1. Create a permanent magnet using ObjHexahedron
2. Compute vector potential A at arbitrary points using rad.Fld()
3. Use RadiaField CoefficientFunction to project A onto NGSolve HCurl space
4. Visualize the vector potential field

Author: Radia Development Team
Date: 2025-12-31
"""
import sys
import os

# Path setup
_script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_script_dir, '..', '..', 'src', 'radia'))

import numpy as np
import radia as rad

# Check NGSolve availability
try:
    from ngsolve import *
    from netgen.occ import Box, Pnt, OCCGeometry
    import radia_ngsolve
    NGSOLVE_AVAILABLE = True
except ImportError as e:
    print('NGSolve not available: %s' % e)
    print('Running Radia-only demo...')
    NGSOLVE_AVAILABLE = False

print('=' * 70)
print('Demo: Vector Potential A from Permanent Magnet')
print('=' * 70)

# =============================================================================
# Step 1: Create permanent magnet
# =============================================================================
print()
print('[Step 1] Creating permanent magnet')
print('-' * 70)

rad.UtiDelAll()
rad.FldUnits('m')  # Use meters (required for NGSolve integration)

# Define hexahedral magnet vertices
# Dimensions: 40mm x 40mm x 60mm centered at origin
dx, dy, dz = 0.02, 0.02, 0.03  # Half-dimensions in meters
vertices = [
    [-dx, -dy, -dz],  # vertex 1
    [ dx, -dy, -dz],  # vertex 2
    [ dx,  dy, -dz],  # vertex 3
    [-dx,  dy, -dz],  # vertex 4
    [-dx, -dy,  dz],  # vertex 5
    [ dx, -dy,  dz],  # vertex 6
    [ dx,  dy,  dz],  # vertex 7
    [-dx,  dy,  dz],  # vertex 8
]

# Magnetization: Br = 1.2 T -> Mr = Br/mu_0 = 954930 A/m
MU_0 = 4 * np.pi * 1e-7
Br = 1.2  # Tesla
Mr = Br / MU_0  # A/m

# Create magnet with z-direction magnetization
magnet = rad.ObjHexahedron(vertices, [0, 0, Mr])

print('  Magnet created: ID = %d' % magnet)
print('  Dimensions: 40mm x 40mm x 60mm')
print('  Magnetization: Mz = %.0f A/m (Br = %.1f T)' % (Mr, Br))

# =============================================================================
# Step 2: Compute A field at sample points using rad.Fld()
# =============================================================================
print()
print('[Step 2] Computing vector potential A at sample points')
print('-' * 70)

# Sample points around the magnet
sample_points = [
    [0.05, 0.0, 0.0],   # +x side
    [0.0, 0.05, 0.0],   # +y side
    [0.0, 0.0, 0.05],   # +z side
    [0.03, 0.03, 0.03], # diagonal
]

print()
print('  %-30s  %-40s' % ('Point (m)', 'A (T*m)'))
print('  ' + '-' * 70)

for pt in sample_points:
    A = rad.Fld(magnet, 'a', pt)
    B = rad.Fld(magnet, 'b', pt)
    print('  [%6.3f, %6.3f, %6.3f]  A=[%10.3e, %10.3e, %10.3e]' % (
        pt[0], pt[1], pt[2], A[0], A[1], A[2]))

# =============================================================================
# Step 3: Use RadiaField with NGSolve (if available)
# =============================================================================
if NGSOLVE_AVAILABLE:
    print()
    print('[Step 3] Creating RadiaField CoefficientFunction for A')
    print('-' * 70)

    # Create vector potential CoefficientFunction
    A_cf = radia_ngsolve.RadiaField(magnet, 'a')
    B_cf = radia_ngsolve.RadiaField(magnet, 'b')

    print('  A_cf = RadiaField(magnet, "a")  # Vector potential')
    print('  B_cf = RadiaField(magnet, "b")  # Magnetic field')

    # =============================================================================
    # Step 4: Create mesh and project A onto HCurl space
    # =============================================================================
    print()
    print('[Step 4] Projecting A onto HCurl finite element space')
    print('-' * 70)

    # Create mesh around the magnet (air region)
    box = Box(Pnt(-0.08, -0.08, -0.08), Pnt(0.08, 0.08, 0.08))
    geo = OCCGeometry(box)
    mesh = Mesh(geo.GenerateMesh(maxh=0.02))

    print('  Mesh: %d elements, %d vertices' % (mesh.ne, mesh.nv))

    # HCurl space for vector potential
    fes_hcurl = HCurl(mesh, order=2)
    gf_A = GridFunction(fes_hcurl)
    gf_A.Set(A_cf)

    print('  HCurl space: %d DOFs' % fes_hcurl.ndof)
    print('  A projected onto GridFunction')

    # Compute curl(A)
    curl_A = curl(gf_A)

    # =============================================================================
    # Step 5: Verify curl(A) ~ B
    # =============================================================================
    print()
    print('[Step 5] Verifying curl(A) ~ B')
    print('-' * 70)

    # Test at a point
    test_pt = [0.05, 0.04, 0.03]
    mip = mesh(test_pt[0], test_pt[1], test_pt[2])

    curl_A_val = [curl_A[i](mip) for i in range(3)]
    curl_A_mag = np.sqrt(sum(v**2 for v in curl_A_val))

    B_direct = rad.Fld(magnet, 'b', test_pt)
    B_mag = np.sqrt(sum(v**2 for v in B_direct))

    print('  Test point: %s m' % test_pt)
    print('  |curl(A)| = %.6e' % curl_A_mag)
    print('  |B|       = %.6e' % B_mag)
    print('  Ratio     = %.4f' % (curl_A_mag / B_mag if B_mag > 0 else 0))

    # =============================================================================
    # Step 6: Export VTK
    # =============================================================================
    print()
    print('[Step 6] Exporting VTK visualization')
    print('-' * 70)

    os.chdir(_script_dir)

    try:
        vtk = VTKOutput(
            mesh,
            coefs=[gf_A, curl_A],
            names=['A_vector_potential', 'curl_A'],
            filename='demo_vector_potential',
            subdivision=2
        )
        vtk.Do()
        print('  [OK] demo_vector_potential.vtu exported')
    except Exception as e:
        print('  [ERROR] VTK export failed: %s' % e)

# =============================================================================
# Summary
# =============================================================================
print()
print('=' * 70)
print('Summary')
print('=' * 70)
print()
print('Vector potential A is available for permanent magnets via:')
print()
print('  1. Direct computation:')
print('     A = rad.Fld(magnet, "a", [x, y, z])')
print()
print('  2. NGSolve CoefficientFunction:')
print('     A_cf = radia_ngsolve.RadiaField(magnet, "a")')
print('     gf_A = GridFunction(HCurl(mesh, order=2))')
print('     gf_A.Set(A_cf)')
print()
print('  3. Curl computation in NGSolve:')
print('     curl_A = curl(gf_A)  # Should equal B')
print()
print('Notes:')
print('  - A field uses FACE INTEGRATION (not dipole approximation)')
print('  - curl(A)/B ratio in Radia internal units is not 1.0')
print('  - In SI units: A_SI = mu_0 * A_Radia, then curl(A_SI) = B')
print('  - Use FldUnits("m") for NGSolve compatibility')
print()
print('=' * 70)

rad.UtiDelAll()
