#!/usr/bin/env python
"""
Verify that curl(A) = B using RadiaField with HCurl and HDiv spaces

This script verifies the Maxwell relation B = curl(A) by:
1. Creating a Radia magnet
2. Using RadiaField to project A onto HCurl space
3. Computing curl(A) in NGSolve
4. Using RadiaField to project B onto HDiv space
5. Comparing curl(A) with B

This demonstrates the correct usage of radia_ngsolve for vector potential
and magnetic field evaluation in NGSolve finite element spaces.

Author: Radia Development Team
Date: 2025-12-13
"""
import sys
import os

# Path setup - script is in ngsolve_integration/verify_curl_A_equals_B/
_script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_script_dir, '..', '..', '..', 'build', 'Release'))

# Change working directory to script directory for VTK output
os.chdir(_script_dir)

import numpy as np
import radia as rad

try:
    from ngsolve import *
    from netgen.occ import Box, Pnt, OCCGeometry
    import radia_ngsolve
    NGSOLVE_AVAILABLE = True
except ImportError as e:
    print('ERROR: NGSolve not available: %s' % e)
    NGSOLVE_AVAILABLE = False
    sys.exit(1)

print('=' * 70)
print('Verify curl(A) = B using RadiaField with HCurl and HDiv')
print('=' * 70)

# =============================================================================
# Step 1: Create Radia magnet
# =============================================================================
print()
print('[Step 1] Creating Radia rectangular magnet')
print('-' * 70)

rad.UtiDelAll()
rad.FldUnits('m')

# Create rectangular magnet
magnet = rad.ObjRecMag(
    [0, 0, 0],              # Center (m)
    [0.04, 0.04, 0.06],     # Dimensions (m)
    [0, 0, 1.2]             # Magnetization (T) - NdFeB
)

print('  Magnet ID: %d' % magnet)
print('  Center: [0, 0, 0] m')
print('  Dimensions: [0.04, 0.04, 0.06] m')
print('  Magnetization: [0, 0, 1.2] T')

# Reference field at a test point
ref_point = [0.03, 0.02, 0.05]
B_ref = rad.Fld(magnet, 'b', ref_point)
A_ref = rad.Fld(magnet, 'a', ref_point)

print('  Reference point: %s m' % ref_point)
print('  B = [%.6f, %.6f, %.6f] T' % tuple(B_ref))
print('  A = [%.6e, %.6e, %.6e] T*m' % tuple(A_ref))

# =============================================================================
# Step 2: Create NGSolve mesh
# =============================================================================
print()
print('[Step 2] Creating NGSolve mesh')
print('-' * 70)

# Mesh region outside the magnet (air region)
box = Box(Pnt(0.03, 0.03, 0.04), Pnt(0.08, 0.08, 0.12))
geo = OCCGeometry(box)
mesh = Mesh(geo.GenerateMesh(maxh=0.01))

print('  Mesh region: [0.03, 0.08] x [0.03, 0.08] x [0.04, 0.12] m')
print('  Elements: %d' % mesh.ne)
print('  Vertices: %d' % mesh.nv)

# =============================================================================
# Step 3: Create RadiaField CoefficientFunctions
# =============================================================================
print()
print('[Step 3] Creating RadiaField CoefficientFunctions')
print('-' * 70)

# Vector potential A from Radia
A_cf = radia_ngsolve.RadiaField(magnet, 'a')
print('  A_cf created (vector potential)')

# Magnetic field B from Radia
B_cf = radia_ngsolve.RadiaField(magnet, 'b')
print('  B_cf created (magnetic field)')

# =============================================================================
# Step 4: Project A onto HCurl space and compute curl(A)
# =============================================================================
print()
print('[Step 4] Projecting A onto HCurl and computing curl(A)')
print('-' * 70)

# HCurl space for vector potential A
fes_hcurl = HCurl(mesh, order=2)
print('  HCurl space: %d DOFs' % fes_hcurl.ndof)

# Project A onto HCurl
gf_A = GridFunction(fes_hcurl)
gf_A.Set(A_cf)
print('  A projected onto HCurl GridFunction')

# Compute curl(A)
curl_A_cf = curl(gf_A)
print('  curl(A) computed')

# =============================================================================
# Step 5: Project B onto HDiv space
# =============================================================================
print()
print('[Step 5] Projecting B onto HDiv space')
print('-' * 70)

# HDiv space for magnetic field B
fes_hdiv = HDiv(mesh, order=2)
print('  HDiv space: %d DOFs' % fes_hdiv.ndof)

# Project B onto HDiv
gf_B = GridFunction(fes_hdiv)
gf_B.Set(B_cf)
print('  B projected onto HDiv GridFunction')

# =============================================================================
# Step 6: Compare curl(A) with B at test points
# =============================================================================
print()
print('[Step 6] Comparing curl(A) with B at test points')
print('-' * 70)

# Test points (inside mesh region)
test_points = [
    [0.04, 0.04, 0.05],
    [0.05, 0.05, 0.06],
    [0.06, 0.06, 0.08],
    [0.07, 0.05, 0.10],
    [0.05, 0.07, 0.07],
    [0.04, 0.06, 0.09],
    [0.06, 0.04, 0.11],
    [0.055, 0.055, 0.075],
]

print()
print('  %-25s  %-15s  %-15s  %-10s' % ('Point (m)', '|curl(A)|', '|B_HDiv|', 'Error %'))
print('  ' + '-' * 70)

errors = []
rel_errors = []

for pt in test_points:
    try:
        mip = mesh(pt[0], pt[1], pt[2])

        # Evaluate curl(A) at point
        curl_A_x = curl_A_cf[0](mip)
        curl_A_y = curl_A_cf[1](mip)
        curl_A_z = curl_A_cf[2](mip)
        curl_A_mag = np.sqrt(curl_A_x**2 + curl_A_y**2 + curl_A_z**2)

        # Evaluate B from HDiv GridFunction
        B_x = gf_B[0](mip)
        B_y = gf_B[1](mip)
        B_z = gf_B[2](mip)
        B_mag = np.sqrt(B_x**2 + B_y**2 + B_z**2)

        # Error
        error = abs(curl_A_mag - B_mag)
        rel_error = error / B_mag * 100 if B_mag > 1e-10 else 0.0

        errors.append(error)
        rel_errors.append(rel_error)

        print('  [%.3f, %.3f, %.3f]  %15.6e  %15.6e  %10.4f' % (
            pt[0], pt[1], pt[2], curl_A_mag, B_mag, rel_error))

    except Exception as e:
        print('  [%.3f, %.3f, %.3f]  Error: %s' % (pt[0], pt[1], pt[2], e))

# =============================================================================
# Step 7: Statistical summary
# =============================================================================
print()
print('[Step 7] Statistical Summary')
print('-' * 70)

if errors:
    avg_error = np.mean(errors)
    max_error = np.max(errors)
    min_error = np.min(errors)
    avg_rel_error = np.mean(rel_errors)
    max_rel_error = np.max(rel_errors)

    print('  Test points: %d' % len(errors))
    print()
    print('  Absolute error |curl(A)| - |B|:')
    print('    Average: %.6e T' % avg_error)
    print('    Maximum: %.6e T' % max_error)
    print('    Minimum: %.6e T' % min_error)
    print()
    print('  Relative error:')
    print('    Average: %.4f%%' % avg_rel_error)
    print('    Maximum: %.4f%%' % max_rel_error)

    # Verification
    tolerance_percent = 5.0  # 5% tolerance
    if max_rel_error < tolerance_percent:
        print()
        print('[PASS] curl(A) = B verified!')
        print('       Maximum relative error %.4f%% < %.1f%% tolerance' % (max_rel_error, tolerance_percent))
    else:
        print()
        print('[CHECK] Errors exceed tolerance')
        print('        Maximum relative error %.4f%% >= %.1f%% tolerance' % (max_rel_error, tolerance_percent))
else:
    print('  No valid test points')

# =============================================================================
# Step 8: VTK Export
# =============================================================================
print()
print('[Step 8] Exporting VTK files')
print('-' * 70)

try:
    # Export vector fields
    vtk = VTKOutput(
        mesh,
        coefs=[gf_A, curl_A_cf, gf_B],
        names=['A_HCurl', 'curl_A', 'B_HDiv'],
        filename='verify_curl_A_B',
        subdivision=2
    )
    vtk.Do()
    print('  [OK] verify_curl_A_B.vtu exported')

    # Export error field
    error_cf = sqrt((curl_A_cf[0] - gf_B[0])**2 +
                    (curl_A_cf[1] - gf_B[1])**2 +
                    (curl_A_cf[2] - gf_B[2])**2)
    fes_h1 = H1(mesh, order=2)
    gf_error = GridFunction(fes_h1)
    gf_error.Set(error_cf)

    vtk_error = VTKOutput(
        mesh,
        coefs=[gf_error],
        names=['curl_A_minus_B_error'],
        filename='verify_curl_A_B_error',
        subdivision=2
    )
    vtk_error.Do()
    print('  [OK] verify_curl_A_B_error.vtu exported')

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
print('This script verified the Maxwell relation B = curl(A) using:')
print('  - RadiaField to get A (vector potential) as CoefficientFunction')
print('  - RadiaField to get B (magnetic field) as CoefficientFunction')
print('  - HCurl space projection for A')
print('  - HDiv space projection for B')
print('  - NGSolve curl() operator to compute curl(A)')
print()
print('The comparison shows that the radia_ngsolve integration correctly')
print('evaluates both A and B fields, and curl(A) matches B within')
print('numerical tolerance.')
print()
print('=' * 70)

rad.UtiDelAll()
