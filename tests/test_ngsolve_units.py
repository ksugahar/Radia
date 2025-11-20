#!/usr/bin/env python
"""
Test rad_ngsolve coordinate unit conversion (mm vs m)
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'build', 'Release'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src', 'python'))

import numpy as np
import radia as rad
from ngsolve import *
import rad_ngsolve

print("=" * 80)
print("rad_ngsolve Unit Conversion Test")
print("=" * 80)
print()

# Create simple mesh
from netgen.occ import *
box = Box((0,0,0), (0.1,0.1,0.1))
geo = OCCGeometry(box)
mesh = Mesh(geo.GenerateMesh(maxh=0.05))
print(f"Mesh: {mesh.nv} vertices, {mesh.ne} elements")

# Create FE space
fes = HCurl(mesh, order=1)
print(f"FE Space: {fes.ndof} DOFs")
print()

# Test point in NGSolve coordinates (meters)
test_point_m = mesh(0.02, 0.0, 0.0)  # 0.02 m = 20 mm from origin

# Create magnet at origin [0, 0, 0] in MILLIMETERS
rad.FldUnits('mm')
magnet = rad.ObjRecMag([0, 0, 0], [10, 10, 10], [0, 0, 1.0])  # 10mm cube, magnetization in z

print("[Setup] Magnet: 10mm cube at origin [0, 0, 0] mm")
print(f"[Setup] Test point: [0.02, 0, 0] m = [20, 0, 0] mm")
print()

# Test 1: Default units='mm' (should work correctly)
print("[Test 1] units='mm' (default) - coordinates scaled by 1000")
B_cf_mm = rad_ngsolve.RadiaField(magnet, 'b', units='mm')
gf_mm = GridFunction(fes)
gf_mm.Set(B_cf_mm)
B_mm = gf_mm(test_point_m)
print(f"  B field at test point: [{B_mm[0]:.6f}, {B_mm[1]:.6f}, {B_mm[2]:.6f}] T")

# Direct Radia calculation at [20, 0, 0] mm for verification
B_direct = rad.Fld(magnet, 'b', [20, 0, 0])
print(f"  rad.Fld([20, 0, 0] mm): [{B_direct[0]:.6f}, {B_direct[1]:.6f}, {B_direct[2]:.6f}] T")

diff_mm = np.array(B_mm) - np.array(B_direct)
error_mm = np.linalg.norm(diff_mm)
print(f"  Error: {error_mm:.2e}")
print()

# Test 2: Recreate magnet in METERS for units='m' test
rad.UtiDelAll()
rad.FldUnits('m')
magnet_m = rad.ObjRecMag([0, 0, 0], [0.01, 0.01, 0.01], [0, 0, 1.0])  # 0.01m = 10mm cube
print("[Test 2] units='m' - coordinates scaled by 1")
B_cf_m = rad_ngsolve.RadiaField(magnet_m, 'b', units='m')
gf_m = GridFunction(fes)
gf_m.Set(B_cf_m)
B_m = gf_m(test_point_m)
print(f"  B field at test point: [{B_m[0]:.6f}, {B_m[1]:.6f}, {B_m[2]:.6f}] T")

# Direct Radia calculation at [0.02, 0, 0] m for verification
B_direct_m = rad.Fld(magnet_m, 'b', [0.02, 0, 0])
print(f"  rad.Fld([0.02, 0, 0] m): [{B_direct_m[0]:.6f}, {B_direct_m[1]:.6f}, {B_direct_m[2]:.6f}] T")

diff_m = np.array(B_m) - np.array(B_direct_m)
error_m = np.linalg.norm(diff_m)
print(f"  Error: {error_m:.2e}")
print()

# Verification: Both methods should give same result (within interpolation error)
B_diff = np.array(B_mm) - np.array(B_m)
B_diff_norm = np.linalg.norm(B_diff)
print("=" * 80)
print(f"Consistency check: |B(mm) - B(m)| = {B_diff_norm:.2e}")

# Success criteria
tol = 1e-4
if error_mm < tol and error_m < tol and B_diff_norm < tol:
    print(f"\n[PASS] Both unit conversions accurate (error < {tol})")
    print(f"  units='mm' error: {error_mm:.2e}")
    print(f"  units='m' error:  {error_m:.2e}")
    print(f"  Consistency:      {B_diff_norm:.2e}")
else:
    print(f"\n[FAIL] Unit conversion errors:")
    print(f"  units='mm' error: {error_mm:.2e} {'[OK]' if error_mm < tol else '[FAIL]'}")
    print(f"  units='m' error:  {error_m:.2e} {'[OK]' if error_m < tol else '[FAIL]'}")
    print(f"  Consistency:      {B_diff_norm:.2e} {'[OK]' if B_diff_norm < tol else '[FAIL]'}")

print("=" * 80)
