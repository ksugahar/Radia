#!/usr/bin/env python
"""
Verify coord_scale_ direction: 2m magnet vs 2000mm magnet should give same field
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
print("Unit Conversion Direction Verification Test")
print("=" * 80)
print()

# Create simple mesh
from netgen.occ import *
box = Box((0,0,0), (0.1,0.1,0.1))
geo = OCCGeometry(box)
mesh = Mesh(geo.GenerateMesh(maxh=0.05))
print(f"Mesh: {mesh.nv} vertices")

# Create FE space
fes = HCurl(mesh, order=1)
print(f"FE Space: {fes.ndof} DOFs")
print()

# Test point in NGSolve coordinates (meters)
test_point_m = mesh(0.05, 0.05, 0.05)  # Center of mesh
print(f"Test point (NGSolve coords): [0.05, 0.05, 0.05] m")
print()

# ============================================================================
# Test 1: Create magnet in METERS (2m cube at origin)
# ============================================================================
print("[Test 1] Magnet in METERS: 2m cube at origin")
rad.UtiDelAll()
rad.FldUnits('m')
magnet_m = rad.ObjRecMag([0, 0, 0], [2, 2, 2], [0, 0, 1.0])  # 2m cube
print("  rad.FldUnits('m')")
print("  Magnet: center=[0, 0, 0] m, size=[2, 2, 2] m, M=[0, 0, 1]")

# Direct Radia calculation at [0.05, 0.05, 0.05] m
B_direct_m = rad.Fld(magnet_m, 'b', [0.05, 0.05, 0.05])
print(f"  rad.Fld([0.05, 0.05, 0.05] m) = [{B_direct_m[0]:.6f}, {B_direct_m[1]:.6f}, {B_direct_m[2]:.6f}] T")

# NGSolve calculation with units='m'
B_cf_m = rad_ngsolve.RadiaField(magnet_m, 'b', units='m')
gf_m = GridFunction(fes)
gf_m.Set(B_cf_m)
B_ngsolve_m = gf_m(test_point_m)
print(f"  rad_ngsolve (units='m') = [{B_ngsolve_m[0]:.6f}, {B_ngsolve_m[1]:.6f}, {B_ngsolve_m[2]:.6f}] T")

error_m = np.linalg.norm(np.array(B_ngsolve_m) - np.array(B_direct_m))
print(f"  Error: {error_m:.2e}")
print()

# ============================================================================
# Test 2: Create magnet in MILLIMETERS (2000mm cube at origin)
# ============================================================================
print("[Test 2] Magnet in MILLIMETERS: 2000mm cube at origin")
rad.UtiDelAll()
rad.FldUnits('mm')
magnet_mm = rad.ObjRecMag([0, 0, 0], [2000, 2000, 2000], [0, 0, 1.0])  # 2000mm = 2m cube
print("  rad.FldUnits('mm')")
print("  Magnet: center=[0, 0, 0] mm, size=[2000, 2000, 2000] mm, M=[0, 0, 1]")

# Direct Radia calculation at [50, 50, 50] mm (= 0.05 m)
B_direct_mm = rad.Fld(magnet_mm, 'b', [50, 50, 50])
print(f"  rad.Fld([50, 50, 50] mm) = [{B_direct_mm[0]:.6f}, {B_direct_mm[1]:.6f}, {B_direct_mm[2]:.6f}] T")

# NGSolve calculation with units='mm'
B_cf_mm = rad_ngsolve.RadiaField(magnet_mm, 'b', units='mm')
gf_mm = GridFunction(fes)
gf_mm.Set(B_cf_mm)
B_ngsolve_mm = gf_mm(test_point_m)
print(f"  rad_ngsolve (units='mm') = [{B_ngsolve_mm[0]:.6f}, {B_ngsolve_mm[1]:.6f}, {B_ngsolve_mm[2]:.6f}] T")

error_mm = np.linalg.norm(np.array(B_ngsolve_mm) - np.array(B_direct_mm))
print(f"  Error: {error_mm:.2e}")
print()

# ============================================================================
# Verification: Both should give IDENTICAL field values
# ============================================================================
print("=" * 80)
print("Verification: 2m magnet vs 2000mm magnet")
print("=" * 80)

diff_direct = np.linalg.norm(np.array(B_direct_m) - np.array(B_direct_mm))
diff_ngsolve = np.linalg.norm(np.array(B_ngsolve_m) - np.array(B_ngsolve_mm))

print(f"rad.Fld() consistency:    |B(2m) - B(2000mm)| = {diff_direct:.2e}")
print(f"rad_ngsolve consistency:  |B(2m) - B(2000mm)| = {diff_ngsolve:.2e}")
print()

# Success criteria
tol = 1e-4
if error_m < tol and error_mm < tol and diff_direct < 1e-10 and diff_ngsolve < tol:
    print(f"[PASS] Unit conversion is correct!")
    print(f"  - units='m' error:     {error_m:.2e} {'[OK]' if error_m < tol else '[FAIL]'}")
    print(f"  - units='mm' error:    {error_mm:.2e} {'[OK]' if error_mm < tol else '[FAIL]'}")
    print(f"  - rad.Fld() matches:   {diff_direct:.2e} {'[OK]' if diff_direct < 1e-10 else '[FAIL]'}")
    print(f"  - rad_ngsolve matches: {diff_ngsolve:.2e} {'[OK]' if diff_ngsolve < tol else '[FAIL]'}")
else:
    print(f"[FAIL] Unit conversion has issues:")
    print(f"  - units='m' error:     {error_m:.2e} {'[OK]' if error_m < tol else '[FAIL]'}")
    print(f"  - units='mm' error:    {error_mm:.2e} {'[OK]' if error_mm < tol else '[FAIL]'}")
    print(f"  - rad.Fld() matches:   {diff_direct:.2e} {'[OK]' if diff_direct < 1e-10 else '[FAIL]'}")
    print(f"  - rad_ngsolve matches: {diff_ngsolve:.2e} {'[OK]' if diff_ngsolve < tol else '[FAIL]'}")

print("=" * 80)

# Additional debug info
print()
print("Debug: Coordinate conversion check")
print(f"  NGSolve point: 0.05 m")
print(f"  Expected Radia point (m):  0.05 m  (with coord_scale_=1.0)")
print(f"  Expected Radia point (mm): 50 mm   (with coord_scale_=1000.0)")
print(f"  Actual coord_scale_ for 'm':  1.0 if units == 'm' else 1000.0")
print(f"  Actual coord_scale_ for 'mm': 1.0 if units == 'm' else 1000.0")
