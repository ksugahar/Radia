import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../build/Release'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/python'))

from ngsolve import *
from netgen.occ import *
import radia as rad
import rad_ngsolve
import numpy as np

print("="*70)
print("Simple GridFunction.Set() Test")
print("="*70)

# Set Radia to use meters (required for NGSolve integration)
rad.FldUnits('m')

# Create simple magnet
magnet = rad.ObjRecMag([0, 0.002, 0], [0.001, 0.001, 0.001], [0, 1, 0])  # meters
print(f"\n[1] Created magnet at (0, 0.002, 0) m, M=[0, 1, 0] T")

# Test point
point_m = [0, 0.002, 0]  # Magnet center in meters
B_radia = rad.Fld(magnet, 'b', point_m)
print(f"\n[2] rad.Fld at magnet center: {B_radia}")

# Create mesh
box = Box((-0.006, -0.006, -0.006), (0.006, 0.006, 0.006))
mesh = Mesh(OCCGeometry(box).GenerateMesh(maxh=0.001))
print(f"\n[9] Created mesh: {mesh.ne} elements")

# Create CoefficientFunction
B_cf = rad_ngsolve.RadiaField(magnet, 'b')
print(f"\n[4] Created CoefficientFunction")

# Test direct CF evaluation
point_m = (0, 0.002, 0)  # Meters (matches Radia units)
mip = mesh(*point_m)
B_cf_direct = B_cf(mip)
print(f"\n[5] Direct CF evaluation: {B_cf_direct}")
print(f"    Expected: {B_radia}")
print(f"    Match: {np.allclose(B_cf_direct, B_radia, rtol=1e-6)}")

# Create GridFunction with HCurl space
fes = HCurl(mesh, order=2)
B_gf = GridFunction(fes)
print(f"\n[6] Created GridFunction with HCurl space ({fes.ndof} DOFs)")

# Set GridFunction from CoefficientFunction
print(f"\n[7] Calling GridFunction.Set(B_cf)...")
B_gf.Set(B_cf)
print(f"    [OK] Set() completed")

# Evaluate GridFunction at same point
B_gf_eval = B_gf(mip)
print(f"\n[8] GridFunction evaluation: {B_gf_eval}")
print(f"    Expected: {B_radia}")
print(f"    Match: {np.allclose(B_gf_eval, B_radia, rtol=1e-6)}")

# Summary
print("\n" + "="*70)
if np.allclose(B_cf_direct, B_radia, rtol=1e-6):
    print("[OK] Direct CF evaluation works correctly")
else:
    print("[FAIL] Direct CF evaluation is wrong!")

if np.allclose(B_gf_eval, B_radia, rtol=1e-6):
    print("[OK] GridFunction.Set() works correctly")
else:
    print("[FAIL] GridFunction.Set() is broken!")
print("="*70)
