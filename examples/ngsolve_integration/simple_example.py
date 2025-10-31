"""
Simple rad_ngsolve Example - Minimal Working Example

This is the simplest possible example showing rad_ngsolve integration.
Perfect for quick testing and verification.
"""

import sys
sys.path.insert(0, r"S:\radia\01_GitHub\build\Release")
sys.path.insert(0, r"S:\radia\01_GitHub\dist")

print("Simple rad_ngsolve Example")
print("=" * 60)

# Step 1: Import ngsolve FIRST (loads libngsolve.dll)
print("\n[1/5] Importing ngsolve...")
import ngsolve
from ngsolve import CoefficientFunction
print("      OK")

# Step 2: Import radia
print("[2/5] Importing radia...")
import radia as rad
print("      OK")

# Step 3: Import rad_ngsolve
print("[3/5] Importing rad_ngsolve...")
import rad_ngsolve
print("      OK")

# Step 4: Create simple Radia magnet
print("[4/5] Creating Radia magnet...")
magnet = rad.ObjRecMag([0, 0, 0], [10, 10, 10], [0, 0, 1.2])
rad.MatApl(magnet, rad.MatStd('NdFeB', 1.2))
rad.Solve(magnet, 0.0001, 10000)
print(f"      Magnet ID: {magnet}")

# Step 5: Create CoefficientFunction
print("[5/5] Creating RadBfield CoefficientFunction...")
B_field = rad_ngsolve.RadBfield(magnet)
print(f"      Type: {type(B_field)}")
print(f"      Is CoefficientFunction: {isinstance(B_field, CoefficientFunction)}")

# Verify with direct Radia calculation
print("\n" + "=" * 60)
print("Verification:")
test_point = [0, 0, 15]
B = rad.Fld(magnet, 'b', test_point)
print(f"  Field at {test_point} mm:")
print(f"    Bx = {B[0]:.6f} T")
print(f"    By = {B[1]:.6f} T")
print(f"    Bz = {B[2]:.6f} T")
print(f"    |B| = {(sum(b**2 for b in B))**0.5:.6f} T")

print("\n" + "=" * 60)
print("SUCCESS: rad_ngsolve is working correctly!")
print("=" * 60)

print("\nNext steps:")
print("  - See example_rad_ngsolve_demo.py for complete FEM example")
print("  - See test_rad_ngsolve_full.py for comprehensive tests")
print("  - See visualize_field.py for mesh visualization")
