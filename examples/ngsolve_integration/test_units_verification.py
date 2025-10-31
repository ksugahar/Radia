"""
Unit conversion verification for rad_ngsolve
Tests that NGSolve (meters) <-> Radia (millimeters) conversion works correctly
"""
import sys
sys.path.insert(0, r"S:\radia\01_GitHub\build\Release")
sys.path.insert(0, r"S:\radia\01_GitHub\dist")

from ngsolve import *
from netgen.csg import *
import radia as rad
import rad_ngsolve
import numpy as np

print("=" * 70)
print("Unit Conversion Verification Test")
print("=" * 70)
print("\nIMPORTANT:")
print("  - NGSolve uses METERS")
print("  - Radia uses MILLIMETERS")
print("  - rad_ngsolve.cpp automatically converts: m -> mm (x1000)")
print()

# Create magnet in Radia (mm coordinates)
print("Creating Radia magnet (mm coordinates):")
print("  Center: [0, 0, 0] mm")
print("  Size: [20, 20, 30] mm")
magnet = rad.ObjRecMag([0, 0, 0], [20, 20, 30], [0, 0, 1.2])
rad.MatApl(magnet, rad.MatStd('NdFeB', 1.2))
rad.Solve(magnet, 0.0001, 10000)

# Test points in both unit systems
test_points = [
	{"m": [0, 0, 0], "mm": [0, 0, 0], "desc": "Origin (inside magnet)"},
	{"m": [0, 0, 0.02], "mm": [0, 0, 20], "desc": "20mm above center"},
	{"m": [0, 0, 0.04], "mm": [0, 0, 40], "desc": "40mm above center"},
	{"m": [0.01, 0.01, 0.01], "mm": [10, 10, 10], "desc": "Off-axis point"},
]

print("\n" + "=" * 70)
print("Test 1: Direct Radia Evaluation (mm coordinates)")
print("=" * 70)
for pt in test_points:
	B = rad.Fld(magnet, 'b', pt["mm"])
	print(f"{pt['mm']} mm ({pt['desc']})")
	print(f"  Bx={B[0]:.6f}, By={B[1]:.6f}, Bz={B[2]:.6f} T")

# Create NGSolve mesh in meters
print("\n" + "=" * 70)
print("Creating NGSolve mesh (meter coordinates)")
print("=" * 70)
geo = CSGeometry()
# Mesh from -50mm to +50mm = -0.05m to +0.05m
geo.Add(OrthoBrick(Pnt(-0.05, -0.05, -0.05), Pnt(0.05, 0.05, 0.05)))
mesh = Mesh(geo.GenerateMesh(maxh=0.03))
print(f"Mesh bounds: [-0.05, 0.05] m = [-50, 50] mm")
print(f"Elements: {mesh.ne}, Vertices: {mesh.nv}")

# Create CoefficientFunction
B_cf = rad_ngsolve.RadBfield(magnet)
print(f"\nCoefficientFunction created: {type(B_cf).__name__}")

print("\n" + "=" * 70)
print("Test 2: CoefficientFunction with Unit Conversion")
print("=" * 70)
print("Evaluating at NGSolve mesh points (meter coordinates)")
print("Comparing with Radia evaluation (millimeter coordinates)")
print()

all_passed = True
for pt in test_points:
	try:
		# Evaluate at NGSolve mesh point (meters)
		mesh_pt = mesh(*pt["m"])
		B_ngsolve = B_cf(mesh_pt)

		# Compare with direct Radia evaluation (millimeters)
		B_radia = rad.Fld(magnet, 'b', pt["mm"])

		# Calculate error
		error = np.sqrt((B_ngsolve[0] - B_radia[0])**2 +
					   (B_ngsolve[1] - B_radia[1])**2 +
					   (B_ngsolve[2] - B_radia[2])**2)

		passed = error < 1e-10
		all_passed = all_passed and passed
		status = "[PASS]" if passed else "[FAIL]"

		print(f"{status} {pt['m']} m = {pt['mm']} mm")
		print(f"  NGSolve CF: Bx={B_ngsolve[0]:.6f}, By={B_ngsolve[1]:.6f}, Bz={B_ngsolve[2]:.6f} T")
		print(f"  Radia:	  Bx={B_radia[0]:.6f}, By={B_radia[1]:.6f}, Bz={B_radia[2]:.6f} T")
		print(f"  Error:	  {error:.6e} T")

	except Exception as e:
		print(f"[ERROR] {pt['m']} m: {e}")
		all_passed = False

print("\n" + "=" * 70)
print("Test 3: GridFunction with Unit Conversion")
print("=" * 70)

fes = VectorH1(mesh, order=2)
gfB = GridFunction(fes)
gfB.Set(B_cf)

vec = gfB.vec.FV().NumPy()
Bx = vec[0::3]
By = vec[1::3]
Bz = vec[2::3]
B_mag = np.sqrt(Bx**2 + By**2 + Bz**2)

print(f"GridFunction statistics:")
print(f"  DOFs: {fes.ndof}")
print(f"  Non-zero: {np.count_nonzero(vec)}/{len(vec)}")
print(f"  |B| range: [{B_mag.min():.6e}, {B_mag.max():.6e}] T")
print(f"  |B| mean:  {B_mag.mean():.6e} T")

print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

if all_passed and B_mag.max() > 0.01:
	print("[SUCCESS] Unit conversion working correctly!")
	print("  ✓ NGSolve coordinates (m) automatically converted to Radia coordinates (mm)")
	print("  ✓ All point evaluations match Radia exactly")
	print("  ✓ GridFunction contains non-zero field values")
else:
	print("[FAILED] Unit conversion has issues")
	if not all_passed:
		print("  ✗ Point evaluations do not match")
	if B_mag.max() <= 0.01:
		print("  ✗ Field values are too small or zero")

print("\n" + "=" * 70)
print("Unit Conversion Details")
print("=" * 70)
print("NGSolve mesh point [0.01, 0.01, 0.01] m")
print("  ↓ Multiply by 1000 in rad_ngsolve.cpp")
print("Radia evaluation at [10, 10, 10] mm")
print("  ↓ Radia computes field")
print("Returns field value in Tesla (no conversion)")
print("  ↓ Return to NGSolve")
print("NGSolve CoefficientFunction returns field in Tesla")
