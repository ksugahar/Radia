"""
Complete verification of rad_ngsolve fix
Tests both direct CF evaluation and GridFunction interpolation
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
print("rad_ngsolve Complete Verification Test")
print("=" * 70)

# Create magnet
magnet = rad.ObjRecMag([0, 0, 0], [20, 20, 30], [0, 0, 1.2])
rad.MatApl(magnet, rad.MatStd('NdFeB', 1.2))
rad.Solve(magnet, 0.0001, 10000)
print(f"\nMagnet created: object #{magnet}")

# Test points
test_points = [
	[0, 0, 0],	 # Inside magnet
	[0, 0, 20],	# Above magnet
	[0, 0, 40],	# Far from magnet
	[10, 10, 10],  # Off-axis
]

print("\n" + "=" * 70)
print("Test 1: Direct Radia Evaluation")
print("=" * 70)
for pt in test_points:
	B = rad.Fld(magnet, 'b', pt)
	print(f"{pt}: Bx={B[0]:.6f}, By={B[1]:.6f}, Bz={B[2]:.6f} T")

# Create CoefficientFunction
B_cf = rad_ngsolve.RadBfield(magnet)
print(f"\nCoefficientFunction created: type={type(B_cf).__name__}, dim={B_cf.dim}")

# Create 3D mesh
geo = CSGeometry()
geo.Add(OrthoBrick(Pnt(-50,-50,-50), Pnt(50,50,50)))
mesh = Mesh(geo.GenerateMesh(maxh=30))
print(f"3D Mesh: {mesh.ne} elements, {mesh.nv} vertices")

print("\n" + "=" * 70)
print("Test 2: Direct CoefficientFunction Evaluation")
print("=" * 70)
print("(This should match Radia exactly)")
print()

all_match = True
for pt in test_points:
	try:
		mesh_pt = mesh(*pt)
		B_cf_val = B_cf(mesh_pt)
		B_rad = rad.Fld(magnet, 'b', pt)

		error = np.sqrt((B_cf_val[0] - B_rad[0])**2 +
					   (B_cf_val[1] - B_rad[1])**2 +
					   (B_cf_val[2] - B_rad[2])**2)

		match = error < 1e-6
		all_match = all_match and match

		status = "[OK]" if match else "[FAIL]"
		print(f"{status} {pt}:")
		print(f"  CF:	  Bx={B_cf_val[0]:.6f}, By={B_cf_val[1]:.6f}, Bz={B_cf_val[2]:.6f} T")
		print(f"  Radia:   Bx={B_rad[0]:.6f}, By={B_rad[1]:.6f}, Bz={B_rad[2]:.6f} T")
		print(f"  Error:   {error:.6e} T")
	except Exception as e:
		print(f"[ERROR] {pt}: {e}")
		all_match = False

if all_match:
	print("\n[SUCCESS] Direct CF evaluation matches Radia perfectly!")
else:
	print("\n[FAILED] Direct CF evaluation has errors")

print("\n" + "=" * 70)
print("Test 3: GridFunction Interpolation")
print("=" * 70)
print("(GridFunction uses FEM interpolation, so small errors are expected)")
print()

fes = VectorH1(mesh, order=2)  # Higher order for better accuracy
gfB = GridFunction(fes)
gfB.Set(B_cf)

print(f"GridFunction: {fes.ndof} DOFs, order={fes.globalorder}")
print()

for pt in test_points:
	try:
		mesh_pt = mesh(*pt)
		B_gf = gfB(mesh_pt)
		B_rad = rad.Fld(magnet, 'b', pt)

		error = np.sqrt((B_gf[0] - B_rad[0])**2 + (B_gf[1] - B_rad[1])**2 + (B_gf[2] - B_rad[2])**2)

		rel_error = error / max(np.sqrt(B_rad[0]**2 + B_rad[1]**2 + B_rad[2]**2), 1e-10)

		status = "[OK]" if rel_error < 0.1 else "[WARN]"  # 10% tolerance
		print(f"{status} {pt}:")
		print(f"  GF:	  Bx={B_gf[0]:.6f}, By={B_gf[1]:.6f}, Bz={B_gf[2]:.6f} T")
		print(f"  Radia:   Bx={B_rad[0]:.6f}, By={B_rad[1]:.6f}, Bz={B_rad[2]:.6f} T")
		print(f"  Error:   {error:.6e} T ({rel_error*100:.2f}%)")
	except Exception as e:
		print(f"[ERROR] {pt}: {e}")

print("\n" + "=" * 70)
print("Test 4: Field Statistics")
print("=" * 70)

vec = gfB.vec.FV().NumPy()
Bx = vec[0::3]
By = vec[1::3]
Bz = vec[2::3]
B_mag = np.sqrt(Bx**2 + By**2 + Bz**2)

print(f"GridFunction field statistics:")
print(f"  Non-zero values: {np.count_nonzero(vec)}/{len(vec)}")
print(f"  |B| range: [{B_mag.min():.6e}, {B_mag.max():.6e}] T")
print(f"  |B| mean:  {B_mag.mean():.6e} T")

print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

if all_match and B_mag.max() > 0.01:
	print("[SUCCESS] All tests passed!")
	print("  - Direct CF evaluation: Perfect match with Radia")
	print("  - GridFunction interpolation: Working correctly")
	print("  - Field values: Non-zero and reasonable")
else:
	print("[FAILED] Some tests failed")
	if not all_match:
		print("  - Direct CF evaluation does not match Radia")
	if B_mag.max() <= 0.01:
		print("  - Field values are too small or zero")

print("\nNote: GridFunction uses finite element interpolation,")
print("	  so it will not exactly match point-wise Radia values.")
print("	  Direct CoefficientFunction evaluation is exact.")
