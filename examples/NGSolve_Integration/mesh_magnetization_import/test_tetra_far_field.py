"""
Test tetrahedral mesh with far-field evaluation points
Implementation analysis:
- Evaluate field far from element surfaces
- Use points >> element size
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../build/Release'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/python'))

from numpy import *
import radia as rad
from netgen.occ import *
from ngsolve import Mesh
from netgen_mesh_import import netgen_mesh_to_radia

print("="*70)
print("Tetrahedral Mesh Test: Far-Field Evaluation")
print("="*70)

rad.FldUnits('m')

# Parameters
cube_size = 0.1
half_size = cube_size / 2
mu_r = 100
mu0 = 4*pi*1e-7
ksi = mu_r - 1

# Create coarse mesh
maxh = cube_size / 3  # ~0.033m elements
cube_geom = Box((-half_size, -half_size, -half_size), (half_size, half_size, half_size))
cube_geom.mat("magnetic")
geo = OCCGeometry(cube_geom)
mesh = Mesh(geo.GenerateMesh(maxh=maxh))

print(f"\nMesh: {mesh.ne} elements, {mesh.nv} vertices")
print(f"  Element size: ~{maxh:.4f} m")

# Import to Radia
print("\nImporting mesh to Radia...")
cube_tet = netgen_mesh_to_radia(
    mesh,
    material={'magnetization': [0, 0, 0]},
    units='m',
    combine=True,
    verbose=False
)

# Apply MatLin
mat = rad.MatLin([ksi, ksi], [0, 0, 1])
rad.MatApl(cube_tet, mat)

# Background field
H_s_tesla = mu0 * 1.0

def uniform_field(pos):
    return [0.0, 0.0, H_s_tesla]

bkg_field = rad.ObjBckgCF(uniform_field)
container_tet = rad.ObjCnt([cube_tet, bkg_field])

# Solve
print("Solving...")
result = rad.Solve(container_tet, 0.0001, 10000)
print(f"  Solve result: {result}")

# Reference: Built-in cube
print("\nReference: Built-in cube...")
cube_builtin = rad.ObjRecMag([float(-half_size), float(-half_size), float(-half_size)],
                              [float(cube_size), float(cube_size), float(cube_size)],
                              [0, 0, 0])
rad.MatApl(cube_builtin, mat)
container_builtin = rad.ObjCnt([cube_builtin, bkg_field])
rad.Solve(container_builtin, 0.0001, 10000)

# Test at FAR-FIELD points
# Rule: Distance >> element_size (use 10x element size as minimum)
min_distance = 10 * maxh  # ~0.33m

print(f"\n" + "="*70)
print("Far-Field Evaluation (distance >> element size)")
print("="*70)
print(f"  Minimum safe distance: {min_distance:.3f} m (~10× element size)")

test_points_far = [
    ([0.5, 0, 0], "x-axis (far)"),
    ([0, 0.5, 0], "y-axis (far)"),
    ([0, 0, 0.5], "z-axis (far)"),
    ([1.0, 0, 0], "x-axis (very far)"),
    ([0, 0, 1.0], "z-axis (very far)"),
]

print(f"\n{'Point [m]':<25} {'Tet Mesh |H|':<20} {'Built-in |H|':<20} {'Error [%]':<15}")
print("-" * 80)

for pt, desc in test_points_far:
    dist = linalg.norm(pt)

    H_tet = rad.Fld(container_tet, 'h', pt)
    H_tet_mag = linalg.norm(H_tet)

    H_ref = rad.Fld(container_builtin, 'h', pt)
    H_ref_mag = linalg.norm(H_ref)

    error = 100 * abs(H_tet_mag - H_ref_mag) / H_ref_mag if H_ref_mag > 1e-10 else 0

    print(f"{desc:<25} {H_tet_mag:>12.6f} A/m  {H_ref_mag:>12.6f} A/m  {error:>10.2f}%")

# Also test NEAR-FIELD for comparison
print(f"\n" + "="*70)
print("Near-Field Evaluation (may have errors)")
print("="*70)

test_points_near = [
    ([0.2, 0, 0], "x-axis (near)"),
    ([0, 0.2, 0], "y-axis (near)"),
    ([0, 0, 0.2], "z-axis (near)"),
]

print(f"\n{'Point [m]':<25} {'Tet Mesh |H|':<20} {'Built-in |H|':<20} {'Error [%]':<15}")
print("-" * 80)

for pt, desc in test_points_near:
    dist = linalg.norm(pt)

    H_tet = rad.Fld(container_tet, 'h', pt)
    H_tet_mag = linalg.norm(H_tet)

    H_ref = rad.Fld(container_builtin, 'h', pt)
    H_ref_mag = linalg.norm(H_ref)

    error = 100 * abs(H_tet_mag - H_ref_mag) / H_ref_mag if H_ref_mag > 1e-10 else 0

    print(f"{desc:<25} {H_tet_mag:>12.6f} A/m  {H_ref_mag:>12.6f} A/m  {error:>10.2f}%")

print("\n" + "="*70)
print("Analysis")
print("="*70)

if isnan(result[0]):
    print("\n[!] Solver did not converge (NaN result)")
    print("    Tetrahedral mesh + MatLin is fundamentally incompatible")
else:
    print("\n[!] Solver converged, but field values may still be incorrect")
    print("    Check far-field vs near-field errors above")

print("\n" + "="*70)
print("Test complete")
print("="*70)
