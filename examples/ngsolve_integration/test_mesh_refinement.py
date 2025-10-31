"""
Test GridFunction.Set() interpolation with mesh refinement
Tests if finer mesh reduces interpolation error
"""
import sys
sys.path.insert(0, r"S:\radia\01_GitHub\build\Release")
sys.path.insert(0, r"S:\radia\01_GitHub\dist")

from ngsolve import *
from netgen.csg import CSGeometry, OrthoBrick, Pnt
import radia as rad
import rad_ngsolve
import numpy as np

print("=" * 70)
print("Mesh Refinement Test - GridFunction.Set() Convergence")
print("=" * 70)

# Create magnet
magnet = rad.ObjRecMag([0, 0, 0], [20, 20, 30], [0, 0, 1.2])
rad.MatApl(magnet, rad.MatStd('NdFeB', 1.2))
rad.Solve(magnet, 0.0001, 10000)

# Create CoefficientFunction
B_cf = rad_ngsolve.RadBfield(magnet)

# Test points
test_point = (0.0, 0.0, 0.0)  # Origin (inside magnet)
pt_mm = [0, 0, 0]
B_radia = rad.Fld(magnet, 'b', pt_mm)

print(f"\nReference (Direct Radia): Bz = {B_radia[2]:.6f} T")
print(f"Test point: {test_point} m = {pt_mm} mm")
print()

# Test different mesh sizes
mesh_sizes = [0.02, 0.01, 0.005, 0.002]
print("=" * 70)
print("Testing different mesh sizes:")
print("=" * 70)

results = []

for maxh in mesh_sizes:
    print(f"\nmaxh = {maxh} m ({maxh*1000} mm)")
    print("-" * 70)

    # Create mesh
    geo = CSGeometry()
    geo.Add(OrthoBrick(Pnt(-0.05, -0.05, -0.05), Pnt(0.05, 0.05, 0.05)))
    mesh = Mesh(geo.GenerateMesh(maxh=maxh))

    print(f"  Mesh: {mesh.ne} elements, {mesh.nv} vertices")

    # Test 1: Direct CoefficientFunction evaluation
    mesh_pt = mesh(*test_point)
    B_cf_val = B_cf(mesh_pt)
    cf_error = abs(B_cf_val[2] - B_radia[2])

    # Test 2: GridFunction interpolation
    fes = VectorH1(mesh, order=2)
    gfB = GridFunction(fes)
    gfB.Set(B_cf)
    B_gf_val = gfB(mesh_pt)
    gf_error = abs(B_gf_val[2] - B_radia[2])
    gf_rel_error = gf_error / B_radia[2] * 100

    print(f"  CoefficientFunction: Bz = {B_cf_val[2]:.6f} T, error = {cf_error:.6e} T")
    print(f"  GridFunction:        Bz = {B_gf_val[2]:.6f} T, error = {gf_error:.6e} T ({gf_rel_error:.2f}%)")

    results.append({
        'maxh': maxh,
        'elements': mesh.ne,
        'vertices': mesh.nv,
        'dofs': fes.ndof,
        'cf_error': cf_error,
        'gf_error': gf_error,
        'gf_rel_error': gf_rel_error
    })

print("\n" + "=" * 70)
print("Summary")
print("=" * 70)
print(f"{'maxh (m)':<12} {'Elements':<12} {'DOFs':<12} {'GF Error (T)':<15} {'Rel Error (%)':<15}")
print("-" * 70)

for r in results:
    print(f"{r['maxh']:<12.4f} {r['elements']:<12} {r['dofs']:<12} {r['gf_error']:<15.6e} {r['gf_rel_error']:<15.2f}")

print("\n" + "=" * 70)
print("Convergence Analysis")
print("=" * 70)

# Check if error is decreasing
decreasing = all(results[i]['gf_error'] > results[i+1]['gf_error'] for i in range(len(results)-1))
if decreasing:
    print("[OK] GridFunction error decreases with mesh refinement")
else:
    print("[WARNING] GridFunction error does not consistently decrease")

# Check if we can get below 1% error
min_error = min(r['gf_rel_error'] for r in results)
if min_error < 1.0:
    print(f"[OK] Achieved < 1% relative error: {min_error:.4f}%")
    for r in results:
        if r['gf_rel_error'] == min_error:
            print(f"    At maxh = {r['maxh']} m ({r['elements']} elements)")
else:
    print(f"[INFO] Minimum relative error: {min_error:.4f}%")
    print(f"      Further refinement may be needed for < 1% error")

print("\n" + "=" * 70)
print("Conclusion")
print("=" * 70)
print("GridFunction.Set() interpolation:")
print(f"  - Initial error (maxh=0.02):  {results[0]['gf_rel_error']:.2f}%")
print(f"  - Final error (maxh={results[-1]['maxh']}): {results[-1]['gf_rel_error']:.2f}%")
print(f"  - Error reduction: {results[0]['gf_rel_error'] / results[-1]['gf_rel_error']:.1f}x")
print()
if min_error < 5.0:
    print("Recommendation: Use GridFunction.Set() with fine mesh (acceptable error)")
else:
    print("Recommendation: Use CoefficientFunction directly (exact, no error)")
