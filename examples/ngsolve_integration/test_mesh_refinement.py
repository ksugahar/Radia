"""
Mesh Refinement Benchmark: GridFunction vs Direct Radia Evaluation

This benchmark verifies that GridFunction interpolation converges to the
exact Radia field calculation as the mesh is refined.

Setup:
- Magnet: 20×20×30 mm rectangular permanent magnet (NdFeB, Br=1.2T)
  - Center: [0, 0, 0] mm
  - Extent: x∈[-10,10], y∈[-10,10], z∈[-15,+15] mm
- Evaluation point: [0, 0, 20] mm (z = 0.02 m)
  - This point is OUTSIDE the magnet, 5 mm above the top surface
- Expected: GridFunction error → 0 as mesh size → 0
"""
import sys
import os

# Add build directory to path (relative)
script_dir = os.path.dirname(os.path.abspath(__file__))
build_dir = os.path.join(script_dir, '..', '..', 'build', 'Release')
sys.path.insert(0, build_dir)

from ngsolve import *
from netgen.csg import CSGeometry, OrthoBrick, Pnt
import radia as rad
import rad_ngsolve
import numpy as np

print("=" * 70)
print("Mesh Refinement Benchmark")
print("GridFunction.Set() Convergence to Direct Radia Evaluation")
print("=" * 70)

# Create magnet: 20×20×30 mm at origin
# Center: [0, 0, 0], Size: [20, 20, 30] mm
# => Magnet extent: x∈[-10,10], y∈[-10,10], z∈[-15,15] mm
magnet = rad.ObjRecMag([0, 0, 0], [20, 20, 30], [0, 0, 1.2])
rad.MatApl(magnet, rad.MatStd('NdFeB', 1.2))
rad.Solve(magnet, 0.0001, 10000)

print("\nMagnet Configuration:")
print("  Material: NdFeB (Br = 1.2 T)")
print("  Center: [0, 0, 0] mm")
print("  Size: [20, 20, 30] mm")
print("  Extent: z ∈ [-15, +15] mm")

# Create CoefficientFunction
B_cf = rad_ngsolve.RadBfield(magnet)

# Evaluation point: OUTSIDE the magnet, 5 mm above the top surface (z=15mm)
test_point = (0.0, 0.0, 0.02)  # m (= 20 mm)
pt_mm = [0, 0, 20]  # mm
B_radia = rad.Fld(magnet, 'b', pt_mm)

print("\nEvaluation Point:")
print(f"  Position: [0, 0, 20] mm = [0, 0, 0.02] m")
print(f"  Location: 5 mm above magnet top surface (z=15mm)")
print(f"  Radia field (reference): Bz = {B_radia[2]:.6f} T")
print()

# Test different mesh sizes
mesh_sizes = [0.02, 0.01, 0.005, 0.002, 0.001]
print("=" * 70)
print("Testing Mesh Refinement:")
print("=" * 70)
print("Strategy: Refine mesh and compare GridFunction to exact Radia value")
print()

results = []

for maxh in mesh_sizes:
	print(f"maxh = {maxh} m ({maxh*1000:.1f} mm)")
	print("-" * 70)

	# Create mesh - centered around evaluation region
	geo = CSGeometry()
	# Mesh domain includes magnet and evaluation point
	# Domain: 30×30×30 mm = x∈[-15,15], y∈[-15,15], z∈[-15,25] mm
	# (extended in +z to include evaluation point at z=20mm)
	geo.Add(OrthoBrick(Pnt(-0.015, -0.015, -0.015), Pnt(0.015, 0.015, 0.025)))
	mesh = Mesh(geo.GenerateMesh(maxh=maxh))

	print(f"  Mesh: {mesh.ne} elements, {mesh.nv} vertices")

	# Test 1: Direct CoefficientFunction evaluation (should be exact)
	mesh_pt = mesh(*test_point)
	B_cf_val = B_cf(mesh_pt)
	cf_error = abs(B_cf_val[2] - B_radia[2])

	# Test 2: GridFunction interpolation (has interpolation error)
	fes = VectorH1(mesh, order=2)
	gfB = GridFunction(fes)
	gfB.Set(B_cf)
	B_gf_val = gfB(mesh_pt)
	gf_error = abs(B_gf_val[2] - B_radia[2])
	gf_rel_error = (gf_error / abs(B_radia[2])) * 100 if B_radia[2] != 0 else 0

	print(f"  CoefficientFunction: Bz = {B_cf_val[2]:.6f} T")
	print(f"    Error vs Radia: {cf_error:.2e} T (should be ~0)")
	print(f"  GridFunction:        Bz = {B_gf_val[2]:.6f} T")
	print(f"    Error vs Radia: {gf_error:.2e} T ({gf_rel_error:.3f}%)")
	print()

	results.append({
		'maxh': maxh,
		'maxh_mm': maxh * 1000,
		'elements': mesh.ne,
		'vertices': mesh.nv,
		'dofs': fes.ndof,
		'cf_Bz': B_cf_val[2],
		'cf_error': cf_error,
		'gf_Bz': B_gf_val[2],
		'gf_error': gf_error,
		'gf_rel_error': gf_rel_error
	})

print("=" * 70)
print("Summary Table")
print("=" * 70)
print(f"{'maxh':<10} {'Elements':<10} {'DOFs':<10} {'GF Bz (T)':<12} {'GF Error':<12} {'Rel Err %':<10}")
print(f"{'(mm)':<10} {'':<10} {'':<10} {'':<12} {'(T)':<12} {'':<10}")
print("-" * 70)

for r in results:
	print(f"{r['maxh_mm']:<10.1f} {r['elements']:<10} {r['dofs']:<10} "
		  f"{r['gf_Bz']:<12.6f} {r['gf_error']:<12.2e} {r['gf_rel_error']:<10.3f}")

print("\n" + "=" * 70)
print("Convergence Analysis")
print("=" * 70)

# Check if error is monotonically decreasing
decreasing = all(results[i]['gf_error'] > results[i+1]['gf_error']
				 for i in range(len(results)-1))

if decreasing:
	print("[OK] GridFunction error decreases monotonically with mesh refinement")
	print("     This confirms proper convergence behavior.")
else:
	print("[WARNING] GridFunction error does not decrease monotonically")
	print("          Check for numerical issues or insufficient resolution")

# Calculate convergence rate (if we have enough data)
if len(results) >= 2:
	print("\nConvergence Rate:")
	for i in range(len(results)-1):
		h_ratio = results[i]['maxh'] / results[i+1]['maxh']
		error_ratio = results[i]['gf_error'] / results[i+1]['gf_error']
		# For order p, error ~ h^p, so error_ratio ~ h_ratio^p
		# p = log(error_ratio) / log(h_ratio)
		if results[i+1]['gf_error'] > 0 and error_ratio > 0:
			conv_rate = np.log(error_ratio) / np.log(h_ratio)
			print(f"  maxh {results[i]['maxh_mm']:.1f}mm → {results[i+1]['maxh_mm']:.1f}mm: "
				  f"rate = {conv_rate:.2f}")

# Error reduction
print(f"\nError Reduction:")
print(f"  Initial (maxh={results[0]['maxh_mm']:.1f}mm): {results[0]['gf_rel_error']:.3f}%")
print(f"  Final   (maxh={results[-1]['maxh_mm']:.1f}mm): {results[-1]['gf_rel_error']:.3f}%")
if results[-1]['gf_error'] > 0:
	reduction = results[0]['gf_error'] / results[-1]['gf_error']
	print(f"  Reduction factor: {reduction:.1f}×")

# Check accuracy targets
print("\n" + "=" * 70)
print("Accuracy Assessment")
print("=" * 70)

min_error = min(r['gf_rel_error'] for r in results)
min_result = [r for r in results if r['gf_rel_error'] == min_error][0]

print(f"Best accuracy achieved: {min_error:.3f}% at maxh={min_result['maxh_mm']:.1f}mm")
print(f"  ({min_result['elements']} elements, {min_result['dofs']} DOFs)")

if min_error < 0.1:
	print("  [OK] Excellent: < 0.1% error")
elif min_error < 1.0:
	print("  [OK] Good: < 1% error")
elif min_error < 5.0:
	print("  [OK] Acceptable: < 5% error")
else:
	print("  [WARNING] Poor: > 5% error - further refinement needed")

print("\n" + "=" * 70)
print("Conclusion")
print("=" * 70)

if decreasing and min_error < 1.0:
	print("[PASS] BENCHMARK PASSED")
	print("  - GridFunction interpolation converges to Radia reference")
	print("  - Error decreases monotonically with mesh refinement")
	print(f"  - Achieved {min_error:.3f}% accuracy")
	print("\nRecommendation:")
	if min_error < 0.1:
		print(f"  Use maxh <= {min_result['maxh_mm']:.1f}mm for high-precision work")
	else:
		print(f"  Use maxh <= {min_result['maxh_mm']:.1f}mm for typical applications")
else:
	print("[WARNING] BENCHMARK INCOMPLETE")
	if not decreasing:
		print("  - Non-monotonic convergence detected")
	if min_error >= 1.0:
		print(f"  - Target accuracy (< 1%) not achieved (got {min_error:.3f}%)")
	print("\nRecommendation:")
	print("  - Use finer mesh (smaller maxh)")
	print("  - Consider higher-order elements (order > 2)")
	print("  - For critical applications, use CoefficientFunction directly")

print("\n" + "=" * 70)
