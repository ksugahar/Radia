import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src/python'))

from ngsolve import *
from netgen.occ import *
import radia as rad
import rad_ngsolve
import numpy as np

print("="*70)
print("Compare ALL Vector FE Spaces for B Projection")
print("="*70)

# Create magnet
rad.UtiDelAll()

# Set Radia to use meters (required for NGSolve integration)
rad.FldUnits('m')

magnet = rad.ObjRecMag([0, 0, 0], [0.04, 0.04, 0.06], [0, 0, 1.2])

def radia_field_with_A(coords):
    x, y, z = coords
    B = rad.Fld(magnet, 'b', [x, y, z])
    A = rad.Fld(magnet, 'a', [x, y, z])
    return {'B': list(B), 'A': list(A)}

bg_field = rad.ObjBckgCF(radia_field_with_A)

# Create mesh
box = Box((0.01, 0.01, 0.02), (0.06, 0.06, 0.08))
mesh = Mesh(OCCGeometry(box).GenerateMesh(maxh=0.010))

print(f"\nMesh: {mesh.ne} elements, {mesh.nv} vertices\n")

B_cf = rad_ngsolve.RadiaField(bg_field, 'b')

# Test all vector spaces
spaces = [
    ("HCurl, order=1", HCurl(mesh, order=1)),
    ("HCurl, order=2", HCurl(mesh, order=2)),
    ("HDiv, order=1", HDiv(mesh, order=1)),
    ("HDiv, order=2", HDiv(mesh, order=2)),
    ("VectorH1, order=1", VectorH1(mesh, order=1)),
    ("VectorH1, order=2", VectorH1(mesh, order=2)),
]

# Test points
test_points = [
    (0.030, 0.020, 0.030),
    (0.030, 0.020, 0.040),
    (0.030, 0.020, 0.050),
    (0.040, 0.040, 0.050),
]

results = []

for name, fes in spaces:
    print(f"{name} (DOFs: {fes.ndof})")
    print("-"*70)
    
    B_gf = GridFunction(fes)
    B_gf.Set(B_cf)
    
    # Point-wise errors
    point_errors = []
    for point in test_points:
        try:
            mip = mesh(*point)
            B_direct = np.array(B_cf(mip))
            B_gf_val = np.array(B_gf(mip))
            
            error = np.linalg.norm(B_gf_val - B_direct)
            B_norm = np.linalg.norm(B_direct)
            rel_error = error / B_norm * 100 if B_norm > 0 else 0
            
            point_errors.append(rel_error)
        except:
            pass
    
    # L2 error
    error_vec = B_gf - B_cf
    error_l2_sq = Integrate(InnerProduct(error_vec, error_vec), mesh)
    error_l2 = np.sqrt(error_l2_sq)
    
    B_l2_sq = Integrate(InnerProduct(B_cf, B_cf), mesh)
    B_l2 = np.sqrt(B_l2_sq)
    
    rel_error_l2 = error_l2 / B_l2 * 100 if B_l2 > 0 else 0
    
    # Sample 100 random points for spatial distribution
    import random
    random.seed(42)
    spatial_errors = []
    
    for _ in range(100):
        x = random.uniform(0.015, 0.055)
        y = random.uniform(0.015, 0.055)
        z = random.uniform(0.025, 0.075)
        
        try:
            mip = mesh(x, y, z)
            B_direct = np.array(B_cf(mip))
            B_gf_val = np.array(B_gf(mip))
            
            error = np.linalg.norm(B_gf_val - B_direct)
            B_norm = np.linalg.norm(B_direct)
            
            if B_norm > 0:
                spatial_errors.append(error / B_norm * 100)
        except:
            pass
    
    mean_point = np.mean(point_errors) if point_errors else 0
    max_point = np.max(point_errors) if point_errors else 0
    
    mean_spatial = np.mean(spatial_errors) if spatial_errors else 0
    std_spatial = np.std(spatial_errors) if spatial_errors else 0
    max_spatial = np.max(spatial_errors) if spatial_errors else 0
    
    print(f"  Point-wise (4 points):    Mean={mean_point:.2f}%, Max={max_point:.2f}%")
    print(f"  Random samples (100):     Mean={mean_spatial:.2f}%, Std={std_spatial:.2f}%, Max={max_spatial:.2f}%")
    print(f"  L2 norm error:            {rel_error_l2:.2f}%")
    print()
    
    results.append({
        'name': name,
        'ndof': fes.ndof,
        'mean_point': mean_point,
        'mean_spatial': mean_spatial,
        'std_spatial': std_spatial,
        'max_spatial': max_spatial,
        'l2_error': rel_error_l2
    })

# Summary
print("="*70)
print("Summary (sorted by L2 error)")
print("="*70)
print(f"{'Space':<20s} {'DOFs':>8s} {'Mean pt%':>9s} {'Mean sp%':>9s} {'Max sp%':>9s} {'L2 %':>8s}")
print("-"*70)

results_sorted = sorted(results, key=lambda x: x['l2_error'])
for r in results_sorted:
    print(f"{r['name']:<20s} {r['ndof']:8d} {r['mean_point']:8.2f}% {r['mean_spatial']:8.2f}% "
          f"{r['max_spatial']:8.2f}% {r['l2_error']:7.2f}%")

best = results_sorted[0]
print(f"\nBest for GridFunction use: {best['name']}")
print(f"  L2 error: {best['l2_error']:.2f}%")
print(f"  Mean spatial error: {best['mean_spatial']:.2f}%")
