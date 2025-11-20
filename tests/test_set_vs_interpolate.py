import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src/python'))

from ngsolve import *
from netgen.occ import *
import radia as rad
import rad_ngsolve
import numpy as np

print("="*70)
print("GridFunction Projection: Set() vs Interpolate() for HDiv")
print("="*70)

# Create magnet
rad.UtiDelAll()
magnet = rad.ObjRecMag([0, 0, 0], [40, 40, 60], [0, 0, 1.2])

def radia_field_with_A(coords):
    x, y, z = coords
    B = rad.Fld(magnet, 'b', [x, y, z])
    A = rad.Fld(magnet, 'a', [x, y, z])
    return {'B': list(B), 'A': list(A)}

bg_field = rad.ObjBckgCF(radia_field_with_A)

# Create mesh
box = Box((0.01, 0.01, 0.02), (0.06, 0.06, 0.08))
mesh = Mesh(OCCGeometry(box).GenerateMesh(maxh=0.008))

print(f"\nMesh: {mesh.ne} elements, {mesh.nv} vertices\n")

# Get CoefficientFunctions
A_cf = rad_ngsolve.RadiaField(bg_field, 'a')
B_cf = rad_ngsolve.RadiaField(bg_field, 'b')

# Test different approaches
test_configs = [
    ("HCurl, order=1, Set", HCurl(mesh, order=1), "set"),
    ("HCurl, order=2, Set", HCurl(mesh, order=2), "set"),
    ("HDiv, order=1, Set", HDiv(mesh, order=1), "set"),
    ("HDiv, order=2, Set", HDiv(mesh, order=2), "set"),
]

# Test points
test_points = [
    (0.030, 0.020, 0.030),
    (0.030, 0.020, 0.040),
    (0.030, 0.020, 0.050),
    (0.040, 0.040, 0.050),
]

results = []

for name, fes, method in test_configs:
    print(f"{name} (DOFs: {fes.ndof})")
    print("-"*70)
    
    B_gf = GridFunction(fes)
    
    if method == "set":
        B_gf.Set(B_cf)
    else:
        # Interpolate not available for HDiv/HCurl
        print("  Interpolate not supported for this space")
        continue
    
    # Point-wise errors
    errors = []
    for point in test_points:
        try:
            mip = mesh(*point)
            B_direct = np.array(B_cf(mip))
            B_gf_val = np.array(B_gf(mip))
            
            error = np.linalg.norm(B_gf_val - B_direct)
            B_norm = np.linalg.norm(B_direct)
            rel_error = error / B_norm * 100 if B_norm > 0 else 0
            
            errors.append(rel_error)
        except:
            pass
    
    if errors:
        mean_err = np.mean(errors)
        max_err = np.max(errors)
        print(f"  Mean pointwise error: {mean_err:.2f}%")
        print(f"  Max pointwise error:  {max_err:.2f}%")
        
        # L2 error
        error_vec = B_gf - B_cf
        error_l2_sq = Integrate(InnerProduct(error_vec, error_vec), mesh)
        error_l2 = np.sqrt(error_l2_sq)
        
        B_l2_sq = Integrate(InnerProduct(B_cf, B_cf), mesh)
        B_l2 = np.sqrt(B_l2_sq)
        
        rel_error_l2 = error_l2 / B_l2 * 100 if B_l2 > 0 else 0
        
        print(f"  L2 norm error:        {rel_error_l2:.2f}%")
        
        results.append({
            'name': name,
            'ndof': fes.ndof,
            'mean_err': mean_err,
            'max_err': max_err,
            'l2_err': rel_error_l2
        })
    
    print()

# Summary
print("="*70)
print("Summary")
print("="*70)
print(f"{'Configuration':<30s} {'DOFs':>8s} {'Mean %':>8s} {'Max %':>8s} {'L2 %':>8s}")
print("-"*70)
for r in results:
    print(f"{r['name']:<30s} {r['ndof']:8d} {r['mean_err']:7.2f}% {r['max_err']:7.2f}% {r['l2_err']:7.2f}%")

# Find best configuration
best_mean = min(results, key=lambda x: x['mean_err'])
best_l2 = min(results, key=lambda x: x['l2_err'])

print(f"\nBest mean pointwise error: {best_mean['name']} ({best_mean['mean_err']:.2f}%)")
print(f"Best L2 error:             {best_l2['name']} ({best_l2['l2_err']:.2f}%)")
