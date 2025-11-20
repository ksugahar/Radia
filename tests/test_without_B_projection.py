import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src/python'))

from ngsolve import *
from netgen.occ import *
import radia as rad
import rad_ngsolve
import numpy as np

print("="*70)
print("Verification: curl(A_gf) vs B_cf (no B projection)")
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

# Test different mesh sizes
mesh_sizes = [0.012, 0.010, 0.008, 0.006]

results = []

for maxh in mesh_sizes:
    print(f"\nMesh size: {maxh*1000:.1f} mm")
    print("-"*70)
    
    # Create mesh
    box = Box((0.01, 0.01, 0.02), (0.06, 0.06, 0.08))
    mesh = Mesh(OCCGeometry(box).GenerateMesh(maxh=maxh))
    
    print(f"  Elements: {mesh.ne}, Vertices: {mesh.nv}")
    
    # A in HCurl
    fes_hcurl = HCurl(mesh, order=2)
    A_cf = rad_ngsolve.RadiaField(bg_field, 'a')
    A_gf = GridFunction(fes_hcurl)
    A_gf.Set(A_cf)
    
    # curl(A)
    curl_A_gf = curl(A_gf)
    
    # B as CoefficientFunction (NO projection!)
    B_cf = rad_ngsolve.RadiaField(bg_field, 'b')
    
    # Sample random points
    n_samples = 50
    errors = []
    
    import random
    random.seed(42)
    
    for _ in range(n_samples):
        x = random.uniform(0.015, 0.055)
        y = random.uniform(0.015, 0.055)
        z = random.uniform(0.025, 0.075)
        
        try:
            mip = mesh(x, y, z)
            curl_A_val = np.array(curl_A_gf(mip))
            B_val = np.array(B_cf(mip))
            
            error = np.linalg.norm(curl_A_val - B_val)
            B_norm = np.linalg.norm(B_val)
            
            if B_norm > 1e-10:
                rel_error = error / B_norm * 100
                errors.append(rel_error)
        except:
            pass
    
    if errors:
        mean_error = np.mean(errors)
        std_error = np.std(errors)
        max_error = np.max(errors)
        min_error = np.min(errors)
        
        print(f"  H(curl) DOFs: {fes_hcurl.ndof}")
        print(f"  Sampled points: {len(errors)}")
        print(f"  Mean error:  {mean_error:.2f}%")
        print(f"  Std dev:     {std_error:.2f}%")
        print(f"  Min error:   {min_error:.2f}%")
        print(f"  Max error:   {max_error:.2f}%")
        
        results.append({
            'maxh': maxh,
            'h_mm': maxh * 1000,
            'ne': mesh.ne,
            'ndof': fes_hcurl.ndof,
            'mean_error': mean_error,
            'max_error': max_error
        })

print("\n" + "="*70)
print("Convergence Summary")
print("="*70)
print(f"{'h (mm)':>8s} {'Elements':>10s} {'H(curl) DOFs':>12s} {'Mean Error':>12s} {'Max Error':>12s}")
print("-"*70)

for r in results:
    print(f"{r['h_mm']:8.1f} {r['ne']:10d} {r['ndof']:12d} {r['mean_error']:11.2f}% {r['max_error']:11.2f}%")

if len(results) >= 2:
    print("\nConvergence trend:")
    for i in range(1, len(results)):
        h1 = results[i-1]['h_mm']
        h2 = results[i]['h_mm']
        e1 = results[i-1]['mean_error']
        e2 = results[i]['mean_error']
        improvement = (e1 - e2) / e1 * 100 if e1 > 0 else 0
        print(f"  h: {h1:.1f}mm → {h2:.1f}mm: error: {e1:.2f}% → {e2:.2f}% "
              f"(improvement: {improvement:.1f}%)")

print("\nConclusion:")
finest = results[-1]
if finest['mean_error'] < 5.0:
    print(f"  [SUCCESS] curl(A) = B verified with {finest['mean_error']:.2f}% mean error")
else:
    print(f"  [INFO] Mean error: {finest['mean_error']:.2f}%")
    print(f"  Further mesh refinement may improve accuracy")
