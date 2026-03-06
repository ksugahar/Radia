import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src/python'))

from ngsolve import *
from netgen.occ import *
import radia as rad
import radia_ngsolve
import numpy as np

print("="*70)
print("HDiv Space: h and p Convergence Study")
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

# h-convergence: different mesh sizes, fixed order=2
print("\n" + "="*70)
print("h-Convergence (order=2, varying mesh size)")
print("="*70)

mesh_sizes = [0.012, 0.010, 0.008, 0.006, 0.005]
h_results = []

for maxh in mesh_sizes:
    box = Box((0.01, 0.01, 0.02), (0.06, 0.06, 0.08))
    mesh = Mesh(OCCGeometry(box).GenerateMesh(maxh=maxh))
    
    fes = HDiv(mesh, order=2)
    B_cf = radia_ngsolve.RadiaField(bg_field, 'b')
    B_gf = GridFunction(fes)
    B_gf.Set(B_cf)
    
    # L2 error
    error_vec = B_gf - B_cf
    error_l2 = sqrt(Integrate(InnerProduct(error_vec, error_vec), mesh))
    B_l2 = sqrt(Integrate(InnerProduct(B_cf, B_cf), mesh))
    rel_error = error_l2 / B_l2 * 100
    
    # Sample errors
    import random
    random.seed(42)
    samples = []
    for _ in range(50):
        x = random.uniform(0.02, 0.05)
        y = random.uniform(0.02, 0.05)
        z = random.uniform(0.03, 0.07)
        try:
            mip = mesh(x, y, z)
            B_direct = np.array(B_cf(mip))
            B_gf_val = np.array(B_gf(mip))
            err = np.linalg.norm(B_gf_val - B_direct) / np.linalg.norm(B_direct) * 100
            samples.append(err)
        except:
            pass
    
    mean_sample = np.mean(samples) if samples else 0
    
    h_results.append({
        'h': maxh * 1000,
        'ne': mesh.ne,
        'ndof': fes.ndof,
        'l2_err': rel_error,
        'mean_sample': mean_sample
    })
    
    print(f"  h={maxh*1000:5.1f}mm, ne={mesh.ne:5d}, DOFs={fes.ndof:6d}: "
          f"L2={rel_error:6.2f}%, Mean={mean_sample:5.2f}%")

# p-convergence: fixed mesh, varying order
print("\n" + "="*70)
print("p-Convergence (h=8mm, varying order)")
print("="*70)

box = Box((0.01, 0.01, 0.02), (0.06, 0.06, 0.08))
mesh = Mesh(OCCGeometry(box).GenerateMesh(maxh=0.008))

orders = [1, 2, 3, 4]
p_results = []

for order in orders:
    fes = HDiv(mesh, order=order)
    B_cf = radia_ngsolve.RadiaField(bg_field, 'b')
    B_gf = GridFunction(fes)
    B_gf.Set(B_cf)
    
    error_vec = B_gf - B_cf
    error_l2 = sqrt(Integrate(InnerProduct(error_vec, error_vec), mesh))
    B_l2 = sqrt(Integrate(InnerProduct(B_cf, B_cf), mesh))
    rel_error = error_l2 / B_l2 * 100
    
    # Sample errors
    import random
    random.seed(42)
    samples = []
    for _ in range(50):
        x = random.uniform(0.02, 0.05)
        y = random.uniform(0.02, 0.05)
        z = random.uniform(0.03, 0.07)
        try:
            mip = mesh(x, y, z)
            B_direct = np.array(B_cf(mip))
            B_gf_val = np.array(B_gf(mip))
            err = np.linalg.norm(B_gf_val - B_direct) / np.linalg.norm(B_direct) * 100
            samples.append(err)
        except:
            pass
    
    mean_sample = np.mean(samples) if samples else 0
    
    p_results.append({
        'order': order,
        'ndof': fes.ndof,
        'l2_err': rel_error,
        'mean_sample': mean_sample
    })
    
    print(f"  order={order}, DOFs={fes.ndof:6d}: L2={rel_error:6.2f}%, Mean={mean_sample:5.2f}%")

# Summary
print("\n" + "="*70)
print("Summary and Recommendations")
print("="*70)

best_h = min(h_results, key=lambda x: x['l2_err'])
best_p = min(p_results, key=lambda x: x['l2_err'])

print(f"\nBest h-refinement: h={best_h['h']:.1f}mm, L2 error={best_h['l2_err']:.2f}%")
print(f"Best p-refinement: order={best_p['order']}, L2 error={best_p['l2_err']:.2f}%")

print("\nFor practical NGSolve use:")
if best_h['l2_err'] < 10:
    print(f"  ✓ Usable with h={best_h['h']:.1f}mm mesh ({best_h['l2_err']:.2f}% error)")
elif best_p['l2_err'] < 10:
    print(f"  ✓ Usable with order={best_p['order']} elements ({best_p['l2_err']:.2f}% error)")
else:
    print(f"  ✗ L2 error still > 10% even with refinement")
    print(f"    Current best: {min(best_h['l2_err'], best_p['l2_err']):.2f}%")
    print(f"    Consider using B_cf directly instead of GridFunction")
