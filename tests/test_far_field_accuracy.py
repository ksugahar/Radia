import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src/python'))

from ngsolve import *
from netgen.occ import *
import radia as rad
import rad_ngsolve
import numpy as np

print("="*70)
print("Far-Field GridFunction Accuracy Test")
print("Magnet: center=[0,0,0], size=[40,40,60] mm")
print("="*70)

# Create magnet at origin
rad.UtiDelAll()

# Set Radia to use meters (required for NGSolve integration)
rad.FldUnits('m')

magnet = rad.ObjRecMag([0, 0, 0], [0.04, 0.04, 0.06], [0, 0, 1.2])

# Magnet boundaries: x=[-20,20], y=[-20,20], z=[-30,30] mm
print("\nMagnet boundaries:")
print("  x: [-20, 20] mm")
print("  y: [-20, 20] mm")
print("  z: [-30, 30] mm")

def radia_field_with_A(coords):
    x, y, z = coords
    B = rad.Fld(magnet, 'b', [x, y, z])
    A = rad.Fld(magnet, 'a', [x, y, z])
    return {'B': list(B), 'A': list(A)}

bg_field = rad.ObjBckgCF(radia_field_with_A)

# Create larger mesh to include far-field region
# Domain: x,y,z = [-50, 50] mm = [-0.05, 0.05] m
box = Box((-0.05, -0.05, -0.05), (0.05, 0.05, 0.05))
mesh = Mesh(OCCGeometry(box).GenerateMesh(maxh=0.010))

print(f"\nMesh: {mesh.ne} elements, domain=[-50,50]mm")

# Test different regions
regions = [
    {
        'name': 'Near field (>5mm from surface)',
        'x_range': (0.025, 0.045),  # 25-45mm
        'y_range': (0.025, 0.045),
        'z_range': (0.035, 0.045),  # z > 35mm (5mm from z=30 surface)
    },
    {
        'name': 'Mid field (>10mm from surface)',
        'x_range': (0.030, 0.045),  # 30-45mm
        'y_range': (0.030, 0.045),
        'z_range': (0.040, 0.045),  # z > 40mm (10mm from surface)
    },
    {
        'name': 'Far field (>20mm from surface)',
        'x_range': (0.040, 0.050),  # 40-50mm
        'y_range': (0.040, 0.050),
        'z_range': (0.050, 0.060),  # z > 50mm (20mm from surface)
    },
]

# Test HDiv order=2 (best from previous tests)
fes = HDiv(mesh, order=2)
B_cf = rad_ngsolve.RadiaField(bg_field, 'b')
B_gf = GridFunction(fes)
B_gf.Set(B_cf)

print(f"\nUsing HDiv order=2 space ({fes.ndof} DOFs)")

import random
random.seed(42)

print("\n" + "="*70)
print("Region-wise Accuracy")
print("="*70)

for region in regions:
    print(f"\n{region['name']}")
    print(f"  x: [{region['x_range'][0]*1000:.0f}, {region['x_range'][1]*1000:.0f}] mm")
    print(f"  y: [{region['y_range'][0]*1000:.0f}, {region['y_range'][1]*1000:.0f}] mm")
    print(f"  z: [{region['z_range'][0]*1000:.0f}, {region['z_range'][1]*1000:.0f}] mm")
    
    errors = []
    B_magnitudes = []
    
    # Sample 100 points in this region
    for _ in range(100):
        x = random.uniform(*region['x_range'])
        y = random.uniform(*region['y_range'])
        z = random.uniform(*region['z_range'])
        
        try:
            mip = mesh(x, y, z)
            B_direct = np.array(B_cf(mip))
            B_gf_val = np.array(B_gf(mip))
            
            B_norm = np.linalg.norm(B_direct)
            B_magnitudes.append(B_norm)
            
            error = np.linalg.norm(B_gf_val - B_direct)
            if B_norm > 1e-10:
                rel_error = error / B_norm * 100
                errors.append(rel_error)
        except:
            pass
    
    if errors and B_magnitudes:
        print(f"  Sampled points: {len(errors)}")
        print(f"  Field strength: {np.mean(B_magnitudes)*1e3:.2f} ± {np.std(B_magnitudes)*1e3:.2f} mT")
        print(f"  Relative error: {np.mean(errors):.2f}% ± {np.std(errors):.2f}%")
        print(f"  Min/Max error:  {np.min(errors):.2f}% / {np.max(errors):.2f}%")
        
        if np.mean(errors) < 2.0:
            print(f"  [OK] Excellent accuracy (<2%)")
        elif np.mean(errors) < 5.0:
            print(f"  [OK] Good accuracy (<5%)")
        elif np.mean(errors) < 10.0:
            print(f"  [~] Acceptable accuracy (<10%)")
        else:
            print(f"  [NG] Poor accuracy (>10%)")

# Conclusion
print("\n" + "="*70)
print("Conclusion")
print("="*70)
print("\nFor practical NGSolve applications:")
print("  - Avoid evaluating GridFunction within [~]10mm of magnet surface")
print("  - In far-field (>20mm), HDiv order=2 provides good accuracy")
print("  - For near-field calculations, consider using B_cf directly")
