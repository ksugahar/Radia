import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src/python'))

from ngsolve import *
from netgen.occ import *
import radia as rad
import radia_ngsolve
import numpy as np

print("="*70)
print("Detailed Analysis: curl(A) vs B")
print("="*70)

# Create magnet
rad.UtiDelAll()

# Set Radia to use meters (required for NGSolve integration)
rad.FldUnits('m')

magnet = rad.ObjRecMag([0, 0, 0], [0.04, 0.04, 0.06], [0, 0, 1.2])

# Background field
def radia_field_with_A(coords):
    x, y, z = coords
    B = rad.Fld(magnet, 'b', [x, y, z])
    A = rad.Fld(magnet, 'a', [x, y, z])
    return {'B': list(B), 'A': list(A)}

bg_field = rad.ObjBckgCF(radia_field_with_A)

# Create mesh
box = Box((0.01, 0.01, 0.02), (0.06, 0.06, 0.08))
mesh = Mesh(OCCGeometry(box).GenerateMesh(maxh=0.010))

print(f"\nMesh: {mesh.ne} elements, {mesh.nv} vertices")

# Finite element spaces
fes_hcurl = HCurl(mesh, order=2)
fes_hdiv = HDiv(mesh, order=2)

# CoefficientFunctions
A_cf = radia_ngsolve.RadiaField(bg_field, 'a')
B_cf = radia_ngsolve.RadiaField(bg_field, 'b')

# Project A to HCurl
A_gf = GridFunction(fes_hcurl)
A_gf.Set(A_cf)

# Compute curl(A) - result is in HDiv
curl_A_gf = curl(A_gf)
print(f"\ncurl(A_gf) type: {type(curl_A_gf)}")

# Project B to HDiv
B_gf = GridFunction(fes_hdiv)
B_gf.Set(B_cf)

# Test multiple points
test_points = [
    (0.030, 0.020, 0.030),
    (0.030, 0.020, 0.040),
    (0.030, 0.020, 0.050),
    (0.040, 0.040, 0.050),
]

print("\n" + "="*70)
print("Point-wise comparison")
print("="*70)
print(f"{'Point (m)':<25s} {'B_cf':<25s} {'curl(A)':<25s} {'B_gf':<25s} {'Error'}")
print("-"*70)

for point in test_points:
    try:
        mip = mesh(*point)
        
        B_direct = np.array(B_cf(mip))
        curl_A_val = np.array(curl_A_gf(mip))
        B_gf_val = np.array(B_gf(mip))
        
        error_curl_vs_direct = np.linalg.norm(curl_A_val - B_direct)
        error_gf_vs_direct = np.linalg.norm(B_gf_val - B_direct)
        
        print(f"{str(point):<25s} "
              f"{np.linalg.norm(B_direct):<8.5f} "
              f"{np.linalg.norm(curl_A_val):<8.5f} "
              f"{np.linalg.norm(B_gf_val):<8.5f} "
              f"{error_curl_vs_direct:.3e}")
    except Exception as e:
        print(f"{str(point):<25s} ERROR: {e}")

# L2 norm comparison
print("\n" + "="*70)
print("L2 Norm Analysis")
print("="*70)

# Method 1: Direct integration of difference
error_vec = curl_A_gf - B_gf
error_l2_1 = sqrt(Integrate(InnerProduct(error_vec, error_vec), mesh))

# Method 2: Separate norms
curl_A_l2 = sqrt(Integrate(InnerProduct(curl_A_gf, curl_A_gf), mesh))
B_gf_l2 = sqrt(Integrate(InnerProduct(B_gf, B_gf), mesh))

print(f"  ||curl(A)||_L2:        {curl_A_l2:.6e}")
print(f"  ||B_gf||_L2:           {B_gf_l2:.6e}")
print(f"  ||curl(A) - B_gf||_L2: {error_l2_1:.6e}")
print(f"  Relative error:        {error_l2_1/B_gf_l2*100:.2f}%")

# Compare curl_A_gf with B_cf directly
print("\n" + "="*70)
print("Comparing curl(A) with B_cf (CoefficientFunction)")
print("="*70)

# Sample points for comparison
n_samples = 20
errors_pointwise = []

import random
random.seed(42)

for _ in range(n_samples):
    # Random point in mesh
    x = random.uniform(0.01, 0.06)
    y = random.uniform(0.01, 0.06)
    z = random.uniform(0.02, 0.08)
    
    try:
        mip = mesh(x, y, z)
        curl_A_val = np.array(curl_A_gf(mip))
        B_cf_val = np.array(B_cf(mip))
        
        error = np.linalg.norm(curl_A_val - B_cf_val)
        rel_error = error / np.linalg.norm(B_cf_val) * 100 if np.linalg.norm(B_cf_val) > 0 else 0
        
        errors_pointwise.append(rel_error)
    except:
        pass

if errors_pointwise:
    print(f"  Mean pointwise error: {np.mean(errors_pointwise):.2f}%")
    print(f"  Std dev:              {np.std(errors_pointwise):.2f}%")
    print(f"  Max error:            {np.max(errors_pointwise):.2f}%")
    print(f"  Min error:            {np.min(errors_pointwise):.2f}%")
