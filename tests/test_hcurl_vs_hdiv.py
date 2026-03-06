import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src/python'))

from ngsolve import *
from netgen.occ import *
import radia as rad
import radia_ngsolve
import numpy as np

print("="*70)
print("Comparison: HCurl vs HDiv for B projection")
print("="*70)

# Create magnet
rad.UtiDelAll()

# Set Radia to use meters (required for NGSolve integration)
rad.FldUnits('m')

magnet = rad.ObjRecMag([0, 0, 0], [0.04, 0.04, 0.06], [0, 0, 1.2])

# Create background field
def radia_field_with_A(coords):
    x, y, z = coords
    B = rad.Fld(magnet, 'b', [x, y, z])
    A = rad.Fld(magnet, 'a', [x, y, z])
    return {'B': list(B), 'A': list(A)}

bg_field = rad.ObjBckgCF(radia_field_with_A)

# Create mesh
box = Box((0.01, 0.01, 0.02), (0.06, 0.06, 0.08))
mesh = Mesh(OCCGeometry(box).GenerateMesh(maxh=0.008))

print(f"\nMesh: {mesh.ne} elements, {mesh.nv} vertices")

# Get CoefficientFunctions
A_cf = radia_ngsolve.RadiaField(bg_field, 'a')
B_cf = radia_ngsolve.RadiaField(bg_field, 'b')

# Test point
test_point = (0.030, 0.020, 0.040)  # meters
mip = mesh(*test_point)

# Direct CF evaluation (baseline)
B_cf_direct = np.array(B_cf(mip))
print(f"\nDirect CF evaluation at {test_point}:")
print(f"  B_cf(mip) = {B_cf_direct}")

# ============================================================
# Test 1: HCurl space for B (mathematically incorrect)
# ============================================================
print("\n" + "="*70)
print("Test 1: B projected to HCurl space")
print("="*70)

fes_hcurl = HCurl(mesh, order=2)
B_gf_hcurl = GridFunction(fes_hcurl)
B_gf_hcurl.Set(B_cf)

B_hcurl_eval = np.array(B_gf_hcurl(mip))
error_hcurl = np.linalg.norm(B_hcurl_eval - B_cf_direct)
rel_error_hcurl = error_hcurl / np.linalg.norm(B_cf_direct) * 100

print(f"  H(curl) DOFs: {fes_hcurl.ndof}")
print(f"  B_gf evaluation: {B_hcurl_eval}")
print(f"  Error: {error_hcurl:.6e} ({rel_error_hcurl:.2f}%)")

# ============================================================
# Test 2: HDiv space for B (mathematically correct)
# ============================================================
print("\n" + "="*70)
print("Test 2: B projected to HDiv space")
print("="*70)

fes_hdiv = HDiv(mesh, order=2)
B_gf_hdiv = GridFunction(fes_hdiv)
B_gf_hdiv.Set(B_cf)

B_hdiv_eval = np.array(B_gf_hdiv(mip))
error_hdiv = np.linalg.norm(B_hdiv_eval - B_cf_direct)
rel_error_hdiv = error_hdiv / np.linalg.norm(B_cf_direct) * 100

print(f"  H(div) DOFs: {fes_hdiv.ndof}")
print(f"  B_gf evaluation: {B_hdiv_eval}")
print(f"  Error: {error_hdiv:.6e} ({rel_error_hdiv:.2f}%)")

# ============================================================
# Comparison
# ============================================================
print("\n" + "="*70)
print("Summary")
print("="*70)
print(f"  Direct CF:      {B_cf_direct}")
print(f"  HCurl result:   {B_hcurl_eval} (error: {rel_error_hcurl:.2f}%)")
print(f"  HDiv result:    {B_hdiv_eval} (error: {rel_error_hdiv:.2f}%)")

if rel_error_hdiv < rel_error_hcurl:
    print(f"\n  HDiv is better: {rel_error_hcurl:.2f}% → {rel_error_hdiv:.2f}%")
elif rel_error_hcurl < rel_error_hdiv:
    print(f"\n  HCurl is better: {rel_error_hdiv:.2f}% → {rel_error_hcurl:.2f}%")
else:
    print(f"\n  Both have similar error: ~{rel_error_hcurl:.2f}%")
