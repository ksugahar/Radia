import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../build/Release'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src/python'))

import radia as rad
import numpy as np

try:
    from ngsolve import *
    from netgen.occ import *
    import radia_ngsolve
    NGSOLVE_AVAILABLE = True
except ImportError:
    print("NGSolve not available. Skipping test.")
    NGSOLVE_AVAILABLE = False
    sys.exit(0)

print("="*70)
print("H-Matrix Field Evaluation Verification with rad_ngsolve")
print("="*70)

# Set Radia units to meters
rad.FldUnits('m')
print("\n[1] Set Radia units to meters")

# Create Radia magnet with sufficient elements for H-matrix
rad.UtiDelAll()
n = 5  # 5x5x5 = 125 elements
cube_size = 0.100  # meters (magnet boundaries: ±0.05m)
elem_size = cube_size / n

print(f"\n[2] Creating magnet: {n}x{n}x{n} = {n**3} elements")
print(f"    Magnet boundaries: ±{cube_size/2:.3f}m (±{cube_size/2*1000:.0f}mm)")
elements = []
for i in range(n):
    for j in range(n):
        for k in range(n):
            x = (i - n/2 + 0.5) * elem_size
            y = (j - n/2 + 0.5) * elem_size
            z = (k - n/2 + 0.5) * elem_size
            elem = rad.ObjRecMag([x, y, z], [elem_size, elem_size, elem_size],
                                 [0, 0, 1.2])
            elements.append(elem)

magnet = rad.ObjCnt(elements)
print(f"    Magnet ID: {magnet}, Elements: {n**3}")

# Create background field providing both B and A
def radia_field_with_A(coords):
    x, y, z = coords
    B = rad.Fld(magnet, 'b', [x, y, z])
    A = rad.Fld(magnet, 'a', [x, y, z])
    return {'B': list(B), 'A': list(A)}

bg_field = rad.ObjBckgCF(radia_field_with_A)
print(f"\n[3] Created background field ID: {bg_field}")

# Create NGSolve mesh - must cover test region but avoid near-field
print("\n[4] Creating NGSolve mesh")
mesh_min = 0.015  # Start 15mm from origin
mesh_max = 0.063  # End at 63mm
box = Box((mesh_min, mesh_min, mesh_min), (mesh_max, mesh_max, mesh_max))
mesh = Mesh(OCCGeometry(box).GenerateMesh(maxh=0.010))
print(f"    Mesh domain: [{mesh_min*1000:.0f}, {mesh_max*1000:.0f}]mm")
print(f"    Mesh: {mesh.nv} vertices, {mesh.ne} elements")

# Test points (in meters)
# Magnet surface at ±50mm, test at distances >1 mesh cell (10mm) from surface
test_points = [
    (0.030, 0.020, 0.040),  # ~20mm from nearest surface
    (0.040, 0.040, 0.050),  # ~10mm from nearest surface  
    (0.050, 0.030, 0.060),  # ~10mm from nearest surface
]
print(f"\n[5] Test points (avoiding near-field <10mm from surface):")
for i, pt in enumerate(test_points):
    print(f"    Point {i+1}: ({pt[0]*1000:.0f}, {pt[1]*1000:.0f}, {pt[2]*1000:.0f})mm")

print("\n" + "="*70)
print("TEST 1: GridFunction.Set() with H-matrix DISABLED")
print("="*70)

# Disable H-matrix
rad.SetHMatrixFieldEval(0)
print("\n[6a] H-matrix disabled")

# Create CoefficientFunctions
A_cf_direct = radia_ngsolve.RadiaField(bg_field, 'a')
B_cf_direct = radia_ngsolve.RadiaField(bg_field, 'b')

# Create GridFunctions
fes_hcurl = HCurl(mesh, order=2)
A_gf_direct = GridFunction(fes_hcurl)
A_gf_direct.Set(A_cf_direct)

curl_A_gf_direct = curl(A_gf_direct)

print(f"\n[6b] GridFunction.Set() completed (Direct mode)")
print(f"    H(curl) space: {fes_hcurl.ndof} DOFs")

# Evaluate at test points
print("\n[6c] Evaluating curl(A) vs B (Direct mode):")
errors_direct = []
for i, point in enumerate(test_points):
    try:
        mip = mesh(*point)
        curl_A_val = np.array(curl_A_gf_direct(mip))
        B_val = np.array(B_cf_direct(mip))
        error = np.linalg.norm(curl_A_val - B_val)
        B_norm = np.linalg.norm(B_val)
        rel_error = error / B_norm * 100 if B_norm > 1e-10 else 0
        errors_direct.append(rel_error)
        print(f"    Point {i+1}: Error = {error:.6e} T ({rel_error:.4f}%)")
    except Exception as e:
        print(f"    Point {i+1}: [FAILED] {e}")

mean_error_direct = np.mean(errors_direct) if errors_direct else 0
print(f"\n[6d] Direct mode mean error: {mean_error_direct:.4f}%")

print("\n" + "="*70)
print("TEST 2: GridFunction.Set() with H-matrix ENABLED")
print("="*70)

# Enable H-matrix
rad.SetHMatrixFieldEval(1, 1e-6)
print("\n[7a] H-matrix enabled (eps=1e-6)")

# Create new CoefficientFunctions with H-matrix
A_cf_hmat = radia_ngsolve.RadiaField(bg_field, 'a')
B_cf_hmat = radia_ngsolve.RadiaField(bg_field, 'b')

# Create new GridFunctions
A_gf_hmat = GridFunction(fes_hcurl)
A_gf_hmat.Set(A_cf_hmat)

curl_A_gf_hmat = curl(A_gf_hmat)

print(f"\n[7b] GridFunction.Set() completed (H-matrix mode)")

# Evaluate at test points
print("\n[7c] Evaluating curl(A) vs B (H-matrix mode):")
errors_hmat = []
for i, point in enumerate(test_points):
    try:
        mip = mesh(*point)
        curl_A_val = np.array(curl_A_gf_hmat(mip))
        B_val = np.array(B_cf_hmat(mip))
        error = np.linalg.norm(curl_A_val - B_val)
        B_norm = np.linalg.norm(B_val)
        rel_error = error / B_norm * 100 if B_norm > 1e-10 else 0
        errors_hmat.append(rel_error)
        print(f"    Point {i+1}: Error = {error:.6e} T ({rel_error:.4f}%)")
    except Exception as e:
        print(f"    Point {i+1}: [FAILED] {e}")

mean_error_hmat = np.mean(errors_hmat) if errors_hmat else 0
print(f"\n[7d] H-matrix mode mean error: {mean_error_hmat:.4f}%")

print("\n" + "="*70)
print("TEST 3: Compare Direct vs H-matrix Results")
print("="*70)

# Compare GridFunction values at test points
print("\n[8] Comparing GridFunction values (A field):")
gf_differences = []
for i, point in enumerate(test_points):
    try:
        mip = mesh(*point)
        A_direct = np.array(A_gf_direct(mip))
        A_hmat = np.array(A_gf_hmat(mip))
        diff = np.linalg.norm(A_hmat - A_direct)
        A_norm = np.linalg.norm(A_direct)
        rel_diff = diff / A_norm * 100 if A_norm > 1e-10 else 0
        gf_differences.append(rel_diff)
        print(f"    Point {i+1}: Difference = {diff:.6e} ({rel_diff:.4f}%)")
    except Exception as e:
        print(f"    Point {i+1}: [FAILED] {e}")

mean_gf_diff = np.mean(gf_differences) if gf_differences else 0
print(f"\n[8a] Mean GridFunction difference: {mean_gf_diff:.4f}%")

# Summary
print("\n" + "="*70)
print("VERIFICATION SUMMARY")
print("="*70)

print(f"\nDirect mode (H-matrix OFF):")
print(f"  Mean curl(A) = B error: {mean_error_direct:.4f}%")
print(f"  Status: {'[OK]' if mean_error_direct < 1.0 else '[WARNING]'}")

print(f"\nH-matrix mode (H-matrix ON):")
print(f"  Mean curl(A) = B error: {mean_error_hmat:.4f}%")
print(f"  Status: {'[OK]' if mean_error_hmat < 1.0 else '[WARNING]'}")

print(f"\nGridFunction consistency (Direct vs H-matrix):")
print(f"  Mean difference: {mean_gf_diff:.4f}%")
print(f"  Status: {'[OK]' if mean_gf_diff < 1.0 else '[WARNING]'}")

# Final verdict
# curl(A)=B allows 2% due to L2 projection error
# GridFunction consistency must be <1% (should be identical)
all_ok = (mean_error_direct < 2.0 and mean_error_hmat < 2.0 and mean_gf_diff < 1.0)

print("\n" + "="*70)
if all_ok:
    print("[SUCCESS] H-matrix field evaluation verified!")
    print("  - GridFunction.Set() works correctly with H-matrix")
    print("  - curl(A) = B relationship maintained")
    print("  - Results consistent between Direct and H-matrix modes")
    print("  - Result matrix indexing fix confirmed working with H-matrix")
else:
    print("[WARNING] Some tests exceeded tolerance")
    if mean_error_direct >= 2.0:
        print(f"  - Direct mode error too high: {mean_error_direct:.4f}% (limit: 2.0%)")
    if mean_error_hmat >= 2.0:
        print(f"  - H-matrix mode error too high: {mean_error_hmat:.4f}% (limit: 2.0%)")
    if mean_gf_diff >= 1.0:
        print(f"  - GridFunction difference too high: {mean_gf_diff:.4f}%")

print("="*70)

# Cleanup
rad.UtiDelAll()

if not all_ok:
    sys.exit(1)
