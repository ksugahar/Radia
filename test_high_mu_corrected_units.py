#!/usr/bin/env python
"""
Test ANALYTICAL tetrahedral method with CORRECTED background field units

Previous tests used rad.ObjBckg([0, 0, H_ext]) where H_ext was in A/m.
This was WRONG - rad.ObjBckg() expects Tesla!

Correct usage: rad.ObjBckg([0, 0, H_ext * mu0]) where mu0 = 4*pi*1e-7
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'build/Release'))

os.environ['RADIA_TETRA_METHOD'] = 'ANALYTICAL'

import numpy as np
import radia as rad

print("="*70)
print("High Permeability Test: CORRECTED Units")
print("="*70)

# Test parameters
cube_size = 0.1  # m
H_ext = 1000.0   # A/m (z-direction)
mu0 = 4*np.pi*1e-7  # H/m
B_ext = H_ext * mu0  # Convert to Tesla (CORRECT)
eval_pt = [0.05, 0.05, 0.2]  # Outside cube

print(f"\nField specification:")
print(f"  H_ext = {H_ext} A/m (desired field)")
print(f"  mu0 = {mu0:.6e} H/m")
print(f"  B_ext = H_ext * mu0 = {B_ext:.6e} T (input to Radia)")

# Test different permeability values
mu_r_values = [100, 1000, 4000, 10000]

results = []

for mu_r in mu_r_values:
    print(f"\n{'='*70}")
    print(f"Testing mu_r = {mu_r}")
    print(f"{'='*70}")

    ksi = mu_r - 1

    # ============================================================
    # Reference: Subdivided Hexahedron (STANDARD method)
    # ============================================================
    print("\n[1] STANDARD Method (Hexahedral 5x5x5 mesh)")

    rad.FldUnits('m')

    cube_hex = rad.ObjRecMag([0, 0, 0], [cube_size, cube_size, cube_size], [0, 0, 0])
    rad.ObjDivMag(cube_hex, [5, 5, 5])

    mat_hex = rad.MatLin([ksi, ksi], [0, 0, 1])
    rad.MatApl(cube_hex, mat_hex)

    # CORRECTED: Use Tesla for background field
    bg_hex = rad.ObjBckg([0, 0, B_ext])
    container_hex = rad.ObjCnt([cube_hex, bg_hex])

    print(f"  Solving with mu_r={mu_r}, ksi={ksi}...")
    result_hex = rad.Solve(container_hex, 0.0001, 10000)
    print(f"  Solve result: {result_hex}")

    try:
        H_hex = rad.Fld(container_hex, 'h', eval_pt)
        H_hex_mag = np.linalg.norm(H_hex)
        print(f"  H = [{H_hex[0]:.4f}, {H_hex[1]:.4f}, {H_hex[2]:.4f}] A/m")
        print(f"  |H| = {H_hex_mag:.4f} A/m")
        H_hex_success = True
    except Exception as e:
        print(f"  ERROR: {e}")
        H_hex = [np.nan, np.nan, np.nan]
        H_hex_mag = np.nan
        H_hex_success = False

    # ============================================================
    # Test: Single Tetrahedron (ANALYTICAL method)
    # ============================================================
    print(f"\n[2] ANALYTICAL Method (Single tetrahedron)")

    rad.FldUnits('m')

    vertices = [
        [0, 0, 0],
        [cube_size, 0, 0],
        [0, cube_size, 0],
        [0, 0, cube_size]
    ]

    TETRA_FACES = [
        [1, 3, 2],
        [1, 2, 4],
        [1, 4, 3],
        [2, 3, 4]
    ]

    tetra = rad.ObjPolyhdr(vertices, TETRA_FACES, [0, 0, 0])
    mat_tetra = rad.MatLin([ksi, ksi], [0, 0, 1])
    rad.MatApl(tetra, mat_tetra)

    # CORRECTED: Use Tesla for background field
    bg_tetra = rad.ObjBckg([0, 0, B_ext])
    container_tetra = rad.ObjCnt([tetra, bg_tetra])

    print(f"  Solving with mu_r={mu_r}, ksi={ksi}...")
    result_tetra = rad.Solve(container_tetra, 0.0001, 10000)
    print(f"  Solve result: {result_tetra}")

    try:
        H_tetra = rad.Fld(container_tetra, 'h', eval_pt)
        H_tetra_mag = np.linalg.norm(H_tetra)
        print(f"  H = [{H_tetra[0]:.4f}, {H_tetra[1]:.4f}, {H_tetra[2]:.4f}] A/m")
        print(f"  |H| = {H_tetra_mag:.4f} A/m")
        H_tetra_success = True
    except Exception as e:
        print(f"  ERROR: {e}")
        H_tetra = [np.nan, np.nan, np.nan]
        H_tetra_mag = np.nan
        H_tetra_success = False

    # ============================================================
    # Comparison
    # ============================================================
    print(f"\n{'='*70}")
    print(f"Comparison for mu_r = {mu_r}")
    print(f"{'='*70}")

    if H_hex_success and H_tetra_success:
        error_mag = abs(H_tetra_mag - H_hex_mag) / H_hex_mag * 100

        print(f"  STANDARD (hex 5x5x5): |H| = {H_hex_mag:.4f} A/m")
        print(f"  ANALYTICAL (tetra):   |H| = {H_tetra_mag:.4f} A/m")
        print(f"\n  Error (magnitude): {error_mag:.2f}%")

        if error_mag < 20:
            status = "PASS"
        elif error_mag < 50:
            status = "WARNING"
        else:
            status = "FAIL"

        print(f"  Status: {status} (target: <20%)")

        results.append({
            'mu_r': mu_r,
            'H_hex_mag': H_hex_mag,
            'H_tetra_mag': H_tetra_mag,
            'error_mag': error_mag,
            'status': status
        })
    else:
        print(f"  FAILED: One or both methods failed")
        results.append({
            'mu_r': mu_r,
            'H_hex_mag': np.nan,
            'H_tetra_mag': np.nan,
            'error_mag': np.nan,
            'status': 'FAILED'
        })

# ============================================================
# Summary
# ============================================================
print("\n" + "="*70)
print("Summary: High Permeability with CORRECTED Units")
print("="*70)

print("\n{:<10} {:<15} {:<15} {:<12} {:<10}".format(
    "mu_r", "Hex |H| (A/m)", "Tetra |H| (A/m)", "Error (%)", "Status"))
print("-"*70)

for r in results:
    if not np.isnan(r['H_hex_mag']) and not np.isnan(r['H_tetra_mag']):
        print("{:<10} {:<15.4f} {:<15.4f} {:<12.2f} {:<10}".format(
            r['mu_r'], r['H_hex_mag'], r['H_tetra_mag'], r['error_mag'], r['status']))
    else:
        print("{:<10} {:<15} {:<15} {:<12} {:<10}".format(
            r['mu_r'], "N/A", "N/A", "N/A", r['status']))

print("\n" + "="*70)
print("Conclusion:")
print("="*70)

pass_count = sum(1 for r in results if r['status'] == 'PASS')
warning_count = sum(1 for r in results if r['status'] == 'WARNING')
fail_count = sum(1 for r in results if r['status'] in ['FAIL', 'FAILED'])

if pass_count == len(results):
    print("  SUCCESS: All tests passed with corrected units!")
    print("  ANALYTICAL method is stable for all tested mu_r values.")
elif pass_count >= len(results) * 0.75:
    print(f"  PARTIAL SUCCESS: {pass_count}/{len(results)} tests passed")
    print("  ANALYTICAL method works for most mu_r values.")
else:
    print(f"  ISSUE PERSISTS: Only {pass_count}/{len(results)} tests passed")
    print("  Units correction did not fully resolve the problem.")

print("\n" + "="*70)
