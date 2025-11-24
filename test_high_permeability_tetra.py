#!/usr/bin/env python
"""
Test ANALYTICAL tetrahedral method with high permeability materials

This script tests the stability and accuracy of the ANALYTICAL method
for tetrahedral meshes with various permeability values:
- mu_r = 100 (baseline, known to work)
- mu_r = 1000 (moderate)
- mu_r = 4000 (high - known to have issues)
- mu_r = 10000 (very high)

Comparison: Single tetrahedron vs subdivided hexahedron (5x5x5)
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'build/Release'))

# Enable ANALYTICAL method
os.environ['RADIA_TETRA_METHOD'] = 'ANALYTICAL'

import numpy as np
import radia as rad

print("="*70)
print("High Permeability Test: ANALYTICAL Tetrahedral Method")
print("="*70)

# Test parameters
cube_size = 0.1  # m
H_ext = 1000.0   # A/m (z-direction)
eval_pt = [0.05, 0.05, 0.2]  # Outside cube

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

    # Create cube with zero magnetization (for MatLin)
    cube_hex = rad.ObjRecMag([0, 0, 0], [cube_size, cube_size, cube_size], [0, 0, 0])

    # Subdivide for better accuracy
    subdivisions = [5, 5, 5]
    rad.ObjDivMag(cube_hex, subdivisions)

    # Apply linear material
    mat_hex = rad.MatLin([ksi, ksi], [0, 0, 1])
    rad.MatApl(cube_hex, mat_hex)

    # Apply background field
    bg_hex = rad.ObjBckg([0, 0, H_ext])
    container_hex = rad.ObjCnt([cube_hex, bg_hex])

    # Solve
    print(f"  Solving with mu_r={mu_r}, ksi={ksi}...")
    result_hex = rad.Solve(container_hex, 0.0001, 10000)
    print(f"  Solve result: {result_hex}")

    # Evaluate field
    try:
        H_hex = rad.Fld(container_hex, 'h', eval_pt)
        H_hex_mag = np.linalg.norm(H_hex)
        print(f"  H = [{H_hex[0]:.4f}, {H_hex[1]:.4f}, {H_hex[2]:.4f}] A/m")
        print(f"  |H| = {H_hex_mag:.4f} A/m")
        H_hex_success = True
    except Exception as e:
        print(f"  ERROR evaluating field: {e}")
        H_hex = [np.nan, np.nan, np.nan]
        H_hex_mag = np.nan
        H_hex_success = False

    # ============================================================
    # Test: Single Tetrahedron (ANALYTICAL method)
    # ============================================================
    print(f"\n[2] ANALYTICAL Method (Single tetrahedron)")

    rad.FldUnits('m')

    # Define single tetrahedron vertices
    vertices = [
        [0, 0, 0],                     # Vertex 1
        [cube_size, 0, 0],              # Vertex 2
        [0, cube_size, 0],              # Vertex 3
        [0, 0, cube_size]               # Vertex 4
    ]

    # Face definitions (1-indexed, counter-clockwise from outside)
    TETRA_FACES = [
        [1, 3, 2],  # Base face (z=0)
        [1, 2, 4],  # Front face (y=0)
        [1, 4, 3],  # Left face (x=0)
        [2, 3, 4]   # Slanted face
    ]

    # Create tetrahedron
    tetra = rad.ObjPolyhdr(vertices, TETRA_FACES, [0, 0, 0])

    # Apply linear material
    mat_tetra = rad.MatLin([ksi, ksi], [0, 0, 1])
    rad.MatApl(tetra, mat_tetra)

    # Apply background field
    bg_tetra = rad.ObjBckg([0, 0, H_ext])
    container_tetra = rad.ObjCnt([tetra, bg_tetra])

    # Solve
    print(f"  Solving with mu_r={mu_r}, ksi={ksi}...")
    result_tetra = rad.Solve(container_tetra, 0.0001, 10000)
    print(f"  Solve result: {result_tetra}")

    # Evaluate field
    try:
        H_tetra = rad.Fld(container_tetra, 'h', eval_pt)
        H_tetra_mag = np.linalg.norm(H_tetra)
        print(f"  H = [{H_tetra[0]:.4f}, {H_tetra[1]:.4f}, {H_tetra[2]:.4f}] A/m")
        print(f"  |H| = {H_tetra_mag:.4f} A/m")
        H_tetra_success = True
    except Exception as e:
        print(f"  ERROR evaluating field: {e}")
        H_tetra = [np.nan, np.nan, np.nan]
        H_tetra_mag = np.nan
        H_tetra_success = False

    # ============================================================
    # Comparison
    # ============================================================
    print(f"\n{'='*70}")
    print(f"Comparison for mu_r = {mu_r}")
    print(f"{'='*70}")

    if H_hex_success and H_tetra_success and not np.isnan(H_hex_mag) and not np.isnan(H_tetra_mag):
        error_mag = abs(H_tetra_mag - H_hex_mag) / H_hex_mag * 100
        error_x = abs(H_tetra[0] - H_hex[0]) / (abs(H_hex[0]) + 1e-10) * 100
        error_y = abs(H_tetra[1] - H_hex[1]) / (abs(H_hex[1]) + 1e-10) * 100
        error_z = abs(H_tetra[2] - H_hex[2]) / (abs(H_hex[2]) + 1e-10) * 100

        print(f"  STANDARD (hex 5x5x5): |H| = {H_hex_mag:.4f} A/m")
        print(f"  ANALYTICAL (tetra):   |H| = {H_tetra_mag:.4f} A/m")
        print(f"\n  Error (magnitude): {error_mag:.2f}%")
        print(f"  Error (components): x={error_x:.2f}%, y={error_y:.2f}%, z={error_z:.2f}%")

        if error_mag < 20:
            status = "PASS"
        elif error_mag < 50:
            status = "WARNING"
        else:
            status = "FAIL"

        print(f"\n  Status: {status} (target: <20%)")

        results.append({
            'mu_r': mu_r,
            'H_hex_mag': H_hex_mag,
            'H_tetra_mag': H_tetra_mag,
            'error_mag': error_mag,
            'status': status
        })
    else:
        print(f"  STANDARD (hex 5x5x5): {'SUCCESS' if H_hex_success else 'FAILED'}")
        print(f"  ANALYTICAL (tetra):   {'SUCCESS' if H_tetra_success else 'FAILED'}")
        print(f"\n  Status: FAILED (computation error)")

        results.append({
            'mu_r': mu_r,
            'H_hex_mag': H_hex_mag if H_hex_success else np.nan,
            'H_tetra_mag': H_tetra_mag if H_tetra_success else np.nan,
            'error_mag': np.nan,
            'status': 'FAILED'
        })

# ============================================================
# Summary Table
# ============================================================
print("\n" + "="*70)
print("Summary: High Permeability Stability Test")
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
            r['mu_r'],
            f"{r['H_hex_mag']:.4f}" if not np.isnan(r['H_hex_mag']) else "N/A",
            f"{r['H_tetra_mag']:.4f}" if not np.isnan(r['H_tetra_mag']) else "N/A",
            "N/A",
            r['status']))

print("\n" + "="*70)
print("Observations:")
print("="*70)

# Count successes
pass_count = sum(1 for r in results if r['status'] == 'PASS')
warning_count = sum(1 for r in results if r['status'] == 'WARNING')
fail_count = sum(1 for r in results if r['status'] in ['FAIL', 'FAILED'])

print(f"  PASS:    {pass_count}/{len(results)} (error < 20%)")
print(f"  WARNING: {warning_count}/{len(results)} (20% <= error < 50%)")
print(f"  FAIL:    {fail_count}/{len(results)} (error >= 50% or computation failed)")

print("\nRecommendations:")
if pass_count == len(results):
    print("  All tests passed! ANALYTICAL method is stable for all tested mu_r values.")
elif fail_count == 0:
    print("  ANALYTICAL method works but accuracy degrades at high mu_r.")
    print("  Consider limiting usage to mu_r < 1000 for best accuracy.")
else:
    max_working_mu = max([r['mu_r'] for r in results if r['status'] in ['PASS', 'WARNING']], default=0)
    print(f"  ANALYTICAL method is stable up to mu_r = {max_working_mu}")
    print(f"  Numerical instability observed for mu_r > {max_working_mu}")
    print("  Possible causes:")
    print("    - Interaction matrix conditioning issues")
    print("    - Magnetization magnitude scaling issues")
    print("    - Solver convergence problems")

print("="*70)
