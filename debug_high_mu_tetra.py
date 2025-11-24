#!/usr/bin/env python
"""
Debug high permeability issues in ANALYTICAL tetrahedral method

This script investigates WHY the ANALYTICAL method fails for mu_r > 100:
1. Check magnetization values inside the tetrahedron
2. Compare solver convergence behavior
3. Investigate interaction matrix conditioning
4. Test with different evaluation points
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'build/Release'))

os.environ['RADIA_TETRA_METHOD'] = 'ANALYTICAL'

import numpy as np
import radia as rad

print("="*70)
print("Debug: High Permeability Tetrahedral Method Issues")
print("="*70)

# Test parameters
cube_size = 0.1  # m
H_ext = 1000.0   # A/m (z-direction)

# Test two cases: working (mu_r=100) and failing (mu_r=1000)
test_cases = [
    {'mu_r': 100, 'label': 'Working case'},
    {'mu_r': 1000, 'label': 'Failing case'}
]

for case in test_cases:
    mu_r = case['mu_r']
    label = case['label']
    ksi = mu_r - 1

    print(f"\n{'='*70}")
    print(f"{label}: mu_r = {mu_r}")
    print(f"{'='*70}")

    # ============================================================
    # Create single tetrahedron
    # ============================================================
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
    mat = rad.MatLin([ksi, ksi], [0, 0, 1])
    rad.MatApl(tetra, mat)

    bg = rad.ObjBckg([0, 0, H_ext])
    container = rad.ObjCnt([tetra, bg])

    # ============================================================
    # Solve and check convergence
    # ============================================================
    print(f"\n[1] Solver Convergence")
    print(f"  Input: mu_r={mu_r}, ksi={ksi}, H_ext={H_ext} A/m")

    result = rad.Solve(container, 0.0001, 10000)
    print(f"  Solve result: {result}")
    print(f"    result[0] = {result[0]:.6e} (final residual)")
    print(f"    result[1] = {result[1]:.6e} (max interaction energy?)")
    print(f"    result[2] = {result[2]:.6e} (scale factor?)")
    print(f"    result[3] = {result[3]:.1f} (iteration count)")

    if result[3] >= 10000:
        print(f"  WARNING: Maximum iterations reached!")
    elif result[0] > 0.0001:
        print(f"  WARNING: Convergence tolerance not met!")

    # ============================================================
    # Check magnetization inside tetrahedron
    # ============================================================
    print(f"\n[2] Magnetization Inside Tetrahedron")

    # Sample points inside the tetrahedron
    # Tetrahedron barycentric coordinates: (0.25, 0.25, 0.25, 0.25)
    center = [cube_size/4, cube_size/4, cube_size/4]

    try:
        M_center = rad.Fld(container, 'm', center)
        M_center_mag = np.linalg.norm(M_center)
        print(f"  Center point [{center[0]:.4f}, {center[1]:.4f}, {center[2]:.4f}]:")
        print(f"    M = [{M_center[0]:.6e}, {M_center[1]:.6e}, {M_center[2]:.6e}] A/m")
        print(f"    |M| = {M_center_mag:.6e} A/m")

        # Expected magnetization: M = ksi * H_total
        # For linear material in uniform field: M_z ~ ksi * H_ext (approximately)
        M_expected = ksi * H_ext
        print(f"    Expected |M| ~ {M_expected:.6e} A/m (ksi * H_ext)")
        M_ratio = M_center_mag / M_expected
        print(f"    Ratio |M|/expected = {M_ratio:.4f}")

        if M_ratio > 10 or M_ratio < 0.1:
            print(f"    WARNING: Magnetization is {M_ratio:.1f}x expected value!")
    except Exception as e:
        print(f"  ERROR evaluating magnetization: {e}")

    # ============================================================
    # Field evaluation at multiple distances
    # ============================================================
    print(f"\n[3] Field at Various Distances")

    # Evaluation points at different distances from tetrahedron surface
    eval_points = [
        ([0.05, 0.05, 0.15], "Near (0.05m from surface)"),
        ([0.05, 0.05, 0.2], "Medium (0.1m from surface)"),
        ([0.05, 0.05, 0.5], "Far (0.4m from surface)"),
        ([0.05, 0.05, 1.0], "Very far (0.9m from surface)")
    ]

    for pt, description in eval_points:
        try:
            H = rad.Fld(container, 'h', pt)
            H_mag = np.linalg.norm(H)
            print(f"  {description}:")
            print(f"    Point: [{pt[0]:.3f}, {pt[1]:.3f}, {pt[2]:.3f}]")
            print(f"    H = [{H[0]:.4e}, {H[1]:.4e}, {H[2]:.4e}] A/m")
            print(f"    |H| = {H_mag:.4e} A/m")

            # Check if field is reasonable (should be ~ H_ext far away)
            if pt[2] > 0.5:  # Far field
                deviation = abs(H_mag - H_ext) / H_ext * 100
                print(f"    Deviation from H_ext: {deviation:.2f}%")
                if deviation > 50:
                    print(f"    WARNING: Far field does not approach H_ext!")

        except Exception as e:
            print(f"  {description}: ERROR - {e}")

    # ============================================================
    # Compare with background field only (no magnet)
    # ============================================================
    print(f"\n[4] Background Field Reference")

    rad.FldUnits('m')
    bg_only = rad.ObjBckg([0, 0, H_ext])

    eval_pt = [0.05, 0.05, 0.2]
    H_bg = rad.Fld(bg_only, 'h', eval_pt)
    H_bg_mag = np.linalg.norm(H_bg)

    print(f"  Background field only (no magnet):")
    print(f"    H = [{H_bg[0]:.4f}, {H_bg[1]:.4f}, {H_bg[2]:.4f}] A/m")
    print(f"    |H| = {H_bg_mag:.4f} A/m")
    print(f"  Expected: [0, 0, {H_ext}] A/m")

    print(f"\n{'='*70}")

# ============================================================
# Summary and Diagnosis
# ============================================================
print(f"\n{'='*70}")
print("Diagnosis Summary")
print(f"{'='*70}")

print("""
Based on the debug output:

1. Solver Convergence:
   - Check if solver reaches maximum iterations
   - Check if final residual is within tolerance
   - Compare result[0] between working and failing cases

2. Magnetization Magnitude:
   - For mu_r=100: M should be ~ 99 * 1000 = 99,000 A/m
   - For mu_r=1000: M should be ~ 999 * 1000 = 999,000 A/m
   - If M is much larger/smaller, interaction matrix may be wrong

3. Field Behavior:
   - Near field: Depends on magnetization distribution
   - Far field: Should approach H_ext = 1000 A/m
   - If far field is wrong, there's a fundamental issue

4. Possible Root Causes:
   a) Interaction matrix scaling issue with high ksi
   b) Magnetization calculation overflow/underflow
   c) Analytical field formula numerical instability
   d) Solver divergence due to poor conditioning

Next steps:
- If M is correct but H is wrong: Issue in field evaluation
- If M is wrong: Issue in interaction matrix or solve
- If solver doesn't converge: Conditioning problem
""")

print("="*70)
