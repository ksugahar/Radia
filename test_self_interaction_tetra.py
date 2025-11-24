#!/usr/bin/env python
"""
Test self-interaction (demagnetization) for tetrahedral elements

The issue with high permeability may be related to incorrect self-interaction
(demagnetization tensor) calculation for tetrahedral elements.

For a tetrahedron in external field H_ext with relative permeability mu_r:
  M = ksi * H_internal
  H_internal = H_ext + H_demag
  H_demag = -N * M  (demagnetization field)

where N is the demagnetization tensor.

If N is wrong, M will be wrong, and the external field will be wrong.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'build/Release'))

os.environ['RADIA_TETRA_METHOD'] = 'ANALYTICAL'

import numpy as np
import radia as rad

print("="*70)
print("Self-Interaction (Demagnetization) Test for Tetrahedra")
print("="*70)

cube_size = 0.1  # m
H_ext = 1000.0   # A/m
mu0 = 4*np.pi*1e-7  # H/m
B_ext = H_ext * mu0  # Tesla

# ============================================================
# Test 1: Hexahedron (Known to work correctly)
# ============================================================
print("\n[1] Reference: Hexahedron (single element, no subdivision)")

rad.FldUnits('m')

cube_hex = rad.ObjRecMag([0, 0, 0], [cube_size, cube_size, cube_size], [0, 0, 0])

mu_r = 100
ksi = mu_r - 1
mat_hex = rad.MatLin([ksi, ksi], [0, 0, 1])
rad.MatApl(cube_hex, mat_hex)

bg_hex = rad.ObjBckg([0, 0, B_ext])
container_hex = rad.ObjCnt([cube_hex, bg_hex])

result_hex = rad.Solve(container_hex, 0.0001, 10000)
print(f"  Solve result: {result_hex}")

# Get magnetization at center
M_hex = rad.Fld(container_hex, 'm', [0, 0, 0])
M_hex_mag = np.linalg.norm(M_hex)

print(f"  Magnetization at center: M = [{M_hex[0]:.2f}, {M_hex[1]:.2f}, {M_hex[2]:.2f}] A/m")
print(f"  |M| = {M_hex_mag:.2f} A/m")

# Theoretical calculation for rectangular prism
# For cube in z-direction field: N_z ~ 1/3
# M_z = ksi * H_internal = ksi * (H_ext - N_z * M_z)
# M_z * (1 + ksi * N_z) = ksi * H_ext
# M_z = ksi * H_ext / (1 + ksi * N_z)

N_z_cube = 1.0/3.0  # Demagnetization factor for cube in z-direction
M_theory_cube = ksi * H_ext / (1 + ksi * N_z_cube)
print(f"  Theoretical M (cube, N_z=1/3): {M_theory_cube:.2f} A/m")
print(f"  Ratio M_actual/M_theory: {M_hex_mag / M_theory_cube:.4f}")

# ============================================================
# Test 2: Single Tetrahedron
# ============================================================
print("\n[2] Test: Single Tetrahedron")

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

bg_tetra = rad.ObjBckg([0, 0, B_ext])
container_tetra = rad.ObjCnt([tetra, bg_tetra])

result_tetra = rad.Solve(container_tetra, 0.0001, 10000)
print(f"  Solve result: {result_tetra}")

# Get magnetization at center
tetra_center = [cube_size/4, cube_size/4, cube_size/4]
M_tetra = rad.Fld(container_tetra, 'm', tetra_center)
M_tetra_mag = np.linalg.norm(M_tetra)

print(f"  Magnetization at center: M = [{M_tetra[0]:.2f}, {M_tetra[1]:.2f}, {M_tetra[2]:.2f}] A/m")
print(f"  |M| = {M_tetra_mag:.2f} A/m")

# For tetrahedron, demagnetization factor is geometry-dependent
# Rough estimate: N_z ~ 0.25 to 0.35 for tetrahedron
N_z_tetra_est = 0.3  # Estimated
M_theory_tetra = ksi * H_ext / (1 + ksi * N_z_tetra_est)
print(f"  Theoretical M (tetra, N_z~0.3): {M_theory_tetra:.2f} A/m")
print(f"  Ratio M_actual/M_theory: {M_tetra_mag / M_theory_tetra:.4f}")

# ============================================================
# Analysis: Extract demagnetization factor from results
# ============================================================
print("\n[3] Analysis: Demagnetization Factor Extraction")

# From M = ksi * (H_ext - N * M)
# M * (1 + ksi * N) = ksi * H_ext
# N = (ksi * H_ext / M - 1) / ksi
# N = H_ext / M - 1 / ksi

if M_hex_mag > 0:
    N_hex_extracted = H_ext / M_hex_mag - 1.0 / ksi
    print(f"  Hexahedron:")
    print(f"    Extracted N_z = {N_hex_extracted:.4f}")
    print(f"    Expected N_z (cube) = 0.3333")
    print(f"    Deviation: {abs(N_hex_extracted - 0.3333) / 0.3333 * 100:.2f}%")

if M_tetra_mag > 0:
    N_tetra_extracted = H_ext / M_tetra_mag - 1.0 / ksi
    print(f"  Tetrahedron:")
    print(f"    Extracted N_z = {N_tetra_extracted:.4f}")
    print(f"    Expected N_z (tetra, estimate) ~ 0.25-0.35")

    if N_tetra_extracted < 0:
        print(f"    ERROR: Negative demagnetization factor!")
        print(f"           This indicates M > ksi * H_ext")
        print(f"           Self-interaction is UNDERESTIMATED")
    elif N_tetra_extracted > 1:
        print(f"    ERROR: Demagnetization factor > 1!")
        print(f"           Self-interaction is OVERESTIMATED")

# ============================================================
# Test 3: Multiple permeability values
# ============================================================
print("\n[4] Demagnetization Factor vs Permeability")

print("\n{:<10} {:<15} {:<15}".format("mu_r", "N_z (hex)", "N_z (tetra)"))
print("-"*40)

for mu_r_test in [10, 100, 1000]:
    ksi_test = mu_r_test - 1

    # Hexahedron
    rad.FldUnits('m')
    cube_hex_test = rad.ObjRecMag([0, 0, 0], [cube_size, cube_size, cube_size], [0, 0, 0])
    mat_hex_test = rad.MatLin([ksi_test, ksi_test], [0, 0, 1])
    rad.MatApl(cube_hex_test, mat_hex_test)
    bg_hex_test = rad.ObjBckg([0, 0, B_ext])
    container_hex_test = rad.ObjCnt([cube_hex_test, bg_hex_test])
    rad.Solve(container_hex_test, 0.0001, 10000)
    M_hex_test = rad.Fld(container_hex_test, 'm', [0, 0, 0])
    M_hex_test_mag = np.linalg.norm(M_hex_test)
    N_hex_test = H_ext / M_hex_test_mag - 1.0 / ksi_test if M_hex_test_mag > 0 else np.nan

    # Tetrahedron
    rad.FldUnits('m')
    tetra_test = rad.ObjPolyhdr(vertices, TETRA_FACES, [0, 0, 0])
    mat_tetra_test = rad.MatLin([ksi_test, ksi_test], [0, 0, 1])
    rad.MatApl(tetra_test, mat_tetra_test)
    bg_tetra_test = rad.ObjBckg([0, 0, B_ext])
    container_tetra_test = rad.ObjCnt([tetra_test, bg_tetra_test])
    rad.Solve(container_tetra_test, 0.0001, 10000)
    M_tetra_test = rad.Fld(container_tetra_test, 'm', tetra_center)
    M_tetra_test_mag = np.linalg.norm(M_tetra_test)
    N_tetra_test = H_ext / M_tetra_test_mag - 1.0 / ksi_test if M_tetra_test_mag > 0 else np.nan

    print(f"{mu_r_test:<10} {N_hex_test:<15.4f} {N_tetra_test:<15.4f}")

# ============================================================
# Summary
# ============================================================
print("\n" + "="*70)
print("Summary")
print("="*70)

print("""
Demagnetization factor N should be:
- Geometry-dependent (NOT permeability-dependent)
- Between 0 and 1 (for ellipsoid/polyhedron)
- For cube: N_z ~ 1/3 = 0.3333
- For tetrahedron: N_z ~ 0.25-0.35 (depends on orientation)

If N varies with mu_r, there's a problem with self-interaction calculation.
If N < 0, magnetization is too high (self-interaction underestimated).
If N > 1, magnetization is too low (self-interaction overestimated).

Expected behavior:
- N should be approximately constant across different mu_r values
- Deviation should be < 10%
""")

print("="*70)
