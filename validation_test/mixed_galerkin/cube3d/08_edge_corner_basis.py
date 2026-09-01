"""08_edge_corner_basis.py

Cube Mixed Galerkin with face + edge + corner envelope enrichment.

The baseline (03_exact_closed_form_Kss.py) uses a single tensor-product
envelope psi_xyz = f(x) f(y) f(z) and hits a 0.33% wall-band floor.
This script enriches the surface block 2 to 7 basis functions:

    face envelopes (3):    psi_x = f(x), psi_y = f(y), psi_z = f(z)
    edge envelopes (3):    psi_xy = f(x)f(y), psi_yz = f(y)f(z), psi_xz = f(x)f(z)
    corner envelope (1):   psi_xyz = f(x)f(y)f(z)

Hypothesis (Paper 1 §III addition): the cuboid Mellin asymptote
    Y(s) = c_0/sqrt(s) + c_1/s + c_2/s^{3/2} + ...
with c_1 ~ edge length and c_2 ~ vertex count requires basis components
of progressively lower codimension.  The tensor-product corner basis
alone reaches c_0 only.  Edge envelopes should add c_1.  Together with
the face envelopes they form a separable approximation of partition-of-
unity on the polyhedron.

Result interpretation:
    enriched error / baseline error -> 1     basis enrichment doesn't help
    enriched error / baseline error -> 0     basis fixes the floor
    in between                               quantifies the contribution
                                             of edge/face vs corner DOFs
"""

from __future__ import annotations

import cmath
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from _references.cube3d_foster import K_SIBC_cube3d, Y_DC_cube3d, Y_mellin_cube3d

SIGMA = 5.8e7
MU = 4 * math.pi * 1e-7
L = 5e-3
MS = MU * SIGMA
Y_DC = Y_DC_cube3d(L, SIGMA)
K_SIBC = K_SIBC_cube3d(L, SIGMA, MU)
V_CUBE = L ** 3


# ---------------------------------------------------------------------------
# 1D primitives (reused from 03_exact_closed_form_Kss.py)
# ---------------------------------------------------------------------------


def _tanh_safe(z):
    if abs(z.real) > 50:
        return complex(1.0, 0.0) if z.real > 0 else complex(-1.0, 0.0)
    return cmath.tanh(z)


def _sech2_safe(z):
    if abs(z.real) > 50:
        return complex(0.0, 0.0)
    return 1.0 / cmath.cosh(z) ** 2


def J_2_exact(t):
    """integral of (df/dx)^2 over [0, L]."""
    Lt2 = L * t / 2
    return t * _tanh_safe(Lt2) - (t ** 2 * L / 2) * _sech2_safe(Lt2)


def J_3_exact(t):
    """integral of f^2 over [0, L]."""
    Lt2 = L * t / 2
    return (L / 2) * _sech2_safe(Lt2) - 3 * _tanh_safe(Lt2) / t + L


def F_0_exact(t):
    """integral of f over [0, L]."""
    return (2 / t) * _tanh_safe(L * t / 2) - L


def phi1_overlap(t):
    """integral of f(x) sin(pi x / L) dx over [0, L]."""
    return 2 * L * math.pi / (math.pi ** 2 + L ** 2 * t ** 2) - 2 * L / math.pi


# ---------------------------------------------------------------------------
# Baseline: rank-1 bulk + 1-DOF corner envelope (= original 03 script)
# ---------------------------------------------------------------------------


def Y_mixed_baseline(s):
    """rank-1 bulk + 1-DOF psi_xyz tensor envelope (the 0.33% floor model)."""
    t = cmath.sqrt(s * MS)
    J_2 = J_2_exact(t)
    J_3 = J_3_exact(t)
    F_0 = F_0_exact(t)
    phi1 = phi1_overlap(t)

    lambda_111 = 3 * math.pi ** 2 / L ** 2
    K_b = (lambda_111 + s * MS) * (L / 2) ** 3
    b_b = (2 * L / math.pi) ** 3

    K_ss = 3 * J_2 * J_3 ** 2 + s * MS * J_3 ** 3
    b_psi = F_0 ** 3
    K_bs = (lambda_111 + s * MS) * phi1 ** 3

    K_mat = np.array([[K_b, K_bs], [K_bs, K_ss]], dtype=complex)
    b_vec = np.array([b_b, b_psi], dtype=complex)
    xi = np.linalg.solve(K_mat, -s * MS * b_vec)
    v_avg = (xi @ b_vec) / V_CUBE
    return Y_DC * (1 + v_avg)


# ---------------------------------------------------------------------------
# Enriched: rank-1 bulk + 7-DOF face/edge/corner basis
# ---------------------------------------------------------------------------
#
# Basis ordering (surface, indices 1..7 of the 8-DOF system):
#   1: psi_x   = f(x)
#   2: psi_y   = f(y)
#   3: psi_z   = f(z)
#   4: psi_xy  = f(x) f(y)
#   5: psi_yz  = f(y) f(z)
#   6: psi_xz  = f(x) f(z)
#   7: psi_xyz = f(x) f(y) f(z)
#
# Index 0 is the bulk phi_{1,1,1} = sin(pi x/L) sin(pi y/L) sin(pi z/L).


def Y_mixed_enriched(s):
    """rank-1 bulk + 7-DOF face/edge/corner enriched surface."""
    t = cmath.sqrt(s * MS)
    J_2 = J_2_exact(t)
    J_3 = J_3_exact(t)
    F_0 = F_0_exact(t)
    phi1 = phi1_overlap(t)

    lambda_111 = 3 * math.pi ** 2 / L ** 2
    K_b = (lambda_111 + s * MS) * (L / 2) ** 3
    b_b = (2 * L / math.pi) ** 3

    # Surface block K_ss (7x7).  All cross terms derived from separable
    # tensor structure of psi_x, psi_xy, psi_xyz factors.
    K_ss = np.zeros((7, 7), dtype=complex)

    # Diagonals
    K_face_diag = (J_2 + s * MS * J_3) * L ** 2          # face envelope
    K_edge_diag = (2 * J_2 * J_3 + s * MS * J_3 ** 2) * L  # edge envelope
    K_corn_diag = 3 * J_2 * J_3 ** 2 + s * MS * J_3 ** 3   # corner envelope
    for i in range(3):
        K_ss[i, i] = K_face_diag
    for i in range(3, 6):
        K_ss[i, i] = K_edge_diag
    K_ss[6, 6] = K_corn_diag

    # Face x Face off-diagonal: gradients orthogonal (each in different axis).
    # Mass cross = ∫ f(x) f(y) dV = F_0 · F_0 · L = F_0^2 L
    K_ff_off = s * MS * F_0 ** 2 * L
    for i in range(3):
        for j in range(i + 1, 3):
            K_ss[i, j] = K_ss[j, i] = K_ff_off

    # Face x Edge: cases depend on whether the face shares an axis with
    # the edge.  Mapping:
    #   psi_x (0) with psi_xy (3): shares x → "stiffness contributes"
    #   psi_x (0) with psi_yz (4): no share → grad orthogonal, mass only
    #   psi_x (0) with psi_xz (5): shares x
    #   psi_y (1) with psi_xy (3): shares y
    #   psi_y (1) with psi_yz (4): shares y
    #   psi_y (1) with psi_xz (5): no share
    #   psi_z (2) with psi_xy (3): no share
    #   psi_z (2) with psi_yz (4): shares z
    #   psi_z (2) with psi_xz (5): shares z
    K_fe_share = (J_2 + s * MS * J_3) * F_0 * L          # share one axis
    K_fe_noshare = s * MS * F_0 ** 3                     # no share (grad zero)
    face_edge_share_table = {
        (0, 3): True,  (0, 4): False, (0, 5): True,
        (1, 3): True,  (1, 4): True,  (1, 5): False,
        (2, 3): False, (2, 4): True,  (2, 5): True,
    }
    for (i, j), shares in face_edge_share_table.items():
        K_ss[i, j] = K_ss[j, i] = K_fe_share if shares else K_fe_noshare

    # Face x Corner (all 3 faces couple to 1 corner via shared axis)
    K_fc = (J_2 + s * MS * J_3) * F_0 ** 2
    for i in range(3):
        K_ss[i, 6] = K_ss[6, i] = K_fc

    # Edge x Edge off-diagonal: each pair shares one axis
    # psi_xy (3) × psi_yz (4): share y
    # psi_xy (3) × psi_xz (5): share x
    # psi_yz (4) × psi_xz (5): share z
    K_ee_off = (J_2 + s * MS * J_3) * F_0 ** 2
    for i in range(3, 6):
        for j in range(i + 1, 6):
            K_ss[i, j] = K_ss[j, i] = K_ee_off

    # Edge x Corner
    K_ec = (2 * J_2 * J_3 + s * MS * J_3 ** 2) * F_0
    for i in range(3, 6):
        K_ss[i, 6] = K_ss[6, i] = K_ec

    # Cross block K_bs (1x7): bulk phi_{1,1,1} = sin·sin·sin overlap with each surface envelope.
    # phi · psi_x integrates over y, z directions: ∫sin(pi y/L) dy = 2L/π.
    K_bs = np.zeros(7, dtype=complex)
    two_L_pi = 2 * L / math.pi
    for i in range(3):
        K_bs[i] = (lambda_111 + s * MS) * phi1 * two_L_pi ** 2
    for i in range(3, 6):
        K_bs[i] = (lambda_111 + s * MS) * phi1 ** 2 * two_L_pi
    K_bs[6] = (lambda_111 + s * MS) * phi1 ** 3

    # b vector (1 bulk + 7 surface) = ∫(·) dV over the cube
    b_psi = np.zeros(7, dtype=complex)
    for i in range(3):
        b_psi[i] = F_0 * L ** 2
    for i in range(3, 6):
        b_psi[i] = F_0 ** 2 * L
    b_psi[6] = F_0 ** 3

    # Assemble 8x8
    K_full = np.zeros((8, 8), dtype=complex)
    K_full[0, 0] = K_b
    K_full[0, 1:] = K_bs
    K_full[1:, 0] = K_bs
    K_full[1:, 1:] = K_ss

    b_full = np.zeros(8, dtype=complex)
    b_full[0] = b_b
    b_full[1:] = b_psi

    xi = np.linalg.solve(K_full, -s * MS * b_full)
    v_avg = (xi @ b_full) / V_CUBE
    return Y_DC * (1 + v_avg)


# ---------------------------------------------------------------------------
# Verification harness
# ---------------------------------------------------------------------------


def main():
    print("=== cube Mixed Galerkin: enrichment of surface DOF ===")
    print(f"L = {L * 1e3} mm, K_SIBC = {K_SIBC:.4e}")
    print(f"baseline = rank-1 bulk + 1-DOF corner envelope (the 0.33% floor model)")
    print(f"enriched = rank-1 bulk + 7-DOF face/edge/corner basis")
    print()
    print(
        f"{'f (Hz)':>10}  {'|Y_Mellin|':>14}  "
        f"{'|Y_baseline|':>14}  {'err_baseline':>12}  "
        f"{'|Y_enriched|':>14}  {'err_enriched':>12}  {'reduction':>10}"
    )

    for f in [1, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9]:
        s = 1j * 2 * math.pi * f
        Y_M = Y_mellin_cube3d(s, L, SIGMA, MU)
        Y_b = Y_mixed_baseline(s)
        Y_e = Y_mixed_enriched(s)
        err_b = abs(Y_b - Y_M) / abs(Y_M) * 100
        err_e = abs(Y_e - Y_M) / abs(Y_M) * 100
        reduction = err_b / err_e if err_e > 0 else float("inf")
        print(
            f"{f:10.2e}  {abs(Y_M):14.6e}  "
            f"{abs(Y_b):14.6e}  {err_b:11.4f}%  "
            f"{abs(Y_e):14.6e}  {err_e:11.6f}%  {reduction:9.1f}x"
        )


if __name__ == "__main__":
    main()
