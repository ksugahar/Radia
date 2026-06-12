"""
2D square mixed Galerkin v1 — DELIBERATELY BROKEN, kept as negative result.

This is the FIRST attempt at extending mixed Galerkin to the 2D square
cross-section (Phase 5, 2026-06-12 morning).  It uses the surface
basis

    psi(x, y; s) = f(x) sin(pi y / L) + sin(pi x / L) f(y)

where f(x) = cosh((x - L/2) t) / cosh(L t / 2) - 1.  This is symmetric
in x <-> y and vanishes on all 4 edges, so it is in V_0.  HOWEVER the
sin(pi y/L) factor is SMOOTH in y — it has no skin-layer structure on
the y faces.  The envelope therefore captures the boundary layer only
along the axis on which the f-factor sits, missing the orthogonal
boundary layer entirely.

Result: wall-band max ~ 89%, SIBC tail ~ 2672%.  Completely useless.

The CORRECTED v2 (`01_corner_envelope.py`) replaces the smooth sin
factor by another f, giving the tensor-product `psi = f(x) f(y)` which
captures BOTH boundary layers AND the corner Wiener-Hopf asymptote.
v2 gives a 0.03 - 0.26% accuracy.

This file is kept under an underscore prefix (`_broken_*`) so it is
clearly marked as not the canonical demonstration.  Its only value is
pedagogical: showing what fails and why is part of the research
record.
"""

from __future__ import annotations

import math
import cmath
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from _references.square2d_foster import (  # noqa: E402
    K_SIBC_square2d,
    Y_DC_square2d,
    Y_exact_square2d,
)

SIGMA = 5.8e7
MU = 4 * math.pi * 1e-7
L = 5e-3
MS = MU * SIGMA
Y_DC = Y_DC_square2d(L, SIGMA)
K_SIBC = K_SIBC_square2d(L, SIGMA, MU)

K0_BB = math.pi**2 / 2
K1_BB = L**2 / 4
B_B = 4 * L**2 / math.pi**2


def f_n(n: int, t):
    return -4 * L**2 * t**2 / (n * math.pi * (n**2 * math.pi**2 + L**2 * t**2))


def F_0(t):
    return (2.0 / t) * cmath.tanh(L * t / 2) - L


def Y_mixed_galerkin_v1(s, N_max: int = 199):
    """The BROKEN envelope: f(x) sin(pi y/L) + sin(pi x/L) f(y)."""
    t = cmath.sqrt(s * MS)

    # Bulk-only contribution from the (1,1) sine mode.
    K_b = K0_BB + s * MS * K1_BB

    # Surface block analytic — see project_mixed_galerkin_schur_dtn for the
    # full derivation.  Each axis term has Fourier coefficients only at
    # (n, 1) (for psi_x) or (1, p) (for psi_y) odd in n or p.
    f1 = f_n(1, t)
    K_ss = (L / 2)**2 * (2 * math.pi**2 / L**2 + s * MS) * (2 * f1)**2
    for n in range(3, N_max + 1, 2):
        fn = f_n(n, t)
        lam_n1 = (n**2 + 1) * math.pi**2 / L**2
        K_ss += 2 * (L / 2)**2 * (lam_n1 + s * MS) * fn**2

    K_bs = (L / 2)**2 * (2 * math.pi**2 / L**2 + s * MS) * (2 * f1)
    b_s = 4 * L * F_0(t) / math.pi

    K_mat = np.array([[K_b, K_bs], [K_bs, K_ss]], dtype=complex)
    b_vec = np.array([B_B, b_s], dtype=complex)
    xi = np.linalg.solve(K_mat, -s * MS * b_vec)
    v_avg = (xi @ b_vec) / L**2
    return Y_DC * (1 + v_avg)


def main():
    print("=== 2D square mixed Galerkin v1 (BROKEN, kept as negative result) ===")
    print("Envelope: psi = f(x) sin(pi y/L) + sin(pi x/L) f(y)")
    print("- The sin smooth factor has no y skin-layer structure.")
    print("- The CORRECT envelope is psi = f(x) f(y) (see 01_corner_envelope.py).")
    print()

    print(f"{'f (Hz)':>10}  {'|Y_exact|':>12}  {'|Y_v1|':>12}  {'rel err':>10}")
    for f in [1.0, 1e3, 1e4, 1e5, 7.94e5, 1e6]:
        s = 1j * 2 * math.pi * f
        Y_e = Y_exact_square2d(s, L, SIGMA, MU)
        Y_m = Y_mixed_galerkin_v1(s)
        err = abs(Y_e - Y_m) / abs(Y_e)
        print(f"{f:10.2e}  {abs(Y_e):12.4e}  {abs(Y_m):12.4e}  {err*100:9.2f}%")


if __name__ == "__main__":
    main()
