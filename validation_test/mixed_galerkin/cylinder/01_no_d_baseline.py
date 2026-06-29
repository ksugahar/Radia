"""
Cylinder mixed Galerkin, rank-1 bulk + 1-DOF surface, NO free parameter.

Phase 2 of the 2026-05-28 -> 2026-06-12 research sprint.

Replaces the IGTE digest's Warburg-terminated CLN (which has a tuned
free parameter d, achieving 17% wall-band error) by an s-DEPENDENT
SIBC envelope used directly as a Galerkin basis function.

The bulk-surface crossover frequency is then determined automatically
by the Galerkin system; no d is needed.

## Setup

    u(r, s) = 1 + v(r, s),   v(a) = 0   (port voltage = 1)
    (-Laplacian + s mu sigma) v = -s mu sigma
    Y(s) = pi a^2 sigma * (1 + <v>)

## Basis

  Bulk (rank-1):   phi_0(r) = (a^2 - r^2) / 4    (-Laplacian phi_0 = 1)
  Surface (1-DOF): psi(r;s) = exp(-(a-r) t) - 1   (t = sqrt(s mu sigma))

ZERO free parameters.  Compare against the canonical
radia.maglev.mixed_galerkin.references cylinder helper.

## Phase 8b corrected result (2026-06-12)

After fixing the asymptote-switch artifact in Y_exact:

    max error anywhere = 0.064%    (was reported as 0.93% under buggy ref)
    wall band max      = 0.039%    (was reported as 0.93%)
    SIBC tail (1e8 Hz) ~ 0.046%

The 1-DOF mixed Galerkin is actually comparable to the sphere case
(0.11% wall band, see ../sphere/01_no_d_baseline.py).  Adding Senior
tower corrections is the route to sub-0.001% (see 02_senior_tower_*.py).
"""

from __future__ import annotations

import math
import cmath

import numpy as np
from scipy.integrate import quad

from radia.maglev.mixed_galerkin.references import (
    K_SIBC_cylinder,
    Y_DC_cylinder,
    Y_exact_cylinder,
)

# Physical parameters: copper, 5 mm-radius cylinder.
SIGMA = 5.8e7
MU = 4 * math.pi * 1e-7
A = 5e-3
MS = MU * SIGMA

Y_DC = Y_DC_cylinder(A, SIGMA)
K_SIBC = K_SIBC_cylinder(A, SIGMA, MU)


# -----------------------------------------------------------------------
# Bulk constants (analytic, since phi_0 is the s=0 Krylov solution).
# -----------------------------------------------------------------------
# phi_0(r) = (a^2 - r^2) / 4,   -Laplacian phi_0 = 1
# grad phi_0 = -r/2 e_r
# Cross-section integrals (2 pi r weight):
#   K0_bb = int |grad phi_0|^2 * 2 pi r dr   = pi a^4 / 8
#   K1_bb = int |phi_0|^2 * 2 pi r dr        = pi a^6 / 48
#   b_b   = int phi_0 * 2 pi r dr            = pi a^4 / 8
K0_BB = math.pi * A**4 / 8
K1_BB = math.pi * A**6 / 48
B_B = math.pi * A**4 / 8


# -----------------------------------------------------------------------
# Surface integrals (numerical quad, with boundary-layer adaptive points).
# Variable u = a - r,  dr = -du,  volume weight 2 pi (a - u).
# -----------------------------------------------------------------------
def _integrate_complex(f, lo, hi, t=None, limit=400):
    """Integrate a complex-valued f over [lo, hi], with boundary-layer
    adaptive breakpoints when |t| is large.
    """
    points = None
    if t is not None and abs(t) > 1.0:
        skin = 1.0 / abs(t)
        if skin < hi - lo:
            points = list(np.geomspace(skin / 100, min(20 * skin, (hi - lo) * 0.99), 30))
            points = [p for p in points if lo < p < hi]
    r_real, _ = quad(lambda x: f(x).real, lo, hi, limit=limit, points=points)
    r_imag, _ = quad(lambda x: f(x).imag, lo, hi, limit=limit, points=points)
    return complex(r_real, r_imag)


def Y_mixed_galerkin(s):
    """rank-1 bulk + 1-DOF surface Galerkin admittance."""
    t = cmath.sqrt(s * MS)

    # Surface block contributions.
    K0_ss = _integrate_complex(
        lambda u: t**2 * cmath.exp(-2 * u * t) * 2 * math.pi * (A - u), 0, A, t=t)
    K1_ss = _integrate_complex(
        lambda u: (cmath.exp(-u * t) - 1)**2 * 2 * math.pi * (A - u), 0, A, t=t)
    # Cross: grad phi_0 . grad psi = (-(A - u)/2)(t e^{-ut})
    K0_bs = _integrate_complex(
        lambda u: (-(A - u) / 2) * (t * cmath.exp(-u * t)) * 2 * math.pi * (A - u), 0, A, t=t)
    K1_bs = _integrate_complex(
        lambda u: ((2 * A * u - u**2) / 4) * (cmath.exp(-u * t) - 1) * 2 * math.pi * (A - u),
        0, A, t=t)
    b_s = _integrate_complex(
        lambda u: (cmath.exp(-u * t) - 1) * 2 * math.pi * (A - u), 0, A, t=t)

    K_mat = np.array(
        [
            [K0_BB + s * MS * K1_BB, K0_bs + s * MS * K1_bs],
            [K0_bs + s * MS * K1_bs, K0_ss + s * MS * K1_ss],
        ],
        dtype=complex,
    )
    b_vec = np.array([B_B, b_s], dtype=complex)

    xi = np.linalg.solve(K_mat, -s * MS * b_vec)
    v_avg = (xi[0] * B_B + xi[1] * b_s) / (math.pi * A**2)
    return Y_DC * (1 + v_avg)


def main():
    print("=== Cylinder mixed Galerkin: rank-1 bulk + 1-DOF surface ===")
    print(f"a = {A*1e3} mm, sigma = {SIGMA:.2e} S/m, mu = {MU:.4e} H/m")
    print(f"Y_DC   = {Y_DC:.4e} S*m")
    print(f"K_SIBC = {K_SIBC:.4e}")
    print()
    print("(reference: radia.maglev.mixed_galerkin.references.Y_exact_cylinder)")
    print()

    f_tests = [1.0, 1e2, 1e3, 1e4, 5e4, 1e5, 1.58e5, 2.51e5, 5e5, 1e6, 1e7, 1e8]
    print(f"{'f (Hz)':>10}  {'|Y_exact|':>12}  {'|Y_mixed|':>12}  {'rel err':>10}")
    for f in f_tests:
        s = 1j * 2 * math.pi * f
        Y_e = Y_exact_cylinder(s, A, SIGMA, MU)
        Y_m = Y_mixed_galerkin(s)
        err = abs(Y_e - Y_m) / abs(Y_e)
        print(f"{f:10.2e}  {abs(Y_e):12.4e}  {abs(Y_m):12.4e}  {err*100:9.4f}%")

    print()
    print("Full sweep (1 Hz to 1e8 Hz, 81 points):")
    fs = np.logspace(0, 8, 81)
    errs = []
    for f in fs:
        s = 1j * 2 * math.pi * f
        Y_e = Y_exact_cylinder(s, A, SIGMA, MU)
        Y_m = Y_mixed_galerkin(s)
        errs.append(abs(Y_e - Y_m) / abs(Y_e))
    errs = np.array(errs)
    wall_mask = (fs > 1e4) & (fs < 1e6)
    print(f"  Max error anywhere:  {errs.max()*100:.4f}% @ f = {fs[errs.argmax()]:.2e} Hz")
    print(f"  Wall-band max:       {errs[wall_mask].max()*100:.4f}%")
    print(f"  DC error (f=1 Hz):   {errs[0]*100:.6f}%")
    print(f"  SIBC tail (f=1e8):   {errs[-1]*100:.4f}%")


if __name__ == "__main__":
    main()
