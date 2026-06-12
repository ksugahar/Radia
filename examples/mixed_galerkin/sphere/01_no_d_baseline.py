"""
Sphere mixed Galerkin, rank-1 bulk + 1-DOF surface, NO free parameter.

Phase 4 of the 2026-05-28 -> 2026-06-12 research sprint.

The natural 3-D extension of the cylinder Phase 2 (../cylinder/01).
Volume element 4 pi r^2 dr; bulk Krylov phi_0 satisfies the spherical
Poisson equation.

## Setup

    u(r, s) = 1 + v(r, s),   v(a) = 0
    (-Laplacian + s mu sigma) v = -s mu sigma
    Y(s) = (4/3) pi a^3 sigma * (1 + <v>)

## Basis

  Bulk (rank-1):   phi_0(r) = (a^2 - r^2) / 6     (-Laplacian phi_0 = 1 in spherical)
  Surface (1-DOF): psi(r;s) = exp(-(a-r) t) - 1   (planar SIBC, same as cylinder)

ZERO free parameters.

## Result

    max error anywhere = 0.14% (at the DC-to-AC transition, ~6 kHz)
    wall band max      = 0.11%
    DC error           = 0.000%
    SIBC tail (1e8 Hz) = 0.000%

The 0.11% wall band is the planar-SIBC envelope intrinsic error: it
misses the leading curvature correction gamma_1 = -1/a (sphere mean
curvature).  Adding gamma_1 in 02_hoibc_gamma1.py drives this down to
0.001% (100x improvement), and since the sphere Senior tower
TERMINATES at gamma_1 (b_k = [0, -1, 0, 0, ...] for the sphere
expansion of Y/K_0), this 2-DOF basis is already optimal.
"""

from __future__ import annotations

import math
import cmath
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import quad

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from _references.sphere_bessel import (  # noqa: E402
    K_SIBC_sphere,
    Y_DC_sphere,
    Y_exact_sphere,
)

SIGMA = 5.8e7
MU = 4 * math.pi * 1e-7
A = 5e-3
MS = MU * SIGMA
V_SPH = (4.0 / 3.0) * math.pi * A**3
Y_DC = Y_DC_sphere(A, SIGMA)
K_SIBC = K_SIBC_sphere(A, SIGMA, MU)

# Bulk constants (analytic, spherical volume weight 4 pi r^2):
#   phi_0(r) = (a^2 - r^2) / 6,  grad phi_0 = -r/3
#   K0_bb = 4 pi a^5 / 45,  K1_bb = 8 pi a^7 / 945,  b_b = 4 pi a^5 / 45
K0_BB = 4 * math.pi * A**5 / 45
K1_BB = 8 * math.pi * A**7 / 945
B_B = 4 * math.pi * A**5 / 45


def _integrate_complex(f, lo, hi, t=None, limit=400):
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
    t = cmath.sqrt(s * MS)
    # Substitution u = a - r, volume weight 4 pi (a - u)^2.
    K0_ss = _integrate_complex(
        lambda u: t**2 * cmath.exp(-2 * u * t) * 4 * math.pi * (A - u)**2, 0, A, t=t)
    K1_ss = _integrate_complex(
        lambda u: (cmath.exp(-u * t) - 1)**2 * 4 * math.pi * (A - u)**2, 0, A, t=t)
    # Cross: grad phi_0 . grad psi = (-(a-u)/3)(t e^{-ut})
    K0_bs = _integrate_complex(
        lambda u: (-(A - u) / 3) * (t * cmath.exp(-u * t)) * 4 * math.pi * (A - u)**2, 0, A, t=t)
    K1_bs = _integrate_complex(
        lambda u: ((2 * A * u - u**2) / 6) * (cmath.exp(-u * t) - 1) * 4 * math.pi * (A - u)**2,
        0, A, t=t)
    b_s = _integrate_complex(
        lambda u: (cmath.exp(-u * t) - 1) * 4 * math.pi * (A - u)**2, 0, A, t=t)

    K_mat = np.array(
        [
            [K0_BB + s * MS * K1_BB, K0_bs + s * MS * K1_bs],
            [K0_bs + s * MS * K1_bs, K0_ss + s * MS * K1_ss],
        ],
        dtype=complex,
    )
    b_vec = np.array([B_B, b_s], dtype=complex)
    xi = np.linalg.solve(K_mat, -s * MS * b_vec)
    v_avg = (xi[0] * B_B + xi[1] * b_s) / V_SPH
    return Y_DC * (1 + v_avg)


def main():
    print("=== Sphere mixed Galerkin: rank-1 bulk + 1-DOF surface ===")
    print(f"a = {A*1e3} mm, sigma = {SIGMA:.2e} S/m, mu = {MU:.4e} H/m")
    print(f"V_sphere = {V_SPH:.4e} m^3")
    print(f"Y_DC     = {Y_DC:.4e} S*m^2")
    print(f"K_SIBC   = {K_SIBC:.4e}")
    print()

    print(f"{'f (Hz)':>10}  {'|Y_exact|':>12}  {'|Y_mixed|':>12}  {'rel err':>10}")
    for f in [1.0, 1e3, 1e4, 5e4, 1e5, 1e6, 1e7, 1e8]:
        s = 1j * 2 * math.pi * f
        Y_e = Y_exact_sphere(s, A, SIGMA, MU)
        Y_m = Y_mixed_galerkin(s)
        err = abs(Y_e - Y_m) / abs(Y_e)
        print(f"{f:10.2e}  {abs(Y_e):12.4e}  {abs(Y_m):12.4e}  {err*100:9.4f}%")

    print()
    fs = np.logspace(0, 8, 81)
    errs = []
    for f in fs:
        s = 1j * 2 * math.pi * f
        Y_e = Y_exact_sphere(s, A, SIGMA, MU)
        errs.append(abs(Y_e - Y_mixed_galerkin(s)) / abs(Y_e))
    errs = np.array(errs)
    wall_mask = (fs > 1e4) & (fs < 1e6)
    print(f"Full sweep summary (1 Hz to 1e8 Hz, 81 points):")
    print(f"  Max error anywhere:  {errs.max()*100:.4f}% @ f = {fs[errs.argmax()]:.2e} Hz")
    print(f"  Wall band max:       {errs[wall_mask].max()*100:.4f}%")


if __name__ == "__main__":
    main()
