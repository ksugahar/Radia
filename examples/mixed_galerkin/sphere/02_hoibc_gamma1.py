"""
Sphere mixed Galerkin with HOIBC gamma_1 = -1/a curvature correction.

Phase 7 of the 2026-05-28 -> 2026-06-12 research sprint.

Adds a SECOND surface basis function reflecting the leading
curvature correction in the Senior tower:

    psi_1(r; s) = exp(-(a-r) t) - 1        (planar SIBC, same as 01)
    psi_2(r; s) = ((a-r) / a) exp(-(a-r) t) (curvature, gamma_1 = -1/a)

The (a-r)/a factor is the leading Taylor expansion term of the
(a/r) geometric factor in the sphere Bessel asymptote
    v(r) ~ -1 + (a/r) exp(-gamma (a-r)).

## Sphere Senior tower terminates at gamma_1

For the sphere, the Senior tower coefficients are b = [0, -1, 0, 0, ...]
(constructive derivation in memory
 project_senior_hoibc_tower_constructive_derivation.md).  Adding gamma_1
captures the entire non-trivial Senior tower; higher-order DOFs would
add nothing.  This is why the 2-DOF sphere result hits machine
precision at deep skin (1e-8 % at f = 10^8 Hz).

## Result

    1-DOF                            2-DOF (gamma_1)
    --------                         ----------------
    wall band max =  0.114%          wall band max =  0.0011% (100x improvement)
    max anywhere  =  0.137%          max anywhere  =  0.0037%
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


def Y_mixed_galerkin(s, with_gamma1: bool):
    """rank-1 bulk + 1-or-2-DOF surface.  with_gamma1=False reproduces 01."""
    t = cmath.sqrt(s * MS)

    K0_ss1 = _integrate_complex(
        lambda u: t**2 * cmath.exp(-2 * u * t) * 4 * math.pi * (A - u)**2, 0, A, t=t)
    K1_ss1 = _integrate_complex(
        lambda u: (cmath.exp(-u * t) - 1)**2 * 4 * math.pi * (A - u)**2, 0, A, t=t)
    K0_bs1 = _integrate_complex(
        lambda u: (-(A - u) / 3) * (t * cmath.exp(-u * t)) * 4 * math.pi * (A - u)**2, 0, A, t=t)
    K1_bs1 = _integrate_complex(
        lambda u: ((2 * A * u - u**2) / 6) * (cmath.exp(-u * t) - 1) * 4 * math.pi * (A - u)**2,
        0, A, t=t)
    b_s1 = _integrate_complex(
        lambda u: (cmath.exp(-u * t) - 1) * 4 * math.pi * (A - u)**2, 0, A, t=t)

    if not with_gamma1:
        K_mat = np.array(
            [
                [K0_BB + s * MS * K1_BB, K0_bs1 + s * MS * K1_bs1],
                [K0_bs1 + s * MS * K1_bs1, K0_ss1 + s * MS * K1_ss1],
            ],
            dtype=complex,
        )
        b_vec = np.array([B_B, b_s1], dtype=complex)
        xi = np.linalg.solve(K_mat, -s * MS * b_vec)
        v_avg = (xi[0] * B_B + xi[1] * b_s1) / V_SPH
        return Y_DC * (1 + v_avg)

    # 2-DOF surface: add psi_2 = (u/a) exp(-u t).
    # psi_2'_r = exp(-ut)/a * (u t - 1)
    K0_ss2 = _integrate_complex(
        lambda u: (cmath.exp(-u * t) / A * (u * t - 1))**2 * 4 * math.pi * (A - u)**2, 0, A, t=t)
    K1_ss2 = _integrate_complex(
        lambda u: ((u / A) * cmath.exp(-u * t))**2 * 4 * math.pi * (A - u)**2, 0, A, t=t)
    K0_ss12 = _integrate_complex(
        lambda u: t * cmath.exp(-u * t) * (cmath.exp(-u * t) / A * (u * t - 1)) *
        4 * math.pi * (A - u)**2, 0, A, t=t)
    K1_ss12 = _integrate_complex(
        lambda u: (cmath.exp(-u * t) - 1) * ((u / A) * cmath.exp(-u * t)) *
        4 * math.pi * (A - u)**2, 0, A, t=t)
    K0_bs2 = _integrate_complex(
        lambda u: (-(A - u) / 3) * ((1 / A) * cmath.exp(-u * t) * (u * t - 1)) *
        4 * math.pi * (A - u)**2, 0, A, t=t)
    K1_bs2 = _integrate_complex(
        lambda u: ((2 * A * u - u**2) / 6) * ((u / A) * cmath.exp(-u * t)) *
        4 * math.pi * (A - u)**2, 0, A, t=t)
    b_s2 = _integrate_complex(
        lambda u: ((u / A) * cmath.exp(-u * t)) * 4 * math.pi * (A - u)**2, 0, A, t=t)

    K11 = K0_ss1 + s * MS * K1_ss1
    K22 = K0_ss2 + s * MS * K1_ss2
    K12 = K0_ss12 + s * MS * K1_ss12
    K_b = K0_BB + s * MS * K1_BB
    K_b1 = K0_bs1 + s * MS * K1_bs1
    K_b2 = K0_bs2 + s * MS * K1_bs2

    K_mat = np.array(
        [[K_b, K_b1, K_b2], [K_b1, K11, K12], [K_b2, K12, K22]], dtype=complex
    )
    b_vec = np.array([B_B, b_s1, b_s2], dtype=complex)
    xi = np.linalg.solve(K_mat, -s * MS * b_vec)
    v_avg = (xi @ b_vec) / V_SPH
    return Y_DC * (1 + v_avg)


def main():
    print("=== Sphere mixed Galerkin: 1-DOF vs 2-DOF (HOIBC gamma_1 = -1/a) ===")
    print(f"a = {A*1e3} mm, sigma = {SIGMA:.2e} S/m, mu = {MU:.4e} H/m")
    print(f"gamma_1 = -1/a = {-1/A:.4e}")
    print()

    print(f"{'f (Hz)':>10}  {'|Y_exact|':>12}  {'err_1':>10}  {'err_2 (gamma_1)':>16}")
    for f in [1.0, 1e3, 1e4, 5e4, 1e5, 1e6, 1e7, 1e8]:
        s = 1j * 2 * math.pi * f
        Y_e = Y_exact_sphere(s, A, SIGMA, MU)
        e1 = abs(Y_e - Y_mixed_galerkin(s, False)) / abs(Y_e)
        e2 = abs(Y_e - Y_mixed_galerkin(s, True)) / abs(Y_e)
        print(f"{f:10.2e}  {abs(Y_e):12.4e}  {e1*100:8.4f}%  {e2*100:13.6f}%")

    print()
    fs = np.logspace(0, 8, 81)
    e1_all, e2_all = [], []
    for f in fs:
        s = 1j * 2 * math.pi * f
        Y_e = Y_exact_sphere(s, A, SIGMA, MU)
        e1_all.append(abs(Y_e - Y_mixed_galerkin(s, False)) / abs(Y_e))
        e2_all.append(abs(Y_e - Y_mixed_galerkin(s, True)) / abs(Y_e))
    e1_all = np.array(e1_all)
    e2_all = np.array(e2_all)
    wall = (fs > 1e4) & (fs < 1e6)
    print(f"Full sweep:")
    print(f"  1-DOF  max anywhere = {e1_all.max()*100:.4f}%,  wall band = {e1_all[wall].max()*100:.4f}%")
    print(f"  2-DOF  max anywhere = {e2_all.max()*100:.4f}%,  wall band = {e2_all[wall].max()*100:.4f}%")
    print(f"  Wall-band improvement factor: {e1_all[wall].max()/e2_all[wall].max():.1f}x")


if __name__ == "__main__":
    main()
