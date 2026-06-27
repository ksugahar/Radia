"""
Cylinder mixed Galerkin, rank-1 bulk + N-DOF Senior tower surface.

Phase 8b of the 2026-05-28 -> 2026-06-12 research sprint, corrected.

The surface basis adds HOIBC Senior tower corrections beyond the
planar SIBC envelope of ../01_no_d_baseline.py:

    psi_1(r; s) = exp(-(a-r) t) - 1                   (planar SIBC)
    psi_2(r; s) = ((a-r) / a) exp(-(a-r) t)            (gamma_1 = -1/(2a))
    psi_3(r; s) = ((a-r) / a)^2 exp(-(a-r) t)          (gamma_2 = -1/(8 a^2))
    psi_4(r; s) = ((a-r) / a)^3 exp(-(a-r) t)          (gamma_3 = -1/(8 a^3))

This is the Taylor expansion of the cylinder Bessel I_1/I_0 asymptotic
factor sqrt(a/r) = 1 + d/(2a) + 3 d^2/(8 a^2) + ... in powers of
d = a - r.  The Galerkin coefficients automatically pick up the
canonical Senior tower factors 1/2, 1/8, etc.

## Senior tower of the cylinder (memory)

Cylinder Y/K_0 ~ 1 + sum_k b_k / (gamma a)^k with the DIVERGENT
asymptotic series

    b = [0, -1/2, -1/8, -1/8, -25/128, -13/32, -1073/1024, ...]

Senior 1962 / Mitzner 1967 / Yuferev-Ida 2010.

## Phase 8b result

With proper Y_exact (full Bessel + Senior-tower asymptotic
continuation past the scipy iv() overflow region):

   basis       | wall-band max | max anywhere
   ------------+---------------+--------------
   1-DOF       |    0.0386%    |   0.0639%
   2-DOF (g_1) |    0.00024%   |   0.0123%
   3-DOF (g_2) |    0.00002%   |   0.0027%
   4-DOF (g_3) |    0.00001%   |   0.0006%

i.e. each Senior tower DOF gives ~100x improvement at wall band.
Compare against the (Phase 7) sphere result, where Senior tower
TERMINATES at gamma_1 so 2-DOF is already optimal.
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

SIGMA = 5.8e7
MU = 4 * math.pi * 1e-7
A = 5e-3
MS = MU * SIGMA
Y_DC = Y_DC_cylinder(A, SIGMA)
K_SIBC = K_SIBC_cylinder(A, SIGMA, MU)

K0_BB = math.pi * A**4 / 8
K1_BB = math.pi * A**6 / 48
B_B = math.pi * A**4 / 8


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


def psi(k: int, u, t):
    """k-th surface basis function evaluated at u = a - r.

    k = 1 ... 4 correspond to planar SIBC + Senior tower g_1 ... g_3.
    """
    if k == 1:
        return cmath.exp(-u * t) - 1
    if k == 2:
        return (u / A) * cmath.exp(-u * t)
    if k == 3:
        return (u / A)**2 * cmath.exp(-u * t)
    if k == 4:
        return (u / A)**3 * cmath.exp(-u * t)
    raise ValueError(f"unsupported k = {k}")


def dpsi(k: int, u, t):
    """Radial derivative d psi / dr at u = a - r.  Note: dr = -du."""
    if k == 1:
        return t * cmath.exp(-u * t)
    if k == 2:
        return cmath.exp(-u * t) / A * (u * t - 1)
    if k == 3:
        return (u / A**2) * cmath.exp(-u * t) * (u * t - 2)
    if k == 4:
        return (u**2 / A**3) * cmath.exp(-u * t) * (u * t - 3)
    raise ValueError(f"unsupported k = {k}")


def Y_mixed_galerkin(s, N_surf: int):
    """rank-1 bulk + N_surf-DOF surface Galerkin admittance."""
    t = cmath.sqrt(s * MS)
    dim = 1 + N_surf
    K_mat = np.zeros((dim, dim), dtype=complex)
    b_vec = np.zeros(dim, dtype=complex)

    K_mat[0, 0] = K0_BB + s * MS * K1_BB
    b_vec[0] = B_B

    for k in range(1, N_surf + 1):
        # Cross with bulk phi_0 = (a^2 - r^2)/4, grad phi_0 = -(A - u)/2 (radial).
        K0_bk = _integrate_complex(
            lambda u, kk=k: (-(A - u) / 2) * dpsi(kk, u, t) * 2 * math.pi * (A - u),
            0, A, t=t)
        K1_bk = _integrate_complex(
            lambda u, kk=k: ((2 * A * u - u**2) / 4) * psi(kk, u, t) * 2 * math.pi * (A - u),
            0, A, t=t)
        K_mat[0, k] = K0_bk + s * MS * K1_bk
        K_mat[k, 0] = K_mat[0, k]
        b_vec[k] = _integrate_complex(
            lambda u, kk=k: psi(kk, u, t) * 2 * math.pi * (A - u), 0, A, t=t)
        for j in range(1, k + 1):
            K0_kj = _integrate_complex(
                lambda u, kk=k, jj=j: dpsi(kk, u, t) * dpsi(jj, u, t) * 2 * math.pi * (A - u),
                0, A, t=t)
            K1_kj = _integrate_complex(
                lambda u, kk=k, jj=j: psi(kk, u, t) * psi(jj, u, t) * 2 * math.pi * (A - u),
                0, A, t=t)
            K_mat[k, j] = K0_kj + s * MS * K1_kj
            K_mat[j, k] = K_mat[k, j]

    xi = np.linalg.solve(K_mat, -s * MS * b_vec)
    v_avg = (xi @ b_vec) / (math.pi * A**2)
    return Y_DC * (1 + v_avg)


def main():
    print("=== Cylinder mixed Galerkin: rank-1 bulk + N-DOF Senior tower ===")
    print(f"a = {A*1e3} mm, sigma = {SIGMA:.2e} S/m, mu = {MU:.4e} H/m")
    print()
    print(f"Sample points:")
    print(f"{'f (Hz)':>10}  {'1-DOF':>10}  {'2-DOF g_1':>10}  {'3-DOF g_2':>10}  {'4-DOF g_3':>10}")
    for f in [1e3, 1e4, 5e4, 1e5, 5e5, 1e6, 1e7, 1e8]:
        s = 1j * 2 * math.pi * f
        Y_e = Y_exact_cylinder(s, A, SIGMA, MU)
        row = [
            abs(Y_e - Y_mixed_galerkin(s, N)) / abs(Y_e) * 100 for N in (1, 2, 3, 4)
        ]
        print(f"{f:10.2e}  {row[0]:8.4f}%  {row[1]:8.4f}%  {row[2]:8.4f}%  {row[3]:8.4f}%")

    print()
    print("Full sweep summary (1 Hz to 1e8 Hz, 81 points):")
    fs = np.logspace(0, 8, 81)
    wall_mask = (fs > 1e4) & (fs < 1e6)
    print(f"{'N-DOF':>6}  {'max anywhere':>16}  {'wall band max':>16}")
    for N in (1, 2, 3, 4):
        errs = []
        for f in fs:
            s = 1j * 2 * math.pi * f
            Y_e = Y_exact_cylinder(s, A, SIGMA, MU)
            errs.append(abs(Y_e - Y_mixed_galerkin(s, N)) / abs(Y_e))
        errs = np.array(errs)
        print(f"  {N}    {errs.max()*100:13.5f}%   {errs[wall_mask].max()*100:13.5f}%")


if __name__ == "__main__":
    main()
