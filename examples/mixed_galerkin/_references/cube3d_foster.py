"""
Reference Y_exact for a solid conducting cube of side L (scalar PDE form).

For the scalar PDE  (-Laplacian + s mu sigma) v = -s mu sigma  on
[0, L]^3 with v = 0 on all 6 faces, the eigenfunction expansion is a
TRIPLE sum over (n, p, q) all odd:

    Y/Y_DC = sum_{n,p,q odd} (512/pi^6) * lambda_{npq} / (lambda_{npq} + s mu sigma)
                                       / (n^2 p^2 q^2)

with lambda_{npq} = (n^2+p^2+q^2) pi^2 / L^2,  Y_DC = L^3 sigma,
and (in the limit s -> infinity) the Mellin asymptote

    Y(s) ~ K_SIBC / sqrt(s) + c_1 / s + c_2 / s^(3/2) + ...
    K_SIBC = 6 L^2 sqrt(sigma/mu)
    c_1    = -16 L (a+b+c)/(pi mu)  = -48 L / (pi mu)  for a = b = c = L
    c_2    = +48 / (pi mu^(3/2) sqrt(sigma))   (dim-independent 90deg-corner term)

(See memory project_cuboid_mellin_asymptote_result for the closed-form
derivation of c_0, c_1, c_2.)

## STATUS: REFERENCE AUDIT PENDING (task #183 as of 2026-06-12)

The Phase 6 mixed Galerkin script (`../cube3d/01_corner_envelope_uncertain.py`)
reports a 20% wall-band peak error.  It was measured against Foster
truncated at N_max = 99, which from the 2D experience (Phase 8c, see
`square2d_foster.py`) is almost certainly under-converged — the 2D
case at N_max = 1999 had 5% bias at the wall-band peak.

Until task #183 is completed, treat the 3D cube result as PROVISIONAL.
The honest expectation, by analogy with 2D, is that most of the 20%
is Foster truncation and the true mixed-Galerkin error is in the
0.1% - 1% range (not yet quantified).

For benchmarking, this module provides:

    Y_foster_cube3d(s, L, sigma, mu, N_max)  - truncated 3D Foster sum
    Y_aitken_cube3d(s, L, sigma, mu)         - Aitken-extrapolated (N=99, 199, 399)
    Y_mellin_cube3d(s, L, sigma, mu)         - Mellin asymptote (3 terms)
    Y_exact_cube3d(s, L, sigma, mu)          - convenience: hybrid Aitken/Mellin

Until the audit (task #183) confirms Aitken extrapolation is reliable
at wall band, the `Y_exact_cube3d` convenience function may be biased.

## Efficient Foster evaluation

The 3D triple sum at N_max = 99 has 50^3 = 125,000 terms; at N_max =
399, 200^3 = 8 million.  We loop over n (outer) and vectorize over
(p, q) per inner block.

## API

    Y_DC_cube3d(L, sigma)              -> real Y_DC (S*m^2)
    K_SIBC_cube3d(L, sigma, mu)        -> real K_SIBC (leading c_0)
    Y_foster_cube3d(s, L, sigma, mu, N_max)
    Y_aitken_cube3d(s, L, sigma, mu)
    Y_mellin_cube3d(s, L, sigma, mu)
    Y_exact_cube3d(s, L, sigma, mu)    [provisional, audit pending]
"""

from __future__ import annotations

import cmath
import math
import numpy as np


def Y_DC_cube3d(L: float, sigma: float) -> float:
    """DC admittance Y_DC = L^3 sigma (S*m^2)."""
    return L**3 * sigma


def K_SIBC_cube3d(L: float, sigma: float, mu: float) -> float:
    """Leading Mellin coefficient K_SIBC = 6 L^2 sqrt(sigma/mu).

    Surface area is 6 L^2, so K_SIBC = (surface_area) * sqrt(sigma/mu).
    """
    return 6.0 * L**2 * math.sqrt(sigma / mu)


def Y_foster_cube3d(s, L: float, sigma: float, mu: float, N_max: int = 99):
    """Truncated 3D Foster sum.  WARNING: convergence is slow at high f."""
    Y_DC = Y_DC_cube3d(L, sigma)
    sMS = s * mu * sigma
    pi_L2 = math.pi**2 / L**2
    odd = np.arange(1, N_max + 1, 2, dtype=np.float64)

    # Loop over n, vectorize over (p, q).  Each n-block is len(odd)^2.
    total = 0.0 + 0j
    pp, qq = np.meshgrid(odd, odd, indexing="ij")
    pq_sq = pp**2 * qq**2
    p2_q2 = pp**2 + qq**2

    for n in odd:
        n2 = n * n
        lam = (n2 + p2_q2) * pi_L2
        denom = lam + sMS
        term = lam / (denom * (n2 * pq_sq))
        total += np.sum(term)
    return Y_DC * (512.0 / math.pi**6) * total


def Y_aitken_cube3d(s, L: float, sigma: float, mu: float):
    """Aitken epsilon-1 extrapolation of (N=99, 199, 399) Foster sums."""
    Y1 = Y_foster_cube3d(s, L, sigma, mu, 99)
    Y2 = Y_foster_cube3d(s, L, sigma, mu, 199)
    Y3 = Y_foster_cube3d(s, L, sigma, mu, 399)
    dY1 = Y2 - Y1
    dY2 = Y3 - Y2
    d2Y = dY2 - dY1
    if abs(d2Y) < abs(Y3) * 1e-13:
        return Y3
    return Y3 - dY2**2 / d2Y


def Y_mellin_cube3d(s, L: float, sigma: float, mu: float):
    """Mellin asymptote Y(s) ~ K_SIBC/sqrt(s) + c_1/s + c_2/s^1.5 .

    Valid at sufficiently high s.  Useful as the "ground truth"
    reference in the deep-SIBC regime where Foster sums become
    impractical even with Aitken.
    """
    K_SIBC = K_SIBC_cube3d(L, sigma, mu)
    c_1 = -48.0 * L / (math.pi * mu)
    c_2 = 48.0 / (math.pi * mu**1.5 * math.sqrt(sigma))
    return K_SIBC / cmath.sqrt(s) + c_1 / s + c_2 / s**1.5


def Y_exact_cube3d(s, L: float, sigma: float, mu: float):
    """Provisional ground-truth reference; uses Aitken Foster at low/mid f
    and Mellin asymptote at high f.

    The transition is at |s| ~ 1e7 rad/s (f ~ 1.6 MHz for the canonical
    L = 5 mm, Cu cube).  Below the transition, Aitken-extrapolated
    Foster sums of (N=99, 199, 399) give ~0.1% accuracy.  Above the
    transition, Foster needs much larger N_max than 399 to converge,
    while the Mellin asymptote (K_SIBC/sqrt(s) + c_1/s + c_2/s^1.5)
    becomes increasingly accurate as |s| grows.

    Note: while task #183 (3D cube Foster reference audit) is open,
    this convenience function should be considered approximate.  In
    the transition region near |s| ~ 1e7 rad/s, both methods may
    disagree by a few percent.
    """
    # Empirically calibrated transition based on the L = 5mm cube self-test:
    # Aitken and Mellin agree to ~0.5% at |s| = 6e6, diverge by 15% at
    # |s| = 6e7.  Choose the threshold conservatively at |s| = 3e7 so
    # that Mellin (which has cleaner asymptotic behavior at high f) is
    # used wherever the Aitken bias starts to matter.
    if abs(s) > 3e7:
        return Y_mellin_cube3d(s, L, sigma, mu)
    return Y_aitken_cube3d(s, L, sigma, mu)


def _self_test() -> None:
    sigma = 5.8e7
    mu = 4 * math.pi * 1e-7
    L = 5e-3

    print("=== cube3d_foster.py self-test ===")
    print(f"L = {L*1e3} mm, sigma = {sigma:.2e} S/m, mu = {mu:.4e} H/m")
    print(f"Y_DC   = {Y_DC_cube3d(L, sigma):.4e} S*m^2")
    print(f"K_SIBC = {K_SIBC_cube3d(L, sigma, mu):.4e}")
    print()

    print("REMINDER: cube3d Foster reference is UNDER AUDIT (task #183).")
    print()

    print("Foster convergence at several f (Aitken-extrapolated last column):")
    print(f"{'f (Hz)':>10}  {'N=99':>14}  {'N=199':>14}  {'N=399':>14}  {'Aitken':>14}")
    for f in [1.0, 1e3, 1e4, 1e5, 1e6, 1e7]:
        s = 1j * 2 * math.pi * f
        Yn1 = abs(Y_foster_cube3d(s, L, sigma, mu, 99))
        Yn2 = abs(Y_foster_cube3d(s, L, sigma, mu, 199))
        Yn3 = abs(Y_foster_cube3d(s, L, sigma, mu, 399))
        Ya = abs(Y_aitken_cube3d(s, L, sigma, mu))
        print(f"{f:10.2e}  {Yn1:14.6e}  {Yn2:14.6e}  {Yn3:14.6e}  {Ya:14.6e}")

    print()
    print("Mellin asymptote vs Aitken at high f:")
    print(f"{'f (Hz)':>10}  {'|Y_Aitken|':>14}  {'|Y_Mellin|':>14}  {'rel diff':>10}")
    for f in [1e6, 1e7, 1e8]:
        s = 1j * 2 * math.pi * f
        Ya = Y_aitken_cube3d(s, L, sigma, mu)
        Ym = Y_mellin_cube3d(s, L, sigma, mu)
        rel = abs(Ya - Ym) / abs(Ya) if abs(Ya) > 0 else float("nan")
        print(f"{f:10.2e}  {abs(Ya):14.6e}  {abs(Ym):14.6e}  {rel*100:9.4f}%")


if __name__ == "__main__":
    _self_test()
