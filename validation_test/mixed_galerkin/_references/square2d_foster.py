"""
Reference Y_exact for an infinite-z conducting bar of square cross-section L x L.

There is NO classical Bessel-style closed form for the 2D rectangular
cross-section.  The natural representation is the Foster eigenmode
expansion (a doubly-infinite series over odd sine modes); for use as a
"ground truth" reference, this module provides:

    1. Y_foster(s, N_max)   - truncated Foster sum
    2. Y_aitken(s)          - Aitken-extrapolated reference (3 Foster
                              sums at N = 4999, 9999, 19999)

The Aitken version is what mixed Galerkin scripts should compare
against.  Plain truncated Foster has slow convergence at high
frequency and produced misleading "wall-band peak" artifacts in
Phase 5 (see ## Phase 8c artifact below).

## Formula

Variational problem on the cross-section [0, L]^2 (infinite z bar):

    (-Laplacian_xy + s mu sigma) v(x,y) = -s mu sigma
    v = 0 on the 4 edges
    Y(s) / length = L^2 sigma * (1 + <v>)

with Y_DC = L^2 sigma (per unit length) and K_SIBC = 4 L sqrt(sigma/mu)
(per unit length; perimeter is 4L).

Sine-mode expansion of v gives the Foster sum:

    Y/Y_DC = sum_{n,p odd} (64 / pi^4) * lambda_{np} / (lambda_{np} + s mu sigma) / (n^2 p^2)

with lambda_{np} = (n^2 + p^2) pi^2 / L^2.  At DC the sum identity
gives Y/Y_DC = 1 (Parseval); at s -> infinity it gives K_SIBC/sqrt(s)
plus corner-Mellin corrections.

## Convergence rate

The series has the structure of a 2D Fourier sum with weights
1/(n^2 p^2) decaying only quadratically.  At a given s the truncation
at N = N_max odd modes per axis leaves residual modes (n, p > N_max)
with lambda_{np} >> s mu sigma; each such mode contributes its full
DC weight (64/pi^4)/(n^2 p^2), so the truncation bias DECREASES with
N_max but the trailing-modes contribution scales like 1/N_max.

Empirically (Phase 8c):

    N_max | bias at wall band peak (f = 7.94e5 Hz)
    -----+---------------------------------
     199 |  10%
     999 |   5%
    4999 |   1%
    9999 |   0.15%
   19999 |   0.07%

So N_max = 19999 gives ~0.07% absolute bias at wall band, and the
Aitken extrapolation of (4999, 9999, 19999) gets down to roughly the
3rd-decimal level — adequate for benchmarking 1-DOF mixed Galerkin
that itself is at the 0.1% level.

## Phase 8c artifact (2026-06-12 - see memory project_mixed_galerkin_schur_dtn)

The original 2D-square mixed Galerkin script (Phase 5, both v1 and
v2) used Foster truncated at N_max ~ 1999, which had ~5% bias at the
wall band peak.  This bias was misattributed to "envelope intrinsic
error" of the mixed Galerkin: in particular the "v2 corner envelope
gives 0.99% wall-band peak" claim was 97% Foster artifact.  The true
v2 mixed Galerkin accuracy against the Aitken-extrapolated reference
is 0.03-0.26% across the spectrum.

## API

    Y_DC_square2d(L, sigma)                     -> real Y_DC per length
    K_SIBC_square2d(L, sigma, mu)               -> real K_SIBC per length
    Y_foster_square2d(s, L, sigma, mu, N_max)   -> complex Y(s) per length
    Y_exact_square2d(s, L, sigma, mu)           -> Aitken-extrapolated reference

## Self-test

Run as `python square2d_foster.py` for a convergence check at DC,
wall band, and high f.
"""

from __future__ import annotations

import cmath
import math
import numpy as np


def Y_DC_square2d(L: float, sigma: float) -> float:
    """DC admittance per unit length: Y_DC = L^2 sigma (S*m)."""
    return L**2 * sigma


def K_SIBC_square2d(L: float, sigma: float, mu: float) -> float:
    """Leading SIBC coefficient per unit length: K_SIBC = 4 L sqrt(sigma/mu).

    Perimeter is 4L, so K_SIBC = (perimeter) * sqrt(sigma/mu).
    """
    return 4.0 * L * math.sqrt(sigma / mu)


def Y_foster_square2d(s, L: float, sigma: float, mu: float, N_max: int = 999):
    """Foster eigenmode sum for the L x L bar, truncated at N_max odd modes per axis.

    Notes
    -----
    Convergence is slow at high f; use Y_exact_square2d() for benchmarking.
    """
    odd_n = np.arange(1, N_max + 1, 2)
    nn, pp = np.meshgrid(odd_n, odd_n, indexing="ij")
    lam = (nn**2 + pp**2) * math.pi**2 / L**2
    term = (64.0 / math.pi**4) * lam / (lam + s * mu * sigma) / (nn**2 * pp**2)
    Y_DC = L**2 * sigma
    return Y_DC * np.sum(term)


def Y_exact_square2d(s, L: float, sigma: float, mu: float):
    """Aitken-extrapolated Foster sum, used as the ground-truth reference.

    Uses Foster at N = (4999, 9999, 19999) and applies the Aitken
    epsilon-1 transformation to estimate the N -> infinity limit.

    For very near-DC (|s mu sigma L^2 / pi^2| < ~1) the truncation
    bias is already tiny and the Aitken step degenerates; in that
    case we fall back on the N=19999 partial sum directly.
    """
    Y1 = Y_foster_square2d(s, L, sigma, mu, 4999)
    Y2 = Y_foster_square2d(s, L, sigma, mu, 9999)
    Y3 = Y_foster_square2d(s, L, sigma, mu, 19999)
    dY1 = Y2 - Y1
    dY2 = Y3 - Y2
    d2Y = dY2 - dY1
    if abs(d2Y) < abs(Y3) * 1e-13:
        return Y3
    return Y3 - dY2**2 / d2Y


def _self_test() -> None:
    sigma = 5.8e7
    mu = 4 * math.pi * 1e-7
    L = 5e-3
    MS = mu * sigma

    print("=== square2d_foster.py self-test ===")
    print(f"L = {L*1e3} mm, sigma = {sigma:.2e} S/m, mu = {mu:.4e} H/m")
    print(f"Y_DC   = {Y_DC_square2d(L, sigma):.4e} S*m (per unit length)")
    print(f"K_SIBC = {K_SIBC_square2d(L, sigma, mu):.4e}")
    print()

    print("Foster convergence at several f (compare to Aitken-extrapolated):")
    print(f"{'f (Hz)':>10}  {'N=999':>14}  {'N=9999':>14}  {'N=19999':>14}  {'Aitken':>14}")
    for f in [1.0, 1e3, 1e4, 1e5, 7.94e5, 1e6, 1e7]:
        s = 1j * 2 * math.pi * f
        Yn1 = abs(Y_foster_square2d(s, L, sigma, mu, 999))
        Yn2 = abs(Y_foster_square2d(s, L, sigma, mu, 9999))
        Yn3 = abs(Y_foster_square2d(s, L, sigma, mu, 19999))
        Ya = abs(Y_exact_square2d(s, L, sigma, mu))
        print(f"{f:10.2e}  {Yn1:14.6e}  {Yn2:14.6e}  {Yn3:14.6e}  {Ya:14.6e}")

    print()
    print("DC limit check (Aitken should converge to Y_DC = L^2 sigma):")
    Y_DC = Y_DC_square2d(L, sigma)
    s_dc = 1j * 2 * math.pi * 1.0  # ~static
    Y_aitken = Y_exact_square2d(s_dc, L, sigma, mu)
    print(f"  Y_DC theory = {Y_DC:.6e}")
    print(f"  Y_Aitken    = {abs(Y_aitken):.6e}")
    print(f"  rel diff    = {abs(abs(Y_aitken) - Y_DC) / Y_DC * 100:.4f}%")


if __name__ == "__main__":
    _self_test()
