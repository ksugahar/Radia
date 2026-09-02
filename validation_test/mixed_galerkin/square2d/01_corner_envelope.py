"""
2D square cross-section mixed Galerkin: rank-1 bulk + corner-aware envelope.

Phase 5 v2 of the 2026-05-28 -> 2026-06-12 research sprint, with Phase
8c corrected reference (Aitken-extrapolated Foster).

The 2-D square is the natural corner-bearing extension of the
cylinder (cylindrical cross-section is corner-free).  An earlier failed
v1 attempt used psi = f(x) sin(pi y/L) + sin(pi x/L) f(y), a smooth
y-factor that lacks the y-direction skin layer.  That dead source was
deleted; the lesson is recorded in memory/mixed_galerkin_examples_prune.md.
The corrected corner-aware envelope is

    psi(x, y; s) = f(x) f(y)

with f(x) = cosh((x - L/2) t) / cosh(L t / 2) - 1 captures BOTH the
x- and y-direction boundary layers AND the right Wiener-Hopf
asymptote v ~ r^2 sin(2 theta) at each 90 degree corner (where
r, theta are local polar coordinates).

## Phase 8c corrected result

The original measurement against truncated Foster (N_max = 1999) gave
"0.99% wall band peak at f = 7.94e5 Hz".  The Aitken-extrapolated
reference revealed that 97% of that "peak" was Foster truncation
artifact; the true mixed Galerkin v2 error is in the 0.03 - 0.26%
range across the spectrum.

    f (Hz)        original report   corrected (Aitken)
    1e3                -            0.070%
    1e4 (true max)     -            0.263%
    5e4                -            0.103%
    1e5                -            0.071%
    7.94e5 (old peak)  0.99%        0.034%   <-- artifact
    1e6                -            0.035%

The remaining residual is the corner Mellin c_2 contribution (90deg
dihedral wedge function, not captured by tensor-product envelope at
leading order).  Adding a wedge basis function would push this down
further.
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

N_PSI = 199  # sine modes for psi expansion per axis


def f_n(n: int, t):
    """sin(n pi x/L) coefficient of f(x; t), n odd."""
    return -4 * L**2 * t**2 / (n * math.pi * (n**2 * math.pi**2 + L**2 * t**2))


def F_0(t):
    """int_0^L f(x; t) dx."""
    return (2.0 / t) * cmath.tanh(L * t / 2) - L


def Y_mixed_galerkin(s):
    """rank-1 bulk (sin sin (1,1) mode) + corner envelope psi = f(x) f(y)."""
    t = cmath.sqrt(s * MS)

    odd_n = np.arange(1, 2 * N_PSI, 2)
    f_arr = np.array([f_n(int(n), t) for n in odd_n], dtype=complex)

    # Surface block: K_ss = sum (lambda + s MS) (L/2)^2 (f_n f_p)^2
    nn, pp = np.meshgrid(odd_n, odd_n, indexing="ij")
    lam_2d = (nn**2 + pp**2) * math.pi**2 / L**2
    d_2d = np.outer(f_arr, f_arr)
    K_ss = (L / 2)**2 * np.sum((lam_2d + s * MS) * d_2d**2)

    b_psi = F_0(t)**2

    # Bulk: phi_{1,1} = sin(pi x/L) sin(pi y/L), only the (1, 1) mode of
    # the surface envelope expansion couples to it.
    lambda_11 = 2 * math.pi**2 / L**2
    K_b = (lambda_11 + s * MS) * (L / 2)**2
    K_bs = (lambda_11 + s * MS) * (L / 2)**2 * f_arr[0]**2
    b_b = (2 * L / math.pi)**2  # = int sin(pi x/L) sin(pi y/L) dx dy

    K_mat = np.array([[K_b, K_bs], [K_bs, K_ss]], dtype=complex)
    b_vec = np.array([b_b, b_psi], dtype=complex)
    xi = np.linalg.solve(K_mat, -s * MS * b_vec)
    v_avg = (xi @ b_vec) / L**2
    return Y_DC * (1 + v_avg)


SWEEP = np.logspace(0, 6, 31)
WALL_BAND_HZ = (1e4, 1e6)
METRIC = "abs(Y_exact - Y_mixed) / abs(Y_exact)"


def summary() -> dict:
    """Run the sweep and return the numbers. main() prints from this."""
    errs = []
    for f in SWEEP:
        s = 1j * 2 * math.pi * f
        Y_e = Y_exact_square2d(s, L, SIGMA, MU)
        errs.append(abs(Y_e - Y_mixed_galerkin(s)) / abs(Y_e))
    errs = np.array(errs)
    wall = (SWEEP > WALL_BAND_HZ[0]) & (SWEEP < WALL_BAND_HZ[1])
    return {
        "case": "square2d_corner_envelope",
        "description": "rank-1 CLN bulk + tensor corner-aware surface envelope",
        "metric": METRIC,
        "geometry": {"L_m": L, "sigma_S_per_m": SIGMA, "mu_H_per_m": MU},
        "sweep": {"f_lo_hz": float(SWEEP[0]), "f_hi_hz": float(SWEEP[-1]),
                  "n_points": int(SWEEP.size),
                  "wall_band_hz": list(WALL_BAND_HZ)},
        "reference": "Aitken-extrapolated Foster sum "
                     "(_references/square2d_foster.py)",
        "max_error_pct": float(errs.max() * 100),
        "max_error_at_hz": float(SWEEP[errs.argmax()]),
        "wall_band_max_error_pct": float(errs[wall].max() * 100),
    }


def main():
    print("=== 2D square cross-section: rank-1 bulk + corner envelope ===")
    print(f"L = {L*1e3} mm, sigma = {SIGMA:.2e} S/m, mu = {MU:.4e} H/m")
    print(f"Y_DC   = {Y_DC:.4e} S*m (per unit length z)")
    print(f"K_SIBC = {K_SIBC:.4e}")
    print()
    print("(reference: Aitken-extrapolated Foster sum; see _references/square2d_foster.py)")
    print()

    print(f"{'f (Hz)':>10}  {'|Y_exact|':>12}  {'|Y_mixed|':>12}  {'rel err':>10}")
    for f in [1.0, 1e3, 1e4, 5e4, 1e5, 7.94e5, 1e6]:
        s = 1j * 2 * math.pi * f
        Y_e = Y_exact_square2d(s, L, SIGMA, MU)
        Y_m = Y_mixed_galerkin(s)
        err = abs(Y_e - Y_m) / abs(Y_e)
        print(f"{f:10.2e}  {abs(Y_e):12.4e}  {abs(Y_m):12.4e}  {err*100:9.4f}%")

    print()
    print(f"Sub-sweep summary (1 Hz to 1e6 Hz, {SWEEP.size} points):")
    r = summary()
    print(f"  Max error anywhere:  {r['max_error_pct']:.4f}% "
          f"@ f = {r['max_error_at_hz']:.2e} Hz")
    print(f"  Wall band max:       {r['wall_band_max_error_pct']:.4f}%")


if __name__ == "__main__":
    main()
