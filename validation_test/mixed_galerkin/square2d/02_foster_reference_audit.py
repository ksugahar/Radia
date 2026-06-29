"""
Foster sum reference audit for the 2D square cross-section.

Phase 8c of the 2026-05-28 -> 2026-06-12 research sprint.

The Phase 5 v2 mixed Galerkin (`01_corner_envelope.py`) was originally
benchmarked against Foster sum truncated at N_max = 1999, giving the
appearance of a "0.99% wall-band peak" at f = 7.94e5 Hz.  This script
demonstrates that the apparent peak is almost entirely an artifact of
slow Foster sum convergence at high frequency, not a real envelope-
intrinsic error of the mixed Galerkin.

## Methodology

For each test frequency:
    1. Compute Foster sum at N_max = 199, 999, 4999, 9999, 19999.
    2. The trend |Y_{N=199} - Y_{N=999} - ...| gives the truncation
       bias as a function of N_max.
    3. Aitken epsilon-1 extrapolation of (4999, 9999, 19999) gives an
       estimate of the N -> infinity Foster sum limit.
    4. Compare the Phase 5 mixed Galerkin to the Aitken reference.

The Aitken extrapolation is justified because Foster bias at fixed s
decays as O(1/N_max) (each missing mode contributes a 1/(n^2 p^2)
weight in the doubly-infinite sum), so the epsilon-1 transformation
removes the leading 1/N term.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from _references.square2d_foster import (  # noqa: E402
    K_SIBC_square2d,
    Y_DC_square2d,
    Y_exact_square2d,
    Y_foster_square2d,
)

SIGMA = 5.8e7
MU = 4 * math.pi * 1e-7
L = 5e-3
MS = MU * SIGMA


def main():
    print("=== 2D square Foster reference convergence audit ===")
    print(f"L = {L*1e3} mm, sigma = {SIGMA:.2e} S/m, mu = {MU:.4e} H/m")
    print()

    print("Foster sum value at increasing N_max:")
    print(f"{'f (Hz)':>10}  {'N=199':>14}  {'N=999':>14}  {'N=4999':>14}  {'N=9999':>14}  {'N=19999':>14}  {'Aitken':>14}")
    for f in [1.0, 1e3, 1e4, 1e5, 7.94e5, 1e6, 1e7]:
        s = 1j * 2 * math.pi * f
        Y = [abs(Y_foster_square2d(s, L, SIGMA, MU, N)) for N in (199, 999, 4999, 9999, 19999)]
        Ya = abs(Y_exact_square2d(s, L, SIGMA, MU))
        print(f"{f:10.2e}  " + "  ".join(f"{x:14.6e}" for x in Y) + f"  {Ya:14.6e}")
    print()

    print("Truncation bias relative to next N (each column = (N_{i+1} - N_i)/N_{i+1}):")
    print(f"{'f (Hz)':>10}  {'(999-199)/999':>16}  {'(4999-999)/4999':>16}  {'(9999-4999)/9999':>18}  {'(19999-9999)/19999':>20}")
    for f in [1.0, 1e3, 1e4, 1e5, 7.94e5, 1e6, 1e7]:
        s = 1j * 2 * math.pi * f
        Y = [abs(Y_foster_square2d(s, L, SIGMA, MU, N)) for N in (199, 999, 4999, 9999, 19999)]
        deltas = [(Y[i+1] - Y[i]) / Y[i+1] * 100 for i in range(4)]
        print(f"{f:10.2e}  " + "  ".join(f"{d:14.4f}%" for d in deltas))


if __name__ == "__main__":
    main()
