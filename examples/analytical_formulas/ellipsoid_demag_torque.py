"""Rotational ellipsoid: demagnetization factor and torque.

Reference: Part 5 (SA-04-44 / RM-04-68, 2004), eq 39-44.

Two scans:
  1. N_z(c/a) over the prolate (c > a) and oblate (c < a) ranges.
     Reproduces Fig. 7 of the PDF.
  2. Torque T_z(alpha) for a fixed prolate (c/a = 5) at chi_r = 999;
     should be sinusoidal with sign convention "alpha -> 0 is
     restoring" (T_z < 0 for alpha in (0, pi/2)).
"""

from __future__ import annotations

import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))

from radia.analytical_formulas import (
    demag_factor_rotational,
    ellipsoid_torque,
)


def main() -> None:
    # ----- N_z scan over c/a -----
    ratios = np.concatenate([
        np.linspace(0.05, 1.0, 30),
        np.linspace(1.0, 20.0, 40),
    ])
    Nz_arr = np.array([demag_factor_rotational(r, 1.0)[2] for r in ratios])

    print("c/a       N_x = N_y     N_z         sum")
    print("-" * 50)
    for r in (0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0):
        Nx, _, Nz = demag_factor_rotational(r, 1.0)
        print(f"{r:6.3f}    {Nx:.6f}     {Nz:.6f}    {Nx+Nx+Nz:.6f}")

    # ----- Torque sweep -----
    a, c = 0.01, 0.05            # prolate, c/a = 5
    chi_r = 999.0
    H0 = 1000.0
    alphas = np.linspace(0.0, math.pi, 91)
    Tz = np.array([
        ellipsoid_torque(H0, alpha, chi_r=chi_r, a=a, c=c) for alpha in alphas
    ])
    print()
    print("Torque scan (c/a = 5 prolate, chi_r = 999, H_0 = 1000 A/m):")
    print(f"  T_z(alpha=0)        = {Tz[0]:.4e} N m  (expect 0)")
    print(f"  T_z(alpha=pi/4)     = {Tz[len(alphas)//4]:.4e} N m  (expect < 0)")
    print(f"  T_z(alpha=pi/2)     = {Tz[len(alphas)//2]:.4e} N m  (expect 0)")
    print(f"  T_z(alpha=3 pi / 4) = {Tz[3*len(alphas)//4]:.4e} N m  (expect > 0)")

    # Optional plot.
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("\n(matplotlib not installed -- skipping plot)")
        return
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    axes[0].plot(ratios, Nz_arr, label="N_z")
    axes[0].plot(ratios, 0.5 * (1.0 - Nz_arr), label="N_x = N_y")
    axes[0].axhline(1.0 / 3.0, color="k", lw=0.5, ls=":")
    axes[0].set_xscale("log")
    axes[0].set_xlabel("c / a")
    axes[0].set_ylabel("Demagnetization factor")
    axes[0].set_title("Demag factor of a rotational ellipsoid (eq 39-41)")
    axes[0].legend()
    axes[0].grid(True, ls=":", alpha=0.5)

    axes[1].plot(np.degrees(alphas), Tz)
    axes[1].set_xlabel("Tilt angle α [deg]")
    axes[1].set_ylabel("T_z  [N·m]")
    axes[1].set_title("Torque (prolate c/a=5, eq 44)")
    axes[1].grid(True, ls=":", alpha=0.5)
    axes[1].axhline(0.0, color="k", lw=0.5)

    fig.tight_layout()
    out = os.path.join(os.path.dirname(__file__), "ellipsoid_demag_torque.png")
    fig.savefig(out, dpi=150)
    print(f"\nSaved: {out}")


if __name__ == "__main__":
    main()
