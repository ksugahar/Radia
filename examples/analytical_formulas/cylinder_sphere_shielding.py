"""Magnetic shielding factor of cylindrical and spherical shells.

Reference: Part 1 (SA-02-28 / RM-02-64, 2002), eq 23-24.

Sweeps over relative permeability mu_r at fixed wall thickness, and
over wall thickness at fixed mu_r. Both shapes are compared on the
same axes.
"""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))

from radia.analytical_formulas import (
    shielding_factor_cylinder,
    shielding_factor_sphere,
)


def main() -> None:
    # -----------------------------------------------------------------
    # 1) Sweep mu_r at fixed thickness ratio t/b = 0.05.
    # -----------------------------------------------------------------
    b = 1.0
    t_over_b = 0.05
    a = b * (1.0 - t_over_b)
    mu_rs = np.geomspace(1.0, 1.0e5, 60)
    S_cyl = np.array([shielding_factor_cylinder(a, b, mu) for mu in mu_rs])
    S_sph = np.array([shielding_factor_sphere(a, b, mu) for mu in mu_rs])

    print("S vs mu_r (t/b = 0.05)")
    print(f"{'mu_r':>10} {'S_cyl':>12} {'S_sph':>12}")
    for k in (0, 10, 20, 30, 40, 50, 59):
        print(f"{mu_rs[k]:>10.2f} {S_cyl[k]:>12.4e} {S_sph[k]:>12.4e}")

    # -----------------------------------------------------------------
    # 2) Sweep wall thickness at fixed mu_r = 1000.
    # -----------------------------------------------------------------
    mu_r = 1000.0
    ts = np.linspace(0.001, 0.5, 50)
    S_cyl_t = np.array([shielding_factor_cylinder(b * (1 - t), b, mu_r) for t in ts])
    S_sph_t = np.array([shielding_factor_sphere(b * (1 - t), b, mu_r) for t in ts])

    print()
    print("S vs t/b (mu_r = 1000)")
    print(f"{'t/b':>8} {'S_cyl':>12} {'S_sph':>12}")
    for k in (0, 10, 20, 30, 40, 49):
        print(f"{ts[k]:>8.4f} {S_cyl_t[k]:>12.4e} {S_sph_t[k]:>12.4e}")

    # -----------------------------------------------------------------
    # Optional plots.
    # -----------------------------------------------------------------
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("\n(matplotlib not installed -- skipping plot)")
        return

    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    axes[0].loglog(mu_rs, S_cyl, label="cylinder", lw=2)
    axes[0].loglog(mu_rs, S_sph, label="sphere", lw=2, ls="--")
    axes[0].set_xlabel("Relative permeability μ_r")
    axes[0].set_ylabel("Shielding factor S")
    axes[0].set_title(f"S vs μ_r  (t/b = {t_over_b:.2f})")
    axes[0].grid(True, ls=":", alpha=0.5, which="both")
    axes[0].legend()

    axes[1].semilogy(ts, S_cyl_t, label="cylinder", lw=2)
    axes[1].semilogy(ts, S_sph_t, label="sphere", lw=2, ls="--")
    axes[1].set_xlabel("Wall thickness ratio  t / b")
    axes[1].set_ylabel("Shielding factor S")
    axes[1].set_title(f"S vs t/b  (μ_r = {mu_r:.0f})")
    axes[1].grid(True, ls=":", alpha=0.5, which="both")
    axes[1].legend()

    fig.tight_layout()
    out = os.path.join(os.path.dirname(__file__), "cylinder_sphere_shielding.png")
    fig.savefig(out, dpi=150)
    print(f"\nSaved: {out}")


if __name__ == "__main__":
    main()
