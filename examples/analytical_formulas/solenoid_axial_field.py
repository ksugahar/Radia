"""Solenoid axial-field profile and Fabri form factor.

Reference: Part 4 (SA-04-1 / RM-04-1, 2004), eq 25-29.
            Montgomery, "Solenoid Magnet Design" (Wiley, 1969).

Plots ``B_z(z)`` along the symmetry axis of a finite rectangular-
section solenoid and reports the central-field Fabri factor
``F(alpha, beta)``. Verifies the long-coil limit
``B_0 -> mu_0 J (a_2 - a_1)``.
"""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))

from radia.analytical_formulas import (
    fabri_F,
    solenoid_axial_field,
    solenoid_central_field,
)
from radia.analytical_formulas.solenoid_central import MU_0


def main() -> None:
    a1, a2 = 0.05, 0.10        # 50 mm inner, 100 mm outer
    J = 1.0e7                  # 10 MA / m**2
    print(f"Geometry: a_1 = {a1*1000:.0f} mm, a_2 = {a2*1000:.0f} mm, "
          f"J = {J:.2e} A/m**2")
    print(f"Long-coil limit B_inf = mu_0 J (a_2 - a_1) = "
          f"{MU_0 * J * (a2 - a1):.4f} T")
    print()
    print("b [m]   alpha   beta    F(alpha, beta)   B_0 [T]   B_0 / B_inf")
    for b in (0.001, 0.01, 0.05, 0.1, 0.5, 5.0):
        alpha = a2 / a1
        beta = b / a1
        F = fabri_F(alpha, beta)
        B0 = solenoid_central_field(a1, a2, b, J)
        ratio = B0 / (MU_0 * J * (a2 - a1))
        print(f"  {b:.3f}   {alpha:.3f}   {beta:6.3f}   {F:.5f}        "
              f"{B0:.4f}   {ratio:.4f}")

    # Axial profile.
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("\n(matplotlib not installed -- skipping plot)")
        return

    fig, ax = plt.subplots(figsize=(7, 4.5))
    z = np.linspace(-0.4, 0.4, 401)
    for b in (0.05, 0.1, 0.2):
        Bz = solenoid_axial_field(z, a1, a2, b, J)
        ax.plot(z, Bz, lw=1.6, label=f"b = {b*1000:.0f} mm")
    ax.axhline(MU_0 * J * (a2 - a1), color="k", ls="--", lw=0.8,
               label=r"$\mu_0 J (a_2 - a_1)$")
    ax.set_xlabel("z [m]")
    ax.set_ylabel(r"$B_z$ [T]")
    ax.set_title("Axial field of a rectangular-section solenoid")
    ax.legend()
    ax.grid(True, ls=":", alpha=0.5)
    fig.tight_layout()
    out = os.path.join(os.path.dirname(__file__), "solenoid_axial_field.png")
    fig.savefig(out, dpi=140)
    print(f"\nSaved: {out}")


if __name__ == "__main__":
    main()
