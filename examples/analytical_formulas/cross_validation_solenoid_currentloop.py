"""Cross-validate the closed-form solenoid axial field against a stacked
set of analytical current loops.

The closed-form `solenoid_axial_field` (Part 4 §4) is the integral of
Biot-Savart over the rectangular cross-section. Numerically the same
quantity is obtained by summing many circular current loops -- here
using `radia.analytical_magnet.CurrentLoop` (Ortner 2022 formulation,
elliptic integrals). Convergence is `O(dr**2 + dz**2)` (midpoint rule).

Both routines accept SI units (m, A, A/m**2, T) so component-wise
compare is meaningful. This is the second cross-validation script,
after `cross_validation_3d_vs_2d.py`.
"""

from __future__ import annotations

import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))

from radia.analytical_formulas import solenoid_axial_field
from radia.analytical_formulas.solenoid_central import MU_0
from radia.analytical_magnet import CurrentLoop


def stacked_loops_axial_Bz(z_obs, a1, a2, b, J, n_r, n_z):
    dr = (a2 - a1) / n_r
    dz = (2.0 * b) / n_z
    Bz = 0.0
    for i in range(n_r):
        r = a1 + (i + 0.5) * dr
        for k in range(n_z):
            z_ring = -b + (k + 0.5) * dz
            dI = J * dr * dz
            loop = CurrentLoop(
                center=[0.0, 0.0, z_ring],
                diameter=2.0 * r,
                current=dI,
                axis="z",
                units="m",
            )
            Bz += loop.get_B([0.0, 0.0, z_obs])[2]
    return Bz


def main() -> None:
    a1, a2, b, J = 0.05, 0.10, 0.10, 1.0e7
    print(f"Solenoid: a1 = {a1*1000:.0f} mm, a2 = {a2*1000:.0f} mm, "
          f"b = {b*1000:.0f} mm, J = {J:.2e} A/m**2")
    print(f"Long-coil limit B_inf = {MU_0 * J * (a2 - a1):.4f} T")
    print()

    # ----- Convergence study at z = 0 -------------------------------
    z_obs = 0.0
    Bz_ref = float(solenoid_axial_field(z_obs, a1, a2, b, J))
    print(f"Closed form B_z(z=0) = {Bz_ref:.6f} T")
    print()
    print(f"{'(n_r, n_z)':>14} {'B_z stacked':>14} {'rel err':>12} {'time [s]':>10}")
    grids = [(5, 20), (10, 40), (20, 80), (40, 160)]
    Bs, errs, times = [], [], []
    for n_r, n_z in grids:
        t0 = time.perf_counter()
        Bz = stacked_loops_axial_Bz(z_obs, a1, a2, b, J, n_r, n_z)
        dt = time.perf_counter() - t0
        rel_err = abs(Bz - Bz_ref) / Bz_ref
        Bs.append(Bz)
        errs.append(rel_err)
        times.append(dt)
        print(f"  ({n_r:3d}, {n_z:3d})    {Bz:.6f}    {rel_err:.3e}    {dt:7.3f}")
    print()

    # Confirm midpoint rule converges as O(h**2): doubling -> ~ 4x reduction.
    if len(errs) >= 3:
        ratio = errs[-2] / errs[-1]
        print(f"Last grid refinement halving error ratio: {ratio:.2f} "
              f"(O(h**2) ideal: 4)")

    # ----- Axial profile compared, fixed grid ------------------------
    print()
    print("Axial profile comparison (n_r=20, n_z=80):")
    print(f"{'z [m]':>8} {'closed':>10} {'stacked':>12} {'rel err':>12}")
    z_grid = np.linspace(-0.2, 0.2, 9)
    Bz_closed = solenoid_axial_field(z_grid, a1, a2, b, J)
    Bz_num = np.array([
        stacked_loops_axial_Bz(z, a1, a2, b, J, 20, 80) for z in z_grid
    ])
    for z, bc, bn in zip(z_grid, Bz_closed, Bz_num):
        rel = abs(bn - bc) / max(abs(bc), 1e-30)
        print(f"  {z:+.3f}   {bc:.5f}    {bn:.5f}     {rel:.3e}")

    # Plot.
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("\n(matplotlib not installed -- skipping plot)")
        return

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    z_dense = np.linspace(-0.3, 0.3, 121)
    axes[0].plot(z_dense, solenoid_axial_field(z_dense, a1, a2, b, J),
                 lw=2, label="closed form")
    axes[0].plot(z_grid, Bz_num, "o", ms=6, label="20 x 80 CurrentLoops")
    axes[0].set_xlabel("z [m]")
    axes[0].set_ylabel(r"$B_z$ [T]")
    axes[0].set_title("Solenoid axial field: closed form vs stacked rings")
    axes[0].legend()
    axes[0].grid(True, ls=":", alpha=0.5)

    nrs = [g[0] for g in grids]
    axes[1].loglog(nrs, errs, "o-", lw=1.7)
    axes[1].set_xlabel("$n_r$ (axial slices = $4 n_r$)")
    axes[1].set_ylabel("relative error at $z = 0$")
    axes[1].set_title("Midpoint-rule convergence (expect ~ 1 / $n_r^2$)")
    axes[1].grid(True, ls=":", alpha=0.5, which="both")

    fig.tight_layout()
    out = os.path.join(os.path.dirname(__file__),
                       "cross_validation_solenoid_currentloop.png")
    fig.savefig(out, dpi=140)
    print(f"\nSaved: {out}")


if __name__ == "__main__":
    main()
