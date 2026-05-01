"""Cross-validate the 2D rectangular bar against the 3D cuboid in the
long-bar limit.

The 2D bar formula (Part 2, eq 2-3) is the ``L_z -> infty`` limit of
the 3D ``CuboidMagnet`` field in :mod:`radia.analytical_magnet`. Both
routines accept ``M`` in A/m and return ``B`` in Tesla, so a direct
component-wise compare is meaningful.

The script sweeps the z-extent ``L_z`` of the 3D cuboid and shows the
component-wise convergence rate of ``B_3d(z = 0) - B_2d``. Convergence
should look like ``(max(a, b) / L_z)**2``.
"""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))

from radia.analytical_formulas import rect_magnet_2d_B
from radia.analytical_formulas.rect_magnet_2d import MU_0
from radia.analytical_magnet import CuboidMagnet


def main() -> None:
    a, b = 0.02, 0.01
    Mx = 1.0e6
    M_3d = [Mx, 0.0, 0.0]

    # Pick observation points along z = 0 outside the cross-section.
    obs_2d = np.array([[0.05, 0.005], [0.10, 0.0], [0.05, -0.02]])

    print(f"2D bar: 2 a x 2 b = {2*a:.3f} x {2*b:.3f} m, M = ({Mx:.1e}, 0) A/m")
    print(f"mu_0 * M_x = {MU_0 * Mx:.4f} T")
    print()

    B_2d_list = [rect_magnet_2d_B(p, a, b, M=(Mx, 0.0)) for p in obs_2d]

    aspects = [3.0, 10.0, 30.0, 100.0, 300.0]
    header = (
        "L_z / max(a,b)" + "".join(f"   |dB|@P{k}" for k in range(len(obs_2d)))
    )
    print(header)
    for L_over_d in aspects:
        Lz = L_over_d * max(a, b)
        cuboid = CuboidMagnet(
            center=[0.0, 0.0, 0.0],
            dimensions=[2 * a, 2 * b, 2 * Lz],
            magnetization=M_3d,
        )
        deltas = []
        for k, (p2, B2) in enumerate(zip(obs_2d, B_2d_list)):
            B3 = np.array(cuboid.get_B([p2[0], p2[1], 0.0]))
            d = float(np.linalg.norm(B3[:2] - B2))
            deltas.append(d)
        ds = "  ".join(f"{d:.4e}" for d in deltas)
        print(f"  {L_over_d:>10.1f}     {ds}")

    # Optional plot.
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("\n(matplotlib not installed -- skipping plot)")
        return
    aspects_plot = np.geomspace(2.0, 1000.0, 30)
    err_per_pt = np.zeros((len(obs_2d), aspects_plot.size))
    for j, asp in enumerate(aspects_plot):
        Lz = asp * max(a, b)
        cuboid = CuboidMagnet(
            center=[0.0, 0.0, 0.0],
            dimensions=[2 * a, 2 * b, 2 * Lz],
            magnetization=M_3d,
        )
        for k, (p2, B2) in enumerate(zip(obs_2d, B_2d_list)):
            B3 = np.array(cuboid.get_B([p2[0], p2[1], 0.0]))
            err_per_pt[k, j] = float(np.linalg.norm(B3[:2] - B2))

    fig, ax = plt.subplots(figsize=(7, 4.5))
    for k in range(len(obs_2d)):
        ax.loglog(
            aspects_plot, err_per_pt[k],
            label=f"obs P{k} = ({obs_2d[k][0]:.2f}, {obs_2d[k][1]:.3f})",
            lw=1.6,
        )
    # 1 / L_z**2 reference line.
    ax.loglog(
        aspects_plot, 0.05 / aspects_plot ** 2, "k--", lw=1.0,
        label="$\\propto 1 / L_z^2$",
    )
    ax.set_xlabel("L_z / max(a, b)")
    ax.set_ylabel("|B_3d(z=0) - B_2d|  [T]")
    ax.set_title("Convergence of 3D cuboid -> 2D bar at z = 0")
    ax.grid(True, ls=":", alpha=0.5, which="both")
    ax.legend()
    fig.tight_layout()
    out = os.path.join(os.path.dirname(__file__), "cross_validation_3d_vs_2d.png")
    fig.savefig(out, dpi=140)
    print(f"\nSaved: {out}")


if __name__ == "__main__":
    main()
