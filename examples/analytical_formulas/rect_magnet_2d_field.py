"""2D field of a uniformly magnetised rectangular bar.

Reference: Part 2 (SA-03-08 / RM-03-08, 2003), eq 2-3.

Plots vector potential contours and B-field streamlines around a 2 a x 2 b
bar magnetised along x. Far-field 2D dipole law is verified along the
horizontal axis.
"""

from __future__ import annotations

import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))

from radia.analytical_formulas import rect_magnet_2d_A, rect_magnet_2d_B
from radia.analytical_formulas.rect_magnet_2d import MU_0


def main() -> None:
    # 4 cm x 2 cm bar, M_x = 1.2 / mu_0 ~ 9.55e5 A/m (typical NdFeB).
    a, b = 0.02, 0.01
    M = (1.2 / MU_0, 0.0)
    print(f"Bar: 2a x 2b = {2*a*1000:.1f} mm x {2*b*1000:.1f} mm, "
          f"M = ({M[0]:.3e}, {M[1]:.3e}) A/m")
    print(f"mu_0 * M_x = {MU_0 * M[0]:.4f} T")

    # Demag-factor check at center.
    B_c = rect_magnet_2d_B([0.0, 0.0], a, b, M=M)
    print(f"B at center: ({B_c[0]:.4f}, {B_c[1]:.4e}) T")
    print(f"  Cross-section demag along x:  N_x = "
          f"{1.0 - B_c[0] / (MU_0 * M[0]):.4f}")

    # Far-field 2D dipole comparison along x-axis.
    rs = np.array([0.05, 0.1, 0.3, 1.0, 5.0])
    print(f"\n  r [m]    |B|_analytic [T]    |B|_dipole-law [T]    rel err")
    m_2D = M[0] * 4.0 * a * b
    for r in rs:
        B = rect_magnet_2d_B(np.array([r, 0.0]), a, b, M=M)
        Bmag = float(np.linalg.norm(B))
        # 2D dipole on its axis (M along x, point on +x axis):
        # |B_x| = mu_0 * m / (2 pi r**2)
        Bd = MU_0 * m_2D / (2.0 * math.pi * r * r)
        print(f"  {r:6.3f}   {Bmag:.4e}        {Bd:.4e}        "
              f"{abs(Bmag - Bd) / Bd:.3e}")

    # Optional 2D plot.
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("\n(matplotlib not installed -- skipping plot)")
        return

    grid = np.linspace(-0.05, 0.05, 240)
    X, Y = np.meshgrid(grid, grid)
    points = np.stack([X.ravel(), Y.ravel()], axis=-1)
    A = rect_magnet_2d_A(points, a, b, M=M).reshape(X.shape)
    B = rect_magnet_2d_B(points, a, b, M=M).reshape(X.shape + (2,))
    Bmag = np.linalg.norm(B, axis=-1)

    fig, ax = plt.subplots(figsize=(7, 6))
    cs = ax.contour(X, Y, A, levels=21, colors="k", linewidths=0.6, alpha=0.7)
    pcm = ax.pcolormesh(X, Y, np.minimum(Bmag, MU_0 * abs(M[0])), shading="gouraud",
                        cmap="viridis", alpha=0.65)
    ax.streamplot(X, Y, B[..., 0], B[..., 1], color="white",
                  density=1.5, linewidth=0.7, arrowsize=0.7)
    rect_x = [-a, +a, +a, -a, -a]
    rect_y = [-b, -b, +b, +b, -b]
    ax.plot(rect_x, rect_y, "r-", lw=1.6)
    ax.set_aspect("equal")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_title("2D rectangular bar (M along x): A_z contours + B streamlines")
    fig.colorbar(pcm, ax=ax, label="|B| [T] (clipped at μ₀|M|)")
    fig.tight_layout()
    out = os.path.join(os.path.dirname(__file__), "rect_magnet_2d_field.png")
    fig.savefig(out, dpi=140)
    print(f"\nSaved: {out}")


if __name__ == "__main__":
    main()
