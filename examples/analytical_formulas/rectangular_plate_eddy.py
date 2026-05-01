"""Eddy current in a thin rectangular plate.

Reference: Part 1 (SA-02-28 / RM-02-64, 2002), eq 26-27.
Source: Weissenburger and Christensen (PPPL-1517, 1979).

Reproduces Figs. 7-8 of the PDF: T_z contours = streamlines, plus J_x
along x = 0.
"""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))

from radia.analytical_formulas import plate_eddy_J, plate_eddy_T


def main() -> None:
    # Square plate as in PDF Fig. 7: 2 a = 2 b = 0.2 m, copper, dB/dt = 1 T/s.
    a = b = 0.1
    sigma = 5.0e7
    Bdot = 1.0
    n_terms = 80
    print(f"Plate: 2 a x 2 b = {2*a:.3f} m x {2*b:.3f} m")
    print(f"sigma = {sigma:.2e} S/m, dB_z/dt = {Bdot} T/s")
    print(f"Series truncation n_terms = {n_terms}")

    # Boundary integrity.
    print()
    print("Boundary check (T_z and J . n should both be small):")
    print(f"  T_z(a, 0)            = "
          f"{plate_eddy_T([a, 0.0], a, b, sigma, Bdot, n_terms=n_terms):.3e}")
    print(f"  T_z(0, b)            = "
          f"{plate_eddy_T([0.0, b], a, b, sigma, Bdot, n_terms=n_terms):.3e}")
    print(f"  J_x(a, 0.05)         = "
          f"{plate_eddy_J([a, 0.05], a, b, sigma, Bdot, n_terms=n_terms)[0]:.3e}")
    print(f"  J_y(0.05, b)         = "
          f"{plate_eddy_J([0.05, b], a, b, sigma, Bdot, n_terms=n_terms)[1]:.3e}")

    # J_x along x = 0 vs y (Fig. 8 of the PDF).
    ys = np.linspace(0.0, b, 21)
    Jx = np.array([
        plate_eddy_J([0.0, y], a, b, sigma, Bdot, n_terms=n_terms)[0]
        for y in ys
    ])
    print()
    print("J_x(0, y) along the symmetry axis x = 0:")
    print(f"{'y [m]':>8} {'J_x [A/m^2]':>16}")
    for y, j in zip(ys[::4], Jx[::4]):
        print(f"{y:>8.3f} {j:>16.4e}")

    # Plots.
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("\n(matplotlib not installed -- skipping plot)")
        return

    nx, ny = 121, 121
    x = np.linspace(-a, +a, nx)
    y = np.linspace(-b, +b, ny)
    X, Y = np.meshgrid(x, y)
    pts = np.stack([X.ravel(), Y.ravel()], axis=-1)
    T = plate_eddy_T(pts, a, b, sigma, Bdot, n_terms=n_terms).reshape(X.shape)
    J = plate_eddy_J(pts, a, b, sigma, Bdot, n_terms=n_terms).reshape(X.shape + (2,))
    Jmag = np.linalg.norm(J, axis=-1)

    fig, axes = plt.subplots(1, 2, figsize=(11, 5))

    cs = axes[0].contour(X, Y, T, levels=15, colors="k", linewidths=0.6)
    pcm = axes[0].pcolormesh(X, Y, Jmag, cmap="magma", shading="gouraud")
    axes[0].set_aspect("equal")
    axes[0].set_xlabel("x [m]")
    axes[0].set_ylabel("y [m]")
    axes[0].set_title("Eddy current streamlines (T_z contours) and |J|")
    fig.colorbar(pcm, ax=axes[0], label="|J| [A/m²]")

    axes[1].plot(ys, Jx, lw=2)
    axes[1].set_xlabel("y [m]   (x = 0)")
    axes[1].set_ylabel("J_x [A/m²]")
    axes[1].set_title("J_x along x = 0 (PDF Fig. 8)")
    axes[1].grid(True, ls=":", alpha=0.5)

    fig.tight_layout()
    out = os.path.join(os.path.dirname(__file__), "rectangular_plate_eddy.png")
    fig.savefig(out, dpi=140)
    print(f"\nSaved: {out}")


if __name__ == "__main__":
    main()
