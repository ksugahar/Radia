"""Fig_1: 2x2 panel contour+streamline plot for 3D dipole H-field.

Ported from Fig_1.m (MATLAB). Loads .mat data files produced by
3D_dipole.py and 3D_dipole_with_Kelvin.py.

Panels:
  (a) Analytical H-field (Kelvin run)
  (b) Finite domain H-field (no Kelvin)
  (c) FEM H-field (Kelvin run)
  (d) Exterior (Kelvin) domain H-field
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat


def main():
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 12.2))
    theta = np.linspace(0, 2 * np.pi, 101)

    panels = [
        ((0, 0), "3D_dipole_with_Kelvin.mat", "analytical", "(a)"),
        ((0, 1), "3D_dipole.mat", "fem", "(b)"),
        ((1, 0), "3D_dipole_with_Kelvin.mat", "fem", "(c)"),
        ((1, 1), "3D_dipole_with_Kelvin.mat", "exterior", "(d)"),
    ]

    for (row, col), matfile, mode, label in panels:
        ax = axes[row, col]
        data = loadmat(matfile)

        if mode == "analytical":
            Hx = data["Hx_analytical"]
            Hz = data["Hz_analytical"]
            xx = data["xx"]
            zz = data["zz"]
        elif mode == "fem":
            Hx = data["Hx"]
            Hz = data["Hz"]
            xx = data["xx"]
            zz = data["zz"]
        elif mode == "exterior":
            xx = data["xx_ext"] - 3
            zz = data["zz_ext"]
            Hx = data["Hx_ext"]
            Hz = data["Hz_ext"]

        H = np.sqrt(Hx**2 + Hz**2)

        ax.contourf(xx, zz, H, levels=np.linspace(0, 1.8, 30))
        ax.streamplot(
            xx[0, :], zz[:, 0], Hx, Hz,
            color="w", linewidth=0.3, density=1.5,
        )

        if mode != "exterior":
            ax.plot(0.5 * np.cos(theta), 0.5 * np.sin(theta), "k-", linewidth=0.5)
            ax.text(0, 0, r"$\mu$=100", ha="center", fontsize=14)

        ax.set_xlim(-1.1, 1.1)
        ax.set_ylim(-1.1, 1.1)
        ax.set_aspect("equal")
        ax.set_xlabel("x (m)")
        ax.set_ylabel("z (m)")
        ax.text(0.75, -1.5, label, fontsize=14)
        ax.grid(True)

    plt.tight_layout()
    plt.savefig("Fig_1.pdf")
    plt.show()


if __name__ == "__main__":
    main()
