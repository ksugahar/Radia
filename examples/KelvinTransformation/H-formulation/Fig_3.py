"""Fig_3: 1x2 panel line plot comparing Hz along x-axis.

Ported from Fig_3.m (MATLAB). Loads .mat data files produced by
3D_dipole.py, 3D_dipole_with_Kelvin.py,
3D_quadrupole.py, 3D_quadrupole_with_Kelvin.py.

Panels:
  (a) Dipole Hz(x) comparison: analytical vs finite domain vs Kelvin
  (b) Quadrupole Hz(x) comparison: analytical vs finite domain vs Kelvin
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat


def main():
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.2))
    row_idx = 110  # corresponds to MATLAB index 111 (0-based)

    # (a) Dipole
    ax = axes[0]

    data_k = loadmat("3D_dipole_with_Kelvin.mat")
    x = data_k["xx"][row_idx, :]
    Hz_ana = data_k["Hz_analytical"][row_idx, :].copy()
    Hz_ana[np.abs(x) >= 0.5] = -Hz_ana[np.abs(x) >= 0.5]
    ax.plot(x, Hz_ana, "k-", linewidth=1.0, label="Analytical")

    data_f = loadmat("3D_dipole.mat")
    x_f = data_f["xx"][row_idx, :]
    Hz_f = data_f["Hz"][row_idx, :]
    ax.plot(x_f, Hz_f, "r-", label="Finite domain")

    Hz_kelvin = data_k["Hz"][row_idx, :]
    ax.plot(x, Hz_kelvin, "c--", linewidth=1.5, label="Kelvin trans.")

    ax.set_xlim(-1.1, 1.1)
    ax.set_ylim(-1.1, 0.0)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("Hz (A/m)")
    ax.text(0.75, -1.35, "(a)", fontsize=14)
    ax.legend(loc="best", frameon=False, fontsize=11)
    ax.grid(True)

    # (b) Quadrupole
    ax = axes[1]

    data_k = loadmat("3D_quadrupole_with_Kelvin.mat")
    x = data_k["xx"][row_idx, :]
    Hz_ana = data_k["Hz_analytical"][row_idx, :]
    ax.plot(x, Hz_ana, "k-", linewidth=1.0, label="Analytical")

    data_f = loadmat("3D_quadrupole.mat")
    x_f = data_f["xx"][row_idx, :]
    Hz_f = data_f["Hz"][row_idx, :]
    ax.plot(x_f, Hz_f, "r-", label="Finite domain")

    Hz_kelvin = data_k["Hz"][row_idx, :]
    ax.plot(x, Hz_kelvin, "c--", linewidth=1.5, label="Kelvin trans.")

    ax.set_xlim(-1.1, 1.1)
    ax.set_ylim(-0.60, 0.60)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("Hz (A/m)")
    ax.text(0.75, -0.90, "(b)", fontsize=14)
    ax.legend(loc="lower right", frameon=False, fontsize=11)
    ax.grid(True)

    plt.tight_layout()
    plt.savefig("Fig_3.pdf")
    plt.show()


if __name__ == "__main__":
    main()
