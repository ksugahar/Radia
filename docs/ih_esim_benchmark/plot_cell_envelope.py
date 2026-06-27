"""Nonlinear cell-solver envelope (validation tier ii).

Sweeps the 1-D radial cell solver across |H_t| in [1, 1e5] A/m for a
representative IH benchmark and overlays the two analytical limits:

  - Low-H limit  : linear-mu Bessel I_0/I_1 at mu_r = mu_r(H -> 0)
                   from the BH curve's initial slope.
  - High-H limit : saturated thin-skin sqrt(omega * rho * mu_0 / 2)
                   where mu -> mu_0.

A smooth monotonic transition through the BH knee is the qualitative
'consistency' check published in the IGTE 2026 ESIM paper as tier (ii)
of the validation hierarchy.

NOTE: this is a *qualitative* check between two analytical limits.
It does NOT validate the cell solver against an independent nonlinear
reference (the Stoll 1974 / Lavers-Biringer 1985 envelopes apply
only to 1-D geometries with uniform H_t; the IH workpiece has a
three-dimensional, spatially non-uniform tangential field).

Usage:
    python docs/ih_esim_benchmark/plot_cell_envelope.py

Outputs:
    cell_envelope.png            (figure)
    cell_envelope_data.json      (numerical sweep data)
"""
import sys
import json
import os
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.special import iv as bessel_iv

from radia_mcp.figure import apply_lab_style, lab_savefig

# Add radia src to path so this works without `pip install`
HERE = Path(__file__).parent.resolve()
SRC = HERE.parent.parent / "src" / "radia"
sys.path.insert(0, str(SRC))

from esim_cell_problem import ESIMFiniteSlabSolver
from em_material import load_bh_file

# Default benchmark configuration matches the IH headline in the
# IGTE 2026 paper (steel cylinder at 50 kHz, sigma typical of
# IH-hardening temperatures).
SIGMA = 2e6                  # S/m (hot steel)
HALF_THICK = 5e-3            # m (cylinder radius)
FREQ = 50e3                  # Hz
MU0 = 4 * np.pi * 1e-7

BH_PATH = HERE.parent.parent / "src" / "radia" / "panels" / "samples" / "em_sample_bh.txt"
OUT_PNG = HERE / "cell_envelope.png"
OUT_JSON = HERE / "cell_envelope_data.json"


def low_H_bessel_Zs(omega, sigma, R, mu_r_linear):
    """Linear-mu Bessel I_0/I_1 closed form for cylindrical conductor."""
    rho = 1.0 / sigma
    delta = np.sqrt(2 * rho / (omega * MU0 * mu_r_linear))
    gamma = (1 + 1j) / delta
    return rho * gamma * bessel_iv(1, gamma * R) / bessel_iv(0, gamma * R)


def main():
    bh = load_bh_file(str(BH_PATH))
    if isinstance(bh, tuple):
        H_arr, B_arr = bh
        bh_curve = list(zip(H_arr.tolist(), B_arr.tolist()))
    else:
        bh_curve = list(bh)

    # mu_r at a low H probe (10 A/m, well below the BH knee)
    H_probe = 10.0
    B_probe = np.interp(H_probe,
                         [r[0] for r in bh_curve],
                         [r[1] for r in bh_curve])
    mu_r_linear = B_probe / (MU0 * H_probe)
    print(f"Linear-mu probe at H={H_probe} A/m: mu_r = {mu_r_linear:.1f}")

    solver = ESIMFiniteSlabSolver(
        geometry="cylinder", half_thickness=HALF_THICK, sigma=SIGMA,
        bh_curve=bh_curve, frequency=FREQ, n_nodes=2000,
    )

    omega = 2 * np.pi * FREQ
    H_sweep = np.logspace(0, 5, 60)
    Z_re, Z_im, Z_abs = [], [], []
    for H_t in H_sweep:
        r = solver.solve(H_t)
        Z = r["Z"]
        Z_re.append(Z.real); Z_im.append(Z.imag); Z_abs.append(abs(Z))
    Z_re = np.array(Z_re); Z_im = np.array(Z_im); Z_abs = np.array(Z_abs)

    Z_low = low_H_bessel_Zs(omega, SIGMA, HALF_THICK, mu_r_linear)
    Z_high = (1 + 1j) * np.sqrt(omega * (1.0 / SIGMA) * MU0 / 2.0)
    print(f"Low-H Bessel limit:  |Z_s| = {abs(Z_low)*1e3:.2f} mOhm")
    print(f"Saturated thin-skin: |Z_s| = {abs(Z_high)*1e3:.4f} mOhm")
    print(f"Cell |Z_s| at H=1:    {Z_abs[0]*1e3:.2f} mOhm")
    print(f"Cell |Z_s| at H=1e5:  {Z_abs[-1]*1e3:.2f} mOhm")

    # IH headline operating range markers (paper Fig. 1 orange band)
    H_t_min, H_t_max = 32.0, 2991.0

    figsize = apply_lab_style(target="paper_single_column", aspect=0.65)
    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    ax.loglog(H_sweep, Z_abs * 1e3, "-", color="C3", linewidth=1.8,
              label=r"cell solver $|Z_s|$")
    ax.axhline(abs(Z_low) * 1e3, color="C0", linestyle="--", linewidth=1.2,
               label=fr"low-$H$ Bessel ($\mu_r={mu_r_linear:.0f}$)")
    ax.axhline(abs(Z_high) * 1e3, color="k", linestyle=":", linewidth=1.2,
               label=r"saturated thin-skin asymptote")
    ax.axvspan(H_t_min, H_t_max, alpha=0.15, color="C1",
               label=r"IH headline $|H_t|$ range")
    ax.set_xlabel(r"$|H_t|$ (A/m)")
    ax.set_ylabel(r"$|Z_s|$ (m$\Omega$)")
    ax.legend(loc="lower left", frameon=False)
    ax.grid(which="both", linestyle=":", alpha=0.5)
    lab_savefig(fig, str(OUT_PNG.with_suffix("")))
    print(f"Saved {OUT_PNG}")

    data = {
        "config": {
            "frequency_hz": FREQ, "sigma_S_per_m": SIGMA,
            "half_thickness_m": HALF_THICK,
            "bh_file": str(BH_PATH), "n_nodes": 2000,
            "mu_r_linear_at_H_10": float(mu_r_linear),
        },
        "asymptotes": {
            "low_H_bessel_Zs_abs": float(abs(Z_low)),
            "saturated_thin_skin_Zs_abs": float(abs(Z_high)),
        },
        "sweep": {
            "H_t_A_per_m": H_sweep.tolist(),
            "Z_s_real": Z_re.tolist(),
            "Z_s_imag": Z_im.tolist(),
            "Z_s_abs": Z_abs.tolist(),
        },
    }
    with open(OUT_JSON, "w") as f:
        json.dump(data, f, indent=2)
    print(f"Saved {OUT_JSON}")


if __name__ == "__main__":
    main()
