"""team28_cln_sweep_full.py -- full 25-point TEAM 28 CLN force-vs-height
sweep + plot, vs the lab full-FEM ground truth.

At every lab height dZ = -7..+17 mm it computes the levitation force from
(a) the direct full split-(K+sN) FEM solve and (b) the 6-stage CLN/Cauer
reduced model, and compares to the lab full-FEM Fz1(dZ).  Saves the data
to JSON and a force-vs-height figure to PNG (both committed next to this
script, per the lab data-persistence policy).

Run:  python team28_cln_sweep_full.py   (~75 s; 25 axisymmetric solves)
"""
import json
import os

import numpy as np

from team28_cln_force import cln_forces, aluminium_z

# Lab full-FEM ground truth Fz1(dZ) [N], NEGATIVE = upward lift,
# from W:\00_CAE\NGSolve\01_菅原\2024_08_TEAM28\50Hz_可動\axisymmetric_mixed.mat
LAB = {
    -7: -6.5794, -6: -5.6713, -5: -4.8770, -4: -4.1833, -3: -3.5780,
    -2: -3.0505, -1: -2.5916, 0: -2.1928, 1: -1.8469, 2: -1.5475,
    3: -1.2887, 4: -1.0655, 5: -0.8736, 6: -0.7090, 7: -0.5683,
    8: -0.4483, 9: -0.3466, 10: -0.2606, 11: -0.1883, 12: -0.1279,
    13: -0.0779, 14: -0.0367, 15: -0.0033, 16: 0.0236, 17: 0.0447,
}
DISK_WEIGHT = 1.055   # N (R=65mm, t=3mm Al, 107.5 g)
NSTAGE = 6
HERE = os.path.dirname(os.path.abspath(__file__))


def equilibrium(dz_mm, force):
    """Linear-interp dZ where lift == -weight (force crosses -DISK_WEIGHT)."""
    tgt = -DISK_WEIGHT
    for i in range(len(dz_mm) - 1):
        if (force[i] - tgt) * (force[i + 1] - tgt) <= 0:
            t = (tgt - force[i]) / (force[i + 1] - force[i])
            return dz_mm[i] + t * (dz_mm[i + 1] - dz_mm[i])
    return None


def run():
    dz_mm = sorted(LAB)
    fz_full, fz_cln, fz_lab = [], [], []
    print(" dZ[mm]  full-FEM   CLN(6)    lab-ref   CLN-vs-lab")
    for dz in dz_mm:
        ff, sf = cln_forces(aluminium_z + dz * 1e-3, max_stage=NSTAGE)
        fc = sf[-1]
        ref = LAB[dz]
        err = abs(fc - ref) / abs(ref) * 100 if ref != 0 else float("nan")
        print(f"  {dz:4d}   {ff:+.4f}  {fc:+.4f}  {ref:+.4f}   {err:6.2f} %")
        fz_full.append(ff); fz_cln.append(fc); fz_lab.append(ref)

    eq_cln = equilibrium(dz_mm, fz_cln)
    eq_lab = equilibrium(dz_mm, fz_lab)
    maxerr = max(abs(c - l) for c, l in zip(fz_cln, fz_lab))
    print(f"\n equilibrium  CLN dZ={eq_cln:.2f} mm   lab dZ={eq_lab:.2f} mm")
    print(f" max |CLN - lab| over the sweep = {maxerr:.4f} N")

    data = {
        "geometry": {"disk_R_mm": 65, "disk_t_mm": 3, "sigma": 3.4e7,
                     "coil1": "960t/+20A r41mm", "coil2": "576t/-20A r87.5mm",
                     "freq_Hz": 50},
        "cln_stages": NSTAGE, "disk_weight_N": DISK_WEIGHT,
        "dZ_mm": dz_mm, "fz_full_N": fz_full, "fz_cln_N": fz_cln,
        "fz_lab_N": fz_lab,
        "equilibrium_dZ_mm": {"cln": eq_cln, "lab": eq_lab},
        "max_abs_cln_minus_lab_N": maxerr,
    }
    with open(os.path.join(HERE, "team28_cln_sweep_results.json"), "w") as f:
        json.dump(data, f, indent=2)
    plot(dz_mm, fz_full, fz_cln, fz_lab, eq_cln)
    return data


def plot(dz_mm, fz_full, fz_cln, fz_lab, eq_cln):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 10, "font.family": "serif",
                         "axes.linewidth": 0.8})
    fig, ax = plt.subplots(figsize=(8 / 2.54, 6 / 2.54))
    dz = np.array(dz_mm)
    ax.plot(dz, np.array(fz_lab), "k-", lw=1.0, label="lab full-FEM")
    ax.plot(dz, np.array(fz_full), "C0--", lw=0.8, label="repo full-FEM")
    ax.plot(dz, np.array(fz_cln), "C3o", ms=3, mfc="none", label="CLN (6 stages)")
    ax.axhline(-1.055, color="0.6", ls=":", lw=0.8)
    ax.text(dz.min(), -1.055, " lift = weight", va="bottom", ha="left",
            fontsize=8, color="0.4")
    if eq_cln is not None:
        ax.axvline(eq_cln, color="0.6", ls=":", lw=0.8)
    ax.set_xlabel(r"disk displacement $\Delta z$ (mm)")
    ax.set_ylabel(r"levitation force $F_z$ (N)")
    ax.legend(frameon=False, fontsize=8)
    ax.tick_params(direction="in", top=True, right=True)
    fig.tight_layout(pad=0.3)
    out = os.path.join(HERE, "team28_cln_force_vs_height.png")
    fig.savefig(out, dpi=300)
    print(f" wrote {out}")


if __name__ == "__main__":
    run()
