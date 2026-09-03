"""team28_cln_sweep_full.py -- full 25-point TEAM 28 CLN force-vs-height
sweep + plot, vs the lab full-FEM ground truth.

At every lab height dZ = -7..+17 mm it computes the levitation force from
(a) the direct full split-(K+sN) FEM solve and (b) the 6-stage CLN/Cauer
reduced model, and compares to the lab full-FEM Fz1(dZ). Saves the durable
JSON under validation_test/maglev/demos/team28 and keeps the public
force-vs-height figure next to this script.

Run:  python team28_cln_sweep_full.py   (~75 s; 25 axisymmetric solves)
"""
import json
import os
import sys

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
sys.path.insert(0, os.path.join(HERE, ".."))
from _validation_output import validation_output  # noqa: E402

# --- Force convention (verified C:/temp/team28_force_convention_check.py) ---
# The reported F_z (fz_full / fz_cln / fz_lab) is the VERBATIM TEAM 28 surface
# integral Re[B_r * J_t] * 2*pi*r, which is EXACTLY 2x the physical
# time-averaged Lorentz force  <f_z> = -(1/2) Re[J_t conj(B_r)]  (measured
# ratio 1.9998; the Im*Im cross term is ~6e-5).  This 2x integral is kept for
# the force column because it matches the lab .mat Fz1(dZ) and the golden, and
# because the CLN-vs-full convergence is convention-independent.  BUT the disk
# floats where the PHYSICAL lift equals its weight, so the levitation
# equilibrium MUST use F_z / 2, not F_z.  (Earlier code balanced the 2x
# integral against the 1x weight -> a spurious dZ=+4.1mm / 14.9mm height; the
# physical equilibrium is ~dZ=+0.2mm / 11.0mm, matching the published 11.5mm.)
PHYS = 0.5                    # physical force = PHYS * (verbatim TEAM 28 integral)
DISK_BOTTOM_DZ0_MM = 10.8     # absolute disk-bottom height (above coil top) at dZ=0
# Published TEAM 28 reference (Karl, Fetzer, Kurz, Lehner, Rucker, Univ.
# Stuttgart): rest z=3.8mm, measured stationary levitation height z=11.5mm
# (laser triangulation, 4-measurement average), i_hat=20A peak, f=50Hz.
PUBLISHED_REST_MM = 3.8
PUBLISHED_LEVITATION_MM = 11.5


def equilibrium(dz_mm, force):
    """Linear-interp dZ where the PHYSICAL lift == weight.

    The disk floats where <f_z> = weight; the reported ``force`` is the 2x
    verbatim integral, so balance PHYS*force against -DISK_WEIGHT.
    """
    tgt = -DISK_WEIGHT
    for i in range(len(dz_mm) - 1):
        a, b = PHYS * force[i], PHYS * force[i + 1]
        if (a - tgt) * (b - tgt) <= 0:
            t = (tgt - a) / (b - a)
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
    z_cln = DISK_BOTTOM_DZ0_MM + eq_cln if eq_cln is not None else None
    pub_err = (abs(z_cln - PUBLISHED_LEVITATION_MM) / PUBLISHED_LEVITATION_MM * 100
               if z_cln is not None else None)
    print(f"\n PHYSICAL levitation equilibrium (lift==weight, F_z/2 used):")
    print(f"   CLN  dZ={eq_cln:.2f} mm  -> absolute disk-bottom z={z_cln:.2f} mm")
    print(f"   published measured steady-state z={PUBLISHED_LEVITATION_MM} mm"
          f"  -> agreement {pub_err:.1f} %")
    print(f" max |CLN - lab(verbatim integral)| over the sweep = {maxerr:.4f} N")

    data = {
        "geometry": {"disk_R_mm": 65, "disk_t_mm": 3, "sigma": 3.4e7,
                     "coil1": "960t/+20A r41mm", "coil2": "576t/-20A r87.5mm",
                     "freq_Hz": 50},
        "cln_stages": NSTAGE, "disk_weight_N": DISK_WEIGHT,
        "force_note": ("fz_*_N are the verbatim TEAM 28 integral Re[B_r J_t]"
                       " = 2x the physical time-averaged Lorentz force;"
                       " equilibrium uses F_z/2 (PHYS) == weight"),
        "dZ_mm": dz_mm, "fz_full_N": fz_full, "fz_cln_N": fz_cln,
        "fz_lab_N": fz_lab,
        "equilibrium_dZ_mm": {"cln": eq_cln, "lab": eq_lab},
        "equilibrium_abs_height_mm": {"cln": z_cln,
                                      "disk_bottom_at_dz0": DISK_BOTTOM_DZ0_MM},
        "published_ref": {"authors": "Karl-Fetzer-Kurz-Lehner-Rucker",
                          "rest_mm": PUBLISHED_REST_MM,
                          "levitation_height_mm": PUBLISHED_LEVITATION_MM,
                          "agreement_percent": pub_err},
        "max_abs_cln_minus_lab_N": maxerr,
    }
    output = validation_output("team28_cln_sweep_results.json", HERE)
    with open(output, "w") as f:
        json.dump(data, f, indent=2, allow_nan=False)
    plot(dz_mm, fz_full, fz_cln, fz_lab, eq_cln, z_cln)
    return data


def plot(dz_mm, fz_full, fz_cln, fz_lab, eq_cln, z_cln):
    """Physical time-averaged levitation force (= PHYS * verbatim integral)
    vs height, so the weight line crosses the curve at the levitation height."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 10, "font.family": "serif",
                         "axes.linewidth": 0.8})
    fig, ax = plt.subplots(figsize=(8 / 2.54, 6 / 2.54))
    dz = np.array(dz_mm)
    ax.plot(dz, PHYS * np.array(fz_lab), "k-", lw=1.0, label="lab full-FEM")
    ax.plot(dz, PHYS * np.array(fz_full), "C0--", lw=0.8, label="repo full-FEM")
    ax.plot(dz, PHYS * np.array(fz_cln), "C3o", ms=3, mfc="none",
            label="CLN (6 stages)")
    ax.axhline(-DISK_WEIGHT, color="0.6", ls=":", lw=0.8)
    ax.text(dz.min(), -DISK_WEIGHT, " lift = weight", va="bottom", ha="left",
            fontsize=8, color="0.4")
    if eq_cln is not None:
        ax.axvline(eq_cln, color="C3", ls=":", lw=0.8)
        ax.text(eq_cln, PHYS * np.array(fz_full).min(),
                f"  z={z_cln:.1f} mm", va="bottom", ha="left",
                fontsize=8, color="C3")
    # published measured steady-state levitation height (z=11.5mm -> dZ)
    dz_pub = PUBLISHED_LEVITATION_MM - DISK_BOTTOM_DZ0_MM
    ax.axvline(dz_pub, color="0.3", ls="--", lw=0.8)
    ax.text(dz_pub, -DISK_WEIGHT * 1.9, "published\n11.5 mm", va="bottom",
            ha="center", fontsize=7, color="0.3")
    ax.set_xlabel(r"disk displacement $\Delta z$ (mm)")
    ax.set_ylabel(r"levitation force $\langle F_z\rangle$ (N)")
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    ax.tick_params(direction="in", top=True, right=True)
    fig.tight_layout(pad=0.3)
    out = os.path.join(HERE, "team28_cln_force_vs_height.png")
    fig.savefig(out, dpi=300)
    print(f" wrote {out}")


if __name__ == "__main__":
    run()
