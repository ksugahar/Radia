"""cube_alpha_sweep_figure.py -- paper figure: Cu cube polarizability alpha(s).

Runs the mixed Galerkin alpha(s) pipeline on a 5 mm Cu cube and renders the
frequency response Re[alpha]/V (reactive flux exclusion) and |Im[alpha]|/V
(eddy-current loss) from 1 Hz to 1 GHz.  Re[alpha]/V rises from 0 (DC, full
penetration) to 1 (high-f perfect-conductor exclusion); the loss |Im| shows
the drag peak in the wall band.

Two forms of the same reduced model are drawn:
  projected  `schur.BoxMixedGalerkin` -- bulk eigenmodes and the surface
             envelope in one Galerkin space, the crossover decided by the
             Schur complement; exact at DC.  The solid curves.
  additive   `alpha.Y_mixed` -- Y_bulk_Foster + K_SIBC/sqrt(s) + c_1/s, the
             composition behind the LTI export; its c_1/s diverges at DC, so
             it is drawn only where a/delta >= 1.  The thin grey curves.

Outputs:
  validation_test/maglev/demos/cube_alpha_sweep_results.json
                                  the swept data of both forms + parameters
  cube_alpha_sweep.pdf / .png     IEEE single-column figure

Figure conventions (lab-wide, self-contained -- no LAB-private helper):
  no in-figure title (goes in the LaTeX caption), Times New Roman,
  ~10 pt on an 8.8 cm (3.5 in) single column, box on, ticks inward.
"""
from __future__ import annotations

import json
import math
import os
import platform
from datetime import datetime, timezone

import numpy as np

from _validation_output import validation_output

SIGMA_CU = 5.8e7
MU0 = 4.0 * math.pi * 1e-7
L = 5e-3                      # cube edge (m)
F_MIN, F_MAX, N_PTS = 1e0, 1e9, 109
N_EIGEN = 80                  # additive form: Foster modes
N_EIGEN_PROJECTED = 40        # projected form: modes sampled on the layer grid
HERE = os.path.dirname(os.path.abspath(__file__))


def run_sweep():
    from netgen.occ import Box, Pnt, OCCGeometry
    from ngsolve import Mesh, TaskManager
    from radia.maglev.mixed_galerkin import alpha as A, cad_edges as CE, schur as SC

    box = Box(Pnt(0, 0, 0), Pnt(L, L, L))
    for f in box.faces:
        f.name = "outer"
    with TaskManager():
        mesh = Mesh(OCCGeometry(box).GenerateMesh(maxh=L / 6))
        lam, tau, g_n, V = A.bulk_foster_via_eigen(mesh, SIGMA_CU, MU0, n_eigen=N_EIGEN)
        mg = SC.BoxMixedGalerkin(mesh, SIGMA_CU, MU0, n_eigen=N_EIGEN_PROJECTED)

    K_SIBC = A.K_SIBC_total(6 * L * L, SIGMA_CU, MU0)
    c1, L_total, n_edges = CE.cad_topology_c1(box, MU0)

    f_arr = np.logspace(math.log10(F_MIN), math.log10(F_MAX), N_PTS)
    re_p, im_p, re_a, im_a, a_over_delta = [], [], [], [], []
    for f in f_arr:
        s = 1j * 2 * math.pi * f
        a_proj = mg.alpha(s)
        re_p.append((a_proj / V).real)
        im_p.append((a_proj / V).imag)
        a_add = A.alpha_from_Y(A.Y_mixed(s, lam, tau, g_n, K_SIBC, c1), V, SIGMA_CU)
        re_a.append((a_add / V).real)
        im_a.append((a_add / V).imag)
        delta = math.sqrt(2.0 / (2 * math.pi * f * MU0 * SIGMA_CU))
        a_over_delta.append((L / 2) / delta)

    return {
        "f_hz": f_arr.tolist(),
        "re_alpha_over_V": re_p,
        "im_alpha_over_V": im_p,
        "additive": {
            "re_alpha_over_V": re_a,
            "im_alpha_over_V": im_a,
            "valid_where": "a_over_delta >= 1 (c_1/s diverges at DC)",
        },
        "a_over_delta": a_over_delta,
        "params": {
            "L_m": L, "sigma": SIGMA_CU, "mu": MU0,
            "V_m3": V, "K_SIBC": K_SIBC, "c1": c1,
            "n_edges": n_edges, "L_total_edges_m": L_total,
            "n_eigen": N_EIGEN, "n_eigen_projected": N_EIGEN_PROJECTED,
            "maxh": L / 6,
            "grid_nodes_per_axis": int(len(mg.nodes[0])),
        },
        "meta": {
            "hostname": platform.node(),
            "script": "cube_alpha_sweep_figure.py",
            "description": ("mixed Galerkin alpha(s) of a 5 mm Cu cube: projected "
                            "(BoxMixedGalerkin, main curves) and additive (Y_mixed)"),
        },
    }


def make_figure(data, stamp_iso):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman"],
        "mathtext.fontset": "stix",
        "font.size": 10,
        "axes.linewidth": 0.8,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.major.size": 4,
        "ytick.major.size": 4,
        "legend.frameon": False,
    })

    f = np.array(data["f_hz"])
    re = np.array(data["re_alpha_over_V"])
    im_abs = np.abs(np.array(data["im_alpha_over_V"]))
    valid = np.array(data["a_over_delta"]) >= 1.0
    re_add = np.where(valid, np.array(data["additive"]["re_alpha_over_V"]), np.nan)
    im_add = np.where(valid, np.abs(np.array(data["additive"]["im_alpha_over_V"])), np.nan)

    # IEEE single column: 3.5 in wide; height ~ 2.4 in.
    fig, ax = plt.subplots(figsize=(3.5, 2.4))
    ax.semilogx(f, re_add, "-", color="0.65", lw=0.8)
    ax.semilogx(f, im_add, "--", color="0.65", lw=0.8, label="additive tail")
    ax.semilogx(f, re, "-", color="#1f3a93", lw=1.4,
                label=r"$\mathrm{Re}\,\alpha/V$ (exclusion)")
    ax.semilogx(f, im_abs, "--", color="#b8001f", lw=1.4,
                label=r"$|\mathrm{Im}\,\alpha/V|$ (loss)")
    ax.axhline(1.0, color="0.5", lw=0.7, ls=":")
    ax.text(1.5e0, 1.02, "PEC limit", fontsize=8, color="0.4")

    ax.set_xlabel(r"frequency $f$ (Hz)")
    ax.set_ylabel(r"$\alpha(s)/V$")
    ax.set_xlim(f[0], f[-1])
    ax.set_ylim(-0.05, 1.15)
    ax.legend(loc="center left", fontsize=8)
    fig.tight_layout(pad=0.3)

    pdf = os.path.join(HERE, "cube_alpha_sweep.pdf")
    png = os.path.join(HERE, "cube_alpha_sweep.png")
    fig.savefig(pdf)
    fig.savefig(png, dpi=300)
    plt.close(fig)
    return pdf, png


def main():
    data = run_sweep()
    # stamp after the (resume-unsafe) clock read, at the very end
    data["meta"]["timestamp"] = datetime.now(timezone.utc).isoformat()
    out_json = validation_output("cube_alpha_sweep_results.json", HERE)
    with open(out_json, "w", encoding="utf-8") as fp:
        json.dump(data, fp, indent=2, allow_nan=False)
    pdf, png = make_figure(data, data["meta"]["timestamp"])

    p = data["params"]
    f = np.array(data["f_hz"])
    i4 = int(np.argmin(np.abs(f - 1e4)))
    print(f"V = {p['V_m3']*1e9:.3f} mm^3, c1 = {p['c1']:.4e} ({p['n_edges']} edges), "
          f"K_SIBC = {p['K_SIBC']:.4e}")
    print(f"projected Re[alpha]/V: {data['re_alpha_over_V'][0]:+.4f} (f={F_MIN:.0e}) "
          f"-> {data['re_alpha_over_V'][-1]:+.4f} (f={F_MAX:.0e}, PEC)")
    print(f"at 10 kHz: projected {data['re_alpha_over_V'][i4]:+.4f}{data['im_alpha_over_V'][i4]:+.4f}j, "
          f"additive {data['additive']['re_alpha_over_V'][i4]:+.4f}"
          f"{data['additive']['im_alpha_over_V'][i4]:+.4f}j")
    print(f"saved: {out_json}")
    print(f"saved: {pdf}")
    print(f"saved: {png}")


if __name__ == "__main__":
    main()
