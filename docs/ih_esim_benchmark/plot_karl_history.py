"""Karl iteration convergence diagnostic.

Reads any ESIM JSON (`calc_inductance.py` / `calc_fem_kelvin.py` /
`calc_fem_coilmesh.py` output) and plots the Karl loop trajectory:

    (a) log dZ vs iteration  - convergence rate (geometric => Lipschitz < 1)
    (b) |Z_s| (and min/max if per-panel) vs iteration  - has the impedance plateaued?
    (c) |H_t|_rms (and per-DOF mean/max if per-panel) vs iteration  - has the driver settled?

Use this to distinguish two failure modes:

  * Karl genuinely diverging (dZ non-monotone, |Z_s| oscillating)
    => lower `--esim-relax` or check BH curve monotonicity.

  * Karl hit `max_iter` on the strict per-DOF `dZ` criterion BUT |Z_s|
    and |H_t|_rms have visibly plateaued
    => not a divergence; the integrated P_wp is trustworthy
       (see docs/esim/IMPLEMENTATION.md section 3.4).  Either raise
       `--esim-tol` after confirming the plateau, or accept the
       max_iter cap as the per-DOF noise floor.

Usage:
    python docs/ih_esim_benchmark/plot_karl_history.py <RESULT.json>
                                                            [<OUTPUT.png>]

If <OUTPUT.png> is omitted the figure is saved next to the JSON as
`<json-stem>_karl_history.png`.
"""
import sys
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from radia_mcp.figure import apply_lab_style, lab_savefig


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(2)
    json_path = Path(sys.argv[1])
    if not json_path.exists():
        print(f"ERROR: {json_path} not found.")
        sys.exit(1)
    out_png = Path(sys.argv[2]) if len(sys.argv) > 2 else \
              json_path.with_name(json_path.stem + "_karl_history.png")

    with open(json_path) as f:
        d = json.load(f)
    hist = d.get("esim_history", [])
    if not hist:
        print(f"ERROR: {json_path.name} has no esim_history (was --impedance-model esim used?)")
        sys.exit(1)

    per_panel = bool(d.get("esim_per_panel", False))
    converged = bool(d.get("esim_converged", False))
    n_iter = int(d.get("esim_iterations", len(hist) - 1))
    P_wp = d.get("P_wp_W", None)

    iters = np.array([h["iteration"] for h in hist])
    # Scalar Karl emits Z_s_abs / H_t_rms / dZ.  Per-panel Karl emits
    # Z_s_abs_mean / H_t_per_dof_mean / dZ_max plus min/max.  Accept either.
    Z_abs = np.array([h.get("Z_s_abs", h.get("Z_s_abs_mean", float("nan")))
                       for h in hist])
    H_rms = np.array([h.get("H_t_rms", h.get("H_t_per_dof_mean", float("nan")))
                       for h in hist])
    dZ_seq = np.array([h.get("dZ", h.get("dZ_max", float("nan")))
                       for h in hist])
    # per-panel auxiliary keys (radia >= 4.55.x emit these)
    Z_min = np.array([h.get("Z_s_abs_min", float("nan")) for h in hist])
    Z_max = np.array([h.get("Z_s_abs_max", float("nan")) for h in hist])
    H_mean = np.array([h.get("H_t_per_dof_mean", float("nan")) for h in hist])
    H_max = np.array([h.get("H_t_per_dof_max", float("nan")) for h in hist])
    dZ_max = np.array([h.get("dZ_max", float("nan")) for h in hist])

    print(f"JSON: {json_path.name}")
    print(f"  per_panel = {per_panel}, converged = {converged}, "
          f"n_iter = {n_iter}, P_wp = {P_wp}")
    print(f"  dZ trajectory : {dZ_seq[1:]}")
    print(f"  |Z_s| start -> end : {Z_abs[0]:.4e} -> {Z_abs[-1]:.4e} Ohm")
    print(f"  H_t_rms start -> end : {H_rms[0]:.2f} -> {H_rms[-1]:.2f} A/m")

    # Use a 3-row stacked layout with a shared x-axis so each quantity
    # gets the full figure width.  A 3-column layout at 16 pt fonts
    # leaves no horizontal room for legends; a stacked layout solves
    # the legend-overlap problem cleanly.
    figsize = apply_lab_style(target="paper_single_column", aspect=1.25)
    fig, axes = plt.subplots(3, 1, figsize=figsize, sharex=True,
                              constrained_layout=True)

    # Row (a): log dZ vs iter
    ax = axes[0]
    mask = iters > 0  # skip iter 0 (relaxation seed dZ == 1.0)
    dZ_label = (r"$dZ_{\max}$ (per-DOF)" if per_panel else
                r"$dZ$ (scalar Karl)")
    ax.semilogy(iters[mask], dZ_seq[mask], "o-", color="C3",
                linewidth=1.6, label=dZ_label)
    tol = d.get("esim_tol", 1e-3)
    ax.axhline(tol, color="k", linestyle="--", linewidth=0.8,
               label=fr"tol = {tol:g}")
    ax.set_ylabel(r"$dZ$")
    ax.legend(loc="lower left", frameon=False)
    ax.grid(which="both", linestyle=":", alpha=0.5)
    ax.text(-0.18, 0.5, "(a)", transform=ax.transAxes,
            ha="right", va="center")

    # Row (b): |Z_s| trajectory
    ax = axes[1]
    ax.plot(iters, Z_abs * 1e3, "o-", color="C0", linewidth=1.6,
            label=r"mean")
    if per_panel and not np.all(np.isnan(Z_min)):
        ax.fill_between(iters, Z_min * 1e3, Z_max * 1e3,
                         color="C0", alpha=0.18,
                         label=r"per-DOF $[\min,\max]$")
    ax.set_ylabel(r"$|Z_s|$ (m$\Omega$)")
    ax.legend(loc="upper right", frameon=False, ncol=2)
    ax.grid(linestyle=":", alpha=0.5)
    ax.text(-0.18, 0.5, "(b)", transform=ax.transAxes,
            ha="right", va="center")

    # Row (c): H_t trajectory
    ax = axes[2]
    ax.plot(iters, H_rms, "o-", color="C2", linewidth=1.6,
            label=r"mean")
    if per_panel and not np.all(np.isnan(H_max)):
        ax.plot(iters, H_max, "^--", color="C4", alpha=0.7,
                label=r"per-DOF max")
    ax.set_xlabel(r"Karl iteration")
    ax.set_ylabel(r"$|H_t|$ (A/m)")
    ax.legend(loc="lower right", frameon=False, ncol=2)
    ax.grid(linestyle=":", alpha=0.5)
    ax.text(-0.18, 0.5, "(c)", transform=ax.transAxes,
            ha="right", va="center")

    lab_savefig(fig, str(out_png.with_suffix("")))
    print(f"Saved {out_png}")


if __name__ == "__main__":
    main()
