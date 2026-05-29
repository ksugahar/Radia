"""Digest 2-panel figure: (a) per-element-vs-uniform P_wp gap CONTOUR
map (narrower) + (b) |H_t| spatial map on the workpiece side wall.

Rendered wide (paper_double_column) and embedded IN-COLUMN at
0.92\\columnwidth of the 1-page digest -- the same footprint as the
old 2-panel sweep figure, so it fits one page.  Panel (a) is made
narrower (width_ratios) to leave room for (b).  NOTE: an in-column
2-panel is scaled down by LaTeX, so on-page fonts are ~font_pt x 0.44;
check_min_font reports the actual visible size.

Inputs (committed; Data Persistence Policy):
  sweep_data/sweep_results.json        -- gap heatmap (32-case sweep)
  sweep_data/I100_f50k_per_panel.json  -- per-DOF |H_t| (100 A, 50 kHz)
  ../../src/radia/panels/samples/ih_bem_sample_p1.vol -- mesh (theta,z)

Output: sweep_heatmap_digest.png  (the digest figure).
"""
import sys
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

from mcp_server_document.graph.tools import (
    apply_lab_style, lab_savefig, check_min_font,
)

HERE = Path(__file__).resolve().parent
SWEEP = HERE / "sweep_data"
VOL = (HERE.parent.parent / "src" / "radia" / "panels" / "samples"
       / "ih_bem_sample_p1.vol")
OUT = HERE / "sweep_heatmap_digest"
MIN_PT = 9.0
EMBED_SCALE = 0.46   # \columnwidth (~8.4 cm) of an 18.2 cm render.


def load_gap(res_path):
    data = json.load(open(res_path))
    runs = data["runs"]
    freqs = sorted({r["frequency_Hz"] for r in runs})
    currents = sorted({r["current_A"] for r in runs})
    P_s = np.full((len(currents), len(freqs)), np.nan)
    P_p = np.full_like(P_s, np.nan)
    for r in runs:
        if r.get("P_wp_W") is None or "error" in r:
            continue
        i = currents.index(r["current_A"])
        j = freqs.index(r["frequency_Hz"])
        (P_p if r["per_panel"] else P_s)[i, j] = r["P_wp_W"]
    gap = (P_p / np.where(P_s > 1e-30, P_s, np.nan) - 1.0) * 100
    return currents, freqs, gap


def load_ht_map(json_path, vol_path):
    from ngsolve import Mesh, BND, TaskManager
    d = json.load(open(json_path))
    H = np.asarray(d["esim_per_panel_H_t"])
    with TaskManager():
        mesh = Mesh(str(vol_path))
        labs = list(mesh.GetBoundaries())
        sibc = {i for i, n in enumerate(labs) if n.lower() == "sibc"}
        bv = sorted({int(v.nr) for el in mesh.Elements(BND)
                     if int(el.index) in sibc for v in el.vertices})
        pts = np.array([list(mesh.ngmesh.Points()[i + 1]) for i in bv])
    x, y, z = pts.T
    r = np.sqrt(x**2 + y**2)
    side = r > 0.95 * r.max()
    return np.degrees(np.arctan2(y[side], x[side])), z[side] * 1e3, H[side]


def main():
    res = SWEEP / "sweep_results.json"
    perp = SWEEP / "I100_f50k_per_panel.json"
    for p in (res, perp, VOL):
        if not p.exists():
            print(f"ERROR: {p} not found.")
            sys.exit(1)
    currents, freqs, gap = load_gap(res)
    th, zs, Ht = load_ht_map(perp, VOL)

    figsize = apply_lab_style(target="paper_double_column", aspect=0.40)
    # Render with inflated fonts so that, after the ~0.44 in-column
    # downscale (0.92\columnwidth of an 18.2 cm render), every glyph
    # still lands >= 9 pt on the page (21 x 0.44 ~= 9.2 pt).  The heatmap
    # (4 labelled columns) gets the wider slot; the |H_t| scatter (3 theta
    # ticks) is happy narrow.
    plt.rcParams.update({
        "font.size": 21, "axes.labelsize": 21,
        "xtick.labelsize": 21, "ytick.labelsize": 21,
    })
    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=figsize, constrained_layout=True,
        gridspec_kw={"width_ratios": [1.45, 1.0]})
    tick_pt = plt.rcParams["ytick.labelsize"]

    # ----- (a) gap CONTOUR map on log-log (I, f) -----
    # 16 sampled (I, f) points; contourf interpolates between them and the
    # actual sampled values are overlaid so the 4x4 sampling stays honest.
    # Evenly-spaced index axes (the sweep points are categorical / unevenly
    # log-spaced -- 100 and 300 A nearly coincide in log), tick-labelled with
    # the actual values; contourf interpolates uniformly between them.
    xf = np.arange(1, len(freqs) + 1, dtype=float)
    yi = np.arange(1, len(currents) + 1, dtype=float)
    Xf, Yi = np.meshgrid(xf, yi)
    norm = TwoSlopeNorm(vmin=min(gap.min(), -1.0), vcenter=0.0,
                        vmax=max(gap.max(), 1.0))
    cf = ax1.contourf(Xf, Yi, gap, levels=np.arange(-50, 31, 5),
                      cmap="RdYlGn_r", norm=norm, extend="both")
    ax1.set_xticks(xf)
    ax1.set_xticklabels([f"{f/1e3:.0f}" for f in freqs])
    ax1.set_yticks(yi)
    ax1.set_yticklabels([f"{I:.0f}" for I in currents])
    ax1.set_xlabel(r"frequency (kHz)")
    ax1.set_ylabel(r"$I_{\mathrm{port}}$ (A)")
    cl = ax1.contour(Xf, Yi, gap, levels=[-40, -20, 0, 20], colors="k",
                     linewidths=0.6)
    ax1.clabel(cl, fmt="%d", fontsize=tick_pt-1)  # >=9 pt after 0.46 downscale
    ax1.scatter(Xf, Yi, c="k", s=6, zorder=5)   # the 16 sampled (I,f) points
    cb1 = fig.colorbar(cf, ax=ax1)
    cb1.set_label(r"gap (\%)")
    cb1.ax.tick_params(labelsize=tick_pt)
    ax1.text(0.5, -0.30, "(a)", transform=ax1.transAxes, ha="center", va="top")

    # ----- (b) |H_t| spatial map on the cylinder side wall -----
    sc = ax2.scatter(th, zs, c=Ht, cmap="magma", s=6, edgecolors="none")
    cb2 = fig.colorbar(sc, ax=ax2)
    cb2.set_label(r"$|H_t|$ (A/m)")
    cb2.ax.tick_params(labelsize=tick_pt)
    ax2.set_xlabel(r"$\theta$ (deg)")
    ax2.set_ylabel(r"$z$ (mm)")
    ax2.set_xticks([-180, 0, 180])
    ax2.text(0.5, -0.30, "(b)", transform=ax2.transAxes, ha="center", va="top")

    fig.canvas.draw()
    bad = check_min_font(fig, min_pt=MIN_PT, embed_scale=EMBED_SCALE)
    if bad:
        vis = min(b["visible_pt"] for b in bad)
        print(f"min-font: {len(bad)} glyph(s) < {MIN_PT} pt at "
              f"embed_scale={EMBED_SCALE}; smallest visible = {vis} pt")
    else:
        print(f"min-font: all text >= {MIN_PT} pt at embed_scale={EMBED_SCALE}")
    lab_savefig(fig, str(OUT))
    print(f"Saved {OUT.with_suffix('.png')}")


if __name__ == "__main__":
    main()
