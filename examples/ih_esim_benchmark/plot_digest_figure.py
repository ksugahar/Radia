"""Digest 2-panel figure: (a) per-element-vs-uniform P_wp gap heatmap
(narrower) + (b) |H_t| filled-contour map on the workpiece side wall.

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
from matplotlib import cm
from matplotlib.colors import Normalize
from scipy.interpolate import griddata

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
    R_mm = float(r[side].mean()) * 1e3
    return (np.degrees(np.arctan2(y[side], x[side])), z[side] * 1e3,
            H[side], R_mm)


def main():
    res = SWEEP / "sweep_results.json"
    perp = SWEEP / "I100_f50k_per_panel.json"
    for p in (res, perp, VOL):
        if not p.exists():
            print(f"ERROR: {p} not found.")
            sys.exit(1)
    currents, freqs, gap = load_gap(res)
    th, zs, Ht, R = load_ht_map(perp, VOL)

    figsize = apply_lab_style(target="paper_double_column", aspect=0.40)
    # Render with inflated fonts so that, after the ~0.44 in-column
    # downscale (0.92\columnwidth of an 18.2 cm render), every glyph
    # still lands >= 9 pt on the page (21 x 0.44 ~= 9.2 pt).  The heatmap
    # (4 labelled columns) gets the wider slot; the |H_t| contour (3 theta
    # ticks) is happy narrow.
    plt.rcParams.update({
        "font.size": 21, "axes.labelsize": 21,
        "xtick.labelsize": 21, "ytick.labelsize": 21,
    })
    fig = plt.figure(figsize=figsize, constrained_layout=True)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.45, 1.0])
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1], projection="3d")
    tick_pt = plt.rcParams["ytick.labelsize"]

    # ----- (a) gap heatmap (narrower), annotated cells -----
    im = ax1.imshow(gap, aspect="auto", origin="lower", cmap="RdYlGn_r",
                    extent=[0.5, len(freqs)+0.5, 0.5, len(currents)+0.5])
    ax1.set_xticks(range(1, len(freqs)+1))
    ax1.set_xticklabels([f"{f/1e3:.0f}" for f in freqs])
    ax1.set_yticks(range(1, len(currents)+1))
    ax1.set_yticklabels([f"{I:.0f}" for I in currents])
    ax1.set_xlabel(r"frequency (kHz)")
    ax1.set_ylabel(r"$I_{\mathrm{port}}$ (A)")
    cb1 = fig.colorbar(im, ax=ax1)
    cb1.set_label(r"gap (\%)")
    cb1.ax.tick_params(labelsize=tick_pt)
    for i in range(len(currents)):
        for j in range(len(freqs)):
            if not np.isnan(gap[i, j]):
                ax1.text(j+1, i+1, f"{gap[i, j]:+.0f}", ha="center",
                         va="center", fontsize=tick_pt,
                         color="white" if abs(gap[i, j]) > 25 else "black")
    ax1.text(0.5, -0.30, "(a)", transform=ax1.transAxes, ha="center", va="top")

    # ----- (b) |H_t| painted on the 3D cylinder side wall -----
    # Interpolate the scattered (theta, z) side-wall DOFs onto a regular
    # grid (theta tiled +-360 so the wrap is seamless), map to the
    # cylinder (x, y, z) and paint |H_t| via plot_surface facecolors.
    th_a = np.concatenate([th - 360, th, th + 360])
    z_a = np.concatenate([zs, zs, zs])
    H_a = np.concatenate([Ht, Ht, Ht])
    TH, ZZ = np.meshgrid(np.linspace(-180, 180, 160),
                         np.linspace(zs.min(), zs.max(), 80))
    Hg = griddata((th_a, z_a), H_a, (TH, ZZ), method="linear")
    Hg = np.where(np.isnan(Hg),
                  griddata((th_a, z_a), H_a, (TH, ZZ), method="nearest"), Hg)
    nrm = Normalize(vmin=Ht.min(), vmax=Ht.max())
    # Lay the cylinder axis along x (horizontal) so it fills the wide,
    # short in-column panel; the hot coil band reads as a bright ring.
    ax2.plot_surface(ZZ, R * np.cos(np.radians(TH)), R * np.sin(np.radians(TH)),
                     facecolors=cm.magma(nrm(Hg)), rstride=1, cstride=1,
                     linewidth=0, antialiased=False, shade=False)
    # Physically-exact proportions (axial span : diameter : diameter);
    # zoom fills the small panel without distorting the cylinder.
    ax2.set_box_aspect((zs.max() - zs.min(), 2 * R, 2 * R), zoom=1.5)
    ax2.view_init(elev=20, azim=-74)
    ax2.grid(False)
    for axp in (ax2.xaxis, ax2.yaxis, ax2.zaxis):
        axp.pane.set_visible(False)
    ax2.set_yticks([])
    ax2.set_zticks([])
    ax2.set_xticks([-10, 0, 10])
    ax2.tick_params(axis="x", pad=-4)
    ax2.set_xlabel(r"$z$ (mm)", labelpad=-9)
    sm = cm.ScalarMappable(norm=nrm, cmap="magma")
    sm.set_array([])
    cb2 = fig.colorbar(sm, ax=ax2, shrink=0.55, pad=0.02)
    cb2.set_label(r"$|H_t|$ (A/m)")
    cb2.ax.tick_params(labelsize=tick_pt)
    ax2.text2D(0.5, -0.12, "(b)", transform=ax2.transAxes, ha="center",
               va="top")

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
