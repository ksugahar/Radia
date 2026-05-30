"""Digest 2-panel figure: (a) per-element-vs-uniform P_wp gap heatmap
(narrower) + (b) |H_t| painted on the 3D workpiece cylinder (axis-equal).

Rendered wide (paper_double_column) and embedded IN-COLUMN at
0.92\\columnwidth of the 1-page digest -- the same footprint as the
old 2-panel sweep figure, so it fits one page.  Panel (a) is made
narrower (width_ratios) to leave room for (b).  NOTE: an in-column
2-panel is scaled down by LaTeX, so on-page fonts are ~font_pt x 0.44;
check_min_font reports the actual visible size.

Inputs (committed; Data Persistence Policy) -- the figure draws ONLY
from these committed files: no .vol mesh, no NGSolve, no ESIM re-solve.
  sweep_data/sweep_results.json        -- gap heatmap (32-case sweep)
  sweep_data/I100_f50k_side_field.json -- self-contained side-wall field
        (theta, z, |H_t|, R) for the 100 A / 50 kHz per-panel case.
        Regenerate from the per-DOF JSON + the (gitignored) .vol mesh
        with:  python plot_digest_figure.py --regen-field

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
SWEEP_DENSE = HERE / "sweep_data_dense"
VOL = (HERE.parent.parent / "src" / "radia" / "panels" / "samples"
       / "ih_bem_sample_p1.vol")
# Committed, self-contained side-wall field (theta, z, |H_t|, R) for the
# 100 A / 50 kHz per-panel case.  The figure draws from THIS file, so it
# needs neither the (gitignored) .vol mesh nor an ESIM re-solve -- it is
# the persisted "solution" the figure is reproduced from (Data
# Persistence Policy).  Regenerate with: python plot_digest_figure.py --regen-field
SIDE_FIELD = SWEEP / "I100_f50k_side_field.json"
# |Z_s| side-wall field for the maximum-gap case (300 A, 10 kHz, -47.5 %).
SIDE_FIELD_ZS = SWEEP / "I300_f10k_Zs_side_field.json"
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


def extract_side_field_from_mesh(json_path, vol_path):
    """Recover side-wall (theta, z, |H_t|, R) from the per-DOF JSON + mesh.

    Needs NGSolve and the (gitignored) .vol mesh -- used ONLY by the
    --regen-field step that writes the committed SIDE_FIELD; the figure
    itself never calls this.
    """
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


def save_side_field(out_path=SIDE_FIELD):
    """One-time: extract the side-wall field and persist it (committed)."""
    perp = SWEEP / "I100_f50k_per_panel.json"
    if not perp.exists() or not VOL.exists():
        print(f"ERROR: need {perp} and {VOL} to regenerate the side field.")
        sys.exit(1)
    th, zs, Ht, R = extract_side_field_from_mesh(perp, VOL)
    json.dump({
        "case": "100 A, 50 kHz, per-panel ESIM",
        "source_mesh": VOL.name,
        "note": ("workpiece side wall (r > 0.95 R_max); theta in deg, "
                 "z in mm, |H_t| in A/m, R in mm.  Self-contained so the "
                 "digest figure draws without the .vol mesh or a re-solve."),
        "R_mm": R,
        "theta_deg": th.tolist(),
        "z_mm": zs.tolist(),
        "H_t": Ht.tolist(),
    }, open(out_path, "w"), indent=1)
    print(f"Saved {out_path}  ({len(th)} side-wall points, R={R:.1f} mm, "
          f"z span={zs.max()-zs.min():.1f} mm)")


def save_side_field_zs(out_path=SIDE_FIELD_ZS):
    """One-time: extract |Z_s| on the side wall for the max-gap case
    (300 A / 10 kHz) and persist it (committed)."""
    perp = SWEEP / "I300_f10k_per_panel.json"
    if not perp.exists() or not VOL.exists():
        print(f"ERROR: need {perp} and {VOL} to regenerate the |Z_s| field.")
        sys.exit(1)
    from ngsolve import Mesh, BND, TaskManager
    d = json.load(open(perp))
    zr = np.asarray(d["esim_per_panel_Z_s_real"])
    zi = np.asarray(d["esim_per_panel_Z_s_imag"])
    Zm = np.sqrt(zr ** 2 + zi ** 2)  # |Z_s| in Ohm
    with TaskManager():
        mesh = Mesh(str(VOL))
        labs = list(mesh.GetBoundaries())
        sibc = {i for i, n in enumerate(labs) if n.lower() == "sibc"}
        bv = sorted({int(v.nr) for el in mesh.Elements(BND)
                     if int(el.index) in sibc for v in el.vertices})
        pts = np.array([list(mesh.ngmesh.Points()[i + 1]) for i in bv])
    x, y, z = pts.T
    r = np.sqrt(x ** 2 + y ** 2)
    side = r > 0.95 * r.max()
    R_mm = float(r[side].mean()) * 1e3
    th = np.degrees(np.arctan2(y[side], x[side]))
    zs_mm = z[side] * 1e3
    Z_side_mOhm = (Zm[side] * 1e3).tolist()  # mOhm
    json.dump({
        "case": "300 A, 10 kHz, per-panel ESIM (max-gap, -47.5 %)",
        "source_mesh": VOL.name,
        "note": ("workpiece side wall (r > 0.95 R_max); theta in deg, "
                 "z in mm, |Z_s| in mOhm, R in mm.  Self-contained so "
                 "the digest figure draws without the .vol mesh or a "
                 "re-solve."),
        "R_mm": R_mm,
        "theta_deg": th.tolist(),
        "z_mm": zs_mm.tolist(),
        "Z_s_mOhm": Z_side_mOhm,
    }, open(out_path, "w"), indent=1)
    print(f"Saved {out_path}  ({len(th)} side-wall points, R={R_mm:.1f} mm, "
          f"|Z_s| range {min(Z_side_mOhm):.2f}--{max(Z_side_mOhm):.2f} mOhm)")


def load_side_field(path=SIDE_FIELD):
    """Read the committed side-wall field -- no NGSolve, no mesh, no re-solve."""
    if not path.exists():
        print(f"ERROR: {path} not found.  Regenerate it once with:\n"
              f"  python {Path(__file__).name} --regen-field")
        sys.exit(1)
    d = json.load(open(path))
    return (np.asarray(d["theta_deg"]), np.asarray(d["z_mm"]),
            np.asarray(d["H_t"]), float(d["R_mm"]))


def load_side_field_zs(path=SIDE_FIELD_ZS):
    """Read the committed |Z_s| side-wall field (max-gap case)."""
    if not path.exists():
        print(f"ERROR: {path} not found.  Regenerate it once with:\n"
              f"  python {Path(__file__).name} --regen-zs-field")
        sys.exit(1)
    d = json.load(open(path))
    return (np.asarray(d["theta_deg"]), np.asarray(d["z_mm"]),
            np.asarray(d["Z_s_mOhm"]), float(d["R_mm"]))


def main():
    if "--regen-field" in sys.argv:
        save_side_field()
        return
    if "--regen-zs-field" in sys.argv:
        save_side_field_zs()
        return
    # Prefer the dense sweep if available (renders as a log-log contour);
    # otherwise fall back to the original 4x4 grid (rendered as annotated
    # heatmap cells).
    dense = SWEEP_DENSE / "sweep_results.json"
    legacy = SWEEP / "sweep_results.json"
    # Prefer the dense sweep only when fully complete (108 runs for the
    # 9 currents x 6 frequencies x 2 methods grid); otherwise fall back
    # to the legacy 4x4 heatmap so the figure stays sensible while the
    # dense sweep is in progress.
    use_dense = False
    if dense.exists():
        try:
            n_runs = len(json.load(open(dense))["runs"])
            if n_runs >= 108:
                use_dense = True
        except Exception:
            pass
    if use_dense:
        res = dense
        gap_mode = "contour"
    elif legacy.exists():
        res = legacy
        gap_mode = "heatmap"
    else:
        print(f"ERROR: neither {dense} nor {legacy} found.")
        sys.exit(1)
    print(f"Using {res} (mode={gap_mode})")
    currents, freqs, gap = load_gap(res)
    # Panel (b): |Z_s| on the side wall at the maximum-gap case
    # (300 A / 10 kHz, -47.5 %).  Reads the committed side-wall field
    # -- no mesh, no re-solve.
    th, zs, field, R = load_side_field_zs()

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

    # ----- (a) gap: contour on log-log (I, f) for the dense sweep,
    # or annotated heatmap cells for the legacy 4x4.
    if gap_mode == "contour":
        from matplotlib.colors import TwoSlopeNorm
        Fg, Ig = np.meshgrid(np.asarray(freqs) / 1e3, np.asarray(currents))
        mask = ~np.isnan(gap)
        vmin = min(np.nanmin(gap), -1.0)
        vmax = max(np.nanmax(gap), 1.0)
        norm = TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax)
        cs = ax1.contourf(Fg, Ig, gap, levels=np.arange(-50, 31, 5),
                          cmap="RdYlGn_r", norm=norm, extend="both")
        cl = ax1.contour(Fg, Ig, gap, levels=[-40, -20, 0, 20], colors="k",
                         linewidths=0.6)
        ax1.clabel(cl, fmt="%d", fontsize=tick_pt - 1)
        ax1.scatter(Fg[mask], Ig[mask], c="k", s=4, zorder=5)
        ax1.set_xscale("log"); ax1.set_yscale("log")
        ax1.set_xticks([10, 50, 100, 500])
        ax1.set_xticklabels(["10", "50", "100", "500"])
        ax1.set_yticks([1, 10, 100, 500])
        ax1.set_yticklabels(["1", "10", "100", "500"])
        cb1 = fig.colorbar(cs, ax=ax1)
    else:
        im = ax1.imshow(gap, aspect="auto", origin="lower", cmap="RdYlGn_r",
                        extent=[0.5, len(freqs)+0.5, 0.5, len(currents)+0.5])
        ax1.set_xticks(range(1, len(freqs)+1))
        ax1.set_xticklabels([f"{f/1e3:.0f}" for f in freqs])
        ax1.set_yticks(range(1, len(currents)+1))
        ax1.set_yticklabels([f"{I:.0f}" for I in currents])
        cb1 = fig.colorbar(im, ax=ax1)
        for i in range(len(currents)):
            for j in range(len(freqs)):
                if not np.isnan(gap[i, j]):
                    ax1.text(j+1, i+1, f"{gap[i, j]:+.0f}", ha="center",
                             va="center", fontsize=tick_pt,
                             color="white" if abs(gap[i, j]) > 25 else "black")
    ax1.set_xlabel(r"frequency (kHz)")
    ax1.set_ylabel(r"$I_{\mathrm{port}}$ (A)")
    cb1.set_label(r"gap (\%)")
    cb1.ax.tick_params(labelsize=tick_pt)
    ax1.text(0.5, -0.30, "(a)", transform=ax1.transAxes, ha="center", va="top")

    # ----- (b) |Z_s| painted on the 3D cylinder side wall (max-gap case)
    # Interpolate the scattered (theta, z) side-wall DOFs onto a regular
    # grid (theta tiled +-360 so the wrap is seamless), map to the
    # cylinder (x, y, z) and paint |Z_s| via plot_surface facecolors.
    th_a = np.concatenate([th - 360, th, th + 360])
    z_a = np.concatenate([zs, zs, zs])
    F_a = np.concatenate([field, field, field])
    TH, ZZ = np.meshgrid(np.linspace(-180, 180, 160),
                         np.linspace(zs.min(), zs.max(), 80))
    Hg = griddata((th_a, z_a), F_a, (TH, ZZ), method="linear")
    Hg = np.where(np.isnan(Hg),
                  griddata((th_a, z_a), F_a, (TH, ZZ), method="nearest"), Hg)
    nrm = Normalize(vmin=field.min(), vmax=field.max())
    # Cylinder axis along z (upright); TRUE axis-equal box so the
    # workpiece keeps its physical proportions (box aspect = data ranges
    # 2R : 2R : axial-span).  The hot coil band reads as a bright ring.
    ax2.plot_surface(R * np.cos(np.radians(TH)), R * np.sin(np.radians(TH)), ZZ,
                     facecolors=cm.magma(nrm(Hg)), rstride=1, cstride=1,
                     linewidth=0, antialiased=False, shade=False)
    ax2.set_box_aspect((2 * R, 2 * R, zs.max() - zs.min()), zoom=1.35)
    ax2.view_init(elev=18, azim=-60)
    ax2.grid(False)
    for axp in (ax2.xaxis, ax2.yaxis, ax2.zaxis):
        axp.pane.set_visible(False)
    # No axis tick labels: the cylinder is shown axis-equal, and its
    # diameter / height are given in the text + caption instead.
    ax2.set_xticks([])
    ax2.set_yticks([])
    ax2.set_zticks([])
    sm = cm.ScalarMappable(norm=nrm, cmap="magma")
    sm.set_array([])
    cb2 = fig.colorbar(sm, ax=ax2, shrink=0.6, pad=0.04)
    cb2.set_label(r"$|Z_s|$ (m$\Omega$)")
    cb2.ax.tick_params(labelsize=tick_pt)
    ax2.text2D(0.5, -0.08, "(b)", transform=ax2.transAxes, ha="center",
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
