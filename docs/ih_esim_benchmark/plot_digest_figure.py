"""Digest 2-panel figure: (a) per-element-vs-uniform P_wp gap heatmap
(narrower) + (b) unwrapped side-wall map of local |Z_s|.

Rendered wide (paper_double_column) and embedded IN-COLUMN at
0.92\\columnwidth of the 1-page digest.  Panel (a) is made
narrower (width_ratios) to leave room for (b).  NOTE: an in-column
2-panel is scaled down by LaTeX, so on-page fonts are ~font_pt x 0.44;
check_min_font reports the actual visible size.

Inputs (committed; Data Persistence Policy) -- the figure draws ONLY
from these committed files: no .vol mesh, no NGSolve, no ESIM re-solve.
  sweep_data_dense/sweep_results.json  -- gap contour (108-case sweep)
  sweep_data_dense/I500_f10k_Zs_side_field.json -- self-contained side-wall
        |Z_s| field for the 500 A / 10 kHz per-DOF case.
        Regenerate from the per-DOF JSON + the (gitignored) .vol mesh
        with:  python plot_digest_figure.py --regen-zs-field

Output: sweep_heatmap_digest.png  (the digest figure).
"""
import sys
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from scipy.interpolate import griddata

from radia_mcp.figure import (
    apply_lab_style, lab_savefig, check_min_font,
)

HERE = Path(__file__).resolve().parent
SWEEP_DENSE = HERE / "sweep_data_dense"
VOL = (HERE.parent.parent / "src" / "radia" / "panels" / "samples"
       / "ih_bem_sample_p1.vol")
# Committed, self-contained side-wall field (theta, z, |H_t|, R) for the
# 100 A / 50 kHz per-panel case.  The figure draws from THIS file, so it
# needs neither the (gitignored) .vol mesh nor an ESIM re-solve -- it is
# the persisted "solution" the figure is reproduced from (Data
# Persistence Policy).  Regenerate with: python plot_digest_figure.py --regen-field
SIDE_FIELD = SWEEP_DENSE / "I100_f50k_side_field.json"
# |Z_s| side-wall field for the dense-grid maximum-gap case
# (500 A, 10 kHz, -49.8 %, ESIM converged).
SIDE_FIELD_ZS = SWEEP_DENSE / "I500_f10k_Zs_side_field.json"
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
    perp = SWEEP_DENSE / "I100_f50k_per_panel.json"
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
    (500 A / 10 kHz, dense-grid, ESIM converged) and persist it
    (committed)."""
    perp = SWEEP_DENSE / "I500_f10k_per_panel.json"
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
        "case": "500 A, 10 kHz, per-panel ESIM (max-gap, -49.75 %)",
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
    dense = SWEEP_DENSE / "sweep_results.json"
    if not dense.exists():
        print(f"ERROR: dense sweep file not found: {dense}")
        sys.exit(1)
    n_runs = len(json.load(open(dense))["runs"])
    if n_runs < 108:
        print(f"ERROR: dense sweep is incomplete: {n_runs}/108 runs in {dense}")
        sys.exit(1)
    res = dense
    print(f"Using {res} (mode=contour)")
    currents, freqs, gap = load_gap(res)
    # Panel (b): unwrapped |Z_s| on the side wall at the maximum-gap case
    # (500 A / 10 kHz, -49.75 %).  Reads the committed side-wall field
    # -- no mesh, no re-solve.
    th, zs, field, _ = load_side_field_zs()

    figsize = apply_lab_style(target="paper_double_column", aspect=0.40)
    # Render with inflated fonts so that, after the ~0.44 in-column
    # downscale (0.92\columnwidth of an 18.2 cm render), every glyph
    # still lands >= 9 pt on the page (21 x 0.44 ~= 9.2 pt).  The heatmap
    # (4 labelled columns) gets the wider slot; the side-wall map needs only
    # three theta ticks and a short height axis.
    plt.rcParams.update({
        "font.size": 21, "axes.labelsize": 21,
        "xtick.labelsize": 21, "ytick.labelsize": 21,
    })
    fig = plt.figure(figsize=figsize, constrained_layout=True)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.45, 1.0])
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    tick_pt = plt.rcParams["ytick.labelsize"]

    # ----- (a) gap: contour on log-log (I, f) for the dense sweep.
    from matplotlib.colors import TwoSlopeNorm
    Fg, Ig = np.meshgrid(np.asarray(freqs) / 1e3, np.asarray(currents))
    mask = ~np.isnan(gap)
    vmin = min(np.nanmin(gap), -1.0)
    vmax = max(np.nanmax(gap), 1.0)
    norm = TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax)
    cs = ax1.contourf(Fg, Ig, gap, levels=np.arange(-50, 31, 5),
                      cmap="RdYlGn_r", norm=norm, extend="both")
    ax1.scatter(Fg[mask], Ig[mask], c="k", s=4, zorder=5)
    ax1.set_xscale("log"); ax1.set_yscale("log")
    ax1.set_xticks([10, 100, 500])
    ax1.set_xticklabels(["10", "100", "500"])
    ax1.set_yticks([1, 10, 100, 500])
    ax1.set_yticklabels(["1", "10", "100", "500"])
    cb1 = fig.colorbar(cs, ax=ax1)
    ax1.set_xlabel(r"frequency (kHz)")
    ax1.set_ylabel(r"$I_{\mathrm{port}}$ (A)")
    cb1.set_label(r"gap (\%)")
    cb1.ax.tick_params(labelsize=tick_pt)
    ax1.text(0.5, -0.30, "(a)", transform=ax1.transAxes, ha="center", va="top")

    # ----- (b) |Z_s| on the unwrapped cylinder side wall (max-gap case)
    # Interpolate the scattered (theta, z) side-wall DOFs onto a regular
    # grid (theta tiled +-360 so the wrap is seamless).  The 2-D unwrapped
    # map makes the local surface variation easier to read at digest scale
    # than a small 3-D cylinder.
    th_a = np.concatenate([th - 360, th, th + 360])
    z_a = np.concatenate([zs, zs, zs])
    F_a = np.concatenate([field, field, field])
    TH, ZZ = np.meshgrid(np.linspace(-180, 180, 160),
                         np.linspace(zs.min(), zs.max(), 80))
    Hg = griddata((th_a, z_a), F_a, (TH, ZZ), method="linear")
    Hg = np.where(np.isnan(Hg),
                  griddata((th_a, z_a), F_a, (TH, ZZ), method="nearest"), Hg)
    nrm = Normalize(vmin=field.min(), vmax=field.max())
    cs2 = ax2.contourf(TH, ZZ, Hg, levels=np.linspace(field.min(), field.max(), 14),
                       cmap="jet", norm=nrm)
    ax2.set_xlim(-180, 180)
    ax2.set_ylim(zs.min(), zs.max())
    ax2.set_xticks([-180, 0, 180])
    ax2.set_yticks([zs.min(), 0, zs.max()])
    ax2.set_yticklabels([f"{zs.min():.1f}", "0", f"{zs.max():.1f}"])
    ax2.set_xlabel(r"azimuth (deg)")
    ax2.set_ylabel(r"$z$ (mm)")
    cb2 = fig.colorbar(cs2, ax=ax2, pad=0.04)
    cb2.set_ticks([5, 10, 15])
    cb2.set_label(r"$|Z_s|$ (m$\Omega$)")
    cb2.ax.tick_params(labelsize=tick_pt)
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
