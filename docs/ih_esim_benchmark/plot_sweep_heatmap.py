"""Per-element-vs-scalar P_wp gap heatmap, (f, I_port) sweep.

Aggregates sweep_results.json from sweep_f_I.py and renders TWO figures,
each sized for its LaTeX embed so every glyph clears 9 pt on the page
(lab rule 6: choose the font size for the EMBEDDED width):

  sweep_heatmap_digest.png  -- single panel: the per-element-vs-scalar
        P_wp gap HEATMAP only.  Single-column (~8.9 cm) in-column
        ``\\begin{figure}`` of the 1-page digest.  A 2-panel + colorbar
        figure cannot carry >= 9 pt fonts inside one 8 cm column (the
        font/canvas ratio needed exceeds the constrained_layout collapse
        threshold), so the digest shows just the headline gap map.

  sweep_heatmap.png         -- two panels: (a) the gap heatmap and (b)
        absolute P_wp vs I curves.  Full width (18.2 cm render) in the
        ``\\begin{figure*}`` of the full paper, embedded at
        0.85\\textwidth (~15.5 cm) where both panels stay >= 9 pt.

Each render self-verifies with check_min_font (no glyph below 9 pt at
the embed scale) and, where a legend exists, check_legend_overlap.

Usage:
  python docs/ih_esim_benchmark/plot_sweep_heatmap.py [SWEEP_DIR]
"""
import sys
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from matplotlib.lines import Line2D
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.cm import ScalarMappable
from radia_mcp.figure import (
    apply_lab_style, lab_savefig, check_legend_overlap, check_min_font,
)

HERE = Path(__file__).resolve().parent
# Data Persistence Policy: read the committed sweep data next to this
# script by default, not a transient C:/temp dir.
DEFAULT_DIR = HERE / "sweep_data_dense"

MIN_PT = 9.0
# Binding embed scale = LaTeX \includegraphics width / figure render width.
# (The digest's combined 2-panel figure is produced by plot_digest_figure.py.)
SCALE_FULL = 0.85     # 0.85\textwidth (15.5 cm) of an 18.2 cm render.


def load_grid(res_path):
    """Load sweep_results.json into (I, f) grids of P_wp."""
    with open(res_path) as f:
        data = json.load(f)
    runs = data["runs"]
    freqs = sorted({r["frequency_Hz"] for r in runs})
    currents = sorted({r["current_A"] for r in runs})
    P_scalar = np.full((len(currents), len(freqs)), np.nan)
    P_per = np.full((len(currents), len(freqs)), np.nan)
    for r in runs:
        if r.get("P_wp_W") is None or "error" in r:
            continue
        i = currents.index(r["current_A"])
        j = freqs.index(r["frequency_Hz"])
        if r["per_panel"]:
            P_per[i, j] = r["P_wp_W"]
        else:
            P_scalar[i, j] = r["P_wp_W"]
    eps = 1e-30
    ratio = P_per / np.where(P_scalar > eps, P_scalar, np.nan)
    gap_pct = (ratio - 1.0) * 100
    return currents, freqs, P_scalar, P_per, ratio, gap_pct


def _draw_heatmap(ax, fig, currents, freqs, gap_pct, *, ann_pt, tick_pt,
                  cb_label=r"$(P_{\mathrm{per}}/P_{\mathrm{uni}} - 1)"
                           r" \times 100$ (\%)"):
    """Gap-% heatmap with annotated cells and a colorbar."""
    im = ax.imshow(gap_pct, aspect="auto", origin="lower", cmap="RdYlGn_r",
                   extent=[0.5, len(freqs)+0.5, 0.5, len(currents)+0.5])
    ax.set_xticks(range(1, len(freqs)+1))
    ax.set_xticklabels([f"{f/1e3:.0f}" for f in freqs])
    ax.set_yticks(range(1, len(currents)+1))
    ax.set_yticklabels([f"{I:.0f}" for I in currents])
    ax.set_xlabel(r"frequency (kHz)")
    ax.set_ylabel(r"$I_{\mathrm{port}}$ (A)")
    cb = fig.colorbar(im, ax=ax)
    cb.set_label(cb_label)
    cb.ax.tick_params(labelsize=tick_pt)
    for i in range(len(currents)):
        for j in range(len(freqs)):
            if not np.isnan(gap_pct[i, j]):
                ax.text(j + 1, i + 1, f"{gap_pct[i, j]:+.0f}",
                        ha="center", va="center", fontsize=ann_pt,
                        color="white" if abs(gap_pct[i, j]) > 25 else "black")
    return im


def _draw_curves(ax, fig, currents, freqs, P_scalar, P_per, *, tick_pt):
    """P_wp vs I, viridis per frequency; method legend above the axes."""
    cols = [plt.cm.viridis(j / max(len(freqs)-1, 1)) for j in range(len(freqs))]
    for j in range(len(freqs)):
        fl = f"{freqs[j]/1e3:.0f}kHz"
        ax.loglog(currents, P_scalar[:, j], "o--", color=cols[j], alpha=0.7,
                  label=f"{fl}-s")
        ax.loglog(currents, P_per[:, j], "s-", color=cols[j], label=fl)
    ax.set_xlabel(r"$I_{\mathrm{port}}$ (A)")
    ax.set_ylabel(r"$P_{\mathrm{wp}}$ (W)")
    cmap_f = ListedColormap(cols)
    norm = BoundaryNorm(np.arange(len(freqs)+1) - 0.5, cmap_f.N)
    sm = ScalarMappable(cmap=cmap_f, norm=norm)
    sm.set_array([])
    cbf = fig.colorbar(sm, ax=ax, ticks=range(len(freqs)), pad=0.02,
                       fraction=0.06)
    cbf.ax.set_yticklabels([f"{f/1e3:.0f}" for f in freqs])
    cbf.set_label(r"$f$ (kHz)")
    cbf.ax.tick_params(labelsize=tick_pt)
    method_handles = [
        Line2D([0], [0], color="0.35", ls="--", marker="o", mfc="none",
               label="uniform"),
        Line2D([0], [0], color="0.35", ls="-", marker="s",
               label="per-element"),
    ]
    ax.legend(handles=method_handles, loc="lower center",
              bbox_to_anchor=(0.5, 1.0), ncol=2, frameon=False,
              handlelength=1.6, columnspacing=1.4, borderaxespad=0.3)
    ax.grid(linestyle=":", alpha=0.5)


def _report(tag, fig, embed_scale, ax_with_legend=None):
    fig.canvas.draw()
    if ax_with_legend is not None:
        ov = check_legend_overlap(ax_with_legend)
        print(f"  [{tag}] legend overlap:",
              "CLEAN" if not ov else
              [(o["label"], f"{o['fraction']*100:.0f}%") for o in ov])
    bad = check_min_font(fig, min_pt=MIN_PT, embed_scale=embed_scale)
    if bad:
        print(f"  [{tag}] !! FONT < {MIN_PT} pt at embed_scale={embed_scale}:")
        for b in bad:
            print(f"      {b['kind']:7s} {b['text']!r}: "
                  f"render {b['render_pt']} -> visible {b['visible_pt']} pt")
    else:
        print(f"  [{tag}] min-font: all text >= {MIN_PT} pt at "
              f"embed_scale={embed_scale} (CLEAN)")


def render_full(out_stem, currents, freqs, P_scalar, P_per, gap_pct):
    """Two-panel (a) heatmap + (b) P_wp curves for the full paper."""
    figsize = apply_lab_style(target="paper_double_column", aspect=0.40)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize,
                                   constrained_layout=True)
    tick_pt = plt.rcParams["ytick.labelsize"]
    _draw_heatmap(ax1, fig, currents, freqs, gap_pct,
                  ann_pt=tick_pt, tick_pt=tick_pt)
    ax1.text(0.5, -0.32, "(a)", transform=ax1.transAxes, ha="center", va="top")
    _draw_curves(ax2, fig, currents, freqs, P_scalar, P_per, tick_pt=tick_pt)
    ax2.text(0.5, -0.32, "(b)", transform=ax2.transAxes, ha="center", va="top")
    _report("full 2-panel", fig, SCALE_FULL, ax_with_legend=ax2)
    lab_savefig(fig, out_stem)
    plt.close(fig)


def main():
    sweep_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_DIR
    res_path = sweep_dir / "sweep_results.json"
    if not res_path.exists():
        print(f"ERROR: {res_path} not found.  Run sweep_f_I.py first.")
        sys.exit(1)
    currents, freqs, P_scalar, P_per, ratio, gap_pct = load_grid(res_path)

    print("Per-element / Scalar ratio (rows=I, cols=f):")
    print(f"  freqs  [kHz]: {[f/1e3 for f in freqs]}")
    print(f"  currents [A]: {currents}")
    for i, I in enumerate(currents):
        print(f"    I={I:6.1f}A: ",
              " ".join(f"{ratio[i,j]:.3f}" for j in range(len(freqs))))

    print("Rendering full 2-panel (full paper, figure*)...")
    render_full(str(HERE / "sweep_heatmap"), currents, freqs,
                P_scalar, P_per, gap_pct)
    print(f"\nSaved {HERE/'sweep_heatmap.png'} (full paper, 2-panel). "
          f"Digest combined figure: see plot_digest_figure.py.")


if __name__ == "__main__":
    main()
