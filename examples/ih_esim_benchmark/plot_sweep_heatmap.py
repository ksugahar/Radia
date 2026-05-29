"""Per-element-vs-scalar P_wp gap heatmap, (f, I_port) sweep.

Aggregates sweep_results.json from sweep_f_I.py and plots a 2D
heatmap of the per-element-vs-scalar P_wp ratio over the (frequency,
port current) operating regime.  Identifies where per-element ESIM
matters most for IH design.

Output:
  sweep_heatmap.png  -- 2-panel figure: (a) P_wp_per/P_wp_scalar
                       ratio heatmap, (b) absolute P_wp per/scalar
                       curves vs I at each frequency.

Usage:
  python examples/ih_esim_benchmark/plot_sweep_heatmap.py [SWEEP_DIR]
"""
import sys
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from matplotlib.lines import Line2D
from mcp_server_document.graph.tools import (
    apply_lab_style, lab_savefig, check_legend_overlap,
)

HERE = Path(__file__).resolve().parent
# Data Persistence Policy: read the committed sweep data next to this
# script by default, not a transient C:/temp dir.
DEFAULT_DIR = HERE / "sweep_data"


def main():
    sweep_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_DIR
    res_path = sweep_dir / "sweep_results.json"
    if not res_path.exists():
        print(f"ERROR: {res_path} not found.  Run sweep_f_I.py first.")
        sys.exit(1)
    with open(res_path) as f:
        data = json.load(f)
    runs = data["runs"]

    # Collect into (I, f, per_panel) -> P_wp.
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

    # Gap ratio.  Avoid div by zero.
    eps = 1e-30
    ratio = P_per / np.where(P_scalar > eps, P_scalar, np.nan)
    gap_pct = (ratio - 1.0) * 100

    print(f"Per-element / Scalar ratio (rows=I, cols=f):")
    print(f"  freqs  [kHz]: {[f/1e3 for f in freqs]}")
    print(f"  currents [A]: {currents}")
    print(f"  ratio:")
    for i, I in enumerate(currents):
        print(f"    I={I:6.1f}A: ", " ".join(f"{ratio[i,j]:.3f}" for j in range(len(freqs))))

    figsize = apply_lab_style(target="paper_double_column", aspect=0.45)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize,
                                    constrained_layout=True)

    # Panel (a): heatmap of gap %
    im = ax1.imshow(gap_pct, aspect="auto", origin="lower",
                     cmap="RdYlGn_r",
                     extent=[0.5, len(freqs)+0.5,
                              0.5, len(currents)+0.5])
    ax1.set_xticks(range(1, len(freqs)+1))
    ax1.set_xticklabels([f"{f/1e3:.0f}" for f in freqs])
    ax1.set_yticks(range(1, len(currents)+1))
    ax1.set_yticklabels([f"{I:.0f}" for I in currents])
    ax1.set_xlabel(r"frequency (kHz)")
    ax1.set_ylabel(r"$I_{\mathrm{port}}$ (A)")
    cb = plt.colorbar(im, ax=ax1)
    cb.set_label(r"$(P_{\mathrm{per}}/P_{\mathrm{scalar}} - 1) \times 100$ (\%)")
    # Annotate cells with the numeric ratio.
    for i in range(len(currents)):
        for j in range(len(freqs)):
            if not np.isnan(gap_pct[i, j]):
                ax1.text(j + 1, i + 1, f"{gap_pct[i, j]:+.0f}",
                          ha="center", va="center",
                          color="white" if abs(gap_pct[i, j]) > 25 else "black")
    ax1.text(0.5, -0.27, "(a)", transform=ax1.transAxes,
              ha="center", va="top")

    # Panel (b): absolute P_wp curves vs I, one line per frequency.
    # Frequency is encoded by viridis colour (named in the caption); the
    # in-figure legend distinguishes only the two METHODS so it stays
    # small enough not to cover the curves at 8 cm embed width.
    for j, f_Hz in enumerate(freqs):
        col = plt.cm.viridis(j / max(len(freqs)-1, 1))
        ax2.loglog(currents, P_scalar[:, j], "o--", color=col, alpha=0.6)
        ax2.loglog(currents, P_per[:, j], "s-", color=col)
    ax2.set_xlabel(r"$I_{\mathrm{port}}$ (A)")
    ax2.set_ylabel(r"$P_{\mathrm{wp}}$ (W)")
    # Two compact legends so all FOUR frequency curves AND both methods are
    # identified without covering the data: rising log-log curves leave the
    # upper-left and lower-right corners empty.
    #  (1) frequency -> viridis colour (4 entries), lower-right.
    freq_handles = [
        Line2D([0], [0], color=plt.cm.viridis(j / max(len(freqs)-1, 1)),
               ls="-", label=fr"{f_Hz/1e3:.0f} kHz")
        for j, f_Hz in enumerate(freqs)
    ]
    leg_freq = ax2.legend(handles=freq_handles, fontsize=7, frameon=False,
                          loc="lower right", handlelength=1.4,
                          labelspacing=0.2, borderaxespad=0.4)
    ax2.add_artist(leg_freq)
    #  (2) method -> line style (2 entries, neutral grey), upper-left.
    method_handles = [
        Line2D([0], [0], color="0.35", ls="--", marker="o", mfc="none",
               label="scalar"),
        Line2D([0], [0], color="0.35", ls="-", marker="s",
               label="per-element"),
    ]
    ax2.legend(handles=method_handles, fontsize=7, frameon=False,
               loc="upper left", handlelength=1.8, labelspacing=0.2)
    fig.canvas.draw()
    # Verify neither legend covers a curve (check both placements).
    _ov = check_legend_overlap(ax2)
    print("panel (b) freq(LR)+method(UL) legends; overlap:",
          "CLEAN" if not _ov else
          [(o["label"], f"{o['fraction']*100:.0f}%") for o in _ov])
    ax2.grid(linestyle=":", alpha=0.5)
    ax2.text(0.5, -0.27, "(b)", transform=ax2.transAxes,
              ha="center", va="top")

    out_png = HERE / "sweep_heatmap.png"
    lab_savefig(fig, str(out_png.with_suffix("")))
    print(f"\nSaved {out_png}")


if __name__ == "__main__":
    main()
