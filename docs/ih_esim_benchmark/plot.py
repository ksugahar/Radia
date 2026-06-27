"""Plot benchmark.py results.json -- P_wp, L, wall-time vs frequency."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def _extract(row, path, field):
    r = row["results"].get(path, {})
    if "error" in r:
        return None
    return r.get(field)


def _coerce(v, default=np.nan):
    try:
        return float(v) if v is not None else default
    except Exception:
        return default


def plot(results_json, out_pdf):
    with open(results_json, encoding="utf-8") as fh:
        d = json.load(fh)
    paths = d["paths"]
    freqs = np.array(d["frequencies_Hz"])

    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    fig.suptitle(
        f"ESIM Karl Benchmark — radia {d.get('radia_version', '?')} "
        f"on {d.get('hostname', '?')}",
        fontsize=11)

    colors = {"inductance": "tab:blue",
              "fem_kelvin": "tab:orange",
              "fem_coilmesh": "tab:green"}
    markers = {"inductance": "o", "fem_kelvin": "s", "fem_coilmesh": "D"}
    label_map = {
        "inductance": "PEEC + scalar BEM-SIBC",
        "fem_kelvin": "PEEC + FEM-HCurl + Kelvin",
        "fem_coilmesh": "Full FEM A-V + Kelvin",
    }

    # ---- P_wp (top-left) ----
    ax = axes[0, 0]
    for p in paths:
        P = []
        for row in d["sweep"]:
            v = _extract(row, p, "P_wp_W")
            if v is None:
                v = _extract(row, p, "P_total_W")
            if v is None:
                v = _extract(row, p, "P_total")
            P.append(_coerce(v))
        ax.loglog(freqs, np.array(P) * 1e3,
                  color=colors[p], marker=markers[p],
                  label=label_map[p], linewidth=1.5)
    ax.set_xlabel("Frequency [Hz]")
    ax.set_ylabel("$P_\\mathrm{wp}$ [mW]")
    ax.set_title("Workpiece dissipation")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)

    # ---- L_total (top-right) ----
    ax = axes[0, 1]
    for p in paths:
        L = []
        for row in d["sweep"]:
            v = _extract(row, p, "L_total_nH")
            L.append(_coerce(v))
        ax.semilogx(freqs, L, color=colors[p], marker=markers[p],
                    label=label_map[p], linewidth=1.5)
    ax.set_xlabel("Frequency [Hz]")
    ax.set_ylabel("$L_\\mathrm{total}$ [nH]")
    ax.set_title("Total inductance")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)

    # ---- wall time (bottom-left) ----
    ax = axes[1, 0]
    for p in paths:
        T = [_coerce(_extract(row, p, "_wall_s")) for row in d["sweep"]]
        ax.semilogx(freqs, T, color=colors[p], marker=markers[p],
                    label=label_map[p], linewidth=1.5)
    ax.set_xlabel("Frequency [Hz]")
    ax.set_ylabel("Wall time [s]")
    ax.set_title("Per-frequency cost")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)

    # ---- Z_s magnitude (bottom-right) ----
    ax = axes[1, 1]
    for p in paths:
        Z = []
        for row in d["sweep"]:
            r = row["results"].get(p, {})
            zr = _coerce(r.get("Z_s_wp_real"))
            zi = _coerce(r.get("Z_s_wp_imag"))
            Z.append(np.hypot(zr, zi))
        ax.loglog(freqs, Z, color=colors[p], marker=markers[p],
                  label=label_map[p], linewidth=1.5)
    ax.set_xlabel("Frequency [Hz]")
    ax.set_ylabel("$|Z_s|$ [$\\Omega$]")
    ax.set_title("Converged surface impedance")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)

    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_pdf.replace(".pdf", ".png"), dpi=150, bbox_inches="tight")
    print(f"wrote {out_pdf} and {out_pdf.replace('.pdf', '.png')}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--results", default=str(Path(__file__).parent / "results.json"))
    p.add_argument("--output", default=str(Path(__file__).parent / "benchmark_plot.pdf"))
    args = p.parse_args()
    if not Path(args.results).is_file():
        sys.exit(f"missing: {args.results}")
    plot(args.results, args.output)


if __name__ == "__main__":
    main()
