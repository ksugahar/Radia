"""
compare_peec_dense_vs_hacapk.py -- Produce scaling comparison tables.

Reads results_bench_peec_dense.json and results_bench_peec_hacapk.json from
this directory and emits:

  - stdout: plain-text table (grouped by N_fil, side-by-side)
  - --md PATH: Markdown table (for the paper draft)
  - --tex PATH: LaTeX booktabs table (for the paper)

The two benchmarks share the helical solenoid geometry (20 turns x 25
seg/turn x nwinc^2 filaments) so rows match on N_fil.  Dense entries
missing beyond the dense benchmark's last case are shown as "--".
"""

import argparse
import json
import os
from typing import Optional


HERE = os.path.dirname(os.path.abspath(__file__))
DENSE_JSON = os.path.join(HERE, "results_bench_peec_dense.json")
HACAPK_JSON = os.path.join(HERE, "results_bench_peec_hacapk.json")


def _load(path: str) -> dict:
    with open(path) as f:
        return json.load(f)


def _by_n_fil(results: list) -> dict:
    return {int(r["n_fil"]): r for r in results if "n_fil" in r}


def _fmt_time(s: Optional[float]) -> str:
    if s is None:
        return "--"
    if s < 1.0:
        return f"{s * 1e3:.1f} ms"
    if s < 60.0:
        return f"{s:.2f} s"
    return f"{s / 60.0:.2f} min"


def _fmt_mb(mb: Optional[float]) -> str:
    if mb is None:
        return "--"
    if mb >= 1024.0:
        return f"{mb / 1024.0:.2f} GB"
    return f"{mb:.0f} MB"


def _fmt_pct(x: Optional[float]) -> str:
    if x is None:
        return "--"
    return f"{x * 100:.1f}%"


def build_rows(dense: dict, hacapk: dict) -> list:
    d_by_n = _by_n_fil(dense["results"])
    h_by_n = _by_n_fil(hacapk["results"])
    all_n = sorted(set(d_by_n.keys()) | set(h_by_n.keys()))

    rows = []
    for n in all_n:
        d = d_by_n.get(n)
        h = h_by_n.get(n)
        row = {
            "n_fil": n,
            # Dense
            "dense_L_mb": d["L_matrix_mb"] if d else None,
            "dense_t_build": d["t_build"] if d else None,
            "dense_t_solve": d["t_solve"] if d else None,
            "dense_peak_mb": d["peak_memory_mb"] if d else None,
            # HACApK
            "hmat_mb": h["hmatrix_mb"] if h else None,
            "compression": h["compression"] if h else None,
            "hacapk_t_build": h["t_build"] if h else None,
            "hacapk_t_solve": h["t_solve"] if h else None,
            "hacapk_iters": h["iterations"] if h else None,
            "hacapk_residual": h["final_residual"] if h else None,
            "hacapk_peak_mb": h["peak_memory_mb"] if h else None,
            # Speedups
            "build_speedup": (d["t_build"] / h["t_build"]) if (d and h) else None,
            "solve_speedup": (d["t_solve"] / h["t_solve"]) if (d and h) else None,
            "mem_ratio": (h["hmatrix_mb"] / d["L_matrix_mb"])
                          if (d and h) else None,
        }
        rows.append(row)
    return rows


# ---------------------------------------------------------------------------
# Plain text
# ---------------------------------------------------------------------------

def print_text(rows: list) -> None:
    header = (
        f"{'N_fil':>6} | "
        f"{'dense L':>8} | {'dense bld':>10} | {'dense slv':>10} | {'dense peak':>10} | "
        f"{'H-mat':>8} | {'cmpr':>5} | {'H bld':>8} | {'H slv':>8} | {'iters':>5} | "
        f"{'res':>7} | {'H peak':>9}"
    )
    sep = "-" * len(header)
    print(sep)
    print(header)
    print(sep)
    for r in rows:
        iters = str(r['hacapk_iters']) if r['hacapk_iters'] is not None else '--'
        res = f"{r['hacapk_residual']:.1e}" if r['hacapk_residual'] is not None else '--'
        print(
            f"{r['n_fil']:>6} | "
            f"{_fmt_mb(r['dense_L_mb']):>8} | {_fmt_time(r['dense_t_build']):>10} | "
            f"{_fmt_time(r['dense_t_solve']):>10} | {_fmt_mb(r['dense_peak_mb']):>10} | "
            f"{_fmt_mb(r['hmat_mb']):>8} | {_fmt_pct(r['compression']):>5} | "
            f"{_fmt_time(r['hacapk_t_build']):>8} | {_fmt_time(r['hacapk_t_solve']):>8} | "
            f"{iters:>5} | {res:>7} | {_fmt_mb(r['hacapk_peak_mb']):>9}"
        )
    print(sep)

    # Speedup summary
    print()
    print("Speedup summary (dense / HACApK):")
    print(f"{'N_fil':>6} | {'build':>10} | {'solve':>10} | {'H-mat / dense L':>18}")
    print("-" * 55)
    for r in rows:
        bs = f"{r['build_speedup']:.1f}x" if r['build_speedup'] else "--"
        ss = f"{r['solve_speedup']:.1f}x" if r['solve_speedup'] else "--"
        mr = _fmt_pct(r['mem_ratio']) if r['mem_ratio'] else "--"
        print(f"{r['n_fil']:>6} | {bs:>10} | {ss:>10} | {mr:>18}")


# ---------------------------------------------------------------------------
# Markdown
# ---------------------------------------------------------------------------

def write_markdown(rows: list, path: str) -> None:
    lines = []
    lines.append("# PEEC: dense LU vs HACApK BiCGSTAB scaling")
    lines.append("")
    lines.append("Helical solenoid, 20 turns x 25 seg/turn, nwinc^2 sub-filaments. "
                 "Freq 10 kHz. BiCGSTAB rtol 1e-6 with Jacobi preconditioner.")
    lines.append("")
    lines.append("| N_fil | Dense L | Dense build | Dense solve | Dense peak mem | "
                 "H-mat | Compression | H-mat build | H-mat solve | Iters | Residual | H-mat peak mem |")
    lines.append("|------:|--------:|------------:|------------:|---------------:|"
                 "------:|------------:|------------:|------------:|------:|---------:|---------------:|")
    for r in rows:
        def T(x): return _fmt_time(x)
        def M(x): return _fmt_mb(x)
        def P(x): return _fmt_pct(x)
        iters = str(r['hacapk_iters']) if r['hacapk_iters'] is not None else '--'
        res = f"{r['hacapk_residual']:.1e}" if r['hacapk_residual'] is not None else '--'
        lines.append(f"| {r['n_fil']} | {M(r['dense_L_mb'])} | {T(r['dense_t_build'])} "
                     f"| {T(r['dense_t_solve'])} | {M(r['dense_peak_mb'])} "
                     f"| {M(r['hmat_mb'])} | {P(r['compression'])} "
                     f"| {T(r['hacapk_t_build'])} | {T(r['hacapk_t_solve'])} "
                     f"| {iters} | {res} | {M(r['hacapk_peak_mb'])} |")
    lines.append("")
    lines.append("## Speedup")
    lines.append("")
    lines.append("| N_fil | Build speedup | Solve speedup | H-mat / dense L memory |")
    lines.append("|------:|--------------:|--------------:|-----------------------:|")
    for r in rows:
        bs = f"{r['build_speedup']:.1f}x" if r['build_speedup'] else "--"
        ss = f"{r['solve_speedup']:.1f}x" if r['solve_speedup'] else "--"
        mr = _fmt_pct(r['mem_ratio']) if r['mem_ratio'] else "--"
        lines.append(f"| {r['n_fil']} | {bs} | {ss} | {mr} |")

    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"[markdown] wrote {path}")


# ---------------------------------------------------------------------------
# LaTeX
# ---------------------------------------------------------------------------

def write_latex(rows: list, path: str) -> None:
    lines = []
    lines.append(r"\begin{table}[t]")
    lines.append(r"\centering")
    lines.append(r"\caption{PEEC filament-level scaling: dense LU vs HACApK BiCGSTAB. "
                 r"Helical solenoid, 20 turns, 25 seg/turn, $n_{\text{winc}}^2$ parallel "
                 r"sub-filaments. 10~kHz. BiCGSTAB rtol $10^{-6}$.}")
    lines.append(r"\label{tab:peec_hacapk_scaling}")
    lines.append(r"\begin{tabular}{rrrrrrrrrrr}")
    lines.append(r"\toprule")
    lines.append(r" & \multicolumn{4}{c}{Dense} & \multicolumn{6}{c}{HACApK} \\")
    lines.append(r"\cmidrule(lr){2-5} \cmidrule(lr){6-11}")
    lines.append(r"$N$ & $L$ & build & solve & peak & H-mat & ratio & build & solve & iter & $\|r\|$ \\")
    lines.append(r"\midrule")
    for r in rows:
        def T(x): return _fmt_time(x) if x is not None else "--"
        def M(x): return _fmt_mb(x) if x is not None else "--"
        def P(x): return _fmt_pct(x) if x is not None else "--"
        iters = str(r['hacapk_iters']) if r['hacapk_iters'] is not None else '--'
        res = f"{r['hacapk_residual']:.1e}" if r['hacapk_residual'] is not None else '--'
        lines.append(
            f"{r['n_fil']} & {M(r['dense_L_mb'])} & {T(r['dense_t_build'])} & "
            f"{T(r['dense_t_solve'])} & {M(r['dense_peak_mb'])} & "
            f"{M(r['hmat_mb'])} & {P(r['compression'])} & "
            f"{T(r['hacapk_t_build'])} & {T(r['hacapk_t_solve'])} & "
            f"{iters} & {res} \\\\"
        )
    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    lines.append(r"\end{table}")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"[latex] wrote {path}")


# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dense", default=DENSE_JSON)
    parser.add_argument("--hacapk", default=HACAPK_JSON)
    parser.add_argument("--md", default=None, help="path to write Markdown")
    parser.add_argument("--tex", default=None, help="path to write LaTeX")
    args = parser.parse_args()

    dense = _load(args.dense)
    hacapk = _load(args.hacapk)

    rows = build_rows(dense, hacapk)
    print_text(rows)

    if args.md:
        write_markdown(rows, args.md)
    if args.tex:
        write_latex(rows, args.tex)


if __name__ == "__main__":
    main()
