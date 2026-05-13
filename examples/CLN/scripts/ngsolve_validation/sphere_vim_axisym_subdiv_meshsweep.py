"""sphere_vim_axisym_subdiv_meshsweep.py — Mesh-convergence study with
self-cell subdivision (ns=2) for proper log-singular quadrature.

For each (N_rho, N_theta) cell count, runs the assembly with ns=2 self-
cell subdivision (proven sufficient on 480 cells), then extracts Cauer
rungs via FP64 eigh + mpmath Hankel-Padé.  Compares against Stoll
analytical 40-stage reference.

This replaces sphere_vim_axisym_sweep.py (which used ns=1 = no
subdivision) as the canonical mesh-convergence study.
"""
from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, str(Path(__file__).parent))

import numpy as np
from scipy.linalg import eigh
import mpmath as mp

from sphere_vim_axisym import (
    R_SPHERE, SIGMA, MU0, B0,
    hankel_pade_cauer,
)
from sphere_vim_axisym_subdiv import assemble_K_M_b_subdiv

mp.mp.dps = 60


def run_one_mesh(N_rho, N_theta, n_quad=2, n_sub_self=2, n_stages=15):
    n_cells = N_rho * N_theta
    print(f"\n--- N_rho={N_rho}, N_theta={N_theta}, n_cells={n_cells}, "
          f"n_sub_self={n_sub_self} ---")
    t0 = time.time()
    K, M, b = assemble_K_M_b_subdiv(N_rho, N_theta,
                                      n_quad=n_quad,
                                      n_sub_self=n_sub_self)
    t_a = time.time() - t0
    print(f"  assembled in {t_a:.1f}s")

    t0 = time.time()
    eigvals, eigvecs = eigh(K, M)
    print(f"  eigh: {time.time() - t0:.2f}s, "
          f"min/max λ = {eigvals.min():.3e}/{eigvals.max():.3e}")

    for k in range(eigvecs.shape[1]):
        v = eigvecs[:, k]
        nrm = math.sqrt(float(v @ M @ v))
        if nrm > 1e-30:
            eigvecs[:, k] = v / nrm

    tau_arr = SIGMA * eigvals
    g_arr = eigvecs.T @ b
    g2 = g_arr * g_arr
    leading_idx = np.argsort(-g2)[0]
    leading_tau_us = float(tau_arr[leading_idx]) * 1e6

    pos = tau_arr > 1e-30
    tau_s = tau_arr[pos]
    g = g_arr[pos]
    g2_pos = g * g
    order = np.argsort(-g2_pos)
    n_use = min(150, np.sum(g2_pos > 1e-32))
    idx = order[:n_use]

    tau_mp = [mp.mpf(float(tau_s[i])) for i in idx]
    g2_mp = [mp.mpf(float(g2_pos[i])) for i in idx]
    n_moments = min(2 * n_stages + 4, 2 * n_use - 2)
    alphas = []
    for k in range(n_moments):
        a = mp.mpf(0)
        for gv, tv in zip(g2_mp, tau_mp):
            a += gv * tv ** (k + 1)
        alphas.append((mp.mpf(-1)) ** k * a)

    p = hankel_pade_cauer(alphas, n_stages)
    rungs = []
    for k in range(len(p) // 2):
        Rv = p[2 * k]
        Linv = p[2 * k + 1]
        if abs(Linv) < mp.mpf(10) ** (-50) or Rv == 0:
            break
        Lv = 1 / Linv
        rungs.append({"k": k, "tau_pair_us": float(Lv / Rv) * 1e6})

    return {
        "N_rho": N_rho, "N_theta": N_theta, "n_cells": n_cells,
        "t_assemble_s": t_a,
        "leading_tau_us": leading_tau_us,
        "cauer_rungs": rungs,
    }


def main():
    here = Path(__file__).parent
    ana = json.load(open(here / "sphere_analytical_cln_highstage.json"))
    stoll_tau1 = float(ana["stoll_tau_us"][0])
    cauer_ref = [float(x) for x in ana["cauer_tau_us"]]

    print("=" * 78)
    print(f"Sphere VIM mesh-convergence with self-cell subdiv (ns=2)")
    print(f"  Stoll Foster τ_1 = {stoll_tau1:.4f} μs")
    print(f"  Stoll Cauer τ_pair[0..7]: "
          f"{[round(cauer_ref[k], 4) for k in range(8)]}")
    print("=" * 78)

    mesh_configs = [
        (20, 24),   # 480
        (30, 36),   # 1080
        (40, 48),   # 1920
        (48, 60),   # 2880
    ]

    results = []
    for N_rho, N_theta in mesh_configs:
        try:
            r = run_one_mesh(N_rho, N_theta, n_quad=2, n_sub_self=2,
                             n_stages=15)
            results.append(r)
            out_path = (here /
                        f"sphere_vim_axisym_subdiv_meshsweep_"
                        f"{N_rho}x{N_theta}.json")
            out_path.write_text(json.dumps(r, indent=2), encoding="utf-8")
        except Exception as e:
            print(f"  FAIL: {type(e).__name__}: {e}")

    print()
    print("=" * 78)
    print(f"Leading Foster τ_1 [μs] vs Stoll = {stoll_tau1:.4f}")
    print("=" * 78)
    print(f"  {'cells':>8s} | {'time_s':>8s} | {'tau_1':>12s} | {'gap [%]':>10s}")
    for r in results:
        g = (r["leading_tau_us"] - stoll_tau1) / stoll_tau1 * 100
        print(f"  {r['n_cells']:>8d} | {r['t_assemble_s']:>8.1f} | "
              f"{r['leading_tau_us']:>12.6f} | {g:>+10.4f}")

    print()
    print("=" * 78)
    print("Cauer τ_pair[k] [μs] (with self-cell ns=2 subdivision)")
    print("=" * 78)
    hdr = [f"{'k':>4s}", f"{'Stoll':>12s}"] + [
        f"{str(r['n_cells']) + 'c':>12s}" for r in results]
    print("  " + " | ".join(hdr))
    print("  " + "-+-".join("-" * 12 for _ in hdr))
    for k in range(10):
        if k >= len(cauer_ref):
            break
        row = [f"{k:>4d}", f"{cauer_ref[k]:>12.5f}"]
        for r in results:
            if k < len(r["cauer_rungs"]):
                row.append(f"{r['cauer_rungs'][k]['tau_pair_us']:>12.5f}")
            else:
                row.append(f"{'--':>12s}")
        print("  " + " | ".join(row))

    print()
    print("Gap to Stoll [%]:")
    print("  " + " | ".join(hdr))
    print("  " + "-+-".join("-" * 12 for _ in hdr))
    for k in range(10):
        if k >= len(cauer_ref):
            break
        row = [f"{k:>4d}", f"{0.0:>+12.4f}"]
        for r in results:
            if k < len(r["cauer_rungs"]):
                g = (r["cauer_rungs"][k]["tau_pair_us"] - cauer_ref[k]) / cauer_ref[k] * 100
                row.append(f"{g:>+12.4f}")
            else:
                row.append(f"{'--':>12s}")
        print("  " + " | ".join(row))

    out = {
        "stoll_tau1_us": stoll_tau1,
        "stoll_cauer_tau_us": cauer_ref,
        "n_sub_self": 2, "n_quad": 2,
        "by_mesh": results,
    }
    (here / "sphere_vim_axisym_subdiv_meshsweep_summary.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8")
    print("\nSaved sphere_vim_axisym_subdiv_meshsweep_summary.json")


if __name__ == "__main__":
    main()
