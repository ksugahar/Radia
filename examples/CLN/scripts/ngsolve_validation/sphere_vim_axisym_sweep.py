"""sphere_vim_axisym_sweep.py — Mesh-convergence sweep of axisym sphere VIM.

For each (N_rho, N_theta) cell count, assembles K, M, b in FP64, eigh,
forms Foster (lambda_n, g_n) at FP64 precision, then computes moments
alpha_k = (-1)^k sum_n g_n^2 * tau_n^(k+1) in mpmath (80 digits) and
Hankel-Pade Cauer extraction.

Compares per-stage tau_pair[k] against the 240-digit Mathematica analytical
reference (`sphere_analytical_cln_highstage.json`).

Output:
  - per-mesh JSON file `sphere_vim_axisym_sweep_NrxNt.json`
  - aggregated comparison table to stdout
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
    assemble_K_M_b,
    hankel_pade_cauer,
    R_SPHERE,
    SIGMA,
    MU0,
    B0,
)

mp.mp.dps = 80


def load_analytical(here: Path):
    p = here / "sphere_analytical_cln_highstage.json"
    if not p.exists():
        raise FileNotFoundError(f"Analytical reference not found: {p}")
    d = json.load(open(p))
    return {
        "stoll_tau_us": [float(x) for x in d["stoll_tau_us"]],
        "stoll_a_n": [float(x) for x in d["stoll_a_n"]],
        "cauer_tau_us": [float(x) for x in d["cauer_tau_us"]],
        "R0": float(d["R0"]),
        "n_reliable": int(d["n_reliable"]),
    }


def run_one_mesh(N_rho, N_theta, n_quad=2, n_stages=15):
    n_cells = N_rho * N_theta
    print(f"\n--- N_rho={N_rho}, N_theta={N_theta}, n_cells={n_cells}, "
          f"n_quad={n_quad} ---")
    t0 = time.time()
    K, M, b = assemble_K_M_b(N_rho, N_theta, n_quad=n_quad)
    t_asm = time.time() - t0
    print(f"  K, M, b assembled in {t_asm:.1f}s")

    t0 = time.time()
    eigvals, eigvecs = eigh(K, M)
    t_eig = time.time() - t0
    print(f"  eigh(K, M): {t_eig:.2f}s, "
          f"min/max λ = {eigvals.min():.3e}/{eigvals.max():.3e}")

    # M-normalize eigenvectors
    for k in range(len(eigvals)):
        v = eigvecs[:, k]
        nrm = math.sqrt(float(v @ M @ v))
        if nrm > 1e-30:
            eigvecs[:, k] = v / nrm

    # Foster: tau_n = sigma * lambda_n, g_n = v_n^T b
    tau_s_arr = SIGMA * eigvals  # in seconds
    g_arr = eigvecs.T @ b

    # Keep only positive tau, sort by |g|^2 descending
    pos = tau_s_arr > 1e-30
    tau_s = tau_s_arr[pos]
    g = g_arr[pos]
    g2 = g * g
    order = np.argsort(-g2)
    n_use = min(150, np.sum(g2 > 1e-32))
    idx = order[:n_use]

    # Moments in mpmath: alpha_k = (-1)^k sum_n g_n^2 * tau_n^(k+1)
    tau_mp = [mp.mpf(float(tau_s[i])) for i in idx]
    g2_mp = [mp.mpf(float(g2[i])) for i in idx]

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
        R_v = p[2 * k]
        Linv = p[2 * k + 1]
        if abs(Linv) < mp.mpf(10) ** (-50) or R_v == 0:
            break
        Lv = 1 / Linv
        tau_p_us = float(Lv / R_v) * 1e6
        rungs.append({
            "k": k,
            "R_2k": float(R_v),
            "L_2k_plus_1": float(Lv),
            "tau_pair_us": tau_p_us,
        })

    leading_tau_us = float(tau_s_arr[order[0] if pos.all() else
                                     np.argsort(-(g_arr * g_arr))[0]]) * 1e6
    # safer: top g² mode
    leading_idx_full = np.argsort(-(g_arr * g_arr))[0]
    leading_tau_us = float(tau_s_arr[leading_idx_full]) * 1e6

    print(f"  Leading τ (top |g|^2): {leading_tau_us:.4f} μs")
    print(f"  Extracted Cauer stages: {len(rungs)}")
    return {
        "N_rho": N_rho, "N_theta": N_theta, "n_cells": n_cells,
        "n_quad": n_quad,
        "t_assemble_s": t_asm, "t_eig_s": t_eig,
        "leading_tau_us": leading_tau_us,
        "cauer_rungs": rungs,
    }


def gap_pct(num, ref):
    return (num - ref) / ref * 100.0


def main():
    here = Path(__file__).parent
    ana = load_analytical(here)
    print("=" * 78)
    print(f"Sphere VIM axisym mesh-convergence sweep vs Stoll analytical")
    print(f"  R={R_SPHERE * 1000:.1f} mm, sigma={SIGMA:.3e} S/m")
    print(f"  Analytical R_0 = {ana['R0']:.6f}")
    print(f"  Analytical Foster τ_1 = {ana['stoll_tau_us'][0]:.6f} μs")
    print(f"  Analytical Cauer τ_pair[0] = {ana['cauer_tau_us'][0]:.6f} μs")
    print("=" * 78)

    # Cell counts: 480, 1080, 1920, 2880 (~12 min budget on this machine)
    mesh_configs = [
        (20, 24),   # 480
        (30, 36),   # 1080 (existing baseline)
        (40, 48),   # 1920
        (48, 60),   # 2880
    ]

    all_results = []
    for N_rho, N_theta in mesh_configs:
        try:
            res = run_one_mesh(N_rho, N_theta, n_quad=2, n_stages=15)
            all_results.append(res)
            out_path = (here /
                        f"sphere_vim_axisym_sweep_{N_rho}x{N_theta}.json")
            out_path.write_text(json.dumps(res, indent=2), encoding="utf-8")
            print(f"  Saved {out_path.name}")
        except Exception as e:
            print(f"  FAIL: {type(e).__name__}: {e}")

    # === Aggregated comparison table ===
    print()
    print("=" * 78)
    print("Cauer rung comparison: VIM τ_pair[k] vs Stoll analytical (μs)")
    print("=" * 78)
    n_show = 10
    hdr = ["k", "Stoll"] + [f"{r['n_cells']}c" for r in all_results]
    print("  " + " | ".join(f"{h:>12s}" for h in hdr))
    print("  " + "-+-".join("-" * 12 for _ in hdr))
    for k in range(n_show):
        if k >= len(ana["cauer_tau_us"]):
            break
        ref = ana["cauer_tau_us"][k]
        row = [f"{k:>12d}", f"{ref:>12.5f}"]
        for r in all_results:
            if k < len(r["cauer_rungs"]):
                row.append(f"{r['cauer_rungs'][k]['tau_pair_us']:>12.5f}")
            else:
                row.append(f"{'--':>12s}")
        print("  " + " | ".join(row))

    print()
    print("Gap to Stoll [%]:")
    print("  " + " | ".join(f"{h:>12s}" for h in hdr))
    print("  " + "-+-".join("-" * 12 for _ in hdr))
    for k in range(n_show):
        if k >= len(ana["cauer_tau_us"]):
            break
        ref = ana["cauer_tau_us"][k]
        row = [f"{k:>12d}", f"{0.0:>12.4f}"]
        for r in all_results:
            if k < len(r["cauer_rungs"]):
                gp = gap_pct(r["cauer_rungs"][k]["tau_pair_us"], ref)
                row.append(f"{gp:>+12.4f}")
            else:
                row.append(f"{'--':>12s}")
        print("  " + " | ".join(row))

    summary = {
        "analytical": ana,
        "mesh_results": all_results,
    }
    (here / "sphere_vim_axisym_sweep_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8")
    print("\nSaved sphere_vim_axisym_sweep_summary.json")


if __name__ == "__main__":
    main()
