"""sphere_vim_axisym_quadsweep.py — Probe quadrature-order sensitivity.

Holds the mesh fixed at N_rho=20, N_theta=24 (480 cells) and varies the
inner Gauss-Legendre quadrature order n_quad in {2, 4, 6, 8}.

If raising n_quad shrinks the Stoll-gap monotonically, the issue is
smooth-Gauss undersampling.  If it saturates with no improvement, the
real bottleneck is the log-type singularity of the axisymmetric Green's
function on the self-cell (m -> 1 in K(m), elliptic-integral kernel),
which requires a polar/Duffy clustering quadrature.

For each n_quad we:
  - Assemble K, M, b in FP64
  - Solve eigh(K, M) for Foster (lambda, g) at FP64
  - Form moments alpha_k = (-1)^k sum_n g^2 tau^(k+1) in mpmath (n_use=150)
  - Hankel-Pade Cauer extraction (15 stages requested)
  - Compare tau_pair[0..5] vs Stoll analytical
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

mp.mp.dps = 60


def load_analytical(here):
    p = here / "sphere_analytical_cln_highstage.json"
    d = json.load(open(p))
    return {
        "stoll_tau_us": [float(x) for x in d["stoll_tau_us"]],
        "cauer_tau_us": [float(x) for x in d["cauer_tau_us"]],
        "R0": float(d["R0"]),
    }


def run_one(N_rho, N_theta, n_quad):
    print(f"\n--- n_quad={n_quad} ---")
    t0 = time.time()
    K, M, b = assemble_K_M_b(N_rho, N_theta, n_quad=n_quad)
    t_a = time.time() - t0
    print(f"  K, M, b assembled in {t_a:.1f}s")

    # Self-cell K_ii statistics (these dominate the diagonal)
    diag_K = np.diag(K)
    print(f"  K_ii: min={diag_K.min():.3e}, max={diag_K.max():.3e}, "
          f"mean={diag_K.mean():.3e}")

    eigvals, eigvecs = eigh(K, M)
    for k in range(eigvecs.shape[1]):
        v = eigvecs[:, k]
        nrm = math.sqrt(float(v @ M @ v))
        if nrm > 1e-30:
            eigvecs[:, k] = v / nrm

    tau_arr = SIGMA * eigvals
    g_arr = eigvecs.T @ b
    g2 = g_arr * g_arr

    # Top |g|^2 mode (Foster leading)
    leading_idx = np.argsort(-g2)[0]
    leading_tau_us = float(tau_arr[leading_idx]) * 1e6

    # Moments + Hankel-Pade
    pos = tau_arr > 1e-30
    tau_s = tau_arr[pos]
    g = g_arr[pos]
    g2_pos = g * g
    order = np.argsort(-g2_pos)
    n_use = min(150, np.sum(g2_pos > 1e-32))
    idx = order[:n_use]

    tau_mp = [mp.mpf(float(tau_s[i])) for i in idx]
    g2_mp = [mp.mpf(float(g2_pos[i])) for i in idx]

    n_stages = 15
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
        tau_p_us = float(Lv / Rv) * 1e6
        rungs.append({"k": k, "tau_pair_us": tau_p_us})

    return {
        "n_quad": n_quad,
        "t_assemble_s": t_a,
        "K_diag_stats": {"min": float(diag_K.min()), "max": float(diag_K.max()),
                          "mean": float(diag_K.mean())},
        "leading_tau_us": leading_tau_us,
        "cauer_rungs": rungs,
    }


def main():
    here = Path(__file__).parent
    ana = load_analytical(here)
    N_rho, N_theta = 20, 24
    print("=" * 78)
    print(f"Sphere VIM quadrature-order sweep, fixed mesh "
          f"{N_rho}x{N_theta}={N_rho*N_theta} cells")
    print(f"  Stoll Foster τ_1 = {ana['stoll_tau_us'][0]:.4f} μs")
    print(f"  Stoll Cauer τ_pair[0] = {ana['cauer_tau_us'][0]:.4f} μs")
    print("=" * 78)

    results = []
    for nq in [2, 3, 4, 5]:
        try:
            r = run_one(N_rho, N_theta, nq)
            results.append(r)
        except Exception as e:
            print(f"  FAIL nq={nq}: {e}")

    # Aggregated table
    print()
    print("=" * 78)
    print("Leading Foster τ_1 [μs] vs Stoll = "
          f"{ana['stoll_tau_us'][0]:.4f}")
    print("=" * 78)
    print(f"  {'n_quad':>8s} | {'time_s':>8s} | {'tau_1':>12s} | {'gap [%]':>10s}")
    for r in results:
        gap = (r["leading_tau_us"] - ana["stoll_tau_us"][0]) / ana["stoll_tau_us"][0] * 100
        print(f"  {r['n_quad']:>8d} | {r['t_assemble_s']:>8.1f} | "
              f"{r['leading_tau_us']:>12.6f} | {gap:>+10.4f}")

    print()
    print("=" * 78)
    print("Cauer τ_pair[k] [μs]")
    print("=" * 78)
    print(f"  {'k':>4s} | {'Stoll':>12s} | " +
          " | ".join(f"{'nq=' + str(r['n_quad']):>12s}" for r in results))
    print(f"  {'-'*4}-+-{'-'*12}-+-" +
          "-+-".join("-" * 12 for _ in results))
    for k in range(8):
        if k >= len(ana["cauer_tau_us"]):
            break
        ref = ana["cauer_tau_us"][k]
        row = [f"{k:>4d}", f"{ref:>12.5f}"]
        for r in results:
            if k < len(r["cauer_rungs"]):
                row.append(f"{r['cauer_rungs'][k]['tau_pair_us']:>12.5f}")
            else:
                row.append(f"{'--':>12s}")
        print("  " + " | ".join(row))

    print()
    print("Gap to Stoll [%]:")
    print(f"  {'k':>4s} | {'Stoll':>12s} | " +
          " | ".join(f"{'nq=' + str(r['n_quad']):>12s}" for r in results))
    print(f"  {'-'*4}-+-{'-'*12}-+-" +
          "-+-".join("-" * 12 for _ in results))
    for k in range(8):
        if k >= len(ana["cauer_tau_us"]):
            break
        ref = ana["cauer_tau_us"][k]
        row = [f"{k:>4d}", f"{0.0:>+12.4f}"]
        for r in results:
            if k < len(r["cauer_rungs"]):
                g = (r["cauer_rungs"][k]["tau_pair_us"] - ref) / ref * 100
                row.append(f"{g:>+12.4f}")
            else:
                row.append(f"{'--':>12s}")
        print("  " + " | ".join(row))

    print()
    print("K diag stats (self-cell K_ii):")
    print(f"  {'n_quad':>8s} | {'K_ii min':>14s} | {'K_ii max':>14s} | {'K_ii mean':>14s}")
    for r in results:
        s = r["K_diag_stats"]
        print(f"  {r['n_quad']:>8d} | {s['min']:>14.4e} | {s['max']:>14.4e} | "
              f"{s['mean']:>14.4e}")

    out_path = here / "sphere_vim_axisym_quadsweep_results.json"
    out = {
        "N_rho": N_rho, "N_theta": N_theta,
        "n_cells": N_rho * N_theta,
        "analytical": ana,
        "by_nq": results,
    }
    out_path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"\nSaved: {out_path.name}")


if __name__ == "__main__":
    main()
