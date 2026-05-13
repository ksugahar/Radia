"""sphere_vim_axisym_RL_extract.py — One-shot extraction of full Cauer R, L
table at a single mesh size with subdiv ns=2.  Used to populate the IGTE
supplement's Cauer rung comparison table for sphere.
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

from sphere_vim_axisym import R_SPHERE, SIGMA, MU0, B0, hankel_pade_cauer
from sphere_vim_axisym_subdiv import assemble_K_M_b_subdiv

mp.mp.dps = 60


def extract(N_rho, N_theta, n_sub_self=2):
    print(f"--- {N_rho}x{N_theta} = {N_rho*N_theta} cells, ns={n_sub_self} ---")
    t0 = time.time()
    K, M, b = assemble_K_M_b_subdiv(N_rho, N_theta, n_quad=2,
                                      n_sub_self=n_sub_self)
    print(f"  assembled in {time.time()-t0:.1f}s")
    eigvals, eigvecs = eigh(K, M)
    for k in range(eigvecs.shape[1]):
        v = eigvecs[:, k]
        nrm = math.sqrt(float(v @ M @ v))
        if nrm > 1e-30:
            eigvecs[:, k] = v / nrm
    tau_arr = SIGMA * eigvals
    g_arr = eigvecs.T @ b
    g2 = g_arr * g_arr
    pos = tau_arr > 1e-30
    tau_pos = tau_arr[pos]
    g2_pos = g_arr[pos] * g_arr[pos]
    order = np.argsort(-g2_pos)
    n_use = min(150, np.sum(g2_pos > 1e-32))
    idx = order[:n_use]

    tau_mp = [mp.mpf(float(tau_pos[i])) for i in idx]
    g2_mp = [mp.mpf(float(g2_pos[i])) for i in idx]
    n_stages = 12
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
        rungs.append({
            "k": k,
            "R_2k": float(Rv),
            "L_2k_plus_1": float(Lv),
            "tau_pair_us": float(Lv / Rv) * 1e6,
        })

    return {
        "n_cells": N_rho * N_theta,
        "N_rho": N_rho, "N_theta": N_theta,
        "n_sub_self": n_sub_self,
        "leading_tau_us": float(tau_arr[np.argsort(-(g_arr*g_arr))[0]]) * 1e6,
        "cauer_rungs": rungs,
    }


def main():
    here = Path(__file__).parent
    ana = json.load(open(here / "sphere_analytical_cln_highstage.json"))
    cauer_ref_R = [float(x) for x in ana["cauer_R"]]
    cauer_ref_L = [float(x) for x in ana["cauer_L"]]
    cauer_ref_tau = [float(x) for x in ana["cauer_tau_us"]]

    # Single 1080-cell extraction for the paper table
    r = extract(30, 36, n_sub_self=2)

    print()
    print("=" * 78)
    print("Cauer rung comparison: Stoll analytical vs VIM 1080 cells (subdiv ns=2)")
    print("=" * 78)
    print(f"  {'k':>3s} | {'Stoll R':>14s} | {'VIM R':>14s} | {'gap R [%]':>10s} | "
          f"{'Stoll L':>14s} | {'VIM L':>14s} | {'gap L [%]':>10s}")
    for k, rung in enumerate(r["cauer_rungs"]):
        if k >= len(cauer_ref_R):
            break
        gR = (rung["R_2k"] - cauer_ref_R[k]) / cauer_ref_R[k] * 100
        gL = (rung["L_2k_plus_1"] - cauer_ref_L[k]) / cauer_ref_L[k] * 100
        print(f"  {k:>3d} | {cauer_ref_R[k]:>14.4e} | {rung['R_2k']:>14.4e} | "
              f"{gR:>+10.3f} | {cauer_ref_L[k]:>14.4e} | "
              f"{rung['L_2k_plus_1']:>14.4e} | {gL:>+10.3f}")

    out_path = here / "sphere_vim_axisym_RL_extract_1080.json"
    out_path.write_text(json.dumps(r, indent=2), encoding="utf-8")
    print(f"\nSaved: {out_path.name}")


if __name__ == "__main__":
    main()
