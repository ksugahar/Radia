"""sphere_vim_axisym_foster_modes.py — Extract top Foster (τ_n, g_n²)
modes from the subdiv ns=2 axisym VIM at one fixed mesh, for paper
table comparison vs Stoll analytical.

Outputs sphere_vim_axisym_foster_modes_<NRxNT>.json with top 8 modes.
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

from sphere_vim_axisym import R_SPHERE, SIGMA, MU0, B0
from sphere_vim_axisym_subdiv import assemble_K_M_b_subdiv


def extract(N_rho, N_theta, n_sub_self=2, top_n=8):
    print(f"\nN_rho={N_rho}, N_theta={N_theta} cells={N_rho*N_theta}")
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
    g_pos = g_arr[pos]
    g2_pos = g_pos * g_pos
    order = np.argsort(-g2_pos)
    sum_g2 = float(g2_pos.sum())
    modes = []
    for r in range(min(top_n, len(order))):
        idx = order[r]
        modes.append({
            "rank": r,
            "tau_us": float(tau_pos[idx]) * 1e6,
            "g2": float(g2_pos[idx]),
            "g2_normalized": float(g2_pos[idx] / sum_g2),
        })
    return {
        "N_rho": N_rho, "N_theta": N_theta,
        "n_cells": N_rho * N_theta,
        "n_sub_self": n_sub_self,
        "sum_g2": sum_g2,
        "top_modes": modes,
    }


def main():
    here = Path(__file__).parent
    ana = json.load(open(here / "sphere_analytical_cln_highstage.json"))
    stoll_tau = [float(x) for x in ana["stoll_tau_us"]]
    stoll_a = [float(x) for x in ana["stoll_a_n"]]

    out_list = []
    # 1080 cells is sufficient and quick (~1 min); 2880 would be slow.
    # We extract at 1080 only for the paper Foster mode table.
    for (Nr, Nt) in [(30, 36)]:
        out_list.append(extract(Nr, Nt, n_sub_self=2, top_n=8))

    print()
    print("=" * 78)
    print("Foster spectrum comparison: Stoll analytical vs VIM (subdiv ns=2)")
    print("=" * 78)
    n_show = 8
    for r in out_list:
        print(f"\n--- VIM {r['n_cells']} cells (subdiv ns=2) vs Stoll ---")
        print(f"  {'n':>3s} | {'Stoll τ_n (μs)':>16s} | {'VIM τ_n (μs)':>14s} | "
              f"{'gap [%]':>9s} | {'Stoll a_n':>11s} | {'VIM g²_norm':>12s} | "
              f"{'ratio':>8s}")
        for n in range(n_show):
            m = r["top_modes"][n] if n < len(r["top_modes"]) else None
            if m is None:
                continue
            st = stoll_tau[n]
            sa = stoll_a[n]
            gap_tau = (m["tau_us"] - st) / st * 100
            ratio = m["g2_normalized"] / sa
            print(f"  {n+1:>3d} | {st:>16.6f} | {m['tau_us']:>14.6f} | "
                  f"{gap_tau:>+9.4f} | {sa:>11.6f} | "
                  f"{m['g2_normalized']:>12.6f} | {ratio:>+8.4f}")

    out_path = here / "sphere_vim_axisym_foster_modes.json"
    out_path.write_text(json.dumps({
        "stoll_tau_us": stoll_tau[:20],
        "stoll_a_n": stoll_a[:20],
        "vim_results": out_list,
    }, indent=2), encoding="utf-8")
    print(f"\nSaved: {out_path.name}")


if __name__ == "__main__":
    main()
