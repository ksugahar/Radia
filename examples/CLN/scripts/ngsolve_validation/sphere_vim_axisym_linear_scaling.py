"""sphere_vim_axisym_linear_scaling.py — Push linear-basis sphere VIM
to the GPU memory limit, study DoF-vs-accuracy convergence.

Linear (bilinear Q1) VIM on (rho, theta) sphere mesh.  Each level
samples Foster spectrum + Cauer τ_pair[k] + top 8 Foster modes.

Output: sphere_vim_axisym_linear_scaling_<NRxNT>.json per case
        sphere_vim_axisym_linear_scaling_summary.json (aggregate)

Memory cost estimate (FP64 on GPU):
  cells  | DoF  | K mem  | eigh CPU time
   2430  | ~2300|  42 MB |  ~5 sec
   4320  | ~4100| 130 MB |  ~30 sec
   9720  | ~9200| 660 MB |  ~5 min
  19440  |~18000|  2.6 GB|  ~40 min   (still on CPU eigh)
  30000  |~28000|  6.3 GB|  ~3 hr     (consider GPU eigh)
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
import mpmath as mp

from sphere_vim_axisym import R_SPHERE, SIGMA, MU0, B0, hankel_pade_cauer
from sphere_vim_axisym_linear import assemble_K_M_b_linear_sphere

mp.mp.dps = 60


def extract_one(N_rho, N_theta, n_quad=3, top_n=10, n_stages=12):
    """Assemble + eigh + Foster + Cauer for a single mesh."""
    print(f"\n--- {N_rho}x{N_theta} = {N_rho*N_theta} cells, n_quad={n_quad} ---")

    t0 = time.time()
    K, M, b, interior_idx, n_vert = assemble_K_M_b_linear_sphere(
        N_rho, N_theta, n_quad=n_quad, use_gpu=True, verbose=False)
    t_K = time.time() - t0
    n_dofs = K.shape[0]
    print(f"  assembled: {t_K:.1f}s, K shape: {K.shape}, "
          f"K mem: {K.nbytes / 1e9:.2f} GB")

    t0 = time.time()
    from scipy.linalg import eigh
    eigvals, eigvecs = eigh(K, M)
    t_eigh = time.time() - t0
    print(f"  eigh: {t_eigh:.1f}s, λ ∈ [{eigvals.min():.3e}, "
          f"{eigvals.max():.3e}]")

    # M-normalize
    for k in range(eigvecs.shape[1]):
        v = eigvecs[:, k]
        nrm = math.sqrt(float(v @ M @ v))
        if nrm > 1e-30:
            eigvecs[:, k] = v / nrm

    tau_arr = SIGMA * eigvals
    g_arr = eigvecs.T @ b
    g2 = g_arr * g_arr

    # Top Foster modes
    pos = tau_arr > 1e-30
    tau_pos = tau_arr[pos]
    g2_pos = g_arr[pos] ** 2
    order = np.argsort(-g2_pos)
    sum_g2 = float(g2_pos.sum())
    top_modes = []
    for r in range(min(top_n, len(order))):
        idx = order[r]
        top_modes.append({
            "rank": r,
            "tau_us": float(tau_pos[idx]) * 1e6,
            "g2": float(g2_pos[idx]),
            "g2_normalized": float(g2_pos[idx] / sum_g2),
        })

    # Cauer extraction (Kameari convention τ^n)
    n_use = min(150, int(np.sum(g2_pos > 1e-32)))
    idx_use = order[:n_use]
    tau_mp = [mp.mpf(float(tau_pos[i])) for i in idx_use]
    g2_mp = [mp.mpf(float(g2_pos[i])) for i in idx_use]
    n_moments = min(2 * n_stages + 4, 2 * n_use - 2)
    alphas = []
    for k in range(n_moments):
        a_val = mp.mpf(0)
        for gv, tv in zip(g2_mp, tau_mp):
            a_val += gv * tv ** k   # KAMEARI convention (matches Stoll wls)
        alphas.append((mp.mpf(-1)) ** k * a_val)
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
        "N_rho": N_rho, "N_theta": N_theta,
        "n_cells": N_rho * N_theta, "n_dofs": n_dofs,
        "K_mem_GB": float(K.nbytes / 1e9),
        "t_assemble_s": t_K, "t_eigh_s": t_eigh,
        "leading_tau_us": float(tau_pos[order[0]]) * 1e6,
        "top_modes": top_modes,
        "cauer_rungs": rungs,
    }


def main():
    here = Path(__file__).parent
    ana = json.load(open(here / "sphere_analytical_cln_highstage.json"))
    stoll_tau = [float(x) for x in ana["stoll_tau_us"]]
    stoll_a = [float(x) for x in ana["stoll_a_n"]]
    stoll_cauer = [float(x) for x in ana["cauer_tau_us"]]

    print("=" * 78)
    print("Sphere VIM linear-basis DoF scaling on GPU (Kameari convention)")
    print(f"  Stoll Foster τ_1 = {stoll_tau[0]:.4f} μs")
    print(f"  Stoll Cauer τ_pair[0] = {stoll_cauer[0]:.4f} μs")
    print("=" * 78)

    # Mesh refinement schedule: roughly double cells each step
    # Stay aware of CPU eigh cost: O(N_dofs^3) — practical limit ~ 20000 DoF
    # CPU intermediate G_full = N_q² × 8 bytes, N_q = N_cells × n_quad².
    # On 16 GB CPU: limit ≈ (G < 8 GB) ⇒ N_q < 31600.
    # n_quad=2: N_cells_max ≈ 31600/4 ≈ 7900
    # n_quad=3: N_cells_max ≈ 31600/9 ≈ 3500 (so 2430 OK, 4320 OOM)
    levels = [
        (45, 54, 3),    # 2430 cells, n_quad=3 (existing baseline)
        (60, 72, 2),    # 4320 cells, n_quad=2 (N_q=17280, G=2.4 GB)
        (75, 90, 2),    # 6750 cells, n_quad=2 (N_q=27000, G=5.9 GB)
        (90, 108, 2),   # 9720 cells, n_quad=2 (N_q=38880, G=12.1 GB — borderline)
    ]

    results = []
    for N_rho, N_theta, n_quad in levels:
        try:
            r = extract_one(N_rho, N_theta, n_quad=n_quad, top_n=10, n_stages=12)
            results.append(r)
            out_path = (here /
                f"sphere_vim_axisym_linear_scaling_"
                f"{N_rho}x{N_theta}.json")
            out_path.write_text(json.dumps(r, indent=2), encoding="utf-8")
            print(f"  Saved {out_path.name}")
        except MemoryError as me:
            print(f"  OUT-OF-MEMORY at {N_rho}x{N_theta}: {me}")
            break
        except Exception as e:
            print(f"  FAIL {N_rho}x{N_theta}: {type(e).__name__}: {e}")

    # Comparison table
    print()
    print("=" * 78)
    print(f"Leading Foster τ_1 vs Stoll = {stoll_tau[0]:.4f} μs")
    print("=" * 78)
    print(f"  {'cells':>7s} | {'DoF':>7s} | {'t_K [s]':>9s} | "
          f"{'t_eigh [s]':>10s} | {'τ_1':>12s} | {'gap [%]':>10s}")
    for r in results:
        g = (r["leading_tau_us"] - stoll_tau[0]) / stoll_tau[0] * 100
        print(f"  {r['n_cells']:>7d} | {r['n_dofs']:>7d} | "
              f"{r['t_assemble_s']:>9.1f} | {r['t_eigh_s']:>10.1f} | "
              f"{r['leading_tau_us']:>12.6f} | {g:>+10.4f}")

    print()
    print("=" * 78)
    print("Cauer τ_pair[k] [μs] mesh convergence (linear basis, Kameari)")
    print("=" * 78)
    hdr = [f"{'k':>3s}", f"{'Stoll':>10s}"] + [
        f"{str(r['n_cells']) + 'c':>10s}" for r in results]
    print("  " + " | ".join(hdr))
    print("  " + "-+-".join("-" * 10 for _ in hdr))
    for k in range(8):
        if k >= len(stoll_cauer):
            break
        row = [f"{k:>3d}", f"{stoll_cauer[k]:>10.4f}"]
        for r in results:
            if k < len(r["cauer_rungs"]):
                row.append(f"{r['cauer_rungs'][k]['tau_pair_us']:>10.4f}")
            else:
                row.append(f"{'--':>10s}")
        print("  " + " | ".join(row))

    print()
    print("Gap to Stoll [%]:")
    print("  " + " | ".join(hdr))
    print("  " + "-+-".join("-" * 10 for _ in hdr))
    for k in range(8):
        if k >= len(stoll_cauer):
            break
        row = [f"{k:>3d}", f"{0.0:>+10.4f}"]
        for r in results:
            if k < len(r["cauer_rungs"]):
                g = (r["cauer_rungs"][k]["tau_pair_us"] - stoll_cauer[k]) / stoll_cauer[k] * 100
                row.append(f"{g:>+10.4f}")
            else:
                row.append(f"{'--':>10s}")
        print("  " + " | ".join(row))

    out = {
        "stoll_tau_us": stoll_tau[:20],
        "stoll_a_n": stoll_a[:20],
        "stoll_cauer_tau_us": stoll_cauer,
        "results": results,
    }
    (here / "sphere_vim_axisym_linear_scaling_summary.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8")
    print("\nSaved sphere_vim_axisym_linear_scaling_summary.json")


if __name__ == "__main__":
    main()
