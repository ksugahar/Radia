"""sphere_vim_axisym_subdiv.py — Self-cell subdivision for VIM K_ii.

Standard tensor Gauss-Legendre on the inner cell badly under-resolves
the log singularity at (rho', theta') = (rho, theta).  This script
replaces the inner integration for self-cells with a uniformly
subdivided quadrature: each self-cell is split into n_sub x n_sub
sub-cells, and Gauss-n_quad applied to each sub-cell.  The total
sample density goes from n_quad^2 to (n_sub n_quad)^2 inside the cell,
clustering many points near the singular outer-point.

Off-diagonal entries (i != j) keep the original n_quad smooth Gauss.

We vary n_sub in {1, 2, 4, 8} and report:
  - Foster leading tau_1 (top |g|^2 mode)
  - Cauer tau_pair[0..5]
  - K_ii diagonal statistics
to test whether better self-cell handling closes the Stoll gap.
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
import scipy.special as sp
from scipy.linalg import eigh
import mpmath as mp

from sphere_vim_axisym import (
    axisym_kernel,
    R_SPHERE, SIGMA, MU0, B0,
    hankel_pade_cauer,
)

mp.mp.dps = 60


def assemble_K_M_b_subdiv(N_rho, N_theta, n_quad=2, n_sub_self=2):
    """Like sphere_vim_axisym.assemble_K_M_b but with subdivided
    inner quadrature for the self-cell.
    """
    n = N_rho * N_theta
    drho = R_SPHERE / N_rho
    dtheta = math.pi / N_theta

    cells = []
    for i in range(N_rho):
        for j in range(N_theta):
            rho_a, rho_b = i * drho, (i + 1) * drho
            th_a, th_b = j * dtheta, (j + 1) * dtheta
            cells.append((rho_a, rho_b, th_a, th_b))

    # Mass matrix (diagonal)
    M = np.zeros((n, n))
    for idx, (ra, rb, ta, tb) in enumerate(cells):
        vol = (2.0 * math.pi * (rb ** 3 - ra ** 3) / 3.0
               * (math.cos(ta) - math.cos(tb)))
        M[idx, idx] = vol

    # b vector
    b_vec = np.zeros(n)
    for idx, (ra, rb, ta, tb) in enumerate(cells):
        rho_int = (rb ** 4 - ra ** 4) / 4.0
        th_int = ((tb - ta) / 2.0
                  - (math.sin(2.0 * tb) - math.sin(2.0 * ta)) / 4.0)
        b_vec[idx] = math.pi * B0 * rho_int * th_int

    # Pre-build n_quad Gauss nodes/weights on [0, 1]
    g_nodes_m1, g_w_m1 = np.polynomial.legendre.leggauss(n_quad)
    g_nodes = 0.5 * (g_nodes_m1 + 1.0)
    g_w = 0.5 * g_w_m1

    def cell_gauss(cell):
        """Return (rho_q, th_q, w_rho, w_th) of n_quad x n_quad Gauss on cell."""
        ra, rb, ta, tb = cell
        rho_q = ra + (rb - ra) * g_nodes
        th_q = ta + (tb - ta) * g_nodes
        w_rho = (rb - ra) * g_w
        w_th = (tb - ta) * g_w
        return rho_q, th_q, w_rho, w_th

    def subdiv_gauss(cell, n_sub):
        """Return (rho_q, th_q, w_rho, w_th) of (n_sub * n_quad) ^ 2 points
        from a uniform n_sub x n_sub subdivision of cell."""
        ra, rb, ta, tb = cell
        dr_sub = (rb - ra) / n_sub
        dt_sub = (tb - ta) / n_sub
        rho_q = np.zeros(n_sub * n_quad)
        th_q = np.zeros(n_sub * n_quad)
        w_rho = np.zeros(n_sub * n_quad)
        w_th = np.zeros(n_sub * n_quad)
        for s in range(n_sub):
            r0 = ra + s * dr_sub
            t0 = ta + s * dt_sub
            for q, (n01, w01) in enumerate(zip(g_nodes, g_w)):
                rho_q[s * n_quad + q] = r0 + dr_sub * n01
                w_rho[s * n_quad + q] = dr_sub * w01
        for s in range(n_sub):
            t0 = ta + s * dt_sub
            for q, (n01, w01) in enumerate(zip(g_nodes, g_w)):
                th_q[s * n_quad + q] = t0 + dt_sub * n01
                w_th[s * n_quad + q] = dt_sub * w01
        return rho_q, th_q, w_rho, w_th

    K = np.zeros((n, n))
    for i, cell_i in enumerate(cells):
        ra_i, rb_i, ta_i, tb_i = cell_i
        rho_qi, th_qi, wrho_i, wth_i = cell_gauss(cell_i)
        r_qi = np.outer(rho_qi, np.sin(th_qi))
        z_qi = np.outer(rho_qi, np.cos(th_qi))
        # jacobian: dr dz = rho dρ dθ ⇒ J = rho_q (independent of θ)
        for j, cell_j in enumerate(cells):
            if i == j:
                rho_qj, th_qj, wrho_j, wth_j = subdiv_gauss(cell_j, n_sub_self)
            else:
                rho_qj, th_qj, wrho_j, wth_j = cell_gauss(cell_j)
            r_qj = np.outer(rho_qj, np.sin(th_qj))
            z_qj = np.outer(rho_qj, np.cos(th_qj))

            kij = 0.0
            for ar in range(rho_qi.shape[0]):
                rho_i = rho_qi[ar]
                for at in range(th_qi.shape[0]):
                    ri = r_qi[ar, at]
                    zi = z_qi[ar, at]
                    if ri < 1e-15:
                        continue
                    wi = rho_i * wrho_i[ar] * wth_i[at]
                    for br in range(rho_qj.shape[0]):
                        rho_j = rho_qj[br]
                        for bt in range(th_qj.shape[0]):
                            rj = r_qj[br, bt]
                            zj = z_qj[br, bt]
                            if rj < 1e-15:
                                continue
                            G = axisym_kernel(ri, zi, rj, zj)
                            wj = rho_j * wrho_j[br] * wth_j[bt]
                            kij += G * (2.0 * math.pi * ri) * wi * wj
            K[i, j] = (MU0 / 2.0) * kij

    K = 0.5 * (K + K.T)
    return K, M, b_vec


def run_one(N_rho, N_theta, n_quad, n_sub_self):
    print(f"\n--- n_quad={n_quad}, n_sub_self={n_sub_self} ---")
    t0 = time.time()
    K, M, b = assemble_K_M_b_subdiv(N_rho, N_theta,
                                      n_quad=n_quad,
                                      n_sub_self=n_sub_self)
    t_a = time.time() - t0
    print(f"  assembled in {t_a:.1f}s")
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
        rungs.append({"k": k, "tau_pair_us": float(Lv / Rv) * 1e6})

    return {
        "n_quad": n_quad, "n_sub_self": n_sub_self,
        "t_assemble_s": t_a,
        "K_diag_mean": float(diag_K.mean()),
        "leading_tau_us": leading_tau_us,
        "cauer_rungs": rungs,
    }


def main():
    here = Path(__file__).parent
    ana = json.load(open(here / "sphere_analytical_cln_highstage.json"))
    stoll_tau1 = float(ana["stoll_tau_us"][0])
    cauer_ref = [float(x) for x in ana["cauer_tau_us"]]

    N_rho, N_theta = 20, 24
    print("=" * 78)
    print(f"Sphere VIM self-cell subdivision sweep, mesh "
          f"{N_rho}x{N_theta}={N_rho*N_theta} cells, n_quad=2")
    print(f"  Stoll Foster τ_1 = {stoll_tau1:.4f} μs")
    print(f"  Stoll Cauer τ_pair[0] = {cauer_ref[0]:.4f} μs")
    print("=" * 78)

    results = []
    for ns in [1, 2, 4, 8]:
        try:
            r = run_one(N_rho, N_theta, n_quad=2, n_sub_self=ns)
            results.append(r)
        except Exception as e:
            print(f"  FAIL n_sub={ns}: {type(e).__name__}: {e}")

    print()
    print("=" * 78)
    print("Leading Foster τ_1 [μs] / K_ii mean / Cauer τ_pair[0]")
    print("=" * 78)
    hdr = f"  {'n_sub':>6s} | {'time_s':>8s} | {'K_ii mean':>14s} | {'τ_1':>10s} | {'gap':>8s} | {'τ_pair[0]':>12s} | {'gap':>8s}"
    print(hdr)
    for r in results:
        g1 = (r["leading_tau_us"] - stoll_tau1) / stoll_tau1 * 100
        rung0 = r["cauer_rungs"][0]["tau_pair_us"] if r["cauer_rungs"] else float('nan')
        g0 = (rung0 - cauer_ref[0]) / cauer_ref[0] * 100
        print(f"  {r['n_sub_self']:>6d} | {r['t_assemble_s']:>8.1f} | "
              f"{r['K_diag_mean']:>14.4e} | {r['leading_tau_us']:>10.4f} | "
              f"{g1:>+8.4f} | {rung0:>12.4f} | {g0:>+8.4f}")

    print()
    print("Cauer τ_pair[k] [μs]:")
    print(f"  {'k':>4s} | {'Stoll':>12s} | " +
          " | ".join(f"{'ns=' + str(r['n_sub_self']):>11s}" for r in results))
    for k in range(8):
        if k >= len(cauer_ref):
            break
        row = [f"{k:>4d}", f"{cauer_ref[k]:>12.5f}"]
        for r in results:
            if k < len(r["cauer_rungs"]):
                row.append(f"{r['cauer_rungs'][k]['tau_pair_us']:>11.4f}")
            else:
                row.append(f"{'--':>11s}")
        print("  " + " | ".join(row))

    print()
    print("Gap to Stoll [%]:")
    print(f"  {'k':>4s} | {'Stoll':>12s} | " +
          " | ".join(f"{'ns=' + str(r['n_sub_self']):>11s}" for r in results))
    for k in range(8):
        if k >= len(cauer_ref):
            break
        row = [f"{k:>4d}", f"{0.0:>+12.4f}"]
        for r in results:
            if k < len(r["cauer_rungs"]):
                g = (r["cauer_rungs"][k]["tau_pair_us"] - cauer_ref[k]) / cauer_ref[k] * 100
                row.append(f"{g:>+11.4f}")
            else:
                row.append(f"{'--':>11s}")
        print("  " + " | ".join(row))

    out = {
        "N_rho": N_rho, "N_theta": N_theta, "n_cells": N_rho * N_theta,
        "stoll_tau1_us": stoll_tau1, "stoll_cauer_tau_us": cauer_ref,
        "by_n_sub": results,
    }
    out_path = here / "sphere_vim_axisym_subdiv_results.json"
    out_path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"\nSaved: {out_path.name}")


if __name__ == "__main__":
    main()
