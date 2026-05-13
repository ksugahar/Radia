"""sphere_vim_axisym_linear_gpu_v2.py — GPU-native linear VIM, v2.

v2 changes:
  - Dense B aggregator (matches local sphere_vim_axisym_linear.py exactly)
  - Explicit GPU memory pool cleanup between iterations
  - Custom ellipe via AGM on GPU (cupy lacks ellipe)
"""
from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
import mpmath as mp
from scipy.linalg import eigh

mp.mp.dps = 60

R_SPHERE = 10e-3
SIGMA = 5.8e7
MU0 = 4 * math.pi * 1e-7
B0 = 1.0


def _gpu_setup():
    import cupy as cp
    from cupyx.scipy.special import ellipk as cp_ellipk

    ellipe_agm = cp.ElementwiseKernel(
        in_params="float64 m",
        out_params="float64 E_out",
        operation=r"""
        double a = 1.0;
        double b = sqrt(1.0 - m);
        double c2 = m;
        double sum_c2 = 0.5 * c2;
        double scale = 1.0;
        for (int n = 0; n < 14; n++) {
            double a_new = 0.5 * (a + b);
            double b_new = sqrt(a * b);
            double diff = 0.5 * (a - b);
            a = a_new;
            b = b_new;
            scale *= 2.0;
            double c2_new = diff * diff;
            sum_c2 += scale * c2_new;
            if (diff < 1e-30) break;
        }
        double K = 1.5707963267948966 / a;
        E_out = K * (1.0 - sum_c2);
        """,
        name="ellipe_agm",
    )
    return cp, cp_ellipk, ellipe_agm


def axisym_kernel_gpu(r, z, r_p, z_p, cp, cp_ellipk, ellipe_agm):
    safe = (r > 1e-15) & (r_p > 1e-15)
    d2 = (r + r_p) ** 2 + (z - z_p) ** 2
    m = cp.where(safe, 4.0 * r * r_p / d2, cp.float64(0.0))
    m = cp.minimum(m, cp.float64(1.0 - 1e-12))
    k = cp.sqrt(m)
    K = cp_ellipk(m)
    E = ellipe_agm(m)
    coef = cp.sqrt(cp.where(safe, r_p / r, cp.float64(1.0))) / math.pi
    inv_k = cp.float64(1.0) / cp.maximum(k, cp.float64(1e-30))
    val = coef * ((2.0 * inv_k - k) * K - 2.0 * inv_k * E)
    return cp.where(safe, val, cp.float64(0.0))


def assemble_K_M_b(N_rho, N_theta, n_quad, verbose=True):
    cp, cp_ellipk, ellipe_agm = _gpu_setup()

    n_vert_rho = N_rho + 1
    n_vert_theta = N_theta + 1
    n_vert = n_vert_rho * n_vert_theta
    drho = R_SPHERE / N_rho
    dtheta = math.pi / N_theta

    interior = np.zeros((n_vert_rho, n_vert_theta), dtype=bool)
    interior[1:, 1:-1] = True
    interior_flat = interior.flatten()
    interior_idx = np.where(interior_flat)[0]
    n_interior = len(interior_idx)

    if verbose:
        print(f"  N={N_rho}x{N_theta} cells={N_rho*N_theta} DoF={n_interior}")

    nodes_m1, w_m1 = np.polynomial.legendre.leggauss(n_quad)
    qp = 0.5 * (nodes_m1 + 1.0)
    qw = 0.5 * w_m1

    # Build q-arrays on CPU first (small)
    rho_v = np.linspace(0, R_SPHERE, n_vert_rho)
    theta_v = np.linspace(0, math.pi, n_vert_theta)

    rho_q_cpu = (rho_v[:-1, None, None, None]
                 + drho * qp[None, None, :, None])
    rho_q_cpu = np.broadcast_to(rho_q_cpu,
                                 (N_rho, N_theta, n_quad, n_quad)).copy()
    theta_q_cpu = (theta_v[None, :-1, None, None]
                   + dtheta * qp[None, None, None, :])
    theta_q_cpu = np.broadcast_to(theta_q_cpu,
                                   (N_rho, N_theta, n_quad, n_quad)).copy()
    r_q_full = rho_q_cpu * np.sin(theta_q_cpu)
    z_q_full = rho_q_cpu * np.cos(theta_q_cpu)
    jac_q_full = (rho_q_cpu * drho * dtheta
                  * qw[None, None, :, None] * qw[None, None, None, :])

    n_q_total = N_rho * N_theta * n_quad * n_quad
    r_flat_cpu = r_q_full.reshape(n_q_total)
    z_flat_cpu = z_q_full.reshape(n_q_total)
    jac_flat_cpu = jac_q_full.reshape(n_q_total)

    # B_dense (CPU): N_00, N_10, N_11, N_01 at each q
    xi = qp
    eta = qp
    N_00 = np.outer(1 - xi, 1 - eta).flatten()
    N_10 = np.outer(xi,     1 - eta).flatten()
    N_11 = np.outer(xi,     eta).flatten()
    N_01 = np.outer(1 - xi, eta).flatten()
    n_qq = n_quad * n_quad

    if verbose:
        t0 = time.time()
    B_dense_cpu = np.zeros((n_q_total, n_vert), dtype=np.float64)
    irho_arr, jtheta_arr, irow_arr, icol_arr = np.meshgrid(
        np.arange(N_rho), np.arange(N_theta),
        np.arange(n_quad), np.arange(n_quad), indexing="ij")
    irho_arr = irho_arr.flatten()
    jtheta_arr = jtheta_arr.flatten()
    irow_arr = irow_arr.flatten()
    icol_arr = icol_arr.flatten()
    q_idx = np.arange(n_q_total)

    v00 = irho_arr * n_vert_theta + jtheta_arr
    v10 = (irho_arr + 1) * n_vert_theta + jtheta_arr
    v11 = (irho_arr + 1) * n_vert_theta + (jtheta_arr + 1)
    v01 = irho_arr * n_vert_theta + (jtheta_arr + 1)

    N_00_v = np.outer(1 - xi, 1 - eta)
    N_10_v = np.outer(xi,     1 - eta)
    N_11_v = np.outer(xi,     eta)
    N_01_v = np.outer(1 - xi, eta)

    B_dense_cpu[q_idx, v00] = N_00_v[irow_arr, icol_arr]
    B_dense_cpu[q_idx, v10] = N_10_v[irow_arr, icol_arr]
    B_dense_cpu[q_idx, v11] = N_11_v[irow_arr, icol_arr]
    B_dense_cpu[q_idx, v01] = N_01_v[irow_arr, icol_arr]
    if verbose:
        print(f"  B_dense ({B_dense_cpu.shape}, "
              f"{B_dense_cpu.nbytes/1e9:.2f} GB CPU): "
              f"{time.time()-t0:.2f}s")

    # Move to GPU
    if verbose:
        t0 = time.time()
    B_gpu = cp.asarray(B_dense_cpu)
    r_flat = cp.asarray(r_flat_cpu)
    z_flat = cp.asarray(z_flat_cpu)
    jac_flat = cp.asarray(jac_flat_cpu)
    del B_dense_cpu, rho_q_cpu, theta_q_cpu, r_q_full, z_q_full, jac_q_full
    if verbose:
        print(f"  GPU transfer: {time.time()-t0:.2f}s, "
              f"B_gpu {B_gpu.nbytes/1e9:.2f} GB")

    # G_full on GPU
    if verbose:
        t0 = time.time()
    R_obs = r_flat[:, None]
    Z_obs = z_flat[:, None]
    R_src = r_flat[None, :]
    Z_src = z_flat[None, :]
    G = axisym_kernel_gpu(R_obs, Z_obs, R_src, Z_src,
                          cp, cp_ellipk, ellipe_agm)
    if verbose:
        print(f"  G_full ({n_q_total}^2): {time.time()-t0:.2f}s, "
              f"{G.nbytes/1e9:.2f} GB")

    # Gw = G * (2 pi r * jac) ⊗ jac
    if verbose:
        t0 = time.time()
    Gw = G * (2.0 * math.pi * r_flat * jac_flat)[:, None] * jac_flat[None, :]
    del G
    cp.get_default_memory_pool().free_all_blocks()

    # K_full = (mu_0/2) B^T Gw B
    Gw_B = Gw @ B_gpu                          # (n_q, n_vert)
    del Gw
    K_full = (MU0 / 2.0) * (B_gpu.T @ Gw_B)
    del Gw_B
    K_full = 0.5 * (K_full + K_full.T)
    if verbose:
        print(f"  K = B^T Gw B: {time.time()-t0:.2f}s, "
              f"{K_full.nbytes/1e9:.2f} GB")

    # M_full = B^T diag(2 pi r * jac) B
    if verbose:
        t0 = time.time()
    weights_obs = 2.0 * math.pi * r_flat * jac_flat
    M_full = B_gpu.T @ (weights_obs[:, None] * B_gpu)
    M_full = 0.5 * (M_full + M_full.T)
    if verbose:
        print(f"  M = B^T D B: {time.time()-t0:.2f}s")

    # b vector
    A_q = (B0 / 2.0) * r_flat
    weights_src = 2.0 * math.pi * r_flat * jac_flat
    b_full = B_gpu.T @ (A_q * weights_src)

    # Restrict + transfer to CPU
    idx_gpu = cp.asarray(interior_idx)
    K_int = K_full[idx_gpu][:, idx_gpu]
    M_int = M_full[idx_gpu][:, idx_gpu]
    b_int = b_full[idx_gpu]

    K_cpu = cp.asnumpy(K_int)
    M_cpu = cp.asnumpy(M_int)
    b_cpu = cp.asnumpy(b_int)

    del K_full, M_full, b_full, K_int, M_int, b_int, B_gpu
    del r_flat, z_flat, jac_flat, weights_obs, weights_src, A_q
    cp.get_default_memory_pool().free_all_blocks()

    return K_cpu, M_cpu, b_cpu, interior_idx, n_vert


def hankel_pade_cauer(alphas, max_stages):
    c = [mp.mpf(a) for a in alphas]
    p = []
    for _ in range(max_stages):
        if len(c) < 2 or abs(c[0]) < mp.mpf(10) ** (-mp.mp.dps + 10):
            break
        n = len(c)
        e = [mp.mpf(0)] * n
        e[0] = 1 / c[0]
        for k in range(1, n):
            e[k] = -sum(c[j + 1] * e[k - j - 1] for j in range(k)
                        if k - j - 1 >= 0) / c[0]
        p.append(e[0])
        c = e[1:]
    return p


def extract(N_rho, N_theta, n_quad, top_n=10, n_stages=12):
    print(f"\n--- {N_rho}x{N_theta} = {N_rho*N_theta} cells, n_quad={n_quad} ---")
    t0 = time.time()
    K, M, b, _, _ = assemble_K_M_b(N_rho, N_theta, n_quad)
    t_K = time.time() - t0
    print(f"  total assemble: {t_K:.1f}s, K shape: {K.shape}")

    t0 = time.time()
    eigvals, eigvecs = eigh(K, M)
    t_eigh = time.time() - t0
    print(f"  eigh CPU: {t_eigh:.1f}s, "
          f"λ ∈ [{eigvals.min():.3e}, {eigvals.max():.3e}]")

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

    n_use = min(150, int(np.sum(g2_pos > 1e-32)))
    idx_use = order[:n_use]
    tau_mp = [mp.mpf(float(tau_pos[i])) for i in idx_use]
    g2_mp = [mp.mpf(float(g2_pos[i])) for i in idx_use]
    n_moments = min(2 * n_stages + 4, 2 * n_use - 2)
    alphas = []
    for k in range(n_moments):
        a_val = mp.mpf(0)
        for gv, tv in zip(g2_mp, tau_mp):
            a_val += gv * tv ** k
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
            "k": k, "R_2k": float(Rv),
            "L_2k_plus_1": float(Lv),
            "tau_pair_us": float(Lv / Rv) * 1e6,
        })

    return {
        "N_rho": N_rho, "N_theta": N_theta,
        "n_cells": N_rho * N_theta, "n_dofs": K.shape[0],
        "t_assemble_s": t_K, "t_eigh_s": t_eigh,
        "leading_tau_us": float(tau_pos[order[0]]) * 1e6,
        "top_modes": top_modes, "cauer_rungs": rungs,
    }


def main():
    here = Path(__file__).parent
    ana = json.load(open(here / "sphere_analytical_cln_highstage.json"))
    stoll_tau = [float(x) for x in ana["stoll_tau_us"]]
    stoll_cauer = [float(x) for x in ana["cauer_tau_us"]]

    print("=" * 78)
    print("Sphere VIM linear-basis GPU scaling v2 (dense B, free GPU mem)")
    print(f"  Stoll Foster τ_1 = {stoll_tau[0]:.6f} μs")
    print(f"  Stoll Cauer τ_pair[0] = {stoll_cauer[0]:.6f} μs")
    print("=" * 78)

    levels = [
        (45, 54, 3),
        (60, 72, 3),
        (75, 90, 3),
        (90, 108, 3),     # 9720 cells, G = 61 GB — too big for A100 40
        (78, 96, 3),
    ]
    # Actually keep budget under 40 GB: cells ≤ 7800 for n_quad=3
    levels = [
        (45, 54, 3),    # 2430 cells, G=3.8GB
        (60, 72, 3),    # 4320 cells, G=12 GB
        (75, 90, 3),    # 6750 cells, G=30 GB
    ]

    results = []
    for N_rho, N_theta, n_quad in levels:
        try:
            r = extract(N_rho, N_theta, n_quad, top_n=10, n_stages=12)
            results.append(r)
            (here / f"sphere_vim_axisym_linear_gpu_v2_{N_rho}x{N_theta}.json").write_text(
                json.dumps(r, indent=2), encoding="utf-8")
        except Exception as e:
            print(f"  FAIL {N_rho}x{N_theta}: {type(e).__name__}: {e}")

    print()
    print("=" * 78)
    print(f"Leading Foster τ_1 vs Stoll = {stoll_tau[0]:.6f} μs")
    print("=" * 78)
    for r in results:
        g = (r["leading_tau_us"] - stoll_tau[0]) / stoll_tau[0] * 100
        print(f"  {r['n_cells']:>7d} cells, {r['n_dofs']:>6d} DoF | "
              f"t_K={r['t_assemble_s']:>5.1f}s eigh={r['t_eigh_s']:>5.1f}s | "
              f"τ_1={r['leading_tau_us']:>10.6f} gap={g:>+8.4f}%")

    print()
    print("Cauer τ_pair[k] gap to Stoll [%]:")
    hdr = [f"{'k':>3s}"] + [f"{str(r['n_cells']) + 'c':>10s}" for r in results]
    print("  " + " | ".join(hdr))
    for k in range(min(8, len(stoll_cauer))):
        row = [f"{k:>3d}"]
        for r in results:
            if k < len(r["cauer_rungs"]):
                gp = (r["cauer_rungs"][k]["tau_pair_us"] - stoll_cauer[k]) / stoll_cauer[k] * 100
                row.append(f"{gp:>+10.4f}")
            else:
                row.append(f"{'--':>10s}")
        print("  " + " | ".join(row))

    out = {"stoll_tau_us": stoll_tau[:20], "stoll_cauer_tau_us": stoll_cauer,
           "results": results}
    (here / "sphere_vim_axisym_linear_gpu_v2_summary.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8")
    print("\nSaved sphere_vim_axisym_linear_gpu_v2_summary.json")


if __name__ == "__main__":
    main()
