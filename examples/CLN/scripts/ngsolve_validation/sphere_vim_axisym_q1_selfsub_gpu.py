"""sphere_vim_axisym_q1_selfsub_gpu.py — Q1 linear sphere VIM with self-cell
subdivision correction. GPU-native via cupy + AGM ellipe.

Idea: standard Gauss quadrature on self-cell (q1, q2 same cell) under-
resolves the log singularity of the axisymmetric kernel at the diagonal.
We correct the self-cell contribution by re-integrating with sub-cell
refinement.

Algorithm:
  1) K_full = (mu_0/2) B^T Gw B   — GPU, includes naive self-cell
  2) For each self-cell A:
       K_naive[A]    = standard n_quad x n_quad x n_quad x n_quad Gauss
       K_subdiv[A]   = same outer, inner uses n_sub x n_sub subdivision
       delta[A]      = K_subdiv[A] - K_naive[A]
       K_full[v_A, w_A] += delta[A]
  3) Eigendecomp K v = lambda M v.

For Q1 (bilinear), each cell has 4 corner DoFs. delta is 4x4 per cell.
"""
from __future__ import annotations

import json
import math
import sys
import time
import gc
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
import mpmath as mp
import scipy.special as sp
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
            a = a_new; b = b_new;
            double c2_new = diff * diff;
            sum_c2 += scale * c2_new;
            scale *= 2.0;
            if (diff < 1e-30) break;
        }
        double K = 1.5707963267948966 / a;
        E_out = K * (1.0 - sum_c2);
        """,
        name="ellipe_agm_q1ss",
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


def axisym_kernel_cpu(r, z, r_p, z_p):
    """CPU axisym kernel via scipy.special.ellipk/ellipe."""
    safe = (r > 1e-15) & (r_p > 1e-15)
    if not np.all(safe):
        return np.zeros_like(r)
    d2 = (r + r_p) ** 2 + (z - z_p) ** 2
    m = 4.0 * r * r_p / d2
    m = np.minimum(m, 1.0 - 1e-12)
    k = np.sqrt(m)
    K = sp.ellipk(m)
    E = sp.ellipe(m)
    coef = np.sqrt(r_p / r) / math.pi
    inv_k = 1.0 / np.maximum(k, 1e-30)
    return coef * ((2.0 * inv_k - k) * K - 2.0 * inv_k * E)


def self_cell_K_q1(rho_a, rho_b, theta_a, theta_b,
                     n_quad_outer, n_quad_inner, n_sub_inner):
    """Compute 4x4 K^A matrix for one Q1 self-cell using:
       - Outer: n_quad_outer x n_quad_outer Gauss on cell.
       - Inner: n_sub_inner x n_sub_inner subdivided cell, each
                 sub-cell with n_quad_inner x n_quad_inner Gauss.
    Includes the mu_0/2 prefactor.
    Returns 4x4 K matrix (node ordering: v00, v10, v11, v01)."""
    drho_cell = rho_b - rho_a
    dtheta_cell = theta_b - theta_a

    # Outer Gauss on [0, 1]
    pts_o, wts_o = np.polynomial.legendre.leggauss(n_quad_outer)
    qp_o = 0.5 * (pts_o + 1.0)
    qw_o = 0.5 * wts_o

    # Inner Gauss on [0, 1]
    pts_i, wts_i = np.polynomial.legendre.leggauss(n_quad_inner)
    qp_i = 0.5 * (pts_i + 1.0)
    qw_i = 0.5 * wts_i

    K_mat = np.zeros((4, 4))
    n_iq = n_sub_inner * n_quad_inner   # total inner samples per dimension
    sub_drho = drho_cell / n_sub_inner
    sub_dtheta = dtheta_cell / n_sub_inner

    # Build all inner samples up front
    xi_in_arr = np.empty(n_iq)
    eta_in_arr = np.empty(n_iq)
    w_in_xi = np.empty(n_iq)
    w_in_eta = np.empty(n_iq)
    for s_xi in range(n_sub_inner):
        for q_xi in range(n_quad_inner):
            idx = s_xi * n_quad_inner + q_xi
            xi_in_arr[idx] = (s_xi + qp_i[q_xi]) / n_sub_inner   # in [0, 1]
            w_in_xi[idx] = qw_i[q_xi] / n_sub_inner               # mass: 1/n_sub
    for s_eta in range(n_sub_inner):
        for q_eta in range(n_quad_inner):
            idx = s_eta * n_quad_inner + q_eta
            eta_in_arr[idx] = (s_eta + qp_i[q_eta]) / n_sub_inner
            w_in_eta[idx] = qw_i[q_eta] / n_sub_inner

    # Inner physical coords
    rho_in = rho_a + drho_cell * xi_in_arr      # (n_iq,)
    theta_in = theta_a + dtheta_cell * eta_in_arr
    r_in = rho_in[:, None] * np.sin(theta_in[None, :])  # outer product over (xi, eta)
    z_in = rho_in[:, None] * np.cos(theta_in[None, :])
    jac_in = (rho_in[:, None] * drho_cell * dtheta_cell
              * w_in_xi[:, None] * w_in_eta[None, :])  # (n_iq, n_iq)

    # Inner Q1 basis
    N_in = np.stack([
        (1 - xi_in_arr[:, None]) * (1 - eta_in_arr[None, :]),  # (xi_i, eta_i)
        xi_in_arr[:, None] * (1 - eta_in_arr[None, :]),
        xi_in_arr[:, None] * eta_in_arr[None, :],
        (1 - xi_in_arr[:, None]) * eta_in_arr[None, :],
    ], axis=0)   # (4, n_iq, n_iq)

    for o_xi_i in range(n_quad_outer):
        xi_o = qp_o[o_xi_i]
        rho_o = rho_a + drho_cell * xi_o
        for o_eta_i in range(n_quad_outer):
            eta_o = qp_o[o_eta_i]
            theta_o = theta_a + dtheta_cell * eta_o
            r_o = rho_o * math.sin(theta_o)
            z_o = rho_o * math.cos(theta_o)
            w_o = qw_o[o_xi_i] * qw_o[o_eta_i]
            jac_o = rho_o * drho_cell * dtheta_cell
            wt_o = w_o * jac_o * 2.0 * math.pi * r_o

            # Q1 basis at outer
            N_o = np.array([
                (1 - xi_o) * (1 - eta_o),
                xi_o * (1 - eta_o),
                xi_o * eta_o,
                (1 - xi_o) * eta_o,
            ])

            # G(P_outer, P_inner) over all inner samples
            G_arr = axisym_kernel_cpu(r_o, z_o, r_in.flatten(),
                                      z_in.flatten()).reshape(r_in.shape)

            # K^A[v, w] += wt_o * sum_inner [ N_o[v] * N_in[w, .] * G * jac_in ]
            # = wt_o * N_o[v] * sum_inner [ N_in[w, .] * G * jac_in ]
            weighted_inner = G_arr * jac_in       # (n_iq, n_iq)
            for w_i in range(4):
                inner_proj = (N_in[w_i] * weighted_inner).sum()
                K_mat[:, w_i] += wt_o * N_o[:] * inner_proj

    return (MU0 / 2.0) * K_mat


def assemble_K_M_b_with_selfsub(N_rho, N_theta, n_quad=3,
                                  n_sub=2, n_quad_inner_sub=None,
                                  verbose=True):
    """Q1 linear sphere VIM + self-cell subdivision correction.
    n_sub is the sub-cell count per dimension for self-cell inner."""
    if n_quad_inner_sub is None:
        n_quad_inner_sub = n_quad

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
        print(f"  Q1+self-sub: N={N_rho}x{N_theta} cells={N_rho*N_theta} "
              f"DoF={n_interior}, n_sub={n_sub}", flush=True)

    # === Step 1: standard GPU K, M, b (naive self-cell) ===
    nodes_m1, w_m1 = np.polynomial.legendre.leggauss(n_quad)
    qp = 0.5 * (nodes_m1 + 1.0)
    qw = 0.5 * w_m1

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

    xi = qp
    eta = qp
    n_qq = n_quad * n_quad

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

    B_gpu = cp.asarray(B_dense_cpu)
    r_flat = cp.asarray(r_flat_cpu)
    z_flat = cp.asarray(z_flat_cpu)
    jac_flat = cp.asarray(jac_flat_cpu)
    del rho_q_cpu, theta_q_cpu, r_q_full, z_q_full, jac_q_full
    gc.collect()
    if verbose:
        print(f"  B + transfer: {time.time()-t0:.2f}s, "
              f"B_gpu {B_gpu.nbytes/1e9:.2f} GB", flush=True)

    t0 = time.time()
    R_obs = r_flat[:, None]; Z_obs = z_flat[:, None]
    R_src = r_flat[None, :]; Z_src = z_flat[None, :]
    G = axisym_kernel_gpu(R_obs, Z_obs, R_src, Z_src,
                          cp, cp_ellipk, ellipe_agm)
    if verbose:
        print(f"  G_full ({n_q_total}^2): {time.time()-t0:.2f}s, "
              f"{G.nbytes/1e9:.2f} GB", flush=True)

    t0 = time.time()
    Gw = G * (2.0 * math.pi * r_flat * jac_flat)[:, None] * jac_flat[None, :]
    del G
    cp.get_default_memory_pool().free_all_blocks()
    Gw_B = Gw @ B_gpu
    del Gw
    cp.get_default_memory_pool().free_all_blocks()
    K_full_gpu = (MU0 / 2.0) * (B_gpu.T @ Gw_B)
    del Gw_B
    K_full_gpu = 0.5 * (K_full_gpu + K_full_gpu.T)
    if verbose:
        print(f"  K = B^T Gw B: {time.time()-t0:.2f}s", flush=True)

    t0 = time.time()
    weights_obs = 2.0 * math.pi * r_flat * jac_flat
    M_full_gpu = B_gpu.T @ (weights_obs[:, None] * B_gpu)
    M_full_gpu = 0.5 * (M_full_gpu + M_full_gpu.T)
    A_q = (B0 / 2.0) * r_flat
    weights_src = 2.0 * math.pi * r_flat * jac_flat
    b_full_gpu = B_gpu.T @ (A_q * weights_src)
    if verbose:
        print(f"  M, b: {time.time()-t0:.2f}s", flush=True)

    K_full = cp.asnumpy(K_full_gpu)
    M_full = cp.asnumpy(M_full_gpu)
    b_full = cp.asnumpy(b_full_gpu)
    del K_full_gpu, M_full_gpu, b_full_gpu, B_gpu
    cp.get_default_memory_pool().free_all_blocks()

    # === Step 2: self-cell correction (CPU) ===
    t0 = time.time()
    for i_cell in range(N_rho):
        rho_a_c = i_cell * drho
        rho_b_c = (i_cell + 1) * drho
        for j_cell in range(N_theta):
            theta_a_c = j_cell * dtheta
            theta_b_c = (j_cell + 1) * dtheta

            K_naive = self_cell_K_q1(rho_a_c, rho_b_c, theta_a_c, theta_b_c,
                                     n_quad_outer=n_quad,
                                     n_quad_inner=n_quad,
                                     n_sub_inner=1)
            K_subdiv = self_cell_K_q1(rho_a_c, rho_b_c, theta_a_c, theta_b_c,
                                      n_quad_outer=n_quad,
                                      n_quad_inner=n_quad_inner_sub,
                                      n_sub_inner=n_sub)
            delta = K_subdiv - K_naive

            v00_c = i_cell * n_vert_theta + j_cell
            v10_c = (i_cell + 1) * n_vert_theta + j_cell
            v11_c = (i_cell + 1) * n_vert_theta + (j_cell + 1)
            v01_c = i_cell * n_vert_theta + (j_cell + 1)
            local_idx = [v00_c, v10_c, v11_c, v01_c]
            for a in range(4):
                for b in range(4):
                    K_full[local_idx[a], local_idx[b]] += delta[a, b]
    if verbose:
        print(f"  Self-cell subdiv correction: {time.time()-t0:.2f}s",
              flush=True)

    K_full = 0.5 * (K_full + K_full.T)

    K_int = K_full[interior_idx][:, interior_idx]
    M_int = M_full[interior_idx][:, interior_idx]
    b_int = b_full[interior_idx]
    return K_int, M_int, b_int, interior_idx, n_vert


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


def extract(N_rho, N_theta, n_quad=3, n_sub=2, n_stages=12, top_n=10):
    print(f"\n--- {N_rho}x{N_theta} n_quad={n_quad} n_sub={n_sub} ---", flush=True)
    t0 = time.time()
    K, M, b, _, _ = assemble_K_M_b_with_selfsub(
        N_rho, N_theta, n_quad=n_quad, n_sub=n_sub)
    t_K = time.time() - t0
    print(f"  total assemble: {t_K:.1f}s, K shape: {K.shape}", flush=True)

    t0 = time.time()
    eigvals, eigvecs = eigh(K, M)
    t_eigh = time.time() - t0
    print(f"  eigh: {t_eigh:.1f}s, λ ∈ [{eigvals.min():.3e}, {eigvals.max():.3e}]",
          flush=True)

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

    top_modes = []
    sum_g2 = float(g2_pos.sum())
    for r in range(min(top_n, len(order))):
        idx = order[r]
        top_modes.append({
            "rank": r, "tau_us": float(tau_pos[idx]) * 1e6,
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

    leading_tau = float(tau_pos[order[0]]) * 1e6
    print(f"  tau_1={leading_tau:.6f} us, "
          f"tau_pair[0]={rungs[0]['tau_pair_us']:.6f} us"
          if rungs else f"  tau_1={leading_tau:.6f} us",
          flush=True)

    return {
        "N_rho": N_rho, "N_theta": N_theta,
        "n_cells": N_rho * N_theta, "n_dofs": K.shape[0],
        "n_quad": n_quad, "n_sub": n_sub,
        "t_assemble_s": t_K, "t_eigh_s": t_eigh,
        "leading_tau_us": leading_tau,
        "top_modes": top_modes, "cauer_rungs": rungs,
    }


def main():
    here = Path(__file__).parent
    ana = json.load(open(here / "sphere_analytical_cln_highstage.json"))
    stoll_tau = [float(x) for x in ana["stoll_tau_us"]]
    stoll_cauer = [float(x) for x in ana["cauer_tau_us"]]

    print("=" * 78, flush=True)
    print("Sphere VIM Q1 + self-cell subdiv correction on GPU", flush=True)
    print(f"  Stoll Foster tau_1 = {stoll_tau[0]:.6f} us", flush=True)
    print(f"  Stoll Cauer tau_pair[0] = {stoll_cauer[0]:.6f} us", flush=True)
    print("=" * 78, flush=True)

    # Test: 2430 cells, n_sub in {1 (= baseline naive), 2, 4}
    levels = [
        (45, 54, 3, 1),   # Q1 baseline (n_sub=1 means no subdiv)
        (45, 54, 3, 2),   # n_sub=2 first try
        (45, 54, 3, 4),   # n_sub=4 stronger
    ]
    results = []
    for N_rho, N_theta, nq, ns in levels:
        try:
            r = extract(N_rho, N_theta, n_quad=nq, n_sub=ns,
                        n_stages=12, top_n=10)
            results.append(r)
            (here / f"sphere_vim_axisym_q1_selfsub_gpu_"
                   f"{N_rho}x{N_theta}_ns{ns}.json").write_text(
                json.dumps(r, indent=2), encoding="utf-8")
        except Exception as e:
            import traceback
            print(f"  FAIL: {type(e).__name__}: {e}", flush=True)
            traceback.print_exc()
        gc.collect()

    print(flush=True)
    print("=" * 78, flush=True)
    print(f"Q1 + self-sub: tau_1 vs Stoll {stoll_tau[0]:.6f} us, "
          f"tau_pair[0] vs Stoll {stoll_cauer[0]:.6f} us", flush=True)
    print("=" * 78, flush=True)
    print(f"  {'cells':>6s} {'n_sub':>5s} | {'t_K':>6s}{'t_eigh':>7s} | "
          f"{'tau_1':>11s}{'gap':>10s} | "
          f"{'tau_pair[0]':>11s}{'gap':>10s}", flush=True)
    for r in results:
        g1 = (r["leading_tau_us"] - stoll_tau[0]) / stoll_tau[0] * 100
        rung0 = r["cauer_rungs"][0]["tau_pair_us"] if r["cauer_rungs"] else float('nan')
        gp0 = (rung0 - stoll_cauer[0]) / stoll_cauer[0] * 100
        print(f"  {r['n_cells']:>6d} {r['n_sub']:>5d} | "
              f"{r['t_assemble_s']:>6.1f}{r['t_eigh_s']:>7.1f} | "
              f"{r['leading_tau_us']:>11.6f}{g1:>+10.4f} | "
              f"{rung0:>11.6f}{gp0:>+10.4f}", flush=True)

    out = {"stoll_tau_us": stoll_tau[:20], "stoll_cauer_tau_us": stoll_cauer,
           "results": results}
    (here / "sphere_vim_axisym_q1_selfsub_gpu_summary.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8")
    print("\nSaved sphere_vim_axisym_q1_selfsub_gpu_summary.json", flush=True)


if __name__ == "__main__":
    main()
