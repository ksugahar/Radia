"""sphere_vim_axisym_q2_gpu.py — Biquadratic Q2 basis sphere VIM on GPU.

Q2 = piecewise biquadratic Lagrange basis on (rho, theta) sphere mesh.
Each cell has 9 nodes: 4 corners + 4 mid-edges + 1 center.
1D Q2 reference shape functions on [0, 1] (nodes at 0, 1/2, 1):
   L_0(t) = (1-t)(1-2t) = 1 - 3t + 2t^2
   L_1(t) = 4 t (1-t)
   L_2(t) = t (2t-1)
2D: N_ij(xi, eta) = L_i(xi) L_j(eta)  for i,j in {0,1,2}.

Strang-Fix theorem: eigenvalue error ~ O(h^{2p}) where p=2 for Q2.
Q1 (p=1) gave gap +0.18% at 2430 cells. Q2 at same h should give
gap ~ (h^2) smaller = ~ 1/(50)^2 = 4e-4 of the Q1 gap, i.e., ~ 7e-5 % gap.
Possibly machine-precision at 2430 cells.

GPU port (cupy + cupyx ellipk + AGM ellipe for cupyx-missing ellipe).
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
        name="ellipe_agm_q2",
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


# 1D Q2 shape functions on [0, 1] (nodes at 0, 1/2, 1)
def _L0(t): return (1.0 - t) * (1.0 - 2.0 * t)
def _L1(t): return 4.0 * t * (1.0 - t)
def _L2(t): return t * (2.0 * t - 1.0)


def assemble_K_M_b_q2(N_rho, N_theta, n_quad=3, verbose=True):
    """Q2 biquadratic basis on (rho, theta) sphere mesh.
    9 nodes per cell. Global vertex grid: (2*N_rho+1) x (2*N_theta+1)."""
    cp, cp_ellipk, ellipe_agm = _gpu_setup()

    # Global node grid (2*N_rho+1) x (2*N_theta+1)
    N_node_rho = 2 * N_rho + 1
    N_node_theta = 2 * N_theta + 1
    n_node = N_node_rho * N_node_theta
    drho_node = R_SPHERE / (2 * N_rho)         # spacing between adjacent nodes
    dtheta_node = math.pi / (2 * N_theta)
    drho_cell = 2 * drho_node                   # cell width
    dtheta_cell = 2 * dtheta_node

    # Boundary: J_phi = 0 where r = rho sin(theta) = 0
    # r=0 at rho=0 (i_node=0) OR theta=0 (j_node=0) OR theta=pi (j_node=N_node_theta-1)
    interior = np.ones((N_node_rho, N_node_theta), dtype=bool)
    interior[0, :] = False                  # rho = 0
    interior[:, 0] = False                  # theta = 0
    interior[:, N_node_theta - 1] = False   # theta = pi
    interior_flat = interior.flatten()
    interior_idx = np.where(interior_flat)[0]
    n_interior = len(interior_idx)

    if verbose:
        print(f"  Q2: N_cells={N_rho}x{N_theta}={N_rho*N_theta}, "
              f"global nodes {N_node_rho}x{N_node_theta}={n_node}, "
              f"interior DoF={n_interior}", flush=True)

    # Quadrature on cell reference [0, 1]^2
    nodes_m1, w_m1 = np.polynomial.legendre.leggauss(n_quad)
    qp = 0.5 * (nodes_m1 + 1.0)
    qw = 0.5 * w_m1

    # Build (N_rho, N_theta, n_quad, n_quad) quadrature arrays
    rho_q_cpu = np.empty((N_rho, N_theta, n_quad, n_quad), dtype=np.float64)
    theta_q_cpu = np.empty_like(rho_q_cpu)
    jac_q_cpu = np.empty_like(rho_q_cpu)

    for i_cell in range(N_rho):
        rho_a = i_cell * drho_cell
        for j_cell in range(N_theta):
            th_a = j_cell * dtheta_cell
            for irow in range(n_quad):
                xi = qp[irow]
                rho_val = rho_a + drho_cell * xi
                for icol in range(n_quad):
                    eta = qp[icol]
                    th_val = th_a + dtheta_cell * eta
                    rho_q_cpu[i_cell, j_cell, irow, icol] = rho_val
                    theta_q_cpu[i_cell, j_cell, irow, icol] = th_val
                    jac_q_cpu[i_cell, j_cell, irow, icol] = (
                        rho_val * drho_cell * dtheta_cell
                        * qw[irow] * qw[icol])

    r_q_full = rho_q_cpu * np.sin(theta_q_cpu)
    z_q_full = rho_q_cpu * np.cos(theta_q_cpu)

    n_q_total = N_rho * N_theta * n_quad * n_quad
    r_flat_cpu = r_q_full.reshape(n_q_total)
    z_flat_cpu = z_q_full.reshape(n_q_total)
    jac_flat_cpu = jac_q_cpu.reshape(n_q_total)

    # --- B aggregator (n_q_total x n_node), 9 nonzeros per row ---
    if verbose:
        t0 = time.time()

    # Shape function values at quadrature: (i_basis_local, j_basis_local) of cell
    # i, j in {0, 1, 2} for (corner_lo, mid, corner_hi)
    xi_q = qp                       # (n_quad,)
    eta_q = qp
    L_xi = np.zeros((3, n_quad))     # L_i(xi_q[k])
    L_xi[0] = _L0(xi_q)
    L_xi[1] = _L1(xi_q)
    L_xi[2] = _L2(xi_q)
    L_eta = np.zeros((3, n_quad))
    L_eta[0] = _L0(eta_q)
    L_eta[1] = _L1(eta_q)
    L_eta[2] = _L2(eta_q)

    # Build dense B for now (will check memory)
    # B[q, node] = N_local at q for the node within q's cell, 0 elsewhere
    B_dense_cpu = np.zeros((n_q_total, n_node), dtype=np.float64)

    # Precompute basis at all (irow, icol) for each (i_basis, j_basis):
    # N_ij_q[i, j, irow, icol] = L_i(xi_q[irow]) * L_j(eta_q[icol])
    N_ij_q = np.einsum("ir,js->ijrs", L_xi, L_eta)  # (3, 3, n_quad, n_quad)

    # Each cell (i_cell, j_cell) has 9 local nodes; global node index:
    # node[(i_cell, j_cell), (i_basis, j_basis)] =
    #    (2 i_cell + i_basis) * N_node_theta + (2 j_cell + j_basis)
    for i_cell in range(N_rho):
        for j_cell in range(N_theta):
            cell_idx = i_cell * N_theta + j_cell
            q_start = cell_idx * n_quad * n_quad
            for i_basis in range(3):
                glob_i = 2 * i_cell + i_basis
                for j_basis in range(3):
                    glob_j = 2 * j_cell + j_basis
                    glob_node = glob_i * N_node_theta + glob_j
                    # Set B[q_start + (irow*n_quad + icol), glob_node] = N_ij_q[i_basis, j_basis, irow, icol]
                    B_block = N_ij_q[i_basis, j_basis]  # (n_quad, n_quad)
                    for irow in range(n_quad):
                        for icol in range(n_quad):
                            q_idx = q_start + irow * n_quad + icol
                            B_dense_cpu[q_idx, glob_node] = B_block[irow, icol]

    if verbose:
        print(f"  B_dense ({B_dense_cpu.shape}, "
              f"{B_dense_cpu.nbytes/1e9:.2f} GB CPU): {time.time()-t0:.2f}s",
              flush=True)

    # Move to GPU
    t0 = time.time()
    B_gpu = cp.asarray(B_dense_cpu)
    r_flat = cp.asarray(r_flat_cpu)
    z_flat = cp.asarray(z_flat_cpu)
    jac_flat = cp.asarray(jac_flat_cpu)
    del B_dense_cpu, rho_q_cpu, theta_q_cpu, r_q_full, z_q_full, jac_q_cpu
    gc.collect()
    if verbose:
        print(f"  GPU transfer: {time.time()-t0:.2f}s, "
              f"B_gpu {B_gpu.nbytes/1e9:.2f} GB", flush=True)

    # G_full
    t0 = time.time()
    R_obs = r_flat[:, None]; Z_obs = z_flat[:, None]
    R_src = r_flat[None, :]; Z_src = z_flat[None, :]
    G = axisym_kernel_gpu(R_obs, Z_obs, R_src, Z_src, cp, cp_ellipk, ellipe_agm)
    if verbose:
        print(f"  G_full ({n_q_total}^2): {time.time()-t0:.2f}s, "
              f"{G.nbytes/1e9:.2f} GB", flush=True)

    # K = (mu_0/2) B^T Gw B with Gw[i,j] = G * 2pi r_i jac_i * jac_j
    t0 = time.time()
    Gw = G * (2.0 * math.pi * r_flat * jac_flat)[:, None] * jac_flat[None, :]
    del G
    cp.get_default_memory_pool().free_all_blocks()
    Gw_B = Gw @ B_gpu
    del Gw
    cp.get_default_memory_pool().free_all_blocks()
    K_full = (MU0 / 2.0) * (B_gpu.T @ Gw_B)
    del Gw_B
    K_full = 0.5 * (K_full + K_full.T)
    if verbose:
        print(f"  K = B^T Gw B: {time.time()-t0:.2f}s", flush=True)

    # M = B^T diag(2 pi r jac) B
    t0 = time.time()
    weights_obs = 2.0 * math.pi * r_flat * jac_flat
    M_full = B_gpu.T @ (weights_obs[:, None] * B_gpu)
    M_full = 0.5 * (M_full + M_full.T)
    if verbose:
        print(f"  M = B^T D B: {time.time()-t0:.2f}s", flush=True)

    # b
    A_q = (B0 / 2.0) * r_flat
    weights_src = 2.0 * math.pi * r_flat * jac_flat
    b_full = B_gpu.T @ (A_q * weights_src)

    # Restrict interior
    idx_gpu = cp.asarray(interior_idx)
    K_int = K_full[idx_gpu][:, idx_gpu]
    M_int = M_full[idx_gpu][:, idx_gpu]
    b_int = b_full[idx_gpu]

    K_cpu = cp.asnumpy(K_int)
    M_cpu = cp.asnumpy(M_int)
    b_cpu = cp.asnumpy(b_int)

    del K_full, M_full, b_full, K_int, M_int, b_int, B_gpu
    del r_flat, z_flat, jac_flat, weights_obs, weights_src, A_q, idx_gpu
    gc.collect()
    cp.get_default_memory_pool().free_all_blocks()

    return K_cpu, M_cpu, b_cpu, interior_idx, n_node


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
    print(f"\n--- Q2 {N_rho}x{N_theta} = {N_rho*N_theta} cells, n_quad={n_quad} ---", flush=True)
    t0 = time.time()
    K, M, b, _, _ = assemble_K_M_b_q2(N_rho, N_theta, n_quad)
    t_K = time.time() - t0
    print(f"  total assemble: {t_K:.1f}s, K shape: {K.shape}", flush=True)

    t0 = time.time()
    eigvals, eigvecs = eigh(K, M)
    t_eigh = time.time() - t0
    print(f"  eigh CPU: {t_eigh:.1f}s, "
          f"λ ∈ [{eigvals.min():.3e}, {eigvals.max():.3e}]", flush=True)

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

    print(f"  leading tau_1 = {tau_pos[order[0]]*1e6:.6f} us", flush=True)
    if rungs:
        print(f"  tau_pair[0] = {rungs[0]['tau_pair_us']:.6f} us", flush=True)

    return {
        "basis": "Q2",
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

    print("=" * 78, flush=True)
    print("Sphere VIM Q2 biquadratic GPU (Strang-Fix h^4 test)", flush=True)
    print(f"  Stoll Foster tau_1 = {stoll_tau[0]:.6f} us", flush=True)
    print(f"  Stoll Cauer tau_pair[0] = {stoll_cauer[0]:.6f} us", flush=True)
    print("=" * 78, flush=True)

    # Q2 doubles DoF density per direction. Start small to test correctness:
    # 2430 cells (45 * 54), DoF ≈ (2*45)*(2*54) = 9720 → close to 9700 interior
    levels = [
        (15, 18, 3),    # 270 cells, DoF ~ 1085  (small test)
        (30, 36, 3),    # 1080 cells, DoF ~ 4495
        (45, 54, 3),    # 2430 cells, DoF ~ 9700
    ]

    results = []
    for N_rho, N_theta, n_quad in levels:
        try:
            r = extract(N_rho, N_theta, n_quad, top_n=10, n_stages=12)
            results.append(r)
            (here / f"sphere_vim_axisym_q2_gpu_{N_rho}x{N_theta}.json").write_text(
                json.dumps(r, indent=2), encoding="utf-8")
        except Exception as e:
            import traceback
            print(f"  FAIL {N_rho}x{N_theta}: {type(e).__name__}: {e}", flush=True)
            traceback.print_exc()
        gc.collect()

    print(flush=True)
    print("=" * 78, flush=True)
    print(f"Q2 Leading Foster tau_1 vs Stoll = {stoll_tau[0]:.6f} us", flush=True)
    print("=" * 78, flush=True)
    for r in results:
        g = (r["leading_tau_us"] - stoll_tau[0]) / stoll_tau[0] * 100
        print(f"  {r['n_cells']:>6d} cells, {r['n_dofs']:>6d} DoF | "
              f"t_K={r['t_assemble_s']:>5.1f}s eigh={r['t_eigh_s']:>5.1f}s | "
              f"tau_1={r['leading_tau_us']:>10.6f} gap={g:>+10.6f}%", flush=True)

    print(flush=True)
    print("Q2 Cauer tau_pair[k] gap to Stoll [%]:", flush=True)
    hdr = [f"{'k':>3s}"] + [f"{str(r['n_cells']) + 'c':>14s}" for r in results]
    print("  " + " | ".join(hdr), flush=True)
    for k in range(min(8, len(stoll_cauer))):
        row = [f"{k:>3d}"]
        for r in results:
            if k < len(r["cauer_rungs"]):
                gp = (r["cauer_rungs"][k]["tau_pair_us"] - stoll_cauer[k]) / stoll_cauer[k] * 100
                row.append(f"{gp:>+14.6f}")
            else:
                row.append(f"{'--':>14s}")
        print("  " + " | ".join(row), flush=True)

    out = {"stoll_tau_us": stoll_tau[:20], "stoll_cauer_tau_us": stoll_cauer,
           "basis": "Q2", "results": results}
    (here / "sphere_vim_axisym_q2_gpu_summary.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8")
    print("\nSaved sphere_vim_axisym_q2_gpu_summary.json", flush=True)


if __name__ == "__main__":
    main()
