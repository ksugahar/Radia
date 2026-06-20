"""sphere_vim_axisym_linear_gpu.py — GPU-native linear-basis sphere VIM.

Eliminates CPU intermediate G_full bottleneck by computing the axisymmetric
Green's function kernel entirely on the GPU via cupy.  Uses cupyx ellipk
plus a custom AGM-based ellipe (cupyx lacks ellipe).

Target machine: A100 40 GB or equivalent, with cupy-cuda12x installed.

Memory budget (FP64 GPU):
  N_q = N_cells * n_quad^2
  G_full = N_q^2 * 8 bytes  (dominant)
  40 GB cap ⇒ N_q < 70710
    n_quad=3: cells < 7853
    n_quad=2: cells < 17677

For each level we report Foster spectrum + Cauer τ_pair[k] + top Foster
modes, comparing to the 240-digit Stoll analytical reference.
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

mp.mp.dps = 60


R_SPHERE = 10e-3
SIGMA = 5.8e7
MU0 = 4 * math.pi * 1e-7
B0 = 1.0


def _import_cupy():
    """Import cupy + cupyx; provide AGM-based ellipe kernel since cupyx lacks it."""
    import cupy as cp
    from cupyx.scipy.special import ellipk as cp_ellipk

    # AGM-based ellipe (custom kernel)
    ellipe_agm = cp.ElementwiseKernel(
        in_params="float64 m",
        out_params="float64 E_out",
        operation=r"""
        double a = 1.0;
        double b = sqrt(1.0 - m);
        double c2 = m;            // c0^2 = m
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
        double K = 1.5707963267948966 / a;   // K(m) = pi/(2a_inf)
        E_out = K * (1.0 - sum_c2);
        """,
        name="ellipe_agm",
    )
    return cp, cp_ellipk, ellipe_agm


def axisym_kernel_gpu(r, z, r_p, z_p, cp, cp_ellipk, ellipe_agm):
    """GPU axisymmetric Green's function (without mu_0/2 prefactor).
    All arguments are cupy arrays. Returns cupy array same shape."""
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


def assemble_K_M_b_linear_gpu(N_rho, N_theta, n_quad=3, verbose=True):
    """GPU bilinear J_phi basis on (rho, theta) sphere mesh."""
    cp, cp_ellipk, ellipe_agm = _import_cupy()

    n_vert_rho = N_rho + 1
    n_vert_theta = N_theta + 1
    n_vert = n_vert_rho * n_vert_theta
    drho = R_SPHERE / N_rho
    dtheta = math.pi / N_theta

    # Interior mask (axis is constrained; sphere surface free)
    interior = np.zeros((n_vert_rho, n_vert_theta), dtype=bool)
    interior[1:, 1:-1] = True
    interior_flat = interior.flatten()
    interior_idx = np.where(interior_flat)[0]
    n_interior = len(interior_idx)

    if verbose:
        print(f"  N_rho={N_rho}, N_theta={N_theta}, cells={N_rho*N_theta}, "
              f"interior DoF={n_interior}")

    # --- Quadrature on reference cell [0,1]^2 ---
    nodes_m1, w_m1 = np.polynomial.legendre.leggauss(n_quad)
    qp_h = 0.5 * (nodes_m1 + 1.0)
    qw_h = 0.5 * w_m1
    qp = cp.asarray(qp_h)
    qw = cp.asarray(qw_h)

    # Build (N_rho, N_theta, n_quad, n_quad) quadrature arrays on GPU
    rho_v = cp.linspace(0, R_SPHERE, n_vert_rho)
    theta_v = cp.linspace(0, math.pi, n_vert_theta)

    rho_q = (rho_v[:-1, None, None, None]
             + drho * qp[None, None, :, None])
    rho_q = cp.broadcast_to(rho_q, (N_rho, N_theta, n_quad, n_quad)).copy()
    theta_q = (theta_v[None, :-1, None, None]
               + dtheta * qp[None, None, None, :])
    theta_q = cp.broadcast_to(theta_q, (N_rho, N_theta, n_quad, n_quad)).copy()
    r_q_full = rho_q * cp.sin(theta_q)
    z_q_full = rho_q * cp.cos(theta_q)
    jac_q_full = (rho_q * drho * dtheta
                  * qw[None, None, :, None] * qw[None, None, None, :])

    n_q_total = N_rho * N_theta * n_quad * n_quad
    r_flat = r_q_full.reshape(n_q_total)
    z_flat = z_q_full.reshape(n_q_total)
    jac_flat = jac_q_full.reshape(n_q_total)

    # --- Bilinear shape values at q-points (CPU, small) ---
    xi = qp_h
    eta = qp_h
    N_00 = np.outer(1 - xi, 1 - eta).flatten()
    N_10 = np.outer(xi,     1 - eta).flatten()
    N_11 = np.outer(xi,     eta).flatten()
    N_01 = np.outer(1 - xi, eta).flatten()
    # Tile to (N_rho * N_theta, n_quad^2) by cell
    n_qq = n_quad * n_quad

    # --- Aggregator B as sparse COO via flat indexing on CPU, then to GPU ---
    if verbose:
        t0 = time.time()
    irho_arr, jtheta_arr = np.meshgrid(
        np.arange(N_rho), np.arange(N_theta), indexing="ij")
    irho_arr = irho_arr.flatten()
    jtheta_arr = jtheta_arr.flatten()

    v00 = irho_arr * n_vert_theta + jtheta_arr
    v10 = (irho_arr + 1) * n_vert_theta + jtheta_arr
    v11 = (irho_arr + 1) * n_vert_theta + (jtheta_arr + 1)
    v01 = irho_arr * n_vert_theta + (jtheta_arr + 1)

    n_cells_tot = N_rho * N_theta
    q_offsets = np.arange(n_cells_tot)[:, None] * n_qq + np.arange(n_qq)[None, :]
    q_offsets_flat = q_offsets.flatten()

    rows = np.repeat(q_offsets_flat, 4)
    cols = np.empty(4 * n_q_total, dtype=np.int64)
    vals = np.empty(4 * n_q_total, dtype=np.float64)

    # COO build (use vectorized assignment)
    cols[0::4] = np.repeat(v00, n_qq)
    cols[1::4] = np.repeat(v10, n_qq)
    cols[2::4] = np.repeat(v11, n_qq)
    cols[3::4] = np.repeat(v01, n_qq)
    vals[0::4] = np.tile(N_00, n_cells_tot)
    vals[1::4] = np.tile(N_10, n_cells_tot)
    vals[2::4] = np.tile(N_11, n_cells_tot)
    vals[3::4] = np.tile(N_01, n_cells_tot)

    # Build sparse B as CSR via scipy (CPU), then move to GPU
    from scipy.sparse import coo_matrix
    import cupyx.scipy.sparse as cpx_sp
    B_coo = coo_matrix((vals, (rows, cols)),
                       shape=(n_q_total, n_vert))
    B_csr_cpu = B_coo.tocsr()
    # Cupy sparse CSR
    B_gpu = cpx_sp.csr_matrix(
        (cp.asarray(B_csr_cpu.data), cp.asarray(B_csr_cpu.indices),
         cp.asarray(B_csr_cpu.indptr)),
        shape=(n_q_total, n_vert))
    if verbose:
        print(f"  B sparse built: nnz={B_csr_cpu.nnz}, "
              f"{time.time()-t0:.2f}s")

    # --- Compute G_full on GPU (dense) ---
    if verbose:
        t0 = time.time()
    R_obs = r_flat[:, None]
    Z_obs = z_flat[:, None]
    R_src = r_flat[None, :]
    Z_src = z_flat[None, :]
    G_full = axisym_kernel_gpu(R_obs, Z_obs, R_src, Z_src,
                                cp, cp_ellipk, ellipe_agm)
    if verbose:
        print(f"  G_full ({n_q_total}^2): {time.time()-t0:.2f}s, "
              f"{G_full.nbytes/1e9:.2f} GB on GPU")

    # Weighted: G_w[qi, qj] = G * (2 pi r_qi) * jac_qi * jac_qj
    if verbose:
        t0 = time.time()
    w_obs = 2.0 * math.pi * r_flat * jac_flat  # (N_q,)
    w_src = jac_flat                             # (N_q,)
    G_full *= w_obs[:, None]
    G_full *= w_src[None, :]
    G_full *= (MU0 / 2.0)
    if verbose:
        print(f"  G weighted: {time.time()-t0:.2f}s")

    # --- K = B^T G_full B  (vertex × vertex) ---
    if verbose:
        t0 = time.time()
    # Step 1: BG = G_full @ B  (n_q × n_vert)
    # Cupy sparse @ dense gives dense (n_q × n_vert)
    # But sparse times dense in cupy: cp.dot ... need explicit
    # B_gpu.T @ G_full @ B_gpu  → (n_vert × n_vert)
    # Use (B^T G_full) first then × B again
    GB = B_gpu.T @ G_full   # n_vert × n_q
    K_full = GB @ B_gpu     # n_vert × n_vert  (cupy dense via sparse @ dense)
    # Convert to dense if not already
    if hasattr(K_full, "toarray"):
        K_full = K_full.toarray()
    if verbose:
        print(f"  K assembled (full): {time.time()-t0:.2f}s, "
              f"K mem {K_full.nbytes/1e9:.2f} GB")

    # Symmetrize
    K_full = 0.5 * (K_full + K_full.T)

    # Restrict to interior DoFs
    K = K_full[cp.asarray(interior_idx)][:, cp.asarray(interior_idx)]

    # --- M (vertex × vertex) diagonal-like ---
    # M[v,w] = int N_v N_w dV  = int N_v N_w 2 pi r dr dz = 2 pi int N_v N_w (rho sin(theta)) rho drho dtheta
    # Build via B aggregator: M = B^T diag(2 pi r * jac) B
    if verbose:
        t0 = time.time()
    M_weight = 2.0 * math.pi * r_flat * jac_flat
    BMb = B_gpu.T.multiply(M_weight[None, :])
    M_full = BMb @ B_gpu
    if hasattr(M_full, "toarray"):
        M_full = M_full.toarray()
    M_full = 0.5 * (M_full + M_full.T)
    M = M_full[cp.asarray(interior_idx)][:, cp.asarray(interior_idx)]
    if verbose:
        print(f"  M assembled: {time.time()-t0:.2f}s, "
              f"M mem {M.nbytes/1e9:.2f} GB")

    # --- b vector ---
    # b_v = int A_imposed N_v dV  with A_imposed = (B0/2) r → b_v = π B0 int r^2 N_v drho dtheta * rho
    # Use B aggregator: b = B^T (π B0 r^2 * jac)
    b_q = math.pi * B0 * r_flat * r_flat * jac_flat
    b_full = B_gpu.T @ b_q
    b = b_full[cp.asarray(interior_idx)]

    # Transfer to CPU numpy for eigh
    K_cpu = cp.asnumpy(K)
    M_cpu = cp.asnumpy(M)
    b_cpu = cp.asnumpy(b)

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


def extract_one(N_rho, N_theta, n_quad=3, top_n=10, n_stages=12):
    n_cells = N_rho * N_theta
    print(f"\n--- {N_rho}x{N_theta} = {n_cells} cells, n_quad={n_quad} ---")
    t0 = time.time()
    K, M, b, _, _ = assemble_K_M_b_linear_gpu(N_rho, N_theta, n_quad=n_quad)
    t_K = time.time() - t0
    print(f"  total assemble: {t_K:.1f}s, K shape: {K.shape}")

    from scipy.linalg import eigh
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
            a_val += gv * tv ** k       # Kameari convention
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
        "n_cells": n_cells, "n_dofs": K.shape[0],
        "t_assemble_s": t_K, "t_eigh_s": t_eigh,
        "leading_tau_us": float(tau_pos[order[0]]) * 1e6,
        "top_modes": top_modes,
        "cauer_rungs": rungs,
    }


def main():
    here = Path(__file__).parent
    ana_path = here / "sphere_analytical_cln_highstage.json"
    ana = json.load(open(ana_path))
    stoll_tau = [float(x) for x in ana["stoll_tau_us"]]
    stoll_a = [float(x) for x in ana["stoll_a_n"]]
    stoll_cauer = [float(x) for x in ana["cauer_tau_us"]]

    print("=" * 78)
    print("Sphere VIM linear-basis DoF scaling on A100 GPU")
    print(f"  Stoll Foster τ_1 = {stoll_tau[0]:.6f} μs")
    print(f"  Stoll Cauer τ_pair[0] = {stoll_cauer[0]:.6f} μs")
    print("=" * 78)

    # n_quad=3 throughout for clean h² convergence. Memory budget per case:
    #   N_q = cells * 9, G = N_q² × 8 bytes
    levels = [
        (45, 54, 3),    # 2430 cells, G = 0.04 GB
        (60, 72, 3),    # 4320 cells, G = 11.3 GB
        (75, 90, 3),    # 6750 cells, G = 29.5 GB (fits 40 GB A100)
        (78, 96, 3),    # 7488 cells, G = 36.3 GB (tight)
    ]

    results = []
    for N_rho, N_theta, n_quad in levels:
        try:
            r = extract_one(N_rho, N_theta, n_quad=n_quad,
                            top_n=10, n_stages=12)
            results.append(r)
            out_path = (here /
                f"sphere_vim_axisym_linear_gpu_{N_rho}x{N_theta}.json")
            out_path.write_text(json.dumps(r, indent=2), encoding="utf-8")
            print(f"  Saved {out_path.name}")
        except Exception as e:
            print(f"  FAIL {N_rho}x{N_theta}: {type(e).__name__}: {e}")

    # Aggregated table
    print()
    print("=" * 78)
    print(f"Leading Foster τ_1 vs Stoll = {stoll_tau[0]:.6f} μs")
    print("=" * 78)
    print(f"  {'cells':>7s} | {'DoF':>7s} | {'t_K [s]':>9s} | "
          f"{'t_eigh':>8s} | {'τ_1':>12s} | {'gap [%]':>10s}")
    for r in results:
        g = (r["leading_tau_us"] - stoll_tau[0]) / stoll_tau[0] * 100
        print(f"  {r['n_cells']:>7d} | {r['n_dofs']:>7d} | "
              f"{r['t_assemble_s']:>9.1f} | {r['t_eigh_s']:>8.1f} | "
              f"{r['leading_tau_us']:>12.6f} | {g:>+10.4f}")

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

    out = {
        "stoll_tau_us": stoll_tau[:20],
        "stoll_a_n": stoll_a[:20],
        "stoll_cauer_tau_us": stoll_cauer,
        "results": results,
    }
    (here / "sphere_vim_axisym_linear_gpu_summary.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8")
    print("\nSaved sphere_vim_axisym_linear_gpu_summary.json")


if __name__ == "__main__":
    main()
