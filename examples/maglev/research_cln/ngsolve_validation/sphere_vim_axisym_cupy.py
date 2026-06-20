"""sphere_vim_axisym_cupy.py — GPU-accelerated axisym VIM for Cu sphere.

Vectorized CuPy port of sphere_vim_axisym.py:
  - (rho, theta) cell grid covering sphere interior
  - Piecewise constant J_phi basis (1 dof per cell)
  - axisym Green's function G(r, z; r', z') with elliptic K(m), E(m)
  - K matrix via einsum on GPU

Reference Stoll analytical leading tau: 694.14 us (Cu sphere R=10mm).

Note: cupyx.scipy.special has ellipk but NOT ellipe. We compute on CPU
(scipy.special.ellipe) for the numerical evaluation (cheap compared to
the 4D quadrature loop), then transfer to GPU for the matmul.

Actually — simpler: do the entire kernel computation on GPU with
manual ellipe via Carlson form, OR keep elliptic eval on CPU since
n_quad scales modestly.
"""
from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
import scipy.special as sp_cpu  # for ellipe (not in cupyx)

import mpmath as mp
mp.mp.dps = 80


R_SPHERE = 10e-3
SIGMA = 5.8e7
MU0 = 4 * math.pi * 1e-7
B0 = 1.0


def axisym_kernel_vec(r, z, r_p, z_p):
    """Vectorized axisym Green's function on numpy arrays.

    Returns the bare bracket without mu0/pi prefactor (matches scalar
    version in cylinder_vim_axisym.py).
    """
    safe = (r > 1e-15) & (r_p > 1e-15)
    d2 = (r + r_p) ** 2 + (z - z_p) ** 2
    m = np.where(safe, 4.0 * r * r_p / d2, 0.0)
    m = np.minimum(m, 1.0 - 1e-12)
    k = np.sqrt(m)
    K = sp_cpu.ellipk(m)
    E = sp_cpu.ellipe(m)
    coef = np.sqrt(np.where(safe, r_p / r, 1.0)) / math.pi
    val = coef * ((2.0 / np.maximum(k, 1e-30) - k) * K
                  - (2.0 / np.maximum(k, 1e-30)) * E)
    return np.where(safe, val, 0.0)


def assemble_K_M_b_vectorized(N_rho, N_theta, n_quad=2, dtype=np.float64,
                                use_gpu=True, verbose=True):
    """Vectorized assembly of K, M, b for axisym sphere VIM.

    Strategy:
      1. Cells: list (rho_a, rho_b, th_a, th_b) for each of N = N_rho * N_theta
      2. Quadrature points: per-cell (n_quad)^2 in (rho, theta), total
         n_q_per_cell = n_quad^2, n_q_total = N * n_quad^2
      3. For each pair (i, j) of cells: K[i,j] = sum_qi sum_qj G(qi, qj) ×
         weight_i × weight_j. Vectorize as matmul:
           G_full[qi, qj] = G(rqi, zqi, rqj, zqj) — large (n_q_total)^2 matrix
           K[i, j] = sum over (qi in cell_i, qj in cell_j)
                   = block sum of G_full × weights

      4. Cleaner: build (n_q_total, N) basis matrix B[q, i] where
         B[q, i] = (2 pi r_q) * weight_i_q  if q ∈ cell i, else 0.
         Then K = B^T G_full B where G is the kernel matrix.

      Memory check: N=1080, n_quad=2 → n_q_total = 4320.
      G_full: 4320^2 × 8 bytes = 150 MB ✓
      For N_rho=60, N_theta=72 → N=4320, n_q_total=17280, G_full=2.4 GB
      For N_rho=80, N_theta=96 → N=7680, n_q_total=30720, G_full=7.5 GB
      For N_rho=120, N_theta=144 → N=17280, n_q_total=69120, G_full=38 GB (no)
    """
    if use_gpu:
        import cupy as xp
    else:
        xp = np
    N = N_rho * N_theta
    drho = R_SPHERE / N_rho
    dtheta = math.pi / N_theta

    # Cell coords (CPU, small arrays)
    rho_a = np.tile(np.arange(N_rho)[:, None], (1, N_theta)) * drho   # (N_rho, N_theta)
    rho_b = rho_a + drho
    th_a = np.tile(np.arange(N_theta)[None, :], (N_rho, 1)) * dtheta
    th_b = th_a + dtheta

    # Mass: M_ii = 2pi/3 * (rho_b^3 - rho_a^3) * (cos th_a - cos th_b)
    M_diag = (2.0 * math.pi / 3.0
              * (rho_b ** 3 - rho_a ** 3)
              * (np.cos(th_a) - np.cos(th_b))).flatten()  # (N,)

    # b_i = pi B0 * (rho_b^4 - rho_a^4)/4 * [(tb - ta)/2 - (sin 2 tb - sin 2 ta)/4]
    b_vec = (math.pi * B0 * (rho_b ** 4 - rho_a ** 4) / 4.0
             * ((th_b - th_a) / 2.0
                - (np.sin(2.0 * th_b) - np.sin(2.0 * th_a)) / 4.0)).flatten()

    # ----- Quadrature points (n_quad x n_quad per cell) ----------------------
    nodes_m1, w_m1 = np.polynomial.legendre.leggauss(n_quad)
    nodes_01 = 0.5 * (nodes_m1 + 1.0)        # in [0, 1]
    w_01 = 0.5 * w_m1

    # Per-cell quadrature: rho_q[N, n_q, n_q] and th_q[N, n_q, n_q]
    # but separable: rho only depends on rho-index, theta on theta-index
    # rho_q[cell_idx, irho] shape (N, n_quad), where cell_idx = rho_idx*N_theta + th_idx
    # Actually cleaner: build (N_rho, N_theta, n_quad) for rho; (N_rho, N_theta, n_quad) for th
    rho_q = (rho_a[..., None]
             + (rho_b - rho_a)[..., None] * nodes_01[None, None, :])   # (N_rho, N_theta, n_quad)
    th_q = (th_a[..., None]
             + (th_b - th_a)[..., None] * nodes_01[None, None, :])      # (N_rho, N_theta, n_quad)
    wrho_q = ((rho_b - rho_a)[..., None] * w_01[None, None, :])          # (N_rho, N_theta, n_quad)
    wth_q = ((th_b - th_a)[..., None] * w_01[None, None, :])             # (N_rho, N_theta, n_quad)

    # Tensor product (irho, ith) → r, z, weight per (cell, quad point combination)
    # cell_idx = rho_idx*N_theta + th_idx
    # For each cell, n_q_per_cell = n_quad^2 quadrature points indexed by (irho, ith)
    # r_qi = rho_qi * sin(th_qi)
    # z_qi = rho_qi * cos(th_qi)
    # jac = rho_qi * wrho_qi * wth_qi  (Jacobian dr dz = rho drho dtheta)
    # 2 pi r factor for K assembly

    # Shape (N_rho, N_theta, n_quad, n_quad) for r and z:
    r_q = rho_q[..., :, None] * np.sin(th_q[..., None, :])        # broadcast to (..., n_quad, n_quad)
    z_q = rho_q[..., :, None] * np.cos(th_q[..., None, :])
    jac_q = (rho_q[..., :, None]                                    # rho × wrho × wth
             * wrho_q[..., :, None]
             * wth_q[..., None, :])

    # Flatten to (n_q_total,) where n_q_total = N_rho * N_theta * n_quad^2
    n_q_total = N_rho * N_theta * n_quad * n_quad
    r_flat = r_q.reshape(n_q_total)
    z_flat = z_q.reshape(n_q_total)
    jac_flat = jac_q.reshape(n_q_total)

    # cell index per quadrature point: cell_idx = (irho * N_theta + ith) for each q
    # When we reshape (N_rho, N_theta, n_quad, n_quad), the cell index varies
    # over the first 2 dims. So q = (irho * N_theta * n_q^2 + ith * n_q^2 + ...)
    # cell_idx_for_q[q] = q // (n_quad^2)
    cell_idx = np.repeat(np.arange(N), n_quad * n_quad)

    if verbose:
        print(f"  N = {N} cells, n_q_total = {n_q_total}")
        mem_K = n_q_total * n_q_total * 8 / 1e9
        print(f"  Kernel matrix G_full: {mem_K:.2f} GB FP64")

    # ----- Build kernel matrix G_full[qi, qj] = G(r_qi, z_qi, r_qj, z_qj) ----
    # G is symmetric in (qi <-> qj) but kernel definition uses sqrt(r'/r) so
    # G(qi, qj) ≠ G(qj, qi) in our formulation. The integration is symmetric
    # under qi <-> qj swap (since (2 pi r) Jacobian on observation side
    # paired with G coefficient sqrt(r'/r)). Let's just compute full matrix.

    if verbose:
        t0 = time.time()
    # Compute on CPU with full vectorization (G uses scipy ellipk, ellipe)
    # Then transfer to GPU for matmul
    R_obs = r_flat[:, None]   # (n_q, 1)
    Z_obs = z_flat[:, None]
    R_src = r_flat[None, :]   # (1, n_q)
    Z_src = z_flat[None, :]
    G_full = axisym_kernel_vec(R_obs, Z_obs, R_src, Z_src)
    if verbose:
        print(f"  G_full computed (CPU scipy.special) in {time.time()-t0:.2f}s")

    # ----- Assemble K[i, j] = sum_{qi in cell_i, qj in cell_j} ...  ---------
    # Aggregator matrix B[q, i]: 1 if q ∈ cell i, else 0. Combined with weights.
    # K_qq = G_full × (jac_qj × 2π r_qi)   -> kij = sum over qi, qj in cells
    # Actually:
    # K[i, j] = (mu0/2) * sum_{qi in i, qj in j} G(qi, qj) × (2 pi r_qi) × jac_qi × jac_qj
    # The (2 pi r_qi) is observation jacobian (toroidal).

    # Weighted G: G_w[qi, qj] = G(qi, qj) × (2 pi r_qi) × jac_qi × jac_qj
    # K[i, j] = (mu0/2) * sum_{qi, qj} G_w[qi, qj] × indicator_i[qi] × indicator_j[qj]

    if use_gpu:
        if verbose:
            t0 = time.time()
        G_gpu = xp.asarray(G_full)
        r_flat_gpu = xp.asarray(r_flat)
        jac_flat_gpu = xp.asarray(jac_flat)

        # Pre-multiply: G[qi, qj] *= (2 pi r_qi) × jac_qi × jac_qj
        Gw = (G_gpu
              * (2.0 * math.pi * r_flat_gpu * jac_flat_gpu)[:, None]
              * jac_flat_gpu[None, :])

        # Sum over (qi, qj) blocks of size n_quad^2. Build assembler matrix
        # A[q, cell] = 1 if q in cell, else 0.
        # K[i, j] = (mu0/2) * A^T Gw A
        n_q_per_cell = n_quad * n_quad
        # A is structured: q ∈ [cell * n_q_per_cell, (cell+1) * n_q_per_cell)
        # So A^T Gw A is a block sum: K[i,j] = sum over n_q_per_cell rows of cell i
        #   and n_q_per_cell cols of cell j of Gw[qi, qj].
        Gw_3d = Gw.reshape(N, n_q_per_cell, N, n_q_per_cell)
        K = Gw_3d.sum(axis=(1, 3))  # (N, N)
        K = (MU0 / 2.0) * K

        if verbose:
            xp.cuda.Stream.null.synchronize()
            print(f"  GPU einsum/sum in {time.time()-t0:.2f}s")
        K = xp.asnumpy(K)
    else:
        if verbose:
            t0 = time.time()
        Gw = (G_full
              * (2.0 * math.pi * r_flat * jac_flat)[:, None]
              * jac_flat[None, :])
        n_q_per_cell = n_quad * n_quad
        Gw_3d = Gw.reshape(N, n_q_per_cell, N, n_q_per_cell)
        K = (MU0 / 2.0) * Gw_3d.sum(axis=(1, 3))
        if verbose:
            print(f"  CPU sum in {time.time()-t0:.2f}s")

    K = 0.5 * (K + K.T)
    M = np.diag(M_diag)
    return K, M, b_vec


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


def main():
    from scipy.linalg import eigh

    cases = [
        (30, 36, "fine (1080 cells, baseline)"),
        (45, 54, "extra (2430 cells)"),
        (60, 72, "very fine (4320 cells)"),
    ]

    out_results = {}
    for N_rho, N_theta, label in cases:
        print(f"\n=== Sphere axisym VIM ({label}) ===")
        print(f"  N_rho={N_rho}, N_theta={N_theta}")

        t0 = time.time()
        K, M, b = assemble_K_M_b_vectorized(N_rho, N_theta, n_quad=2,
                                              use_gpu=True, verbose=True)
        t_K = time.time() - t0
        print(f"  Total K, M, b: {t_K:.2f}s")

        eigvals, eigvecs = eigh(K, M)
        for k in range(len(eigvals)):
            v = eigvecs[:, k]
            nrm = math.sqrt(float(v @ M @ v))
            if nrm > 1e-30:
                eigvecs[:, k] = v / nrm

        tau_us = np.where(eigvals > 1e-30, 1e6 * SIGMA * eigvals, 0.0)
        g = eigvecs.T @ b
        g2 = g * g
        order = np.argsort(-g2)
        leading_foster_us = float(tau_us[order[0]])
        print(f"  Foster top tau = {leading_foster_us:.4f} us")

        n_use = min(120, int(np.sum(g2 > 1e-32)))
        sorted_idx = order[:n_use]
        tau_seconds = [mp.mpf(float(tau_us[i])) * mp.mpf("1e-6")
                       for i in sorted_idx]
        g2_mp = [mp.mpf(float(g2[i])) for i in sorted_idx]
        n_moments = min(40, 2 * n_use - 4)
        alphas = []
        for n in range(n_moments):
            a_val = mp.mpf(0)
            for gv, tv in zip(g2_mp, tau_seconds):
                a_val += gv * tv ** n
            alphas.append((mp.mpf(-1)) ** n * a_val)

        cauer_p = hankel_pade_cauer(alphas, 12)
        rungs = []
        print(f"  Kameari Cauer-I rungs:")
        print(f"  k | R_2k       | L_2k+1     | tau_pair us")
        for k in range(min(6, len(cauer_p) // 2)):
            Rv = cauer_p[2 * k]
            Linv = cauer_p[2 * k + 1]
            if abs(Linv) < mp.mpf(10) ** (-50):
                break
            Lv = 1 / Linv
            tau_p = float(Lv / Rv) * 1e6
            print(f"  {k} | {float(Rv):>10.4e} | {float(Lv):>10.4e} | {tau_p:.4f}")
            rungs.append({"k": k, "R_2k": float(Rv),
                          "L_2k_plus_1": float(Lv), "tau_pair_us": tau_p})

        out_results[label] = {
            "N_rho": N_rho, "N_theta": N_theta, "n_dofs": N_rho * N_theta,
            "K_assembly_s": t_K,
            "leading_foster_us": leading_foster_us,
            "Cauer_rungs_Kameari": rungs,
        }

    out_path = Path(__file__).parent / "sphere_vim_axisym_cupy_results.json"
    out_path.write_text(json.dumps(out_results, indent=2), encoding="utf-8")
    print(f"\nSaved: {out_path}")


if __name__ == "__main__":
    main()
