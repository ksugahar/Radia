"""cylinder_vim_axisym_cupy.py — GPU-accelerated axisym VIM for Cu disk.

Vectorized port of cylinder_vim_axisym.py:
  - rectangular (r, z) cell grid, Nr × Nz cells
  - Piecewise constant J_phi basis (1 dof per cell)
  - axisym Green's function G with elliptic K(m), E(m) (scipy CPU fallback)
  - K matrix via GPU einsum

Reference (BEM Kameari): leading tau ~ 209 us. axifemm Q1 Kameari: 207 us.
"""
from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
import scipy.special as sp_cpu
from scipy.linalg import eigh
import mpmath as mp
mp.mp.dps = 80


R_DISK = 10e-3
T_DISK = 2e-3
SIGMA = 5.8e7
MU0 = 4 * math.pi * 1e-7
B0 = 1.0


def axisym_kernel_vec(r, z, r_p, z_p):
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


def assemble_K_M_b_cylinder(Nr, Nz, n_quad=2, use_gpu=True, verbose=True):
    """Vectorized cylinder axisym VIM K, M, b assembly."""
    if use_gpu:
        import cupy as xp
    else:
        xp = np

    N = Nr * Nz
    dr = R_DISK / Nr
    dz = T_DISK / Nz

    # Cell coords (CPU)
    r_a = np.tile(np.arange(Nr)[:, None], (1, Nz)) * dr   # (Nr, Nz)
    r_b = r_a + dr
    z_a = np.tile(np.arange(Nz)[None, :], (Nr, 1)) * dz - T_DISK / 2.0
    z_b = z_a + dz

    # Mass M_ii = pi * (r_b² - r_a²) * (z_b - z_a)  (toroidal volume)
    M_diag = (math.pi * (r_b ** 2 - r_a ** 2) * (z_b - z_a)).flatten()  # (N,)

    # b_i = (B0/2) * (2 pi/3) * (r_b³ - r_a³) * (z_b - z_a)
    b_vec = ((B0 / 2.0) * (2.0 * math.pi / 3.0)
             * (r_b ** 3 - r_a ** 3) * (z_b - z_a)).flatten()

    # Quadrature points (n_quad × n_quad per cell)
    nodes_m1, w_m1 = np.polynomial.legendre.leggauss(n_quad)
    nodes_01 = 0.5 * (nodes_m1 + 1.0)
    w_01 = 0.5 * w_m1

    r_q = (r_a[..., None]
            + (r_b - r_a)[..., None] * nodes_01[None, None, :])  # (Nr, Nz, n_quad)
    z_q = (z_a[..., None]
            + (z_b - z_a)[..., None] * nodes_01[None, None, :])
    wr_q = ((r_b - r_a)[..., None] * w_01[None, None, :])
    wz_q = ((z_b - z_a)[..., None] * w_01[None, None, :])

    # Tensor product to get (Nr, Nz, n_quad, n_quad) per cell
    # In rectangular (r, z), point at (irow, izw) in cell (i, j):
    #   r = r_q[i, j, irow], z = z_q[i, j, izw]
    r_full = np.broadcast_to(r_q[..., :, None],
                              (Nr, Nz, n_quad, n_quad))   # (..., n_quad, n_quad)
    z_full = np.broadcast_to(z_q[..., None, :],
                              (Nr, Nz, n_quad, n_quad))
    jac_full = (wr_q[..., :, None] * wz_q[..., None, :])    # (..., n_quad, n_quad)

    n_q_total = Nr * Nz * n_quad * n_quad
    r_flat = np.ascontiguousarray(r_full).reshape(n_q_total)
    z_flat = np.ascontiguousarray(z_full).reshape(n_q_total)
    jac_flat = np.ascontiguousarray(jac_full).reshape(n_q_total)

    if verbose:
        print(f"  N = {N} cells, n_q_total = {n_q_total}")
        mem_K = n_q_total * n_q_total * 8 / 1e9
        print(f"  Kernel matrix G_full: {mem_K:.2f} GB FP64")

    # G_full[qi, qj] = G(r_qi, z_qi, r_qj, z_qj)
    if verbose:
        t0 = time.time()
    R_obs = r_flat[:, None]
    Z_obs = z_flat[:, None]
    R_src = r_flat[None, :]
    Z_src = z_flat[None, :]
    G_full = axisym_kernel_vec(R_obs, Z_obs, R_src, Z_src)
    if verbose:
        print(f"  G_full computed in {time.time()-t0:.2f}s")

    # Weighted G: G_w[qi, qj] = G * (2 pi r_qi) * jac_qi * jac_qj
    if use_gpu:
        if verbose:
            t0 = time.time()
        G_gpu = xp.asarray(G_full)
        r_flat_gpu = xp.asarray(r_flat)
        jac_flat_gpu = xp.asarray(jac_flat)
        Gw = (G_gpu
              * (2.0 * math.pi * r_flat_gpu * jac_flat_gpu)[:, None]
              * jac_flat_gpu[None, :])
        n_q_per_cell = n_quad * n_quad
        Gw_3d = Gw.reshape(N, n_q_per_cell, N, n_q_per_cell)
        K = (MU0 / 2.0) * Gw_3d.sum(axis=(1, 3))
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
    return K, np.diag(M_diag), b_vec


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
    cases = [
        (48, 12, "576 cells (baseline)"),
        (96, 24, "2304 cells (4x)"),
        (144, 36, "5184 cells (9x)"),
    ]

    out_results = {}
    for Nr, Nz, label in cases:
        print(f"\n=== Cylinder axisym VIM ({label}) ===")
        print(f"  Nr={Nr}, Nz={Nz}")

        t0 = time.time()
        K, M, b = assemble_K_M_b_cylinder(Nr, Nz, n_quad=2,
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
        leading_foster = float(tau_us[order[0]])
        print(f"  Foster top tau = {leading_foster:.4f} us")

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
            "Nr": Nr, "Nz": Nz, "n_dofs": Nr * Nz,
            "K_assembly_s": t_K,
            "leading_foster_us": leading_foster,
            "Cauer_rungs_Kameari": rungs,
        }

    out_path = Path(__file__).parent / "cylinder_vim_axisym_cupy_results.json"
    out_path.write_text(json.dumps(out_results, indent=2), encoding="utf-8")
    print(f"\nSaved: {out_path}")


if __name__ == "__main__":
    main()
