"""cylinder_vim_axisym_linear.py — Cylinder axisym VIM with LINEAR (bilinear) J_phi basis.

DOFs: 1 per vertex of (r, z) grid.
Shape functions: bilinear in (xi, eta) on each cell.
Boundary conditions: J_phi = 0 only at r=0 (axis). At r=R, z=±t/2 the
J_phi is FREE (J_n is automatic since we model only J_phi component).

Reference (BEM Kameari): leading tau ~ 209 us, NGSolve 3D ~ 208.34 us.
Expected accuracy: 0.1-0.2% leading rung gap with linear basis + n_quad=3.
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


def assemble_K_M_b_linear_cylinder(Nr, Nz, n_quad=3,
                                     use_gpu=True, verbose=True):
    """Bilinear J_phi basis on (r, z) cylinder mesh."""
    if use_gpu:
        import cupy as xp
    else:
        xp = np

    n_vert_r = Nr + 1
    n_vert_z = Nz + 1
    n_vert = n_vert_r * n_vert_z

    dr = R_DISK / Nr
    dz = T_DISK / Nz

    # ----- Boundary: only r = 0 (axis) -----
    interior = np.ones((n_vert_r, n_vert_z), dtype=bool)
    interior[0, :] = False   # r=0 axis vertices J_phi = 0
    interior_flat = interior.flatten()
    interior_idx = np.where(interior_flat)[0]
    n_interior = len(interior_idx)

    if verbose:
        print(f"  Nr={Nr}, Nz={Nz}")
        print(f"  Total vertices: {n_vert}, interior: {n_interior}")
        print(f"  Cells: {Nr * Nz}")

    # ----- Quadrature -----
    nodes_m1, w_m1 = np.polynomial.legendre.leggauss(n_quad)
    qp = 0.5 * (nodes_m1 + 1.0)
    qw = 0.5 * w_m1

    # Vertex coords
    r_v = np.linspace(0, R_DISK, n_vert_r)
    z_v = np.linspace(-T_DISK / 2, T_DISK / 2, n_vert_z)

    # Quadrature points per cell (Nr, Nz, n_quad, n_quad)
    r_q_full = (r_v[:-1, None, None, None]                 # (Nr, 1, 1, 1) lower edge
                + dr * qp[None, None, :, None])             # (1, 1, n_quad, 1) offset
    z_q_full = (z_v[None, :-1, None, None]                 # (1, Nz, 1, 1)
                + dz * qp[None, None, None, :])             # (1, 1, 1, n_quad)
    r_q_full = np.broadcast_to(r_q_full, (Nr, Nz, n_quad, n_quad)).copy()
    z_q_full = np.broadcast_to(z_q_full, (Nr, Nz, n_quad, n_quad)).copy()
    # Jacobian: dr dz factor per cell = dr * dz × GL weights
    jac_q_full = (dr * dz
                  * qw[None, None, :, None] * qw[None, None, None, :]
                  * np.ones((Nr, Nz, 1, 1)))

    # Shape functions
    xi = qp
    eta = qp
    N_00 = np.outer(1 - xi, 1 - eta)
    N_10 = np.outer(xi,     1 - eta)
    N_11 = np.outer(xi,     eta)
    N_01 = np.outer(1 - xi, eta)

    n_q_total = Nr * Nz * n_quad * n_quad
    r_flat = r_q_full.reshape(n_q_total)
    z_flat = z_q_full.reshape(n_q_total)
    jac_flat = jac_q_full.reshape(n_q_total)

    # Aggregator B[q, v]
    if verbose:
        t0 = time.time()
    B_dense = np.zeros((n_q_total, n_vert), dtype=np.float64)

    ir_arr, iz_arr, irow_arr, icol_arr = np.meshgrid(
        np.arange(Nr), np.arange(Nz),
        np.arange(n_quad), np.arange(n_quad), indexing="ij")
    ir_arr = ir_arr.flatten()
    iz_arr = iz_arr.flatten()
    irow_arr = irow_arr.flatten()
    icol_arr = icol_arr.flatten()
    q_idx = np.arange(n_q_total)

    v00 = ir_arr * n_vert_z + iz_arr
    v10 = (ir_arr + 1) * n_vert_z + iz_arr
    v11 = (ir_arr + 1) * n_vert_z + (iz_arr + 1)
    v01 = ir_arr * n_vert_z + (iz_arr + 1)

    B_dense[q_idx, v00] = N_00[irow_arr, icol_arr]
    B_dense[q_idx, v10] = N_10[irow_arr, icol_arr]
    B_dense[q_idx, v11] = N_11[irow_arr, icol_arr]
    B_dense[q_idx, v01] = N_01[irow_arr, icol_arr]
    if verbose:
        print(f"  Aggregator B in {time.time()-t0:.2f}s, shape "
              f"({n_q_total}, {n_vert})")

    # G_full
    if verbose:
        t0 = time.time()
    R_obs = r_flat[:, None]
    Z_obs = z_flat[:, None]
    R_src = r_flat[None, :]
    Z_src = z_flat[None, :]
    G_full = axisym_kernel_vec(R_obs, Z_obs, R_src, Z_src)
    if verbose:
        print(f"  G_full in {time.time()-t0:.2f}s, "
              f"size {G_full.nbytes/1e9:.2f} GB")

    if verbose:
        t0 = time.time()
    Gw = (G_full
          * (2.0 * math.pi * r_flat * jac_flat)[:, None]
          * jac_flat[None, :])

    if use_gpu:
        Gw_gpu = xp.asarray(Gw)
        B_gpu = xp.asarray(B_dense)
        Gw_B = Gw_gpu @ B_gpu
        K_full = (MU0 / 2.0) * (B_gpu.T @ Gw_B)
        K_full = xp.asnumpy(K_full)
    else:
        Gw_B = Gw @ B_dense
        K_full = (MU0 / 2.0) * (B_dense.T @ Gw_B)

    if verbose:
        print(f"  K = B^T Gw B in {time.time()-t0:.2f}s")

    K_full = 0.5 * (K_full + K_full.T)

    # M = B^T diag(2 pi r * jac) B
    weights_obs = 2.0 * math.pi * r_flat * jac_flat
    M_full = B_dense.T @ (weights_obs[:, None] * B_dense)

    # b_v = sum_q B[q, v] * A_q * (2pi r_q * jac_q) with A_q = (B0/2) r_q
    A_q = (B0 / 2.0) * r_flat
    b_full = B_dense.T @ (A_q * weights_obs)

    K_int = K_full[interior_idx[:, None], interior_idx[None, :]]
    M_int = M_full[interior_idx[:, None], interior_idx[None, :]]
    b_int = b_full[interior_idx]
    return K_int, M_int, b_int, n_vert, n_interior


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
        (48, 12, 3, "576 cells, n_quad=3"),
        (96, 24, 3, "2304 cells, n_quad=3"),
    ]

    out_results = {}
    for Nr, Nz, n_quad, label in cases:
        print(f"\n=== Cylinder axisym VIM linear basis ({label}) ===")
        t0 = time.time()
        K, M, b, n_vert, n_int = assemble_K_M_b_linear_cylinder(
            Nr, Nz, n_quad=n_quad, use_gpu=True, verbose=True)
        t_K = time.time() - t0
        print(f"  Total assembly: {t_K:.2f}s, K shape: {K.shape}")

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

        n_use = min(150, int(np.sum(g2 > 1e-32)))
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
            "Nr": Nr, "Nz": Nz, "n_quad": n_quad,
            "n_vert": n_vert, "n_interior": n_int,
            "K_assembly_s": t_K,
            "leading_foster_us": leading_foster,
            "Cauer_rungs_Kameari": rungs,
        }

    out_path = Path(__file__).parent / "cylinder_vim_axisym_linear_results.json"
    out_path.write_text(json.dumps(out_results, indent=2), encoding="utf-8")
    print(f"\nSaved: {out_path}")


if __name__ == "__main__":
    main()
