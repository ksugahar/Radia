"""sphere_vim_axisym_linear.py — Sphere axisym VIM with LINEAR (bilinear) J_phi basis.

DOFs: 1 per vertex of (rho, theta) grid.
Shape functions: bilinear in (xi, eta) on each cell.
Boundary conditions: J_phi = 0 on sphere surface (rho = R) and axis (rho = 0,
or theta = 0 / theta = pi where r = rho sin(theta) = 0). Eliminate boundary
vertices, eigenproblem on interior only.

Expected: O(h^4) convergence in eigenvalues vs O(h^2) for piecewise constant.
At Nrho=30, Ntheta=36 (1080 cells, 1147 vertices, ~1006 interior), expect
leading rung gap ~0.01-0.1% to Stoll, vs 0.88% for constant basis.
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


R_SPHERE = 10e-3
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


def assemble_K_M_b_linear_sphere(N_rho, N_theta, n_quad=2,
                                   use_gpu=True, verbose=True):
    """Bilinear J_phi basis on (rho, theta) sphere mesh."""
    if use_gpu:
        import cupy as xp
    else:
        xp = np

    # Vertices: (N_rho+1) x (N_theta+1) grid
    n_vert_rho = N_rho + 1
    n_vert_theta = N_theta + 1
    n_vert = n_vert_rho * n_vert_theta

    # Vertex coordinates
    rho_v = np.linspace(0, R_SPHERE, n_vert_rho)
    theta_v = np.linspace(0, math.pi, n_vert_theta)
    # Vertex indexing: v = i * n_vert_theta + j  (i = rho idx, j = theta idx)

    # Cells: N_rho x N_theta
    drho = R_SPHERE / N_rho
    dtheta = math.pi / N_theta

    # ----- Boundary identification -----
    # J_phi = 0 only where r = rho sin(theta) = 0, i.e. on the z-axis:
    #   rho = 0 (i = 0), or theta = 0 (j = 0), or theta = pi (j = N_theta).
    # At rho = R (sphere surface), J_phi is FREE — surface normal is radial,
    # J_phi is azimuthal, so J_n = J_rho = 0 automatically (we don't model J_rho).
    interior = np.zeros((n_vert_rho, n_vert_theta), dtype=bool)
    interior[1:, 1:-1] = True   # exclude i=0, j=0, j=N_theta. Include i=N_rho
    interior_flat = interior.flatten()
    interior_idx = np.where(interior_flat)[0]
    n_interior = len(interior_idx)

    if verbose:
        print(f"  N_rho={N_rho}, N_theta={N_theta}")
        print(f"  Total vertices: {n_vert}, interior: {n_interior}")
        print(f"  Cells: {N_rho * N_theta}")

    # ----- Quadrature on reference cell [0,1]^2 -----
    nodes_m1, w_m1 = np.polynomial.legendre.leggauss(n_quad)
    qp = 0.5 * (nodes_m1 + 1.0)        # quadrature points in [0, 1]
    qw = 0.5 * w_m1                     # weights summing to 1

    # All q-points: (N_rho, N_theta, n_quad, n_quad)
    # Local (xi, eta) for q-point (irow, icol) in cell (i_rho, j_theta):
    #   xi = qp[irow] (along rho direction)
    #   eta = qp[icol] (along theta direction)
    rho_q = (rho_v[:-1, None, None, None]
              + drho * qp[None, None, :, None])  # (N_rho, 1, n_quad, 1)
    theta_q = (theta_v[None, :-1, None, None]
                + dtheta * qp[None, None, None, :])  # (1, N_theta, 1, n_quad)
    rho_q = np.broadcast_to(rho_q, (N_rho, N_theta, n_quad, n_quad)).copy()
    theta_q = np.broadcast_to(theta_q, (N_rho, N_theta, n_quad, n_quad)).copy()
    # r = rho sin(theta), z = rho cos(theta)
    r_q_full = rho_q * np.sin(theta_q)
    z_q_full = rho_q * np.cos(theta_q)
    # Jacobian dr dz = rho drho dtheta -> include rho factor
    # Cell jacobian (from xi,eta in [0,1]^2 to physical drho dtheta) is drho * dtheta
    jac_q_full = (rho_q                                       # rho factor
                  * drho * dtheta
                  * qw[None, None, :, None] * qw[None, None, None, :])  # GL weights

    # Shape functions at quadrature points (xi=qp[irow], eta=qp[icol])
    # N_00 = (1 - xi) (1 - eta)
    # N_10 = xi (1 - eta)
    # N_11 = xi eta
    # N_01 = (1 - xi) eta
    # Indexing: irow=xi index, icol=eta index, vertex_local in {00, 10, 11, 01}
    xi = qp                                        # (n_quad,)
    eta = qp
    N_00 = np.outer(1 - xi, 1 - eta)               # (n_quad, n_quad) [irow, icol]
    N_10 = np.outer(xi,     1 - eta)
    N_11 = np.outer(xi,     eta)
    N_01 = np.outer(1 - xi, eta)

    # Flatten q-points: order (i_rho, j_theta, irow, icol)
    n_q_total = N_rho * N_theta * n_quad * n_quad
    r_flat = r_q_full.reshape(n_q_total)
    z_flat = z_q_full.reshape(n_q_total)
    jac_flat = jac_q_full.reshape(n_q_total)

    # ----- Aggregator B[q, v] -----
    # For q at (i_rho, j_theta, irow, icol), 4 contributing vertices:
    #   v_00 = i_rho * n_vert_theta + j_theta
    #   v_10 = (i_rho+1) * n_vert_theta + j_theta
    #   v_11 = (i_rho+1) * n_vert_theta + (j_theta+1)
    #   v_01 = i_rho * n_vert_theta + (j_theta+1)
    # B[q, v_00] = N_00[irow, icol], etc.
    if verbose:
        t0 = time.time()
    # Build sparse aggregator (dense for now since n_vert is moderate)
    # Use COO format then convert to CSR or dense
    B_dense = np.zeros((n_q_total, n_vert), dtype=np.float64)
    # Iterate over q-points and assign
    # Vectorized: build row, col, val arrays
    irho_arr, jtheta_arr, irow_arr, icol_arr = np.meshgrid(
        np.arange(N_rho), np.arange(N_theta),
        np.arange(n_quad), np.arange(n_quad), indexing="ij")
    irho_arr = irho_arr.flatten()
    jtheta_arr = jtheta_arr.flatten()
    irow_arr = irow_arr.flatten()
    icol_arr = icol_arr.flatten()

    q_idx = np.arange(n_q_total)

    # Vertex indices per q
    v00 = irho_arr * n_vert_theta + jtheta_arr
    v10 = (irho_arr + 1) * n_vert_theta + jtheta_arr
    v11 = (irho_arr + 1) * n_vert_theta + (jtheta_arr + 1)
    v01 = irho_arr * n_vert_theta + (jtheta_arr + 1)

    # Shape function values per q
    N00v = N_00[irow_arr, icol_arr]
    N10v = N_10[irow_arr, icol_arr]
    N11v = N_11[irow_arr, icol_arr]
    N01v = N_01[irow_arr, icol_arr]

    B_dense[q_idx, v00] = N00v
    # Note: cells share vertices, so multiple q-points assign to same (q, v)?
    # No — each q is in exactly ONE cell, so B[q, v00] for that q is N_00 from
    # that cell only. v00 of cell c is ALSO v11 of cell (c shifted up-left).
    # But for THIS q (in cell c), B[q, v_global] = shape function from cell c.
    # The vertex DOF is shared across cells through other q-points.
    B_dense[q_idx, v10] = N10v
    B_dense[q_idx, v11] = N11v
    B_dense[q_idx, v01] = N01v

    if verbose:
        print(f"  Aggregator B built in {time.time()-t0:.2f}s, "
              f"shape ({n_q_total}, {n_vert})")

    # ----- Compute G_full (dense kernel matrix on quadrature points) -----
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
    if verbose:
        t0 = time.time()
    Gw = (G_full
          * (2.0 * math.pi * r_flat * jac_flat)[:, None]
          * jac_flat[None, :])

    # K_full[v, w] = (mu0/2) * B[:, v]^T @ Gw @ B[:, w]
    # = (mu0/2) * (Gw @ B)[v, w] * B[v]... no wait
    # K = (mu0/2) * B^T @ Gw @ B  (n_vert, n_vert)
    if use_gpu:
        Gw_gpu = xp.asarray(Gw)
        B_gpu = xp.asarray(B_dense)
        # K = (mu0/2) B^T Gw B
        Gw_B = Gw_gpu @ B_gpu                  # (n_q, n_vert)
        K_full = (MU0 / 2.0) * (B_gpu.T @ Gw_B)  # (n_vert, n_vert)
        K_full = xp.asnumpy(K_full)
    else:
        Gw_B = Gw @ B_dense
        K_full = (MU0 / 2.0) * (B_dense.T @ Gw_B)

    if verbose:
        print(f"  K = B^T Gw B in {time.time()-t0:.2f}s")

    # K may have small asymmetry due to G non-symmetric formulation
    K_full = 0.5 * (K_full + K_full.T)

    # ----- Mass matrix: M[v, w] = ∫ N_v N_w * 2 pi r dV (sparse, 4x4 per cell) -----
    if verbose:
        t0 = time.time()
    M_full = np.zeros((n_vert, n_vert), dtype=np.float64)
    # M_local 4x4: ∫ N_a N_b * 2 pi r dr dz over cell
    # = ∫ N_a N_b * 2 pi (rho sin theta) * rho drho dtheta
    # For each cell, evaluate via quadrature using B values
    # M_local[a, b] = sum_q (2 pi r_q) * jac_q * N_a(q) * N_b(q)  (within cell)
    # Equivalently: B[:, v]^T diag(2 pi r * jac) B[:, w] — but only for
    # (v, w) sharing a cell.
    # Easiest: M = B^T diag(2 pi r * jac) B (DENSE n_vert x n_vert)
    weights_obs = 2.0 * math.pi * r_flat * jac_flat
    M_full = B_dense.T @ (weights_obs[:, None] * B_dense)  # (n_vert, n_vert)
    if verbose:
        print(f"  M = B^T diag(2 pi r jac) B in {time.time()-t0:.2f}s")

    # ----- b vector: b_v = ∫ A_imposed * N_v * 2 pi r dV -----
    # A_imposed_phi = (B0 / 2) * r
    # b_v = ∫ (B0/2) r * N_v * 2 pi r dr dz
    #     = pi B0 ∫ r^2 N_v dr dz
    # Computed via quadrature
    if verbose:
        t0 = time.time()
    A_q = (B0 / 2.0) * r_flat                          # A_imposed at q
    weights_src = 2.0 * math.pi * r_flat * jac_flat    # observation weight
    b_full = B_dense.T @ (A_q * weights_src)
    if verbose:
        print(f"  b vector in {time.time()-t0:.2f}s")

    # ----- Restrict to interior DOFs -----
    K_int = K_full[interior_idx[:, None], interior_idx[None, :]]
    M_int = M_full[interior_idx[:, None], interior_idx[None, :]]
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


def main():
    cases = [
        (45, 54, 3, "2430 cells, n_quad=3 (baseline)"),
        (60, 72, 3, "4320 cells, n_quad=3 (push)"),
    ]

    out_results = {}
    for N_rho, N_theta, n_quad, label in cases:
        print(f"\n=== Sphere axisym VIM linear basis ({label}) ===")
        t0 = time.time()
        K, M, b, interior_idx, n_vert = assemble_K_M_b_linear_sphere(
            N_rho, N_theta, n_quad=n_quad, use_gpu=True, verbose=True)
        t_K = time.time() - t0
        n_int = len(interior_idx)
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
            "N_rho": N_rho, "N_theta": N_theta, "n_interior": n_int,
            "K_assembly_s": t_K,
            "leading_foster_us": leading_foster,
            "Cauer_rungs_Kameari": rungs,
        }

    out_path = (Path(__file__).parent
                / "sphere_vim_axisym_linear_results.json")
    out_path.write_text(json.dumps(out_results, indent=2), encoding="utf-8")
    print(f"\nSaved: {out_path}")


if __name__ == "__main__":
    main()
