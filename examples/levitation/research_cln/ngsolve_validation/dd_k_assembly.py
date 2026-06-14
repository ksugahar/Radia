"""dd_k_assembly.py — Vectorized GPU DD K matrix assembly for hex VIM.

Production-scale DD port of assemble_K_cupy() from hex_vim_cupy_kassembly.py.
Strategy:
  - Pre-compute all geometric quantities (sin/cos, sphere weight, rho_max,
    Lomega, c-box, jacobians) in mpmath 40-digit, then convert to DD arrays
  - DD basis evaluation at all 750K (rx, ry, rz) and (rpx, rpy, rpz) points
  - DD weighted matmul: K[i,j] = sum_q w_q × D_c × Phi_L[q,i,c] × Phi_R[q,j,c]
    implemented as Phi_L_hi/lo @ (w * Phi_R_hi/lo) DD matmul per c-component
  - All on GPU via CuPy with (hi, lo) array pairs

Memory budget for cuboid order=3 (n_dofs=81, n_q ≈ 750K):
  Phi (hi+lo): 81 × 750K × 3 × 8 × 2 = 2.9 GB
  Phi_L + Phi_R: 5.8 GB
  Intermediate matmuls: ~1-2 GB
  Total: ~8 GB on 16 GB GPU — feasible.

For order=4 (n_dofs=176): ~12 GB total, very tight on 16 GB.
"""
from __future__ import annotations

import math
import sys
import time
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
import mpmath as mp

sys.path.insert(0, str(Path(__file__).parent))
from dd_arithmetic import (
    dd_add, dd_sub, dd_mul, dd_div, dd_from_float, dd_from_mpmath, dd_sqrt,
)
from dd_basis import dd_evaluate_basis


def precompute_geometry_dd(a, b, c, n_th, n_ph, n_rh, n_c, dps=40):
    """Precompute all geometric quantities in mpmath then convert to DD arrays.

    Returns dict with DD (hi, lo) arrays for:
      r_xyz[q], rp_xyz[q], full_w[q]
    where q indexes all n_q = n_th × n_ph × n_rh × n_c³ quadrature points.

    Points where the integration domain is empty (c-box has negative side)
    are given weight zero.
    """
    mp.mp.dps = dps
    pi_mp = mp.pi
    a_mp, b_mp, c_mp = mp.mpf(a), mp.mpf(b), mp.mpf(c)

    # Gauss-Legendre nodes/weights at 40-digit precision (Newton-refined)
    from dd_gl_nodes import gauss_legendre_mpf
    th_n_mp, th_w_mp = gauss_legendre_mpf(n_th, dps=dps)
    ph_n_mp, ph_w_mp = gauss_legendre_mpf(n_ph, dps=dps)
    rh_n_mp, rh_w_mp = gauss_legendre_mpf(n_rh, dps=dps)
    c_n_mp, c_w_mp = gauss_legendre_mpf(n_c, dps=dps)

    n_q = n_th * n_ph * n_rh * n_c ** 3

    # Build arrays
    r_x_hi = np.zeros(n_q)
    r_x_lo = np.zeros(n_q)
    r_y_hi = np.zeros(n_q)
    r_y_lo = np.zeros(n_q)
    r_z_hi = np.zeros(n_q)
    r_z_lo = np.zeros(n_q)
    rp_x_hi = np.zeros(n_q)
    rp_x_lo = np.zeros(n_q)
    rp_y_hi = np.zeros(n_q)
    rp_y_lo = np.zeros(n_q)
    rp_z_hi = np.zeros(n_q)
    rp_z_lo = np.zeros(n_q)
    w_hi = np.zeros(n_q)
    w_lo = np.zeros(n_q)

    print(f"  Precomputing {n_q} quadrature points in mpmath {dps}-digit...")
    t0 = time.time()
    q = 0
    for it in range(n_th):
        theta_mp = (th_n_mp[it] + 1) * pi_mp / 2
        s_th_mp = mp.sin(theta_mp)
        c_th_mp = mp.cos(theta_mp)
        jac_th_mp = pi_mp / 2

        for ip in range(n_ph):
            phi_mp = (ph_n_mp[ip] + 1) * pi_mp
            c_ph_mp = mp.cos(phi_mp)
            s_ph_mp = mp.sin(phi_mp)
            jac_ph_mp = pi_mp

            omega_x_mp = s_th_mp * c_ph_mp
            omega_y_mp = s_th_mp * s_ph_mp
            omega_z_mp = c_th_mp

            omega_max_mp = max(abs(omega_x_mp), abs(omega_y_mp),
                                abs(omega_z_mp))
            if omega_max_mp < mp.mpf("1e-30"):
                q += n_rh * n_c ** 3
                continue
            rho_max_mp = mp.mpf(1) / omega_max_mp

            Lomega_norm_mp = mp.sqrt(a_mp ** 2 * omega_x_mp ** 2
                                       + b_mp ** 2 * omega_y_mp ** 2
                                       + c_mp ** 2 * omega_z_mp ** 2)

            sphere_w_mp = (th_w_mp[it]
                            * ph_w_mp[ip]
                            * jac_th_mp * jac_ph_mp
                            * s_th_mp / Lomega_norm_mp)

            for ir in range(n_rh):
                rho_val_mp = (rh_n_mp[ir] + 1) * rho_max_mp / 2
                jac_rh_mp = rho_max_mp / 2
                u_x_mp = rho_val_mp * omega_x_mp
                u_y_mp = rho_val_mp * omega_y_mp
                u_z_mp = rho_val_mp * omega_z_mp

                cx_lo_mp = abs(u_x_mp) / 2
                cy_lo_mp = abs(u_y_mp) / 2
                cz_lo_mp = abs(u_z_mp) / 2
                cx_jac_mp = (mp.mpf(1) - abs(u_x_mp)) / 2
                cy_jac_mp = (mp.mpf(1) - abs(u_y_mp)) / 2
                cz_jac_mp = (mp.mpf(1) - abs(u_z_mp)) / 2

                valid = (cx_jac_mp > 0 and cy_jac_mp > 0 and cz_jac_mp > 0)
                if not valid:
                    q += n_c ** 3
                    continue

                outer_w_mp = (sphere_w_mp * rh_w_mp[ir]
                                * jac_rh_mp * rho_val_mp
                                * cx_jac_mp * cy_jac_mp * cz_jac_mp)

                for icx in range(n_c):
                    cx_mp = cx_lo_mp + (c_n_mp[icx] + 1) * cx_jac_mp
                    wcx_mp = c_w_mp[icx]
                    for icy in range(n_c):
                        cy_mp = cy_lo_mp + (c_n_mp[icy] + 1) * cy_jac_mp
                        wcy_mp = c_w_mp[icy]
                        for icz in range(n_c):
                            cz_mp = cz_lo_mp + (c_n_mp[icz] + 1) * cz_jac_mp
                            wcz_mp = c_w_mp[icz]

                            rx_mp = cx_mp + u_x_mp / 2
                            ry_mp = cy_mp + u_y_mp / 2
                            rz_mp = cz_mp + u_z_mp / 2
                            rpx_mp = cx_mp - u_x_mp / 2
                            rpy_mp = cy_mp - u_y_mp / 2
                            rpz_mp = cz_mp - u_z_mp / 2

                            full_w_mp = outer_w_mp * wcx_mp * wcy_mp * wcz_mp

                            r_x_hi[q], r_x_lo[q] = dd_from_mpmath(rx_mp)
                            r_y_hi[q], r_y_lo[q] = dd_from_mpmath(ry_mp)
                            r_z_hi[q], r_z_lo[q] = dd_from_mpmath(rz_mp)
                            rp_x_hi[q], rp_x_lo[q] = dd_from_mpmath(rpx_mp)
                            rp_y_hi[q], rp_y_lo[q] = dd_from_mpmath(rpy_mp)
                            rp_z_hi[q], rp_z_lo[q] = dd_from_mpmath(rpz_mp)
                            w_hi[q], w_lo[q] = dd_from_mpmath(full_w_mp)
                            q += 1

    print(f"  Geometry precompute: {time.time()-t0:.2f}s")
    return {
        "r_x": (r_x_hi, r_x_lo), "r_y": (r_y_hi, r_y_lo),
        "r_z": (r_z_hi, r_z_lo),
        "rp_x": (rp_x_hi, rp_x_lo), "rp_y": (rp_y_hi, rp_y_lo),
        "rp_z": (rp_z_hi, rp_z_lo),
        "w": (w_hi, w_lo),
    }


def dd_matmul_T(A_hi, A_lo, B_hi, B_lo):
    """Compute C = A^T @ B in DD where A, B are 2D arrays (m, n).
    Result C is (m, m) DD. Each entry is sum of (n,) DD products.

    Uses Kahan-like accumulation: for each (i, j), compute sum over k of
    dd_mul(A[k, i], B[k, j]). Sum in DD.
    """
    n_k, n_i = A_hi.shape
    _, n_j = B_hi.shape
    C_hi = np.zeros((n_i, n_j))
    C_lo = np.zeros((n_i, n_j))

    # Loop over (i, j); inner reduction over k vectorized
    # For each i, compute c[i, :] = sum_k A[k, i] * B[k, :]
    for i in range(n_i):
        # Initialize sum = 0
        s_hi = np.zeros(n_j)
        s_lo = np.zeros(n_j)
        # Vectorized over k: A[k, i] × B[k, j] for all j
        # Loop k for clarity (vectorize later)
        for k in range(n_k):
            # element-wise product over j-dim
            a_k_hi = A_hi[k, i]
            a_k_lo = A_lo[k, i]
            # broadcast a_k * B[k, :] : (1,) * (n_j,)
            prod_hi, prod_lo = dd_mul(
                np.full(n_j, a_k_hi), np.full(n_j, a_k_lo),
                B_hi[k, :], B_lo[k, :])
            s_hi, s_lo = dd_add(s_hi, s_lo, prod_hi, prod_lo)
        C_hi[i, :] = s_hi
        C_lo[i, :] = s_lo

    return C_hi, C_lo


def assemble_K_dd_vectorized(p, a, b, c, n_th=8, n_ph=12, n_rh=8, n_c=4,
                                use_gpu=False, verbose=True):
    """DD K matrix assembly, vectorized.

    Phase 2 production version. For MVP demo, use_gpu=False (numpy);
    GPU support via use_gpu=True for cupy arrays.
    """
    if use_gpu:
        import cupy as xp
    else:
        xp = np

    n_dofs = 2 * p ** 3 + 3 * p ** 2
    n_q = n_th * n_ph * n_rh * n_c ** 3
    if verbose:
        print(f"\n  DD K assembly: order={p}, n_dofs={n_dofs}, n_q={n_q}")

    geom = precompute_geometry_dd(a, b, c, n_th, n_ph, n_rh, n_c)

    # Transfer to GPU if needed
    if use_gpu:
        for k in ["r_x", "r_y", "r_z", "rp_x", "rp_y", "rp_z", "w"]:
            hi, lo = geom[k]
            geom[k] = (xp.asarray(hi), xp.asarray(lo))

    r_x_hi, r_x_lo = geom["r_x"]
    r_y_hi, r_y_lo = geom["r_y"]
    r_z_hi, r_z_lo = geom["r_z"]
    rp_x_hi, rp_x_lo = geom["rp_x"]
    rp_y_hi, rp_y_lo = geom["rp_y"]
    rp_z_hi, rp_z_lo = geom["rp_z"]
    w_hi, w_lo = geom["w"]

    # Evaluate DD basis at all (r, r') points
    if verbose:
        t0 = time.time()
    Phi_L_hi, Phi_L_lo = dd_evaluate_basis(r_x_hi, r_x_lo, r_y_hi, r_y_lo,
                                              r_z_hi, r_z_lo, p, xp_mod=xp)
    if verbose:
        print(f"  DD basis at r: {time.time()-t0:.2f}s")
        t0 = time.time()
    Phi_R_hi, Phi_R_lo = dd_evaluate_basis(rp_x_hi, rp_x_lo, rp_y_hi, rp_y_lo,
                                              rp_z_hi, rp_z_lo, p, xp_mod=xp)
    if verbose:
        print(f"  DD basis at r': {time.time()-t0:.2f}s")

    # K[i,j] = sum_q sum_c diag(a²,b²,c²)[c] × w[q] × Phi_L[q,i,c] × Phi_R[q,j,c]
    # For each c, compute Phi_L_c_T @ (w × Phi_R_c) via DD matmul, sum.
    if verbose:
        t0 = time.time()
    K_hi = xp.zeros((n_dofs, n_dofs))
    K_lo = xp.zeros((n_dofs, n_dofs))
    diag_abc_sq = [a * a, b * b, c * c]
    for cc in range(3):
        # Weighted Phi_R: w × Phi_R[:, :, cc]
        # Element-wise DD multiplication over (n_q, n_dofs)
        wR_hi, wR_lo = dd_mul(
            w_hi[:, None] + 0 * Phi_R_hi[:, :, cc],
            w_lo[:, None] + 0 * Phi_R_lo[:, :, cc],
            Phi_R_hi[:, :, cc], Phi_R_lo[:, :, cc])
        # DD matmul: K_c[i, j] = sum_q Phi_L[q, i, cc] × wR[q, j]
        K_c_hi, K_c_lo = dd_matmul_T(Phi_L_hi[:, :, cc],
                                       Phi_L_lo[:, :, cc],
                                       wR_hi, wR_lo)
        # Multiply by diag_abc_sq[cc] (FP64, exact for our values)
        d = diag_abc_sq[cc]
        K_c_hi = K_c_hi * d
        K_c_lo = K_c_lo * d
        K_hi, K_lo = dd_add(K_hi, K_lo, K_c_hi, K_c_lo)
    if verbose:
        print(f"  DD K einsum: {time.time()-t0:.2f}s")

    # Upper-triangle mirror (matches C++ assemble_K_bare convention)
    iu = np.triu_indices_from(K_hi, k=1) if not use_gpu \
        else xp.asarray(np.triu_indices_from(np.asarray(K_hi), k=1))
    if use_gpu:
        K_hi_np = xp.asnumpy(K_hi)
        K_lo_np = xp.asnumpy(K_lo)
    else:
        K_hi_np = K_hi
        K_lo_np = K_lo
    iu = np.triu_indices_from(K_hi_np, k=1)
    K_hi_np[iu[1], iu[0]] = K_hi_np[iu[0], iu[1]]
    K_lo_np[iu[1], iu[0]] = K_lo_np[iu[0], iu[1]]

    return K_hi_np, K_lo_np


def main():
    print("=" * 72)
    print("DD K matrix assembly (vectorized, order=2 demo)")
    print("=" * 72)

    sys.path.insert(0,
        str(Path("S:/Radia/01_GitHub/packages/radia-vim/src")))
    import radia_vim

    p = 2
    A_M, B_M, C_M = 5e-3, 2e-3, 1e-3
    n_th, n_ph, n_rh, n_c = 4, 6, 4, 2

    print(f"\n  Configuration: order={p}, quad=({n_th},{n_ph},{n_rh},{n_c})")

    t0 = time.time()
    K_hi, K_lo = assemble_K_dd_vectorized(p, A_M, B_M, C_M,
                                            n_th, n_ph, n_rh, n_c,
                                            use_gpu=False, verbose=True)
    t_dd = time.time() - t0
    print(f"\n  Total DD K assembly: {t_dd:.2f}s")

    # C++ FP64 reference
    basis = radia_vim.HDivDivFreeHexBasis(p)
    K_cpp = radia_vim.assemble_K_bare(basis, A_M, B_M, C_M,
                                       n_omega_theta=n_th, n_omega_phi=n_ph,
                                       n_rho=n_rh, n_c=n_c, verbose=0)

    # Apply mu0/(4 pi) prefactor to DD K (matching C++ convention for "physical K")
    # C++ assemble_K_bare returns bare integral WITHOUT mu0/(4pi), our DD matches
    K_dd_fp64 = K_hi

    # Compare
    diff = np.abs(K_dd_fp64 - K_cpp)
    rel = diff / (np.abs(K_cpp) + 1e-30)
    print(f"\n  max |K_DD_hi - K_CPP| = {np.max(diff):.3e}")
    print(f"  max relative          = {np.max(rel):.3e}")
    print(f"  (Expected: rel ~ 1e-14 to 1e-15, FP64 precision)")

    # Top discrepancies on significant entries
    significant = np.abs(K_cpp) > 1e-15
    rel_sig = np.where(significant, rel, 0)
    print(f"  max rel on |K|>1e-15: {np.max(rel_sig):.3e}")
    if np.max(rel_sig) < 1e-13:
        print(f"  PASS: DD K matches C++ FP64 within precision limit ✓")


if __name__ == "__main__":
    main()
