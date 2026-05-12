"""dd_M_b_assembly.py — DD mass matrix M and excitation vector b for hex VIM.

Production-scale DD port of the M assembly + b assembly in
hex_vim_gpu_order_sweep.py.

Mass matrix (HDiv Piola, diag Jacobian J = diag(a, b, c)):
  M_phys[i, j] = ∫ phi_phys_i · phi_phys_j dV_phys
              = (a/(bc)) ∫ φ_x_i × φ_x_j dV_ref
                + (b/(ac)) ∫ φ_y_i × φ_y_j dV_ref
                + (c/(ab)) ∫ φ_z_i × φ_z_j dV_ref

Excitation vector for uniform B_z applied (A_ext = (-y/2, x/2, 0) × B_0):
  b_i = sum_q W[q] × (Ax_phys × a × Phi_ref[q, i, 0]
                    + Ay_phys × b × Phi_ref[q, i, 1])
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
    dd_add, dd_sub, dd_mul, dd_div, dd_from_float, dd_from_mpmath,
)
from dd_basis import dd_evaluate_basis


B0 = 1.0


def assemble_M_b_dd_gpu(p, a, b, c, n_gauss=8, use_gpu=False, verbose=True):
    """DD M matrix and b vector assembly.

    Returns:
      M_hi, M_lo: (n_dofs, n_dofs) DD mass matrix
      b_hi, b_lo: (n_dofs,) DD excitation vector
    """
    if use_gpu:
        import cupy as xp
    else:
        xp = np

    mp.mp.dps = 40

    # Gauss-Legendre nodes/weights on [0, 1] (from [-1, 1] mapping)
    nodes_m1, w_m1 = np.polynomial.legendre.leggauss(n_gauss)
    nodes_f = 0.5 * (nodes_m1 + 1.0)
    weights_f = 0.5 * w_m1

    # Refine via mpmath for DD precision (or just use FP64 — error is small)
    nodes_hi = nodes_f.astype(np.float64)
    nodes_lo = np.zeros(n_gauss, dtype=np.float64)
    weights_hi = weights_f.astype(np.float64)
    weights_lo = np.zeros(n_gauss, dtype=np.float64)

    # Tensor-product grid for (xi, eta, zeta) on [0,1]^3
    XI, ETA, ZETA = np.meshgrid(nodes_hi, nodes_hi, nodes_hi, indexing="ij")
    xi_flat_hi = XI.flatten()
    eta_flat_hi = ETA.flatten()
    zeta_flat_hi = ZETA.flatten()
    xi_flat_lo = np.zeros_like(xi_flat_hi)
    eta_flat_lo = np.zeros_like(xi_flat_hi)
    zeta_flat_lo = np.zeros_like(xi_flat_hi)
    n_q = xi_flat_hi.size

    W3 = (weights_hi[:, None, None] * weights_hi[None, :, None]
          * weights_hi[None, None, :]).flatten()

    n_dofs = 2 * p ** 3 + 3 * p ** 2

    if verbose:
        print(f"  DD M, b assembly: order={p}, n_dofs={n_dofs}, n_q={n_q}")
        t0 = time.time()

    # Evaluate DD basis at reference quad points
    if use_gpu:
        xi_arr_hi = xp.asarray(xi_flat_hi)
        xi_arr_lo = xp.asarray(xi_flat_lo)
        eta_arr_hi = xp.asarray(eta_flat_hi)
        eta_arr_lo = xp.asarray(eta_flat_lo)
        zeta_arr_hi = xp.asarray(zeta_flat_hi)
        zeta_arr_lo = xp.asarray(zeta_flat_lo)
    else:
        xi_arr_hi = xi_flat_hi
        xi_arr_lo = xi_flat_lo
        eta_arr_hi = eta_flat_hi
        eta_arr_lo = eta_flat_lo
        zeta_arr_hi = zeta_flat_hi
        zeta_arr_lo = zeta_flat_lo

    Phi_hi, Phi_lo = dd_evaluate_basis(xi_arr_hi, xi_arr_lo,
                                          eta_arr_hi, eta_arr_lo,
                                          zeta_arr_hi, zeta_arr_lo,
                                          p, xp_mod=xp)
    if use_gpu:
        xp.cuda.Stream.null.synchronize()
    if verbose:
        print(f"  DD Phi at reference grid: {time.time()-t0:.2f}s")

    # ----- M matrix -----
    # M_ij = (a/(bc)) sum_q W[q] × Phi[q,i,0] × Phi[q,j,0]
    #      + (b/(ac)) sum_q W[q] × Phi[q,i,1] × Phi[q,j,1]
    #      + (c/(ab)) sum_q W[q] × Phi[q,i,2] × Phi[q,j,2]
    # Vectorized: for each c, weight Phi[:,:,c] by sqrt(W[q] × coef_c),
    # then matmul.
    diag_M = [a / (b * c), b / (a * c), c / (a * b)]

    if verbose:
        t0 = time.time()
    W_arr = xp.asarray(W3) if use_gpu else W3
    M_hi = xp.zeros((n_dofs, n_dofs))
    M_lo = xp.zeros((n_dofs, n_dofs))
    for cc in range(3):
        d = diag_M[cc]
        # Phi_w[q, i] = sqrt(W[q] × d) × Phi[q, i, cc]
        # Actually simpler: M_c = Phi^T @ (W × Phi) per component
        # Then sum with d coefficients
        w_phi_hi, w_phi_lo = dd_mul(
            (W_arr[:, None] + 0 * Phi_hi[:, :, cc]),
            xp.zeros_like(Phi_hi[:, :, cc]),
            Phi_hi[:, :, cc], Phi_lo[:, :, cc])
        # Now M_c = Phi[:,:,cc]^T @ w_phi : (n_dofs, n_dofs)
        # Use DD matmul (small: n_q=512, n_dofs=81)
        from dd_k_assembly_gpu import dd_matmul_T_chunked
        M_c_hi, M_c_lo = dd_matmul_T_chunked(
            Phi_hi[:, :, cc], Phi_lo[:, :, cc],
            w_phi_hi, w_phi_lo, xp, chunk_size=128)
        # Scale by d
        M_c_hi = M_c_hi * d
        M_c_lo = M_c_lo * d
        M_hi, M_lo = dd_add(M_hi, M_lo, M_c_hi, M_c_lo)
    if verbose:
        print(f"  DD M matrix: {time.time()-t0:.2f}s")

    # ----- b vector -----
    # A_ext_phys = B0 × (-y_phys/2, x_phys/2, 0) on centered cube
    # x_phys = a × (xi - 0.5), y_phys = b × (eta - 0.5)
    # b_i = sum_q W[q] × (Ax × a × Phi[q,i,0] + Ay × b × Phi[q,i,1])
    if verbose:
        t0 = time.time()
    # Compute A_ext at quad points in DD
    # xi - 0.5 (DD): use dd_sub
    xi_minus_half_hi, xi_minus_half_lo = dd_sub(
        xi_arr_hi, xi_arr_lo,
        xp.zeros_like(xi_arr_hi) + 0.5, xp.zeros_like(xi_arr_hi))
    eta_minus_half_hi, eta_minus_half_lo = dd_sub(
        eta_arr_hi, eta_arr_lo,
        xp.zeros_like(eta_arr_hi) + 0.5, xp.zeros_like(eta_arr_hi))
    # x_phys, y_phys (DD)
    a_arr_hi = xp.zeros_like(xi_arr_hi) + float(a)
    a_arr_lo = xp.zeros_like(xi_arr_hi)
    b_arr_hi = xp.zeros_like(xi_arr_hi) + float(b)
    b_arr_lo = xp.zeros_like(xi_arr_hi)
    x_phys_hi, x_phys_lo = dd_mul(a_arr_hi, a_arr_lo,
                                    xi_minus_half_hi, xi_minus_half_lo)
    y_phys_hi, y_phys_lo = dd_mul(b_arr_hi, b_arr_lo,
                                    eta_minus_half_hi, eta_minus_half_lo)
    # Ax = -B0/2 × y_phys; Ay = B0/2 × x_phys
    half_B0_hi = xp.zeros_like(xi_arr_hi) + (B0 / 2.0)
    half_B0_lo = xp.zeros_like(xi_arr_hi)
    Ay_hi, Ay_lo = dd_mul(half_B0_hi, half_B0_lo, x_phys_hi, x_phys_lo)
    Ax_neg_hi, Ax_neg_lo = dd_mul(half_B0_hi, half_B0_lo,
                                    y_phys_hi, y_phys_lo)
    Ax_hi = -Ax_neg_hi
    Ax_lo = -Ax_neg_lo

    # wax = W × Ax × a, way = W × Ay × b
    # Build wax_q for each q
    wax_temp_hi, wax_temp_lo = dd_mul(W_arr, xp.zeros_like(W_arr),
                                        Ax_hi, Ax_lo)
    wax_hi, wax_lo = dd_mul(wax_temp_hi, wax_temp_lo,
                              a_arr_hi, a_arr_lo)
    way_temp_hi, way_temp_lo = dd_mul(W_arr, xp.zeros_like(W_arr),
                                        Ay_hi, Ay_lo)
    way_hi, way_lo = dd_mul(way_temp_hi, way_temp_lo,
                              b_arr_hi, b_arr_lo)

    # b[i] = sum_q (wax[q] × Phi[q,i,0] + way[q] × Phi[q,i,1])
    # DD-sum over q axis
    # First compute prod[q, i] for c=0 and c=1, sum over q
    b_hi = xp.zeros(n_dofs)
    b_lo = xp.zeros(n_dofs)
    for cc, w_c_hi, w_c_lo in [(0, wax_hi, wax_lo), (1, way_hi, way_lo)]:
        # prod = w_c × Phi[:,:,cc] : (n_q, n_dofs)
        prod_hi, prod_lo = dd_mul(
            w_c_hi[:, None] + 0 * Phi_hi[:, :, cc],
            xp.zeros_like(Phi_hi[:, :, cc]),
            Phi_hi[:, :, cc], Phi_lo[:, :, cc])
        # Sum over q axis (DD pairwise sum)
        from dd_k_assembly_gpu import dd_pairwise_sum_axis0
        s_hi, s_lo = dd_pairwise_sum_axis0(prod_hi, prod_lo, xp)
        b_hi, b_lo = dd_add(b_hi, b_lo, s_hi, s_lo)
    if verbose:
        print(f"  DD b vector: {time.time()-t0:.2f}s")

    if use_gpu:
        M_hi_np = xp.asnumpy(M_hi)
        M_lo_np = xp.asnumpy(M_lo)
        b_hi_np = xp.asnumpy(b_hi)
        b_lo_np = xp.asnumpy(b_lo)
    else:
        M_hi_np, M_lo_np = M_hi, M_lo
        b_hi_np, b_lo_np = b_hi, b_lo

    return M_hi_np, M_lo_np, b_hi_np, b_lo_np


def main():
    print("=" * 72)
    print("DD M, b assembly — cuboid 5×2×1 order=3")
    print("=" * 72)

    sys.path.insert(0,
        str(Path("S:/Radia/01_GitHub/packages/radia-vim/src")))
    import radia_vim
    sys.path.insert(0, str(Path(__file__).parent))
    # Legacy FP64 reference (in legacy_fp64/ since 2026-05-12)
    sys.path.insert(0, str(Path(__file__).parent / "legacy_fp64"))
    from hex_vim_gpu_order_sweep import (
        assemble_mass_matrix_python, assemble_K_cupy, B0,
    )
    from hex_vim_cupy import evaluate_basis

    p = 3
    A_M, B_M, C_M = 5e-3, 2e-3, 1e-3
    basis = radia_vim.HDivDivFreeHexBasis(p)

    # DD assembly
    t0 = time.time()
    M_hi, M_lo, b_hi, b_lo = assemble_M_b_dd_gpu(
        p, A_M, B_M, C_M, n_gauss=8, use_gpu=True, verbose=True)
    print(f"  Total DD M+b: {time.time()-t0:.2f}s")

    # FP64 reference (existing code)
    M_fp64 = assemble_mass_matrix_python(basis, A_M, B_M, C_M, n_gauss=8)

    # Build b in FP64 reference
    n_gauss = 8
    nodes_m1, w_m1 = np.polynomial.legendre.leggauss(n_gauss)
    nodes = 0.5 * (nodes_m1 + 1.0)
    weights = 0.5 * w_m1
    XI, ETA, ZETA = np.meshgrid(nodes, nodes, nodes, indexing="ij")
    xi_f = XI.flatten()
    eta_f = ETA.flatten()
    zeta_f = ZETA.flatten()
    Phi_ref = evaluate_basis(xi_f, eta_f, zeta_f, p, xp_mod=np)
    x_phys = A_M * (xi_f - 0.5)
    y_phys = B_M * (eta_f - 0.5)
    Ax_phys = B0 * (-y_phys / 2)
    Ay_phys = B0 * (x_phys / 2)
    W3 = (weights[:, None, None] * weights[None, :, None]
          * weights[None, None, :]).flatten()
    b_fp64 = (Phi_ref[..., 0].T @ (W3 * Ax_phys * A_M)) \
             + (Phi_ref[..., 1].T @ (W3 * Ay_phys * B_M))

    # Compare
    M_diff = np.abs(M_hi - M_fp64)
    M_sig = np.abs(M_fp64) > 1e-15
    M_rel = M_diff / (np.abs(M_fp64) + 1e-30)
    M_max_rel = float(np.max(np.where(M_sig, M_rel, 0)))
    print(f"\n  M_hi vs M_FP64:")
    print(f"    max abs diff = {np.max(M_diff):.3e}")
    print(f"    max rel sig  = {M_max_rel:.3e}")
    print(f"    max |M_lo|   = {np.max(np.abs(M_lo)):.3e}")

    b_diff = np.abs(b_hi - b_fp64)
    b_sig = np.abs(b_fp64) > 1e-15
    b_rel = b_diff / (np.abs(b_fp64) + 1e-30)
    b_max_rel = float(np.max(np.where(b_sig, b_rel, 0)))
    print(f"\n  b_hi vs b_FP64:")
    print(f"    max abs diff = {np.max(b_diff):.3e}")
    print(f"    max rel sig  = {b_max_rel:.3e}")
    print(f"    max |b_lo|   = {np.max(np.abs(b_lo)):.3e}")

    np.savez(Path(__file__).parent / "dd_Mb_cuboid521_order3.npz",
             M_hi=M_hi, M_lo=M_lo, b_hi=b_hi, b_lo=b_lo,
             p=p, a=A_M, b=B_M, c=C_M)
    print(f"\n  Saved: dd_Mb_cuboid521_order3.npz")


if __name__ == "__main__":
    main()
