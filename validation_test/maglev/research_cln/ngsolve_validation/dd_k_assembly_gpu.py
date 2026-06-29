"""dd_k_assembly_gpu.py — GPU-vectorized DD K matrix assembly.

Production-scale version. Uses CuPy for DD basis evaluation and chunked
DD matmul to keep memory under 16 GB on Quadro RTX 5000.

Strategy:
  - Geometry precompute (mpmath sin/cos, omega, c-box, r/r', weights):
    on CPU, then transfer DD arrays to GPU
  - DD basis at (r, r') points: GPU vectorized (works since dd_basis uses
    elementwise (hi, lo) ops)
  - DD matmul K = Phi_L^T @ (w × Phi_R) per c-component, chunked over q-axis
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
from dd_k_assembly import precompute_geometry_dd


def dd_pairwise_sum_axis0(arr_hi, arr_lo, xp):
    """DD pairwise sum along axis 0.

    Pairs adjacent elements and sums in DD. Repeats until one element remains.
    More numerically stable than sequential accumulation.
    """
    n = arr_hi.shape[0]
    while n > 1:
        half = n // 2
        # Sum even pairs
        s_hi, s_lo = dd_add(arr_hi[:half], arr_lo[:half],
                             arr_hi[half:2 * half], arr_lo[half:2 * half])
        if 2 * half < n:
            # Append leftover
            new_hi = xp.concatenate([s_hi, arr_hi[2 * half:2 * half + 1]], axis=0)
            new_lo = xp.concatenate([s_lo, arr_lo[2 * half:2 * half + 1]], axis=0)
        else:
            new_hi = s_hi
            new_lo = s_lo
        arr_hi = new_hi
        arr_lo = new_lo
        n = arr_hi.shape[0]
    return arr_hi[0], arr_lo[0]


def dd_matmul_T_chunked(A_hi, A_lo, B_hi, B_lo, xp, chunk_size=2048,
                          verbose=False):
    """DD matrix multiplication C = A^T @ B, chunked over k-axis.

    A shape (n_k, n_i), B shape (n_k, n_j).
    Result C shape (n_i, n_j).

    Each chunk: outer-product per k, pairwise sum over chunk, add to C.
    """
    n_k, n_i = A_hi.shape
    _, n_j = B_hi.shape
    C_hi = xp.zeros((n_i, n_j))
    C_lo = xp.zeros((n_i, n_j))

    n_chunks = (n_k + chunk_size - 1) // chunk_size
    if verbose:
        print(f"    DD matmul: ({n_k}, {n_i}) @ ({n_k}, {n_j}) "
              f"-> ({n_i}, {n_j}), {n_chunks} chunks of {chunk_size}")

    for chunk_idx in range(n_chunks):
        k_start = chunk_idx * chunk_size
        k_end = min(k_start + chunk_size, n_k)

        # A_chunk: shape (chunk, n_i, 1)
        # B_chunk: shape (chunk, 1, n_j)
        # prod: shape (chunk, n_i, n_j)
        prod_hi, prod_lo = dd_mul(
            A_hi[k_start:k_end, :, None], A_lo[k_start:k_end, :, None],
            B_hi[k_start:k_end, None, :], B_lo[k_start:k_end, None, :])

        # Pairwise sum along chunk axis (axis 0)
        chunk_sum_hi, chunk_sum_lo = dd_pairwise_sum_axis0(prod_hi, prod_lo, xp)

        # Add to C
        C_hi, C_lo = dd_add(C_hi, C_lo, chunk_sum_hi, chunk_sum_lo)

    return C_hi, C_lo


def assemble_K_dd_gpu(p, a, b, c, n_th=12, n_ph=24, n_rh=12, n_c=6,
                       use_gpu=True, chunk_size=2048, verbose=True):
    """GPU DD K matrix assembly (production scale)."""
    if use_gpu:
        import cupy as xp
    else:
        xp = np

    n_dofs = 2 * p ** 3 + 3 * p ** 2
    n_q = n_th * n_ph * n_rh * n_c ** 3
    if verbose:
        print(f"\n=== DD K assembly (GPU) ===")
        print(f"  order={p}, n_dofs={n_dofs}")
        print(f"  quad: th={n_th}, ph={n_ph}, rh={n_rh}, c={n_c}, n_q={n_q}")
        mem_phi = n_q * n_dofs * 3 * 8 * 2 / 1e9
        print(f"  Phi memory (hi+lo, single side): {mem_phi:.2f} GB")

    if verbose:
        t0 = time.time()
    geom = precompute_geometry_dd(a, b, c, n_th, n_ph, n_rh, n_c)
    if verbose:
        print(f"  Geometry precompute (CPU mpmath): {time.time()-t0:.2f}s")

    # Transfer geometry to GPU
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

    if verbose:
        t0 = time.time()
    Phi_L_hi, Phi_L_lo = dd_evaluate_basis(r_x_hi, r_x_lo, r_y_hi, r_y_lo,
                                              r_z_hi, r_z_lo, p, xp_mod=xp)
    if use_gpu:
        xp.cuda.Stream.null.synchronize()
    if verbose:
        print(f"  DD Phi_L (basis at r): {time.time()-t0:.2f}s")
        t0 = time.time()
    Phi_R_hi, Phi_R_lo = dd_evaluate_basis(rp_x_hi, rp_x_lo, rp_y_hi, rp_y_lo,
                                              rp_z_hi, rp_z_lo, p, xp_mod=xp)
    if use_gpu:
        xp.cuda.Stream.null.synchronize()
    if verbose:
        print(f"  DD Phi_R (basis at r'): {time.time()-t0:.2f}s")

    # K[i,j] = sum_c diag_c × sum_q w_q × Phi_L[q,i,c] × Phi_R[q,j,c]
    # For each c: K_c = Phi_L[:,:,c]^T @ (w × Phi_R[:,:,c]) via DD matmul
    if verbose:
        t0 = time.time()
    K_hi = xp.zeros((n_dofs, n_dofs))
    K_lo = xp.zeros((n_dofs, n_dofs))
    diag_abc_sq = [a * a, b * b, c * c]

    for cc in range(3):
        if verbose:
            print(f"  c-component {cc}:")
        # wR = w × Phi_R[:, :, cc] in DD, shape (n_q, n_dofs)
        w_bcast_hi = w_hi[:, None] + 0 * Phi_R_hi[:, :, cc]
        w_bcast_lo = w_lo[:, None] + 0 * Phi_R_lo[:, :, cc]
        wR_hi, wR_lo = dd_mul(w_bcast_hi, w_bcast_lo,
                                Phi_R_hi[:, :, cc], Phi_R_lo[:, :, cc])

        # DD matmul K_c = Phi_L[:,:,cc]^T @ wR
        K_c_hi, K_c_lo = dd_matmul_T_chunked(
            Phi_L_hi[:, :, cc], Phi_L_lo[:, :, cc],
            wR_hi, wR_lo, xp, chunk_size=chunk_size, verbose=verbose)

        # K += diag_c × K_c (FP64 scaling; for a, b, c in meters this is exact)
        d = diag_abc_sq[cc]
        K_c_hi = K_c_hi * d
        K_c_lo = K_c_lo * d
        K_hi, K_lo = dd_add(K_hi, K_lo, K_c_hi, K_c_lo)
        if use_gpu:
            xp.cuda.Stream.null.synchronize()
    if verbose:
        print(f"  DD matmul total: {time.time()-t0:.2f}s")

    # Upper-triangle mirror (C++ convention)
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
    print("DD K matrix assembly GPU version — cuboid 5×2×1 order=3")
    print("=" * 72)

    sys.path.insert(0,
        str(Path("S:/Radia/01_GitHub/src/ext/radia_vim/src")))
    import radia_vim

    # First test: reduced quad to verify pipeline
    p = 3
    A_M, B_M, C_M = 5e-3, 2e-3, 1e-3
    n_th, n_ph, n_rh, n_c = 8, 12, 8, 4
    print(f"\n  Config: order={p}, quad=({n_th},{n_ph},{n_rh},{n_c})")

    t0 = time.time()
    K_hi, K_lo = assemble_K_dd_gpu(p, A_M, B_M, C_M, n_th, n_ph, n_rh, n_c,
                                     use_gpu=True, chunk_size=2048,
                                     verbose=True)
    t_dd = time.time() - t0
    print(f"\n  Total DD K assembly: {t_dd:.2f}s")

    # C++ FP64 reference
    basis = radia_vim.HDivDivFreeHexBasis(p)
    K_cpp = radia_vim.assemble_K_bare(basis, A_M, B_M, C_M,
                                       n_omega_theta=n_th, n_omega_phi=n_ph,
                                       n_rho=n_rh, n_c=n_c, verbose=0)

    diff = np.abs(K_hi - K_cpp)
    rel = diff / (np.abs(K_cpp) + 1e-30)
    significant = np.abs(K_cpp) > 1e-12
    max_rel_sig = float(np.max(np.where(significant, rel, 0)))
    print(f"\n  max |K_DD - K_CPP|        = {np.max(diff):.3e}")
    print(f"  max rel diff (|K|>1e-12) = {max_rel_sig:.3e}")
    print(f"  (Expected: rel ~ 1e-14 to 1e-15)")
    print(f"\n  DD lo part statistics:")
    print(f"    max |K_lo|         = {np.max(np.abs(K_lo)):.3e}")
    print(f"    max |K_lo/K_hi|    = "
          f"{float(np.max(np.abs(K_lo) / (np.abs(K_hi) + 1e-30))):.3e}")

    # Save for further pipeline (mpmath eigh + Cauer)
    np.savez(Path(__file__).parent / "dd_K_cuboid521_order3.npz",
             K_hi=K_hi, K_lo=K_lo, p=p, a=A_M, b=B_M, c=C_M,
             n_th=n_th, n_ph=n_ph, n_rh=n_rh, n_c=n_c)
    print(f"\n  Saved: dd_K_cuboid521_order3.npz")


if __name__ == "__main__":
    main()
