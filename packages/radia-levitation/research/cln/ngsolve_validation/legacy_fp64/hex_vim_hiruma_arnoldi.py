"""hex_vim_hiruma_arnoldi.py — Hiruma 3-term Lanczos + Arnoldi reorthogonalization
for hex VIM Cauer extraction.

Inspired by Nagamine's FreeFEM++ rectangle/arnoldi.edp which extends the
direct Sugahara-Kameari Cauer iteration with explicit Gram-Schmidt
reorthogonalization against ALL previous stage vectors. The key claim:
Arnoldi reorthog maintains Krylov M-orthogonality at high stages, extending
the reliable stage count beyond what plain Lanczos achieves.

Hiruma 3-term Lanczos (no reorthog) applied to a generic (K, M, b) system:

    w[1] = K^{-1} b
    lam[1] = w[1]^T K w[1]       # = 1/R_0
    w[2] = -w[1] / lam[1]
    for n = 1..N:
        lam[2n]   = w[2n]^T M w[2n]      # = L_{2n-1}
        Cw  = M w[2n]
        Aw  = -K^{-1} Cw
        w[2n+1] = w[2n-1] - Aw/lam[2n]
        lam[2n+1] = w[2n+1]^T K w[2n+1]  # = 1/R_{2n}
        w[2n+2] = w[2n] - w[2n+1]/lam[2n+1]

Arnoldi extension: after computing w[2n+1] and w[2n+2], explicitly
reorthogonalize against ALL previous stage vectors with the appropriate
inner product (K for odd vectors, M for even vectors).

Outputs Cauer rungs in Hiruma (Krylov-Padé impedance) convention.
"""
from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
import scipy.linalg as la

sys.path.insert(0, str(Path(__file__).parent))
from hex_vim_cupy_kassembly import assemble_K_cupy
from hex_vim_cupy import evaluate_basis

sys.path.insert(0, str(Path("S:/Radia/01_GitHub/packages/radia-vim/src")))
import radia_vim


# Cuboid 5x2x1 mm Cu
A_M = 5e-3
B_M = 2e-3
C_M = 1e-3
SIGMA = 5.8e7
MU0 = 4 * math.pi * 1e-7
B0 = 1.0


def assemble_mass_matrix_python(basis, a, b, c, n_gauss=8):
    n = basis.n_dofs
    nodes_m1, weights_m1 = np.polynomial.legendre.leggauss(n_gauss)
    nodes = 0.5 * (nodes_m1 + 1.0)
    weights = 0.5 * weights_m1

    Phi = np.zeros((n, 3, n_gauss, n_gauss, n_gauss))
    for ix, x in enumerate(nodes):
        for iy, y in enumerate(nodes):
            for iz, z in enumerate(nodes):
                for d in range(n):
                    v = basis.evaluate(d, float(x), float(y), float(z))
                    Phi[d, 0, ix, iy, iz] = v[0]
                    Phi[d, 1, ix, iy, iz] = v[1]
                    Phi[d, 2, ix, iy, iz] = v[2]

    W3 = weights[:, None, None] * weights[None, :, None] * weights[None, None, :]
    Mref = np.zeros((3, n, n))
    for comp in range(3):
        Pcomp = Phi[:, comp]
        Pw = Pcomp * W3[None, :, :, :]
        Mref[comp] = np.tensordot(Pcomp.reshape(n, -1),
                                  Pw.reshape(n, -1), axes=([1], [1]))

    Mphys = (a / (b * c)) * Mref[0] + (b / (a * c)) * Mref[1] + (c / (a * b)) * Mref[2]
    return Mphys


def assemble_b_vector(basis, a, b, c, p, n_gauss=8):
    """b_i = sum_q W3[q] * (Ax_phys * a * Phi_ref[q,i,0] + Ay_phys * b * Phi_ref[q,i,1])"""
    nodes_m1, w_m1 = np.polynomial.legendre.leggauss(n_gauss)
    nodes = 0.5 * (nodes_m1 + 1.0)
    weights = 0.5 * w_m1
    XI, ETA, ZETA = np.meshgrid(nodes, nodes, nodes, indexing="ij")
    xi_flat = XI.flatten()
    eta_flat = ETA.flatten()
    zeta_flat = ZETA.flatten()
    Phi_ref = evaluate_basis(xi_flat, eta_flat, zeta_flat, p, xp_mod=np)
    x_phys = a * (xi_flat - 0.5)
    y_phys = b * (eta_flat - 0.5)
    Ax_phys = B0 * (-y_phys / 2)
    Ay_phys = B0 * (x_phys / 2)
    W3 = (weights[:, None, None] * weights[None, :, None]
          * weights[None, None, :]).flatten()
    wax = W3 * Ax_phys * a
    way = W3 * Ay_phys * b
    return (Phi_ref[..., 0].T @ wax) + (Phi_ref[..., 1].T @ way)


def hiruma_3term_lanczos(K, M, b, N_stages, use_arnoldi=False, verbose=True):
    """Hiruma 3-term Lanczos with optional Arnoldi reorthogonalization.

    Returns: list of stage dicts with R_2k, L_2k+1, tau_pair_us, orth_drift.
    """
    K_lu = la.lu_factor(K)

    n_dofs = len(b)
    w = [None] * (2 * N_stages + 3)
    lam = [None] * (2 * N_stages + 3)

    w[1] = la.lu_solve(K_lu, b)
    lam[1] = float(w[1] @ K @ w[1])         # = 1/R_0
    w[2] = -w[1] / lam[1]

    R0 = 1.0 / lam[1]
    if verbose:
        prefix = "[Arnoldi] " if use_arnoldi else "[Lanczos] "
        print(f"  {prefix}R_0 = {R0:.4e}")
        print(f"  {prefix}{'k':>3} {'R_2k':>13} {'L_2k+1':>13} "
              f"{'tau_pair us':>13} {'orth drift':>12}")

    out = []
    for n in range(1, N_stages + 1):
        # lam[2n] = w[2n]^T M w[2n] = L_{2n-1}
        lam[2 * n] = float(w[2 * n] @ M @ w[2 * n])

        Cw = M @ w[2 * n]
        Aw = -la.lu_solve(K_lu, Cw)
        w[2 * n + 1] = w[2 * n - 1].copy() - Aw / lam[2 * n]

        # ----- Arnoldi reorthogonalization (odd vectors, K-inner product) -----
        if use_arnoldi:
            for k in range(1, n):
                prod = float(w[2 * n + 1] @ K @ w[2 * k + 1])
                # Subtract projection: w_new -= w_k * (prod / lam[2k+1])
                w[2 * n + 1] -= w[2 * k + 1] * prod / lam[2 * k + 1]

        lam[2 * n + 1] = float(w[2 * n + 1] @ K @ w[2 * n + 1])

        # Compute orthogonality drift (max cross-inner-product w/ previous odd vecs)
        drift = 0.0
        for kk in range(1, n):
            cross = float(w[2 * kk + 1] @ K @ w[2 * n + 1])
            drift = max(drift, abs(cross) / max(abs(lam[2 * n + 1]), 1e-30))

        # Cauer rung values
        R_2k = 1.0 / lam[2 * n - 1]
        L_2k_plus_1 = lam[2 * n]
        tau_pair_us = (L_2k_plus_1 / R_2k * 1e6
                       if R_2k > 0 and L_2k_plus_1 > 0 else None)

        w[2 * n + 2] = w[2 * n].copy() - w[2 * n + 1] / lam[2 * n + 1]

        # Arnoldi reorthog for even vectors (M-inner product)
        if use_arnoldi:
            for k in range(1, n + 1):
                prod = float(w[2 * n + 2] @ M @ w[2 * k])
                w[2 * n + 2] -= w[2 * k] * prod / lam[2 * k]

        ts = f"{tau_pair_us:.4f}" if tau_pair_us else "N/A"
        if verbose:
            k_idx = n - 1
            print(f"  {prefix}{k_idx:>3} {R_2k:>13.6e} "
                  f"{L_2k_plus_1:>13.6e} {ts:>13} {drift:>12.3e}")
        out.append({"k": n - 1, "R_2k": R_2k, "L_2k_plus_1": L_2k_plus_1,
                    "tau_pair_us": tau_pair_us, "orth_drift": drift})

    return out, R0


def main():
    p = 4
    n_th, n_ph, n_rh, n_c = 12, 24, 12, 6
    print("=" * 72)
    print(f"Hex VIM cuboid 5x2x1 order={p}: Hiruma vs Arnoldi 3-term Lanczos")
    print("=" * 72)

    # Build K, M, b
    print(f"\nAssembling K, M, b for order={p}...")
    t0 = time.time()
    K_bare = assemble_K_cupy(p, A_M, B_M, C_M, n_th, n_ph, n_rh, n_c,
                              use_gpu=True, verbose=False)
    K_phys = K_bare * (MU0 / (4 * math.pi))
    print(f"  K assembled in {time.time()-t0:.2f}s, shape {K_phys.shape}")

    basis = radia_vim.HDivDivFreeHexBasis(p)
    M = assemble_mass_matrix_python(basis, A_M, B_M, C_M, n_gauss=8)
    b = assemble_b_vector(basis, A_M, B_M, C_M, p, n_gauss=8)
    print(f"  M, b assembled")

    # axifemm convention: K (no sigma), M (with sigma), b (with sigma).
    # tau in seconds directly from lam (no extra sigma factor needed).
    M_sigma = SIGMA * M
    b_sigma = SIGMA * b

    N_stages = 10

    # Plain Lanczos
    print(f"\n--- Plain Hiruma 3-term Lanczos (no reorthog) ---")
    out_plain, R0 = hiruma_3term_lanczos(K_phys, M_sigma, b_sigma, N_stages,
                                          use_arnoldi=False, verbose=True)

    # Arnoldi (full reorthog)
    print(f"\n--- Hiruma 3-term + Arnoldi reorthogonalization ---")
    out_arnoldi, R0 = hiruma_3term_lanczos(K_phys, M_sigma, b_sigma, N_stages,
                                            use_arnoldi=True, verbose=True)

    # Compare
    print(f"\n--- Comparison (tau_pair us) ---")
    print(f"  {'k':>3} {'Plain Lanczos':>15} {'Arnoldi':>15} "
          f"{'diff %':>10}")
    for k in range(len(out_plain)):
        tau_p = out_plain[k]["tau_pair_us"]
        tau_a = out_arnoldi[k]["tau_pair_us"]
        if tau_p and tau_a:
            diff = abs(tau_p - tau_a) / abs(tau_a) * 100
            print(f"  {k:>3} {tau_p:>15.6f} {tau_a:>15.6f} {diff:>10.4f}")

    # Reference: QD-Padé Kameari (from existing JSON)
    print(f"\n--- Reference: QD-Padé Kameari τ_pair (existing GPU sweep) ---")
    ref_json = (Path(__file__).parent / "hex_vim_gpu_cuboid521_results.json")
    if ref_json.exists():
        ref = json.loads(ref_json.read_text(encoding="utf-8"))
        for key in ref:
            if "order_4" in key:
                print(f"  {key}:")
                for r in ref[key]["Cauer_rungs_Kameari"][:6]:
                    print(f"    k={r['k']}: tau_pair = {r['tau_pair_us']:.4f}")

    out_path = (Path(__file__).parent
                / "hex_vim_hiruma_arnoldi_cuboid521_results.json")
    out_path.write_text(json.dumps({
        "order": p,
        "n_dofs": basis.n_dofs,
        "lanczos_plain": out_plain,
        "arnoldi": out_arnoldi,
        "R_0": R0,
    }, indent=2), encoding="utf-8")
    print(f"\nSaved: {out_path}")


if __name__ == "__main__":
    main()
