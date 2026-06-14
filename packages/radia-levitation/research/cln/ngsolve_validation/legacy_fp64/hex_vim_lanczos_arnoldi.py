"""hex_vim_lanczos_arnoldi.py — VIM-native Lanczos+Arnoldi Cauer extraction.

Apply Lanczos process to T = sigma * K^{-1} M with M-inner product,
starting from u_0 = M^{-1} b normalized. Each stage produces (alpha_n, beta_n)
Jacobi matrix entries. With Arnoldi reorthogonalization, maintains M-orthogonality
of Krylov vectors against ALL previous stages (FP64 precision).

The Jacobi tridiagonal J = tridiag(alpha, sqrt(beta)) has eigenvalues that
are a compressed Foster spectrum {tau_k} (in seconds, since T eigenvalues
= sigma × lambda_K-M = tau). Weights g_k^2 = b_0 * (first eigenvector component)^2.

This gives an INDEPENDENT path from K, M, b to Foster spectrum, vs scipy.eigh
on the full N×N problem. Plus Arnoldi maintains numerical orthogonality.

Then we apply QD-Pade on the reconstructed spectrum to extract Cauer rungs,
and compare with the direct QD-Pade Kameari pipeline.

Inspired by Nagamine's rep/rectangle/arnoldi.edp.
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
from scipy.linalg import eigh_tridiagonal
import mpmath as mp
mp.mp.dps = 80

sys.path.insert(0, str(Path(__file__).parent))
from hex_vim_cupy_kassembly import assemble_K_cupy
from hex_vim_cupy import evaluate_basis

sys.path.insert(0, str(Path("S:/Radia/01_GitHub/packages/radia-vim/src")))
import radia_vim


# Cuboid 5x2x1 mm Cu
A_M, B_M, C_M = 5e-3, 2e-3, 1e-3
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
    return (a / (b * c)) * Mref[0] + (b / (a * c)) * Mref[1] + (c / (a * b)) * Mref[2]


def assemble_b_vector(basis, a, b, c, p, n_gauss=8):
    nodes_m1, w_m1 = np.polynomial.legendre.leggauss(n_gauss)
    nodes = 0.5 * (nodes_m1 + 1.0)
    weights = 0.5 * w_m1
    XI, ETA, ZETA = np.meshgrid(nodes, nodes, nodes, indexing="ij")
    Phi_ref = evaluate_basis(XI.flatten(), ETA.flatten(), ZETA.flatten(), p,
                              xp_mod=np)
    x_phys = a * (XI.flatten() - 0.5)
    y_phys = b * (ETA.flatten() - 0.5)
    Ax = B0 * (-y_phys / 2)
    Ay = B0 * (x_phys / 2)
    W3 = (weights[:, None, None] * weights[None, :, None]
          * weights[None, None, :]).flatten()
    return (Phi_ref[..., 0].T @ (W3 * Ax * a)) + (Phi_ref[..., 1].T @ (W3 * Ay * b))


def lanczos_arnoldi_M(K, M, b, N_stages, use_arnoldi=True, verbose=True):
    """Lanczos on T = sigma * K^{-1} M with M-inner product.

    u_0 = M^{-1} b / ||M^{-1} b||_M
    alpha_n = u_n^T M T u_n
    beta_n  = ||residual||_M after subtracting projections
    Cross-validate: eigenvalues of J = tridiag(alpha, sqrt(beta)) ≈ tau_k (seconds)
    Weights g_k^2 = b_M_norm^2 × (eigvec_first_component)^2

    Returns: alpha, beta, b_M_norm_sq (= b^T M^{-1} b)
    """
    K_lu = la.lu_factor(K)
    M_lu = la.lu_factor(M)

    u = la.lu_solve(M_lu, b)        # M^{-1} b
    b_M_norm_sq = float(b @ u)       # = b^T M^{-1} b
    u_normalized = u / math.sqrt(b_M_norm_sq)

    V = [u_normalized]
    alpha = []
    beta = []
    drift_log = []

    for n in range(N_stages):
        # T u_n = sigma * M^{-1} K u_n (eigenvalues = tau_k = sigma * lambda)
        KV = K @ V[n]
        Tv = SIGMA * la.lu_solve(M_lu, KV)

        a_n = float(V[n] @ M @ Tv)
        alpha.append(a_n)

        if n == 0:
            w = Tv - a_n * V[n]
        else:
            w = Tv - a_n * V[n] - beta[n - 1] * V[n - 1]

        # Arnoldi reorthogonalization (with M-inner product)
        if use_arnoldi:
            for k in range(n + 1):
                proj = float(V[k] @ M @ w)
                w -= V[k] * proj
            # Double reorthog for safety (Gram-Schmidt twice)
            for k in range(n + 1):
                proj = float(V[k] @ M @ w)
                w -= V[k] * proj

        b_n_squared = float(w @ M @ w)
        if b_n_squared <= 0:
            if verbose:
                print(f"  Stage {n}: breakdown beta²={b_n_squared:.3e}, truncating")
            break
        b_n = math.sqrt(b_n_squared)
        beta.append(b_n)

        # Measure orthogonality drift (without Arnoldi case)
        max_cross = 0.0
        for k in range(n + 1):
            cross = abs(float(V[k] @ M @ (w / b_n)))
            max_cross = max(max_cross, cross)
        drift_log.append(max_cross)

        V.append(w / b_n)

    return alpha, beta, b_M_norm_sq, drift_log


def cauer_via_jacobi_eig(alpha, beta, b_norm_sq, n_extract=12):
    """Build Jacobi matrix from (alpha, beta), compute eigenvalues = compressed
    Foster spectrum tau_k, weights g_k^2 = b_norm_sq × (eigvec[0])^2,
    then apply QD-Padé Kameari moment formula."""
    N = len(alpha)
    diag = np.array(alpha)
    offdiag = np.array(beta[:N - 1])  # length N-1

    eigvals, eigvecs = eigh_tridiagonal(diag, offdiag)
    weights = b_norm_sq * eigvecs[0, :] ** 2

    # Sort eigenvalues descending (largest tau first)
    order = np.argsort(-eigvals)
    tau_us = eigvals[order] * 1e6     # tau in microseconds (since alpha in seconds)
    g2 = weights[order]

    # QD-Padé Kameari moments
    tau_mp = [mp.mpf(float(t)) * mp.mpf("1e-6") for t in tau_us]
    g2_mp = [mp.mpf(float(g)) for g in g2]
    n_moments = min(40, 2 * N - 4)
    alphas_mp = []
    for n in range(n_moments):
        a_val = mp.mpf(0)
        for gv, tv in zip(g2_mp, tau_mp):
            a_val = a_val + gv * (tv ** n)
        alphas_mp.append((mp.mpf(-1)) ** n * a_val)

    c = [mp.mpf(a) for a in alphas_mp]
    cauer_p = []
    for _ in range(n_extract):
        if len(c) < 2 or abs(c[0]) < mp.mpf(10) ** (-mp.mp.dps + 10):
            break
        n_c = len(c)
        e = [mp.mpf(0)] * n_c
        e[0] = 1 / c[0]
        for k in range(1, n_c):
            e[k] = -sum(c[j + 1] * e[k - j - 1] for j in range(k)
                        if k - j - 1 >= 0) / c[0]
        cauer_p.append(e[0])
        c = e[1:]

    rungs = []
    for k in range(min(6, len(cauer_p) // 2)):
        Rv = cauer_p[2 * k]
        Linv = cauer_p[2 * k + 1]
        if abs(Linv) < mp.mpf(10) ** (-50):
            break
        Lv = 1 / Linv
        tau_p = float(Lv / Rv) * 1e6
        rungs.append({"k": k, "R_2k": float(Rv),
                      "L_2k_plus_1": float(Lv), "tau_pair_us": tau_p})
    return rungs, tau_us[:N], g2[:N]


def main():
    p = 4
    n_th, n_ph, n_rh, n_c = 12, 24, 12, 6
    print("=" * 72)
    print(f"Hex VIM cuboid 5x2x1 order={p}: VIM-direct Lanczos+Arnoldi vs Foster+QD")
    print("=" * 72)

    print(f"\nAssembling K, M, b...")
    t0 = time.time()
    K_bare = assemble_K_cupy(p, A_M, B_M, C_M, n_th, n_ph, n_rh, n_c,
                              use_gpu=True, verbose=False)
    K_phys = K_bare * (MU0 / (4 * math.pi))
    basis = radia_vim.HDivDivFreeHexBasis(p)
    M = assemble_mass_matrix_python(basis, A_M, B_M, C_M, n_gauss=8)
    b = assemble_b_vector(basis, A_M, B_M, C_M, p, n_gauss=8)
    print(f"  done in {time.time()-t0:.2f}s, n_dofs={K_phys.shape[0]}")

    # Reference: scipy.eigh + QD-Padé (existing baseline)
    print(f"\n--- Reference: scipy.eigh full Foster + QD-Padé ---")
    eigvals_ref, eigvecs_ref = la.eigh(K_phys, M)
    for k in range(len(eigvals_ref)):
        v = eigvecs_ref[:, k]
        nrm = math.sqrt(float(v @ M @ v))
        if nrm > 1e-30:
            eigvecs_ref[:, k] = v / nrm
    tau_us_ref = np.where(eigvals_ref > 1e-30, 1e6 * SIGMA * eigvals_ref, 0.0)
    g_ref = eigvecs_ref.T @ b
    g2_ref = g_ref * g_ref
    order = np.argsort(-g2_ref)
    print(f"  Top 5 Foster modes (full eigh):")
    for r in range(5):
        idx = order[r]
        print(f"    rank={r}: tau={tau_us_ref[idx]:.4f} us, g2={g2_ref[idx]:.4e}")

    # Plain Lanczos (no reorthog)
    print(f"\n--- Plain Lanczos (no reorthog) ---")
    a_pl, b_pl, bn_pl, drift_pl = lanczos_arnoldi_M(
        K_phys, M, b, N_stages=12, use_arnoldi=False, verbose=True)
    print(f"  n_stages computed: {len(a_pl)}")
    print(f"  alpha[0..4]: {[f'{a*1e6:.4e}' for a in a_pl[:5]]}")
    print(f"  beta[0..3]:  {[f'{b:.4e}' for b in b_pl[:4]]}")
    print(f"  M-orth drift max: {[f'{d:.2e}' for d in drift_pl[:6]]}")
    rungs_pl, _, _ = cauer_via_jacobi_eig(a_pl, b_pl, bn_pl, n_extract=12)
    print(f"  Cauer rungs (via Jacobi → QD-Padé):")
    for r in rungs_pl[:6]:
        print(f"    k={r['k']}: tau_pair = {r['tau_pair_us']:.4f} us")

    # Arnoldi (full reorthog)
    print(f"\n--- Lanczos + Arnoldi reorthog ---")
    a_a, b_a, bn_a, drift_a = lanczos_arnoldi_M(
        K_phys, M, b, N_stages=12, use_arnoldi=True, verbose=True)
    print(f"  n_stages computed: {len(a_a)}")
    print(f"  alpha[0..4]: {[f'{a*1e6:.4e}' for a in a_a[:5]]}")
    print(f"  beta[0..3]:  {[f'{b:.4e}' for b in b_a[:4]]}")
    print(f"  M-orth drift max: {[f'{d:.2e}' for d in drift_a[:6]]}")
    rungs_a, tau_arn, g2_arn = cauer_via_jacobi_eig(a_a, b_a, bn_a, n_extract=12)
    print(f"  Cauer rungs (via Jacobi → QD-Padé):")
    for r in rungs_a[:6]:
        print(f"    k={r['k']}: tau_pair = {r['tau_pair_us']:.4f} us")

    print(f"\n  Lanczos-compressed top Foster modes:")
    for k in range(min(5, len(tau_arn))):
        print(f"    rank={k}: tau={float(tau_arn[k]):.4f} us, "
              f"g2={float(g2_arn[k]):.4e}")

    out_path = (Path(__file__).parent
                / "hex_vim_lanczos_arnoldi_cuboid521_results.json")
    out_path.write_text(json.dumps({
        "order": p, "n_dofs": K_phys.shape[0],
        "plain_lanczos": {
            "alpha_seconds": a_pl, "beta": b_pl,
            "b_M_norm_sq": bn_pl, "drift": drift_pl,
            "cauer_rungs": rungs_pl,
        },
        "arnoldi": {
            "alpha_seconds": a_a, "beta": b_a,
            "b_M_norm_sq": bn_a, "drift": drift_a,
            "cauer_rungs": rungs_a,
        },
    }, indent=2), encoding="utf-8")
    print(f"\nSaved: {out_path}")


if __name__ == "__main__":
    main()
