"""sphere_vim_axisym_ir.py — Sphere VIM Foster moments via direct
Krylov iteration in mpmath (no FP64 eigh, no Cholesky).

Convention (matches sphere_vim_axisym.py):
  K (dense)    : axisymmetric Green's function kernel matrix (mu_0 / 2) ∬ G 2πr dV dV'
  M (diagonal) : per-cell volume V_i
  K v = λ M v with τ_n = σ * λ_n
  Foster:  f(s) = b^T (M + s σ K)^-1 b = Σ_n g_n^2 / (1 + s τ_n)

Moments (Taylor expansion of f around s=0):
  f(s) = M^-1 (I + s σ M^-1 K)^-1 acted on b, with the (-s)^k expansion:

      α_k_physical = σ^k * b^T (M^-1 K)^k M^-1 b      (k = 0, 1, 2, ...)

Cauer extraction uses the existing alpha convention
      alphas[k] = (-1)^k * α_physical[k+1]

Algorithm:
  u_0 = M^-1 b              (M diagonal: u_0[i] = b[i] / V[i],  exact in mpmath)
  for k = 1, 2, ...:
      u_k = M^-1 (K u_{k-1})       (K @ u in mpmath, then divide by V_i)
  α_physical[k] = σ^k * b · u_k    (mpmath dot product)
  alphas_hankel[k] = (-1)^k * α_physical[k+1]   (k = 0..n_moments-1)

The only mpmath-cost step is the dense K @ u matvec; this is O(n^2) per
iteration with no Cholesky required.  K's near-zero eigenvalues (the
sweep showed min λ ~ 1e-15 from spurious modes) are irrelevant: we never
solve K^-1 anything; we only multiply by K.

This pipeline computes Foster moments to mpmath precision from FP64 K, M, b
inputs.  Hankel-Padé Cauer extraction then yields R, L rungs whose
accuracy is limited only by the underlying FE discretisation (and not by
FP64 eigh precision as in sphere_vim_axisym.py).
"""
from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, str(Path(__file__).parent))

import numpy as np
import mpmath as mp

from sphere_vim_axisym import (
    assemble_K_M_b,
    hankel_pade_cauer,
    R_SPHERE,
    SIGMA,
    MU0,
    B0,
)

mp.mp.dps = 60   # working precision


def _to_mp_vec(v_fp64):
    return mp.matrix([mp.mpf(float(x)) for x in v_fp64])


def _to_mp_mat(A_fp64):
    n = A_fp64.shape[0]
    Mm = mp.matrix(n, n)
    for i in range(n):
        for j in range(n):
            Mm[i, j] = mp.mpf(float(A_fp64[i, j]))
    return Mm


def _mp_dot(u, v):
    """Column-vector dot product."""
    n = u.rows
    s = mp.mpf(0)
    for i in range(n):
        s += u[i, 0] * v[i, 0]
    return s


def compute_moments_krylov(K, M_diag, b, n_alphas, verbose=True):
    """Compute Foster moments α_physical[k] = σ^k * b^T (M^-1 K)^k M^-1 b
    for k = 0..n_alphas-1, using mpmath dense matvec.

    Args:
      K       : (n, n) FP64 dense kernel matrix (= VIM K)
      M_diag  : (n,)   FP64 diagonal of M (per-cell volume)
      b       : (n,)   FP64 source vector
      n_alphas: number of moments to compute (we need k=0..n_alphas-1)
      verbose : print progress

    Returns:
      list of mpmath alphas (α_physical[k])
    """
    n = K.shape[0]
    sigma_mp = mp.mpf(float(SIGMA))

    if verbose:
        print(f"  Lifting K (n={n}) to mpmath (dps={mp.mp.dps})...")
    t0 = time.time()
    K_mp = _to_mp_mat(K)
    if verbose:
        print(f"    K lift: {time.time() - t0:.1f}s")
    t0 = time.time()
    b_mp = _to_mp_vec(b)
    M_inv = [mp.mpf(1) / mp.mpf(float(M_diag[i])) for i in range(n)]
    if verbose:
        print(f"    b, M^-1 lift: {time.time() - t0:.2f}s")

    # u_0 = M^-1 b
    u_mp = mp.matrix(n, 1)
    for i in range(n):
        u_mp[i, 0] = b_mp[i, 0] * M_inv[i]

    alphas = []
    sigma_pow = mp.mpf(1)  # σ^k, k=0..
    for k in range(n_alphas):
        # α_physical[k] = σ^k * b · u_k
        bdotu = _mp_dot(b_mp, u_mp)
        alpha_k = sigma_pow * bdotu
        alphas.append(alpha_k)
        if verbose and (k < 8 or k % 5 == 0):
            print(f"    α_phys[{k}] = {mp.nstr(alpha_k, 10)}")

        if k == n_alphas - 1:
            break

        # u_{k+1} = M^-1 K u_k
        t0 = time.time()
        Ku = K_mp * u_mp
        for i in range(n):
            u_mp[i, 0] = Ku[i, 0] * M_inv[i]
        sigma_pow *= sigma_mp
        if verbose and k < 3:
            print(f"      step k={k+1}: matvec+scale {time.time() - t0:.2f}s")

    return alphas


def load_analytical(here: Path):
    p = here / "sphere_analytical_cln_highstage.json"
    d = json.load(open(p))
    return {
        "stoll_tau_us": [float(x) for x in d["stoll_tau_us"]],
        "cauer_tau_us": [float(x) for x in d["cauer_tau_us"]],
        "R0": float(d["R0"]),
        "n_reliable": int(d["n_reliable"]),
    }


def main():
    here = Path(__file__).parent
    ana = load_analytical(here)
    print("=" * 78)
    print(f"Sphere VIM axisym + Krylov-mpmath Foster moments")
    print(f"  Stoll Cauer τ_pair[0..7] (μs): "
          f"{[round(x, 4) for x in ana['cauer_tau_us'][:8]]}")
    print(f"  mpmath dps = {mp.mp.dps}")
    print("=" * 78)

    N_rho, N_theta = 20, 24       # 480 cells (validation)
    print(f"\nAssembling K, M, b "
          f"(N_rho={N_rho}, N_theta={N_theta}, n_cells={N_rho*N_theta})...")
    t0 = time.time()
    K, M, b = assemble_K_M_b(N_rho, N_theta, n_quad=2)
    print(f"  done in {time.time() - t0:.1f}s")
    M_diag = np.diag(M)
    print(f"  M diagonal? max|M - diag(M)| = "
          f"{np.max(np.abs(M - np.diag(M_diag))):.3e}")

    # Sanity check vs eigh-baseline α_physical[1] (= Σ g² τ)
    print("\nSanity reference (FP64 eigh):")
    from scipy.linalg import eigh
    evals, evecs = eigh(K, M)
    for k in range(evecs.shape[1]):
        nrm = math.sqrt(float(evecs[:, k] @ M @ evecs[:, k]))
        if nrm > 1e-30:
            evecs[:, k] = evecs[:, k] / nrm
    tau = SIGMA * evals
    pos = tau > 1e-30
    g = evecs.T @ b
    g2 = g * g
    alpha_phys_1_baseline = float((g2[pos] * tau[pos]).sum())
    print(f"  Σ g² τ (FP64 eigh) = {alpha_phys_1_baseline:.6e}")

    n_alphas = 30
    print(f"\nKrylov-mpmath moments α_physical[0..{n_alphas-1}]...")
    alphas_phys = compute_moments_krylov(K, M_diag, b, n_alphas)
    print(f"  α_physical[1] (Krylov) = {mp.nstr(alphas_phys[1], 10)}  "
          f"(eigh baseline {alpha_phys_1_baseline:.6e})")
    rel_diff = abs(float(alphas_phys[1]) - alpha_phys_1_baseline) / alpha_phys_1_baseline
    print(f"  rel diff = {rel_diff:.3e}")

    # Convert physical alphas to Hankel-Padé alphas convention:
    #   alphas_hankel[k] = (-1)^k * α_physical[k+1]   for k = 0..n_alphas-2
    n_h = n_alphas - 1
    alphas_hankel = [(mp.mpf(-1)) ** k * alphas_phys[k + 1] for k in range(n_h)]

    print(f"\nHankel-Padé Cauer extraction "
          f"(mpmath dps={mp.mp.dps}, n_moments={n_h})...")
    p = hankel_pade_cauer(alphas_hankel, 15)
    rungs = []
    for k in range(len(p) // 2):
        Rv = p[2 * k]
        Linv = p[2 * k + 1]
        if abs(Linv) < mp.mpf(10) ** (-50) or Rv == 0:
            break
        Lv = 1 / Linv
        tau_p_us = float(Lv / Rv) * 1e6
        rungs.append({"k": k, "R_2k": float(Rv),
                      "L_2k_plus_1": float(Lv),
                      "tau_pair_us": tau_p_us})

    print(f"\nCauer rungs (Krylov + mpmath moments, FE on {N_rho * N_theta} cells):")
    print(f"  k | τ_pair (VIM)   | Stoll τ_pair    | gap [%]")
    print(f"  " + "-" * 60)
    for k, rung in enumerate(rungs):
        if k < len(ana["cauer_tau_us"]):
            ref = ana["cauer_tau_us"][k]
            gap = (rung["tau_pair_us"] - ref) / ref * 100.0
            print(f"  {k} | {rung['tau_pair_us']:>14.6f} | "
                  f"{ref:>14.6f} | {gap:>+8.4f}")
        else:
            print(f"  {k} | {rung['tau_pair_us']:>14.6f} | "
                  f"{'--':>14s} | {'--':>8s}")

    out_path = here / "sphere_vim_axisym_ir_results.json"
    out = {
        "method": "Sphere axisym VIM + Krylov-mpmath Foster moments + Hankel-Padé",
        "n_cells": N_rho * N_theta,
        "N_rho": N_rho, "N_theta": N_theta,
        "mp_dps": mp.mp.dps,
        "n_alphas_physical": n_alphas,
        "alphas_physical_str": [mp.nstr(a, 20) for a in alphas_phys],
        "cauer_rungs": rungs,
    }
    out_path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"\nSaved: {out_path.name}")


if __name__ == "__main__":
    main()
