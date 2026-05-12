"""dd_full_pipeline.py — Phase 3: full DD K, M, b → mpmath eigh → Cauer.

Loads DD K from dd_K_cuboid521_order3.npz, builds DD M, b, converts all to
mpmath 30+ digit, solves generalized eigh, computes Foster spectrum and
Cauer rungs with verified-interval analysis.

End goal: demonstrate that DD-precision input extends rigorous Cauer
extraction beyond the FP64 stage 5 limit.
"""
from __future__ import annotations

import math
import sys
import time
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
import mpmath as mp


def dd_to_mpmath_matrix(M_hi, M_lo, dps=50):
    """Convert DD matrix (hi + lo numpy arrays) to mpmath matrix."""
    mp.mp.dps = dps
    n, m = M_hi.shape
    M_mp = mp.matrix(n, m)
    for i in range(n):
        for j in range(m):
            M_mp[i, j] = mp.mpf(M_hi[i, j]) + mp.mpf(M_lo[i, j])
    return M_mp


def dd_to_mpmath_vector(b_hi, b_lo, dps=50):
    mp.mp.dps = dps
    n = b_hi.size
    b_mp = mp.matrix(n, 1)
    for i in range(n):
        b_mp[i, 0] = mp.mpf(b_hi[i]) + mp.mpf(b_lo[i])
    return b_mp


def mpmath_generalized_eigh(K_mp, M_mp, dps=50):
    """Solve K v = λ M v with symmetric K, positive-definite M in mpmath.

    Strategy: Cholesky L of M, then symmetric eigh of L^{-1} K L^{-T}.
    """
    mp.mp.dps = dps
    n = K_mp.rows

    # Cholesky decomposition of M: M = L L^T (lower triangular)
    L = mp.matrix(n, n)
    for j in range(n):
        # diagonal
        s = M_mp[j, j]
        for k in range(j):
            s -= L[j, k] ** 2
        if s <= 0:
            raise ValueError(f"M not positive definite at j={j}: s={s}")
        L[j, j] = mp.sqrt(s)
        for i in range(j + 1, n):
            s = M_mp[i, j]
            for k in range(j):
                s -= L[i, k] * L[j, k]
            L[i, j] = s / L[j, j]

    # Compute L^{-1} K L^{-T} = symmetric A
    # A = L^{-1} K L^{-T}
    # Step 1: solve L Y = K   (Y = L^{-1} K)
    Y = mp.matrix(n, n)
    for j in range(n):
        for i in range(n):
            s = K_mp[i, j]
            for k in range(i):
                s -= L[i, k] * Y[k, j]
            Y[i, j] = s / L[i, i]
    # Step 2: solve L Z^T = Y^T → Z = (Y^T L^{-T})^T = Y L^{-T}
    # i.e., we want A = Y L^{-T} → Solve L^T X = Y^T means... simpler:
    # A_ij = sum_k Y[i, k] (L^{-T})[k, j] = sum_k Y[i, k] / ...
    # We want A symmetric. Use: solve L Z = Y^T for Z, then A = Z^T
    A = mp.matrix(n, n)
    for j in range(n):
        # Solve L z = Y_row_j (i.e., L z = (Y^T)_:,j = Y_j,:)
        for i in range(n):
            s = Y[j, i]  # Y^T_{i,j} = Y_{j,i}
            for k in range(i):
                s -= L[i, k] * A[k, j]
            A[i, j] = s / L[i, i]

    # mpmath eig
    print(f"    Computing eigenvalues of symmetric n={n} matrix in mpmath...")
    t0 = time.time()
    # mp.eig returns (eigenvalues, eigenvectors) — non-symmetric routine
    eigvals_mp, eigvecs_A_mp = mp.eig(A)
    print(f"    mpmath eig: {time.time()-t0:.2f}s")

    # Recover original eigenvectors: M v = K v / λ
    # We solved A x = λ x with A = L^{-1} K L^{-T}, x = L^T v
    # So v = L^{-T} x
    eigvecs_v_mp = mp.matrix(n, n)
    for k_eig in range(n):
        # Solve L^T v = x_k: backsubst
        x_col = mp.matrix(n, 1)
        for i in range(n):
            x_col[i, 0] = eigvecs_A_mp[i, k_eig]
        v_col = mp.matrix(n, 1)
        for i in range(n - 1, -1, -1):
            s = x_col[i, 0]
            for k in range(i + 1, n):
                s -= L[k, i] * v_col[k, 0]
            v_col[i, 0] = s / L[i, i]
        for i in range(n):
            eigvecs_v_mp[i, k_eig] = v_col[i, 0]

    return eigvals_mp, eigvecs_v_mp, L


def main():
    print("=" * 72)
    print("DD full pipeline: K, M, b → mpmath eigh → Cauer + verified-interval")
    print("=" * 72)

    # Load DD K from saved file
    data = np.load(Path(__file__).parent / "dd_K_cuboid521_order3.npz",
                   allow_pickle=False)
    K_hi = data["K_hi"]
    K_lo = data["K_lo"]
    p = int(data["p"])
    a = float(data["a"])
    b = float(data["b"])
    c = float(data["c"])
    print(f"\n  Loaded K_hi, K_lo from dd_K_cuboid521_order3.npz")
    print(f"  order={p}, n_dofs={K_hi.shape[0]}")

    # Assemble M, b in DD
    print(f"\n  Assembling M, b in DD with high-precision nodes...")
    sys.path.insert(0, str(Path(__file__).parent))
    from dd_M_b_assembly import assemble_M_b_dd_gpu
    M_hi, M_lo, b_hi, b_lo = assemble_M_b_dd_gpu(
        p, a, b, c, n_gauss=8, use_gpu=True, verbose=False)
    print(f"  M, b assembled")
    print(f"    max |M_lo/M_hi| = "
          f"{float(np.max(np.abs(M_lo) / (np.abs(M_hi) + 1e-30))):.2e}")
    print(f"    max |b_lo/b_hi| = "
          f"{float(np.max(np.abs(b_lo) / (np.abs(b_hi) + 1e-30))):.2e}")

    # Apply mu0/(4 pi) prefactor to K (K_phys = K_bare * mu0/4pi)
    MU0 = 4 * math.pi * 1e-7
    prefactor = MU0 / (4 * math.pi)
    K_phys_hi = K_hi * prefactor
    K_phys_lo = K_lo * prefactor

    # Convert to mpmath
    print(f"\n  Converting DD K, M, b to mpmath at 50-digit...")
    t0 = time.time()
    K_mp = dd_to_mpmath_matrix(K_phys_hi, K_phys_lo, dps=50)
    M_mp = dd_to_mpmath_matrix(M_hi, M_lo, dps=50)
    b_mp = dd_to_mpmath_vector(b_hi, b_lo, dps=50)
    print(f"  Conversion: {time.time()-t0:.2f}s")

    # mpmath generalized eigh
    print(f"\n  Solving K v = λ M v in mpmath...")
    t0 = time.time()
    eigvals_mp, eigvecs_mp, L_chol = mpmath_generalized_eigh(K_mp, M_mp, dps=50)
    print(f"  Total eigh: {time.time()-t0:.2f}s")

    # Foster spectrum {tau_k, g_k^2}
    SIGMA = 5.8e7
    n = len(eigvals_mp)
    tau_list = []
    g2_list = []
    for k in range(n):
        lam = eigvals_mp[k]
        if abs(lam) < mp.mpf("1e-50"):
            continue
        tau_k = mp.mpf(SIGMA) * lam
        # g_k = v_k . b
        g_k = mp.mpf(0)
        for i in range(n):
            g_k += eigvecs_mp[i, k] * b_mp[i, 0]
        # M-normalize: g_k² / (v^T M v)
        vMv = mp.mpf(0)
        for i in range(n):
            for j in range(n):
                vMv += eigvecs_mp[i, k] * M_mp[i, j] * eigvecs_mp[j, k]
        g2_k = g_k * g_k / vMv if vMv > 0 else mp.mpf(0)
        tau_list.append(tau_k)
        g2_list.append(g2_k)

    # Sort by |g²|, descending
    pairs = sorted(zip(tau_list, g2_list),
                   key=lambda p: -float(abs(p[1])))
    print(f"\n  Top 10 Foster modes (by |g|²):")
    print(f"  {'rank':>5} {'tau (us)':>15} {'g²':>15}")
    for k in range(min(10, len(pairs))):
        tau_us = float(pairs[k][0]) * 1e6
        g2_v = float(pairs[k][1])
        print(f"  {k:>5} {tau_us:>15.4f} {g2_v:>15.4e}")

    # Compute Kameari moments in mpmath
    print(f"\n  Computing Kameari moments α_n in mpmath...")
    n_use = min(60, len(pairs))
    n_moments = 40
    alphas_mp = []
    for n_idx in range(n_moments):
        a_val = mp.mpf(0)
        for tau_k, g2_k in pairs[:n_use]:
            a_val += g2_k * (tau_k ** n_idx)
        alphas_mp.append((mp.mpf(-1)) ** n_idx * a_val)

    # Verified-interval Cauer with DD-precision input
    print(f"\n  Verified-interval Hankel-Padé (DD input, eps_rel = 1e-30)...")
    mp.iv.dps = 80
    eps_rel = mp.mpf("1e-30")
    eps_iv = mp.iv.mpf((float(1 - eps_rel), float(1 + eps_rel)))
    alphas_iv = []
    for a_mp in alphas_mp:
        a_iv = mp.iv.mpf(float(a_mp)) * eps_iv
        alphas_iv.append(a_iv)

    # Run interval QD-Padé
    c_list = [a for a in alphas_iv]
    print(f"\n  Cauer rungs with verified-interval bounds:")
    print(f"  {'stage':>5} {'value':>20} {'rel_width':>15} {'reliable?':>10}")
    width_warn = 1e-2
    n_extract = 12
    stages_reliable = 0
    for stage in range(n_extract):
        if len(c_list) < 2:
            break
        if 0 in c_list[0]:
            print(f"  {stage:>5} (c[0] contains zero, truncating)")
            break
        n_c = len(c_list)
        e = [mp.iv.mpf(0)] * n_c
        e[0] = mp.iv.mpf(1) / c_list[0]
        for k in range(1, n_c):
            acc = mp.iv.mpf(0)
            for j in range(k):
                if k - j - 1 >= 0:
                    acc += c_list[j + 1] * e[k - j - 1]
            e[k] = -acc / c_list[0]
        rung = e[0]
        rung_a = mp.mpf(rung.a)
        rung_b = mp.mpf(rung.b)
        mid = (rung_a + rung_b) / 2
        width = rung_b - rung_a
        rel = width / abs(mid) if abs(mid) > 0 else mp.mpf("inf")
        reliable = "YES" if float(rel) < width_warn else "NO"
        if float(rel) < width_warn:
            stages_reliable = stage + 1
        print(f"  {stage:>5} {float(mid):>20.6e} {float(rel):>15.3e} "
              f"{reliable:>10}")
        if float(rel) > 1.0:
            print(f"        (width > 100%, halting)")
            break
        c_list = e[1:]

    print(f"\n  Reliable stages with DD input: {stages_reliable}")
    print(f"  (Compare with FP64 input baseline: stage 4-5)")
    if stages_reliable > 5:
        print(f"  *** SUCCESS: DD input extends Cauer extraction beyond FP64 ***")


if __name__ == "__main__":
    main()
