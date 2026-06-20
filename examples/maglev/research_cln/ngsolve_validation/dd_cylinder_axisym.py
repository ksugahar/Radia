"""dd_cylinder_axisym.py — DD axisymmetric VIM for Cu cylinder.

Port of cylinder_vim_axisym_cupy.py to DD precision. Uses mpmath elliptic
integrals (K(m), E(m)) per pair for DD-precision Green's function.

CPU-bound (mpmath calls ~1 ms per pair). For Nr=24 Nz=6 (144 cells) with
n_quad=2 → 576 q-points → 332K pairs → ~6 min. Scales to ~90 min at full
576 cells.
"""
from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
import mpmath as mp

sys.path.insert(0, str(Path(__file__).parent))
from dd_arithmetic import dd_add, dd_mul, dd_from_mpmath
from dd_axisym_kernel import axisym_kernel_dd_mpmath


R_DISK = 10e-3
T_DISK = 2e-3
SIGMA = 5.8e7
MU0 = 4 * math.pi * 1e-7
B0 = 1.0


def assemble_K_M_b_cylinder_dd(Nr, Nz, n_quad=2, verbose=True):
    """DD K, M, b for cylinder axisym VIM. Returns (K_hi, K_lo, M, b) numpy."""
    mp.mp.dps = 40
    N = Nr * Nz
    dr = R_DISK / Nr
    dz = T_DISK / Nz

    r_a = np.tile(np.arange(Nr)[:, None], (1, Nz)) * dr
    r_b = r_a + dr
    z_a = np.tile(np.arange(Nz)[None, :], (Nr, 1)) * dz - T_DISK / 2.0
    z_b = z_a + dz

    M_diag = (math.pi * (r_b ** 2 - r_a ** 2) * (z_b - z_a)).flatten()
    b_vec = ((B0 / 2.0) * (2.0 * math.pi / 3.0)
             * (r_b ** 3 - r_a ** 3) * (z_b - z_a)).flatten()

    # Gauss-Legendre nodes/weights at high precision
    from dd_gl_nodes import gauss_legendre_mpf
    nodes_mp, weights_mp = gauss_legendre_mpf(n_quad, dps=40)
    nodes_01_mp = [(n + 1) / 2 for n in nodes_mp]
    w_01_mp = [w / 2 for w in weights_mp]

    # Cell quadrature points (DD)
    n_q_total = Nr * Nz * n_quad * n_quad
    r_flat_hi = np.zeros(n_q_total)
    r_flat_lo = np.zeros(n_q_total)
    z_flat_hi = np.zeros(n_q_total)
    z_flat_lo = np.zeros(n_q_total)
    jac_flat_hi = np.zeros(n_q_total)
    jac_flat_lo = np.zeros(n_q_total)
    cell_of_q = np.zeros(n_q_total, dtype=int)

    q = 0
    dr_mp = mp.mpf(dr)
    dz_mp = mp.mpf(dz)
    for i in range(Nr):
        for j in range(Nz):
            cell_idx = i * Nz + j
            r_a_mp = mp.mpf(r_a[i, j])
            z_a_mp = mp.mpf(z_a[i, j])
            for qr in range(n_quad):
                r_q_mp = r_a_mp + dr_mp * nodes_01_mp[qr]
                wr_mp = dr_mp * w_01_mp[qr]
                for qz in range(n_quad):
                    z_q_mp = z_a_mp + dz_mp * nodes_01_mp[qz]
                    wz_mp = dz_mp * w_01_mp[qz]
                    r_flat_hi[q], r_flat_lo[q] = dd_from_mpmath(r_q_mp)
                    z_flat_hi[q], z_flat_lo[q] = dd_from_mpmath(z_q_mp)
                    jac_mp = wr_mp * wz_mp
                    jac_flat_hi[q], jac_flat_lo[q] = dd_from_mpmath(jac_mp)
                    cell_of_q[q] = cell_idx
                    q += 1

    if verbose:
        print(f"  N={N} cells, n_q_total={n_q_total}, "
              f"n_pairs={n_q_total ** 2}")
        print(f"  Building DD G matrix via mpmath elliptic...")
        t0 = time.time()

    # K assembly: K[i, j] = (mu0/2) sum_qi sum_qj G(qi, qj) × (2 pi r_qi)
    #                       × jac_qi × jac_qj × indicator(qi in cell i)
    #                       × indicator(qj in cell j)
    # Compute G_full once (n_q × n_q), weight, then aggregate per cell pair.
    K_hi = np.zeros((N, N))
    K_lo = np.zeros((N, N))
    progress = 0
    progress_step = max(1, n_q_total // 20)
    for qi in range(n_q_total):
        if verbose and qi % progress_step == 0:
            print(f"    qi={qi}/{n_q_total} ({100*qi/n_q_total:.0f}%), "
                  f"elapsed {time.time()-t0:.0f}s")
        ci = cell_of_q[qi]
        # weight on observation side: 2 pi r_qi × jac_qi (DD)
        two_pi_mp = mp.mpf(2) * mp.pi
        r_qi_mp = mp.mpf(r_flat_hi[qi]) + mp.mpf(r_flat_lo[qi])
        jac_qi_mp = mp.mpf(jac_flat_hi[qi]) + mp.mpf(jac_flat_lo[qi])
        w_obs_mp = two_pi_mp * r_qi_mp * jac_qi_mp

        for qj in range(n_q_total):
            cj = cell_of_q[qj]
            # G in DD
            g_hi, g_lo = axisym_kernel_dd_mpmath(
                np.array(r_flat_hi[qi]), np.array(r_flat_lo[qi]),
                np.array(z_flat_hi[qi]), np.array(z_flat_lo[qi]),
                np.array(r_flat_hi[qj]), np.array(r_flat_lo[qj]),
                np.array(z_flat_hi[qj]), np.array(z_flat_lo[qj]),
                dps=40)
            # Multiply by w_obs * jac_qj (mpmath for precision)
            g_mp = mp.mpf(float(g_hi)) + mp.mpf(float(g_lo))
            jac_qj_mp = mp.mpf(jac_flat_hi[qj]) + mp.mpf(jac_flat_lo[qj])
            contrib_mp = g_mp * w_obs_mp * jac_qj_mp
            contrib_hi, contrib_lo = dd_from_mpmath(contrib_mp)
            # Add to K[ci, cj] in DD
            K_hi_old = K_hi[ci, cj]
            K_lo_old = K_lo[ci, cj]
            K_hi_new, K_lo_new = dd_add(
                np.array(K_hi_old), np.array(K_lo_old),
                np.array(contrib_hi), np.array(contrib_lo))
            K_hi[ci, cj] = float(K_hi_new)
            K_lo[ci, cj] = float(K_lo_new)

    # Multiply by (mu0/2) — exact in FP64 since mu0 = 4 pi × 1e-7
    mu0_over_2 = MU0 / 2.0
    K_hi = K_hi * mu0_over_2
    K_lo = K_lo * mu0_over_2

    if verbose:
        print(f"  K assembly: {time.time()-t0:.1f}s")

    return K_hi, K_lo, np.diag(M_diag), b_vec


def main():
    print("=" * 72)
    print("DD axisym VIM cylinder Cu disk R=10mm t=2mm")
    print("=" * 72)

    # Reduced cells for demo (production: Nr=48, Nz=12)
    Nr, Nz = 24, 6
    print(f"\n  Nr={Nr}, Nz={Nz}, total {Nr*Nz} cells")

    t0 = time.time()
    K_hi, K_lo, M, b = assemble_K_M_b_cylinder_dd(Nr, Nz, n_quad=2,
                                                    verbose=True)
    print(f"  Total assembly: {time.time()-t0:.1f}s")

    # K_lo/K_hi ratio
    rel_lo = np.max(np.abs(K_lo) / (np.abs(K_hi) + 1e-30))
    print(f"  K_lo/K_hi: {rel_lo:.2e}")

    # Convert to mpmath, eigh, Cauer
    print(f"\n  Converting to mpmath + eigh...")
    from dd_full_pipeline import dd_to_mpmath_matrix, mpmath_generalized_eigh
    K_mp = dd_to_mpmath_matrix(K_hi, K_lo, dps=50)
    # M and b are FP64 exact (polynomial) — convert directly
    mp.mp.dps = 50
    N = K_hi.shape[0]
    M_mp = mp.matrix(N, N)
    for i in range(N):
        for j in range(N):
            M_mp[i, j] = mp.mpf(float(M[i, j]))
    b_mp = mp.matrix(N, 1)
    for i in range(N):
        b_mp[i, 0] = mp.mpf(float(b[i]))

    t0 = time.time()
    eigvals_mp, eigvecs_mp, L_chol = mpmath_generalized_eigh(K_mp, M_mp, dps=50)
    print(f"  mpmath eigh: {time.time()-t0:.1f}s")

    # Foster spectrum: τ_k = σ × λ_k (cylinder convention)
    tau_list = []
    g2_list = []
    for k in range(N):
        lam = eigvals_mp[k]
        if abs(lam) < mp.mpf("1e-50"):
            continue
        tau_k = mp.mpf(SIGMA) * lam
        g_k = mp.mpf(0)
        for i in range(N):
            g_k += eigvecs_mp[i, k] * b_mp[i, 0]
        vMv = mp.mpf(0)
        for i in range(N):
            for j in range(N):
                vMv += eigvecs_mp[i, k] * M_mp[i, j] * eigvecs_mp[j, k]
        g2_k = g_k * g_k / vMv if vMv > 0 else mp.mpf(0)
        tau_list.append(tau_k)
        g2_list.append(g2_k)

    pairs = sorted(zip(tau_list, g2_list), key=lambda p: -float(abs(p[1])))
    print(f"\n  Top 5 Foster:")
    for k in range(5):
        print(f"    rank={k}: tau={float(pairs[k][0])*1e6:.4f} us, "
              f"g²={float(pairs[k][1]):.4e}")

    # Kameari moments + verified-interval QD-Padé
    n_use = min(60, len(pairs))
    n_moments = 40
    alphas_mp = []
    for n_idx in range(n_moments):
        a_val = mp.mpf(0)
        for tau_k, g2_k in pairs[:n_use]:
            a_val += g2_k * (tau_k ** n_idx)
        alphas_mp.append((mp.mpf(-1)) ** n_idx * a_val)

    mp.iv.dps = 80
    eps_rel = mp.mpf("1e-30")
    eps_iv = mp.iv.mpf((float(1 - eps_rel), float(1 + eps_rel)))
    alphas_iv = [mp.iv.mpf(float(a)) * eps_iv for a in alphas_mp]
    c_list = list(alphas_iv)
    rungs_raw = []
    for stage in range(20):
        if len(c_list) < 2 or 0 in c_list[0]:
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
        rel_width = float((rung_b - rung_a) / abs(mid)) if abs(mid) > 0 else float("inf")
        rungs_raw.append({"stage": stage, "mid": float(mid),
                          "rel_width": rel_width,
                          "reliable": rel_width < 1e-2})
        if rel_width > 1.0:
            break
        c_list = e[1:]

    print(f"\n  Cauer rungs (DD axisym):")
    print(f"  {'k':>3} {'R_2k':>15} {'L_2k+1':>15} {'tau_pair us':>13}")
    cauer_rungs = []
    for k in range(len(rungs_raw) // 2):
        R_raw = rungs_raw[2 * k]
        Linv_raw = rungs_raw[2 * k + 1]
        if abs(Linv_raw["mid"]) < 1e-50:
            break
        L_mid = 1.0 / Linv_raw["mid"]
        tau_p = L_mid / R_raw["mid"] * 1e6
        print(f"  {k:>3} {R_raw['mid']:>15.4e} {L_mid:>15.4e} {tau_p:>13.4f}")
        cauer_rungs.append({"k": k, "R_2k": R_raw["mid"], "L_2k_plus_1": L_mid,
                            "tau_pair_us": tau_p,
                            "reliable": R_raw["reliable"] and Linv_raw["reliable"]})

    out_path = (Path(__file__).parent
                / f"dd_cylinder_axisym_Nr{Nr}_Nz{Nz}_results.json")
    out_path.write_text(json.dumps({
        "shape": "cylinder", "Nr": Nr, "Nz": Nz, "n_quad": 2,
        "Cauer_rungs": cauer_rungs,
    }, indent=2), encoding="utf-8")
    print(f"\nSaved: {out_path}")


if __name__ == "__main__":
    main()
