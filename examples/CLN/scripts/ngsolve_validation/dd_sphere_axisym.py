"""dd_sphere_axisym.py — DD axisymmetric VIM for Cu sphere.

Port of sphere_vim_axisym.py to DD via mpmath elliptic integrals.
(ρ, θ) spherical grid, piecewise constant J_phi basis.

For Nρ=15 Nθ=18 = 270 cells × n_quad=2 → 1080 q-points → 1.2M pairs × 1 ms
= ~20 min CPU.
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


R_SPHERE = 10e-3
SIGMA = 5.8e7
MU0 = 4 * math.pi * 1e-7
B0 = 1.0


def assemble_K_M_b_sphere_dd(N_rho, N_theta, n_quad=2, verbose=True):
    mp.mp.dps = 40
    N = N_rho * N_theta
    drho = R_SPHERE / N_rho
    dtheta = math.pi / N_theta

    # Cell coords for M, b (FP64 exact closed-form)
    rho_a = np.tile(np.arange(N_rho)[:, None], (1, N_theta)) * drho
    rho_b = rho_a + drho
    th_a = np.tile(np.arange(N_theta)[None, :], (N_rho, 1)) * dtheta
    th_b = th_a + dtheta

    M_diag = ((2.0 * math.pi / 3.0)
              * (rho_b ** 3 - rho_a ** 3)
              * (np.cos(th_a) - np.cos(th_b))).flatten()
    b_vec = (math.pi * B0 * (rho_b ** 4 - rho_a ** 4) / 4.0
             * ((th_b - th_a) / 2.0
                - (np.sin(2.0 * th_b) - np.sin(2.0 * th_a)) / 4.0)).flatten()

    # Quadrature points: (rho_q, theta_q) per cell, n_quad^2 each
    from dd_gl_nodes import gauss_legendre_mpf
    nodes_mp, weights_mp = gauss_legendre_mpf(n_quad, dps=40)
    nodes_01 = [(n + 1) / 2 for n in nodes_mp]
    w_01 = [w / 2 for w in weights_mp]

    n_q_total = N_rho * N_theta * n_quad * n_quad
    r_flat_hi = np.zeros(n_q_total)
    r_flat_lo = np.zeros(n_q_total)
    z_flat_hi = np.zeros(n_q_total)
    z_flat_lo = np.zeros(n_q_total)
    jac_flat_hi = np.zeros(n_q_total)
    jac_flat_lo = np.zeros(n_q_total)
    cell_of_q = np.zeros(n_q_total, dtype=int)

    drho_mp = mp.mpf(drho)
    dtheta_mp = mp.mpf(dtheta)

    q = 0
    for i in range(N_rho):
        for j in range(N_theta):
            cell_idx = i * N_theta + j
            rho_a_mp = mp.mpf(rho_a[i, j])
            th_a_mp = mp.mpf(th_a[i, j])
            for qrho in range(n_quad):
                rho_q_mp = rho_a_mp + drho_mp * nodes_01[qrho]
                wrho_mp = drho_mp * w_01[qrho]
                for qth in range(n_quad):
                    th_q_mp = th_a_mp + dtheta_mp * nodes_01[qth]
                    wth_mp = dtheta_mp * w_01[qth]
                    # r = rho sin(theta), z = rho cos(theta)
                    r_q_mp = rho_q_mp * mp.sin(th_q_mp)
                    z_q_mp = rho_q_mp * mp.cos(th_q_mp)
                    # Jacobian dr dz = rho drho dtheta
                    jac_q_mp = rho_q_mp * wrho_mp * wth_mp
                    r_flat_hi[q], r_flat_lo[q] = dd_from_mpmath(r_q_mp)
                    z_flat_hi[q], z_flat_lo[q] = dd_from_mpmath(z_q_mp)
                    jac_flat_hi[q], jac_flat_lo[q] = dd_from_mpmath(jac_q_mp)
                    cell_of_q[q] = cell_idx
                    q += 1

    if verbose:
        print(f"  N={N} cells, n_q_total={n_q_total}, "
              f"n_pairs={n_q_total ** 2}")
        print(f"  Computing G(qi, qj) via mpmath elliptic...")
        t0 = time.time()

    two_pi_mp = mp.mpf(2) * mp.pi
    K_hi = np.zeros((N, N))
    K_lo = np.zeros((N, N))
    progress_step = max(1, n_q_total // 20)
    for qi in range(n_q_total):
        if verbose and qi % progress_step == 0:
            print(f"    qi={qi}/{n_q_total} ({100*qi/n_q_total:.0f}%), "
                  f"elapsed {time.time()-t0:.0f}s")
        if r_flat_hi[qi] < 1e-15:
            continue
        ci = cell_of_q[qi]
        r_qi_mp = mp.mpf(r_flat_hi[qi]) + mp.mpf(r_flat_lo[qi])
        jac_qi_mp = mp.mpf(jac_flat_hi[qi]) + mp.mpf(jac_flat_lo[qi])
        w_obs_mp = two_pi_mp * r_qi_mp * jac_qi_mp

        for qj in range(n_q_total):
            if r_flat_hi[qj] < 1e-15:
                continue
            cj = cell_of_q[qj]
            g_hi, g_lo = axisym_kernel_dd_mpmath(
                np.array(r_flat_hi[qi]), np.array(r_flat_lo[qi]),
                np.array(z_flat_hi[qi]), np.array(z_flat_lo[qi]),
                np.array(r_flat_hi[qj]), np.array(r_flat_lo[qj]),
                np.array(z_flat_hi[qj]), np.array(z_flat_lo[qj]),
                dps=40)
            g_mp = mp.mpf(float(g_hi)) + mp.mpf(float(g_lo))
            jac_qj_mp = mp.mpf(jac_flat_hi[qj]) + mp.mpf(jac_flat_lo[qj])
            contrib_mp = g_mp * w_obs_mp * jac_qj_mp
            contrib_hi, contrib_lo = dd_from_mpmath(contrib_mp)
            K_hi_new, K_lo_new = dd_add(
                np.array(K_hi[ci, cj]), np.array(K_lo[ci, cj]),
                np.array(contrib_hi), np.array(contrib_lo))
            K_hi[ci, cj] = float(K_hi_new)
            K_lo[ci, cj] = float(K_lo_new)

    mu0_over_2 = MU0 / 2.0
    K_hi = K_hi * mu0_over_2
    K_lo = K_lo * mu0_over_2

    if verbose:
        print(f"  K assembly: {time.time()-t0:.1f}s")

    return K_hi, K_lo, np.diag(M_diag), b_vec


def main():
    print("=" * 72)
    print("DD axisym VIM sphere Cu R=10mm")
    print("=" * 72)

    # Reduced for demo (~20 min)
    N_rho, N_theta = 15, 18
    print(f"\n  N_rho={N_rho}, N_theta={N_theta}, total {N_rho*N_theta} cells")

    t0 = time.time()
    K_hi, K_lo, M, b = assemble_K_M_b_sphere_dd(N_rho, N_theta, n_quad=2,
                                                  verbose=True)
    print(f"  Total assembly: {time.time()-t0:.1f}s")
    print(f"  K_lo/K_hi: "
          f"{float(np.max(np.abs(K_lo)/(np.abs(K_hi)+1e-30))):.2e}")

    # mpmath eigh + Cauer
    from dd_full_pipeline import dd_to_mpmath_matrix, mpmath_generalized_eigh
    K_mp = dd_to_mpmath_matrix(K_hi, K_lo, dps=50)
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

    tau_list = []; g2_list = []
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
        tau_list.append(tau_k); g2_list.append(g2_k)
    pairs = sorted(zip(tau_list, g2_list), key=lambda p: -float(abs(p[1])))
    print(f"\n  Top 5 Foster:")
    for k in range(5):
        print(f"    rank={k}: tau={float(pairs[k][0])*1e6:.4f} us, "
              f"g²={float(pairs[k][1]):.4e}")

    # Kameari + verified-interval
    n_use = min(80, len(pairs))
    n_moments = 40
    alphas_mp = []
    for n_idx in range(n_moments):
        a_val = mp.mpf(0)
        for tau_k, g2_k in pairs[:n_use]:
            a_val += g2_k * (tau_k ** n_idx)
        alphas_mp.append((mp.mpf(-1)) ** n_idx * a_val)

    mp.iv.dps = 80
    eps_iv = mp.iv.mpf((1 - 1e-30, 1 + 1e-30))
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
        mid = (mp.mpf(rung.a) + mp.mpf(rung.b)) / 2
        rel_w = float((mp.mpf(rung.b) - mp.mpf(rung.a)) / abs(mid)) if abs(mid) > 0 else float("inf")
        rungs_raw.append({"stage": stage, "mid": float(mid),
                          "rel_width": rel_w})
        if rel_w > 1.0:
            break
        c_list = e[1:]

    print(f"\n  Cauer rungs (DD sphere axisym):")
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
                            "tau_pair_us": tau_p})

    out_path = Path(__file__).parent / f"dd_sphere_axisym_Nrho{N_rho}_Nth{N_theta}_results.json"
    out_path.write_text(json.dumps({
        "shape": "sphere", "N_rho": N_rho, "N_theta": N_theta,
        "Cauer_rungs": cauer_rungs}, indent=2), encoding="utf-8")
    print(f"\nSaved: {out_path}")


if __name__ == "__main__":
    main()
