"""dd_sphere_axisym_mp.py — DD axisym VIM sphere, multiprocessing K assembly.

Speedup over serial dd_sphere_axisym.py via Pool of N workers (default cpu_count).
K matrix is partitioned by qi chunks; each worker accumulates a local DD K_hi/K_lo
sub-sum, then main reduces via dd_add (preserves DD precision).

For production scale:
  - Nρ=20, Nθ=27 (540 cells, 2160 q-points, 4.67M pairs): ~20 min K (8 cores) + ~2 h eigh
  - Nρ=30, Nθ=36 (1080 cells, 4320 q-points, 18.7M pairs): ~80 min K + ~15 h eigh
"""
from __future__ import annotations

import argparse
import json
import math
import multiprocessing
import os
import sys
import time
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
import mpmath as mp

sys.path.insert(0, str(Path(__file__).parent))
from dd_arithmetic import dd_add, dd_from_mpmath
from dd_axisym_kernel import axisym_kernel_dd_mpmath


R_SPHERE = 10e-3
SIGMA = 5.8e7
MU0 = 4 * math.pi * 1e-7
B0 = 1.0


# Worker globals (set by initializer)
_W = {}


def _worker_init(r_hi, r_lo, z_hi, z_lo, jac_hi, jac_lo, cell_of_q, N, dps):
    _W['r_hi'] = r_hi
    _W['r_lo'] = r_lo
    _W['z_hi'] = z_hi
    _W['z_lo'] = z_lo
    _W['jac_hi'] = jac_hi
    _W['jac_lo'] = jac_lo
    _W['cell_of_q'] = cell_of_q
    _W['N'] = N
    _W['dps'] = dps
    mp.mp.dps = dps


def _compute_chunk(qi_range):
    """Worker: accumulate K contributions for qi in qi_range."""
    r_hi = _W['r_hi']
    r_lo = _W['r_lo']
    z_hi = _W['z_hi']
    z_lo = _W['z_lo']
    jac_hi = _W['jac_hi']
    jac_lo = _W['jac_lo']
    cell_of_q = _W['cell_of_q']
    N = _W['N']
    dps = _W['dps']
    mp.mp.dps = dps

    n_q = len(r_hi)
    K_hi = np.zeros((N, N))
    K_lo = np.zeros((N, N))
    two_pi_mp = mp.mpf(2) * mp.pi

    for qi in qi_range:
        if r_hi[qi] < 1e-15:
            continue
        ci = int(cell_of_q[qi])
        r_qi_mp = mp.mpf(float(r_hi[qi])) + mp.mpf(float(r_lo[qi]))
        jac_qi_mp = mp.mpf(float(jac_hi[qi])) + mp.mpf(float(jac_lo[qi]))
        w_obs_mp = two_pi_mp * r_qi_mp * jac_qi_mp

        for qj in range(n_q):
            if r_hi[qj] < 1e-15:
                continue
            cj = int(cell_of_q[qj])
            g_hi, g_lo = axisym_kernel_dd_mpmath(
                np.array(r_hi[qi]), np.array(r_lo[qi]),
                np.array(z_hi[qi]), np.array(z_lo[qi]),
                np.array(r_hi[qj]), np.array(r_lo[qj]),
                np.array(z_hi[qj]), np.array(z_lo[qj]),
                dps=dps)
            g_mp = mp.mpf(float(g_hi)) + mp.mpf(float(g_lo))
            jac_qj_mp = mp.mpf(float(jac_hi[qj])) + mp.mpf(float(jac_lo[qj]))
            contrib_mp = g_mp * w_obs_mp * jac_qj_mp
            contrib_hi, contrib_lo = dd_from_mpmath(contrib_mp)
            new_hi, new_lo = dd_add(
                np.array(K_hi[ci, cj]), np.array(K_lo[ci, cj]),
                np.array(contrib_hi), np.array(contrib_lo))
            K_hi[ci, cj] = float(new_hi)
            K_lo[ci, cj] = float(new_lo)

    return K_hi, K_lo


def assemble_K_M_b_sphere_dd_mp(N_rho, N_theta, n_quad=2, n_workers=None,
                                  dps=40, verbose=True):
    if n_workers is None:
        n_workers = max(1, multiprocessing.cpu_count() - 1)
    mp.mp.dps = dps
    N = N_rho * N_theta
    drho = R_SPHERE / N_rho
    dtheta = math.pi / N_theta

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

    from dd_gl_nodes import gauss_legendre_mpf
    nodes_mp, weights_mp = gauss_legendre_mpf(n_quad, dps=dps)
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
                    r_q_mp = rho_q_mp * mp.sin(th_q_mp)
                    z_q_mp = rho_q_mp * mp.cos(th_q_mp)
                    jac_q_mp = rho_q_mp * wrho_mp * wth_mp
                    r_flat_hi[q], r_flat_lo[q] = dd_from_mpmath(r_q_mp)
                    z_flat_hi[q], z_flat_lo[q] = dd_from_mpmath(z_q_mp)
                    jac_flat_hi[q], jac_flat_lo[q] = dd_from_mpmath(jac_q_mp)
                    cell_of_q[q] = cell_idx
                    q += 1

    if verbose:
        print(f"  N={N} cells, n_q_total={n_q_total}, "
              f"n_pairs={n_q_total ** 2}")
        print(f"  Parallelizing over {n_workers} workers...")
        t0 = time.time()

    # Partition qi indices into chunks (one chunk per worker, balanced).
    chunks = [list(range(i, n_q_total, n_workers)) for i in range(n_workers)]
    # Round-robin gives better load balance than contiguous when r_hi < 1e-15
    # cells are uneven.

    init_args = (r_flat_hi, r_flat_lo, z_flat_hi, z_flat_lo,
                 jac_flat_hi, jac_flat_lo, cell_of_q, N, dps)

    K_hi = np.zeros((N, N))
    K_lo = np.zeros((N, N))

    with multiprocessing.Pool(processes=n_workers,
                              initializer=_worker_init,
                              initargs=init_args) as pool:
        # imap_unordered for progress feedback
        for w_idx, (K_h, K_l) in enumerate(
                pool.imap_unordered(_compute_chunk, chunks)):
            K_hi, K_lo = dd_add(K_hi, K_lo, K_h, K_l)
            if verbose:
                print(f"    worker {w_idx + 1}/{n_workers} done "
                      f"({time.time() - t0:.1f}s)")

    mu0_over_2 = MU0 / 2.0
    K_hi = K_hi * mu0_over_2
    K_lo = K_lo * mu0_over_2

    if verbose:
        print(f"  K assembly (MP): {time.time() - t0:.1f}s")

    return K_hi, K_lo, np.diag(M_diag), b_vec


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--Nrho", type=int, default=15)
    parser.add_argument("--Ntheta", type=int, default=18)
    parser.add_argument("--nquad", type=int, default=2)
    parser.add_argument("--workers", type=int, default=None,
                        help="default: cpu_count() - 1")
    parser.add_argument("--eigh_dps", type=int, default=35,
                        help="dps for mpmath eigh (default 35, was 50)")
    args = parser.parse_args()

    print("=" * 72)
    print("DD axisym VIM sphere (multiprocessing)")
    print("=" * 72)
    print(f"\n  N_rho={args.Nrho}, N_theta={args.Ntheta}, "
          f"total {args.Nrho * args.Ntheta} cells")
    print(f"  Workers: {args.workers or (multiprocessing.cpu_count() - 1)}, "
          f"n_quad={args.nquad}, eigh_dps={args.eigh_dps}")

    t0 = time.time()
    K_hi, K_lo, M, b = assemble_K_M_b_sphere_dd_mp(
        args.Nrho, args.Ntheta, n_quad=args.nquad,
        n_workers=args.workers, dps=40, verbose=True)
    print(f"  Total assembly: {time.time() - t0:.1f}s")
    print(f"  K_lo/K_hi: "
          f"{float(np.max(np.abs(K_lo) / (np.abs(K_hi) + 1e-30))):.2e}")

    from dd_full_pipeline import dd_to_mpmath_matrix, mpmath_generalized_eigh
    K_mp = dd_to_mpmath_matrix(K_hi, K_lo, dps=args.eigh_dps)
    mp.mp.dps = args.eigh_dps
    N = K_hi.shape[0]
    M_mp = mp.matrix(N, N)
    for i in range(N):
        for j in range(N):
            M_mp[i, j] = mp.mpf(float(M[i, j]))
    b_mp = mp.matrix(N, 1)
    for i in range(N):
        b_mp[i, 0] = mp.mpf(float(b[i]))

    t0 = time.time()
    eigvals_mp, eigvecs_mp, L_chol = mpmath_generalized_eigh(
        K_mp, M_mp, dps=args.eigh_dps)
    print(f"  mpmath eigh: {time.time() - t0:.1f}s")

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
        print(f"    rank={k}: tau={float(pairs[k][0]) * 1e6:.4f} us, "
              f"g²={float(pairs[k][1]):.4e}")

    n_use = min(120, len(pairs))
    n_moments = 60
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
    for stage in range(40):
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
        rel_w = (float((mp.mpf(rung.b) - mp.mpf(rung.a)) / abs(mid))
                 if abs(mid) > 0 else float("inf"))
        rungs_raw.append({"stage": stage, "mid": float(mid),
                          "rel_width": rel_w})
        if rel_w > 1.0:
            break
        c_list = e[1:]

    print(f"\n  Cauer rungs (DD sphere axisym, MP):")
    print(f"  {'k':>3} {'R_2k':>15} {'L_2k+1':>15} "
          f"{'tau_pair us':>13} {'rel_w':>10}")
    cauer_rungs = []
    for k in range(len(rungs_raw) // 2):
        R_raw = rungs_raw[2 * k]
        Linv_raw = rungs_raw[2 * k + 1]
        if abs(Linv_raw["mid"]) < 1e-50:
            break
        L_mid = 1.0 / Linv_raw["mid"]
        tau_p = L_mid / R_raw["mid"] * 1e6
        rw_max = max(R_raw["rel_width"], Linv_raw["rel_width"])
        print(f"  {k:>3} {R_raw['mid']:>15.4e} {L_mid:>15.4e} "
              f"{tau_p:>13.4f} {rw_max:>10.2e}")
        cauer_rungs.append({"k": k, "R_2k": R_raw["mid"],
                            "L_2k_plus_1": L_mid,
                            "tau_pair_us": tau_p,
                            "R_2k_rel_width": R_raw["rel_width"],
                            "L_inv_rel_width": Linv_raw["rel_width"]})

    out_path = (Path(__file__).parent
                / f"dd_sphere_axisym_mp_Nrho{args.Nrho}_Nth{args.Ntheta}_"
                f"nq{args.nquad}_results.json")
    out_path.write_text(json.dumps({
        "shape": "sphere", "N_rho": args.Nrho, "N_theta": args.Ntheta,
        "n_quad": args.nquad, "workers": args.workers,
        "eigh_dps": args.eigh_dps,
        "Cauer_rungs": cauer_rungs}, indent=2), encoding="utf-8")
    print(f"\nSaved: {out_path}")


if __name__ == "__main__":
    main()
