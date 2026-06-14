"""dd_production_A1.py — DD Cauer extraction for A1 square prism (17.72×17.72×2 mm).

Same as dd_production_cuboid521.py but with A1 geometry (V=628 mm³, matched
to cylinder R=10mm t=2mm volume).
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
from dd_k_assembly_gpu import assemble_K_dd_gpu
from dd_M_b_assembly import assemble_M_b_dd_gpu
from dd_full_pipeline import (
    dd_to_mpmath_matrix, dd_to_mpmath_vector, mpmath_generalized_eigh,
)


SIGMA = 5.8e7
MU0 = 4 * math.pi * 1e-7


def main():
    print("=" * 72)
    print("Production DD Cauer extraction: A1 (17.72×17.72×2 mm)")
    print("=" * 72)

    p = 3
    A_M, B_M, C_M = 17.72e-3, 17.72e-3, 2e-3
    n_th, n_ph, n_rh, n_c = 12, 24, 12, 6
    print(f"\n  order={p}, n_dofs={2*p**3+3*p**2}")
    print(f"  Production quad: ({n_th}, {n_ph}, {n_rh}, {n_c})")

    print(f"\n--- DD K assembly (GPU) ---")
    t0 = time.time()
    K_bare_hi, K_bare_lo = assemble_K_dd_gpu(p, A_M, B_M, C_M,
                                               n_th, n_ph, n_rh, n_c,
                                               use_gpu=True, verbose=True)
    print(f"  Total: {time.time()-t0:.1f}s")

    pref = MU0 / (4 * math.pi)
    K_hi = K_bare_hi * pref
    K_lo = K_bare_lo * pref

    print(f"\n--- DD M, b assembly ---")
    t0 = time.time()
    M_hi, M_lo, b_hi, b_lo = assemble_M_b_dd_gpu(p, A_M, B_M, C_M,
                                                    n_gauss=8, use_gpu=True,
                                                    verbose=False)
    print(f"  Total: {time.time()-t0:.1f}s")
    print(f"  K_lo/K_hi: {float(np.max(np.abs(K_lo)/(np.abs(K_hi)+1e-30))):.2e}")
    print(f"  M_lo/M_hi: {float(np.max(np.abs(M_lo)/(np.abs(M_hi)+1e-30))):.2e}")

    print(f"\n--- mpmath eigh ---")
    t0 = time.time()
    K_mp = dd_to_mpmath_matrix(K_hi, K_lo, dps=50)
    M_mp = dd_to_mpmath_matrix(M_hi, M_lo, dps=50)
    b_mp = dd_to_mpmath_vector(b_hi, b_lo, dps=50)
    eigvals_mp, eigvecs_mp, L_chol = mpmath_generalized_eigh(K_mp, M_mp, dps=50)
    print(f"  Total: {time.time()-t0:.1f}s")

    # Foster spectrum
    n = len(eigvals_mp)
    tau_list = []
    g2_list = []
    for k in range(n):
        lam = eigvals_mp[k]
        if abs(lam) < mp.mpf("1e-50"):
            continue
        tau_k = mp.mpf(SIGMA) * lam
        g_k = mp.mpf(0)
        for i in range(n):
            g_k += eigvecs_mp[i, k] * b_mp[i, 0]
        vMv = mp.mpf(0)
        for i in range(n):
            for j in range(n):
                vMv += eigvecs_mp[i, k] * M_mp[i, j] * eigvecs_mp[j, k]
        g2_k = g_k * g_k / vMv if vMv > 0 else mp.mpf(0)
        tau_list.append(tau_k)
        g2_list.append(g2_k)

    pairs = sorted(zip(tau_list, g2_list), key=lambda p: -float(abs(p[1])))
    print(f"\n  Top 5 Foster modes:")
    for k in range(5):
        print(f"    rank={k}: tau={float(pairs[k][0])*1e6:.4f} us, "
              f"g²={float(pairs[k][1]):.4e}")

    # Kameari moments + verified-interval Hankel-Padé
    n_use = min(60, len(pairs))
    n_moments = 40
    alphas_mp = []
    for n_idx in range(n_moments):
        a_val = mp.mpf(0)
        for tau_k, g2_k in pairs[:n_use]:
            a_val += g2_k * (tau_k ** n_idx)
        alphas_mp.append((mp.mpf(-1)) ** n_idx * a_val)

    print(f"\n--- Verified-interval Hankel-Padé (DD input) ---")
    mp.iv.dps = 80
    eps_rel = mp.mpf("1e-30")
    eps_iv = mp.iv.mpf((float(1 - eps_rel), float(1 + eps_rel)))
    alphas_iv = [mp.iv.mpf(float(a)) * eps_iv for a in alphas_mp]

    c_list = [a for a in alphas_iv]
    rungs_raw = []
    width_warn = 1e-2
    for stage in range(20):
        if len(c_list) < 2:
            break
        if 0 in c_list[0]:
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
        rel = float(width / abs(mid)) if abs(mid) > 0 else float("inf")
        rungs_raw.append({"stage": stage, "mid": float(mid),
                          "low": float(rung_a), "high": float(rung_b),
                          "rel_width": rel,
                          "reliable": rel < width_warn})
        if rel > 1.0:
            break
        c_list = e[1:]

    print(f"\n  Cauer rungs (DD input, verified-interval):")
    print(f"  {'k':>3} {'R_2k':>15} {'L_2k+1':>15} {'tau_pair us':>13} "
          f"{'reliable?':>10}")
    cauer_rungs = []
    for k in range(len(rungs_raw) // 2):
        R_raw = rungs_raw[2 * k]
        Linv_raw = rungs_raw[2 * k + 1]
        if abs(Linv_raw["mid"]) < 1e-50:
            break
        L_mid = 1.0 / Linv_raw["mid"]
        tau_p = L_mid / R_raw["mid"] * 1e6
        reliable = R_raw["reliable"] and Linv_raw["reliable"]
        print(f"  {k:>3} {R_raw['mid']:>15.4e} {L_mid:>15.4e} "
              f"{tau_p:>13.4f} {str(reliable):>10}")
        cauer_rungs.append({
            "k": k, "R_2k": R_raw["mid"], "L_2k_plus_1": L_mid,
            "tau_pair_us": tau_p,
            "R_2k_rel_width": R_raw["rel_width"],
            "L_2k_plus_1_inv_rel_width": Linv_raw["rel_width"],
            "reliable": reliable,
        })

    out_path = (Path(__file__).parent
                / "dd_production_A1_order3_results.json")
    out_path.write_text(json.dumps({
        "shape": "A1_17.72x17.72x2",
        "p": p, "n_dofs": K_hi.shape[0],
        "quad": [n_th, n_ph, n_rh, n_c],
        "top_foster": [
            {"rank": k, "tau_us": float(pairs[k][0]) * 1e6,
             "g2": float(pairs[k][1])}
            for k in range(min(20, len(pairs)))
        ],
        "Cauer_rungs_verified": cauer_rungs,
    }, indent=2), encoding="utf-8")
    print(f"\nSaved: {out_path}")


if __name__ == "__main__":
    main()
