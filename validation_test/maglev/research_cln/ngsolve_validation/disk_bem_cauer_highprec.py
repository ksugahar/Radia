"""disk_bem_cauer_highprec.py -- Cauer-I CFE on the high-precision BEM moments.

Reads `bem_disk_axisym_cauer_highprec.json` (n_moments=40, WP=80) and applies
classical Cauer extraction at mpmath dps=120 to push the negative-tau onset
beyond the standard run's k=6 limit.

Output: bem_disk_axisym_cauer_highprec_python_results.json (same schema as
        the standard pipeline output, just longer p_cauer/R/L arrays).

Track B-ii: the standard 20-moment pipeline produces tau_pair[k=6] = -43 us
            (negative, Cauer extraction has gone unstable). This script tests
            whether higher precision pushes the breakdown further out.
"""
from __future__ import annotations

import json
from pathlib import Path

import mpmath as mp

mp.mp.dps = 120

REPO = Path(__file__).resolve().parents[0]
BEM_JSON = REPO / "bem_disk_axisym_cauer_highprec.json"


def invert_taylor(c):
    n = len(c)
    e = [mp.mpf(0)] * n
    e[0] = 1 / c[0]
    for k in range(1, n):
        sum_term = mp.fsum(c[j] * e[k - j] for j in range(1, k + 1))
        e[k] = -sum_term / c[0]
    return e


def cauer_extract(c_taylor, max_stages):
    c = list(c_taylor)
    p_list = []
    for _ in range(max_stages):
        if len(c) < 2 or abs(c[0]) < mp.mpf("1e-1000"):
            break
        e = invert_taylor(c)
        p_list.append(e[0])
        c = e[1:]
    return p_list


def main():
    print("=" * 76)
    print("Track B-ii: BEM Cauer-I extraction at HIGH PRECISION")
    print(f"            mpmath dps={mp.mp.dps}, BEM moments WP=80, n=40")
    print("=" * 76)
    print()

    bem = json.loads(BEM_JSON.read_text(encoding="utf-8"))
    print(f"BEM: {bem['method']}")
    p = bem["params"]
    print(f"  ne = {p['ne']}, n_modes = {p['n_modes']}, "
          f"n_moments = {p['n_moments']}, WP = {p.get('wp_digits', '?')}")

    # alphas are stored as 30-digit strings (from N[..., 30]); load to mpf.
    alphas = [mp.mpf(repr(a)) for a in bem["alphas"]]
    bem_taus = [mp.mpf(repr(t)) for t in bem["bem_tau_us"]]

    print(f"  BEM Foster tau_1..tau_8 (us) = "
          f"{[float(t) for t in bem_taus[:8]]}")
    print(f"  alpha_0 = {mp.nstr(alphas[0], 6)}")
    print(f"  alpha_1 = {mp.nstr(alphas[1], 6)}")
    print(f"  alpha_2 = {mp.nstr(alphas[2], 6)}")
    print()

    tau1_moment = -alphas[1] / alphas[0]
    print(f"Leading-mode estimate tau_1 = -alpha_1/alpha_0 = "
          f"{mp.nstr(tau1_moment * 1e6, 12)} us")
    print(f"  vs BEM Foster tau_1: {mp.nstr(bem_taus[0], 12)} us")
    print(f"  rel err: {mp.nstr(abs((tau1_moment - bem_taus[0]/mp.mpf('1e6'))/(bem_taus[0]/mp.mpf('1e6'))), 6)}")
    print()

    print(f"Cauer extraction on {len(alphas)}-term Taylor series at dps={mp.mp.dps}...")
    p_cauer = cauer_extract(alphas, len(alphas) - 1)
    print(f"  Extracted {len(p_cauer)} p-coefficients")
    print()

    n_pairs = len(p_cauer) // 2
    bem_R_2k, bem_L_2k_plus_1, bem_tau_pair = [], [], []
    print(f"{'k':>3} {'R_{2k}':>22} {'L_{2k+1}':>22} {'tau_pair (us)':>22}  flag")
    print("-" * 88)
    for k in range(n_pairs):
        R   = p_cauer[2*k]
        Lin = p_cauer[2*k + 1]
        L   = 1 / Lin if abs(Lin) > mp.mpf("1e-1000") else mp.mpf("nan")
        tau = L / R if abs(R) > mp.mpf("1e-1000") else mp.mpf("nan")
        bem_R_2k.append(R)
        bem_L_2k_plus_1.append(L)
        bem_tau_pair.append(tau)
        flag = ""
        tau_us = float(tau * 1e6) if tau == tau else float("nan")
        if tau_us < 0:
            flag = "  <-- NEGATIVE"
        elif k > 0:
            prev_us = float(bem_tau_pair[k-1] * 1e6)
            if prev_us > 0 and tau_us > prev_us:
                flag = "  <-- INCREASE (decay should be monotone)"
        print(f"{k:>3} {mp.nstr(R, 12):>22} {mp.nstr(L, 12):>22} "
              f"{mp.nstr(tau * 1e6, 12):>22}{flag}")
    print()

    n_neg = sum(1 for t in bem_tau_pair if t < 0)
    first_neg = next((k for k, t in enumerate(bem_tau_pair) if t < 0), None)
    print(f"Negative-tau stages: {n_neg}")
    print(f"First negative at k = {first_neg}" if first_neg is not None
          else "No negative-tau in extracted range")
    print()

    out = {
        "method": ("BEM Foster -> Cauer ladder, HIGH PRECISION variant "
                   f"(WP=80, mpmath dps={mp.mp.dps}, n_moments=40)"),
        "n_moments": len(alphas),
        "n_pairs": n_pairs,
        "first_negative_k": first_neg,
        "p_cauer":         [mp.nstr(p, 30) for p in p_cauer],
        "R_2k":            [mp.nstr(r, 30) for r in bem_R_2k],
        "L_2k_plus_1":     [mp.nstr(l, 30) for l in bem_L_2k_plus_1],
        "tau_pair_us":     [mp.nstr(t * 1e6, 30) for t in bem_tau_pair],
    }
    out_path = REPO / "bem_disk_axisym_cauer_highprec_python_results.json"
    out_path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
