"""BEM-Foster Cauer extraction with high moment count (50-digit mpmath).

Reuses the existing BEM eigenvalues + Foster amplitudes from
`bem_disk_axisym_cauer.json` (Mathematica output, 50 modes), but recomputes
the alpha_n moments and the Cauer ladder at arbitrary precision and arbitrary
moment count -- without rebuilding K_sym (the slow step).

Goal: demonstrate BEM-Foster *graceful degradation* by extending the Cauer
ladder past the 9-stage limit set by NMOMENTS=20 in the .wls file. Comparison
endpoint with FE-Hiruma in `breakdown_study.py`.

Output: scripts/bem_cauer_high_moments_results.json
"""
from __future__ import annotations

import json
from pathlib import Path

import mpmath as mp

mp.mp.dps = 60

HERE = Path(__file__).resolve().parent
BEM_JSON = Path("W:/30_CauerLadderNetwork/2026_04_01_長方形CLN/ngsolve_validation/"
                "bem_disk_axisym_cauer.json")


def compute_alphas_hp(eigvals_lam, gampl, sigma, n_moments):
    """alpha_n = (-1)^n * Sum_k g_k^2 * tau_k^(n+1),  tau_k = sigma * lambda_k.

    All arithmetic at mp.dps precision. eigvals_lam are floats from JSON
    (machine precision ~16 digits), but we lift them to mpf so that
    tau_k^(n+1) for large n does not underflow IEEE doubles.
    """
    lam = [mp.mpf(repr(l)) for l in eigvals_lam]
    g   = [mp.mpf(repr(x)) for x in gampl]
    sig = mp.mpf(repr(sigma))
    tau = [sig * l for l in lam]

    alphas = []
    for n in range(n_moments):
        s = mp.mpf(0)
        for k in range(len(tau)):
            s += g[k]**2 * tau[k]**(n + 1)
        if n % 2 == 1:
            s = -s
        alphas.append(s)
    return alphas


def invert_taylor(c):
    n = len(c)
    e = [mp.mpf(0)] * n
    e[0] = 1 / c[0]
    for k in range(1, n):
        e[k] = -mp.fsum(c[j] * e[k - j] for j in range(1, k + 1)) / c[0]
    return e


def cauer_extract(c_taylor, max_stages):
    c = list(c_taylor)
    p_list = []
    for _ in range(max_stages):
        if len(c) < 2 or abs(c[0]) < mp.mpf("1e-300"):
            break
        e = invert_taylor(c)
        p_list.append(e[0])
        c = e[1:]
    return p_list


def main():
    bem = json.loads(BEM_JSON.read_text(encoding="utf-8"))
    sigma = bem["params"]["sigma"]
    eigvals = bem["bem_eigvals_lambda"]
    gampl = bem["foster_amplitudes_g"]
    print(f"BEM JSON: {BEM_JSON.name}")
    print(f"  ne = {bem['params']['ne']}, n_modes = {len(eigvals)}, "
          f"sigma = {sigma}")
    print(f"  mp.dps = {mp.mp.dps}")
    print()

    n_moments = 60
    print(f"Recomputing alpha_n at {mp.mp.dps}-digit precision, "
          f"n_moments = {n_moments} ...")
    alphas = compute_alphas_hp(eigvals, gampl, sigma, n_moments)
    print(f"  alpha_0  = {mp.nstr(alphas[0], 10)}")
    print(f"  alpha_1  = {mp.nstr(alphas[1], 10)}")
    print(f"  alpha_19 = {mp.nstr(alphas[19], 10)}")
    print(f"  alpha_59 = {mp.nstr(alphas[-1], 10)}")
    print()

    print("Cauer extraction ...")
    p_cauer = cauer_extract(alphas, n_moments - 1)
    n_pairs = len(p_cauer) // 2
    print(f"  Extracted {len(p_cauer)} p-coefficients => {n_pairs} pairs")
    print()

    bem_R_2k, bem_L_2k_plus_1, bem_tau_pair_us = [], [], []
    print(f"{'k':>3} {'R_{2k}':>16} {'L_{2k+1}':>16} {'tau_pair (us)':>18}")
    for k in range(n_pairs):
        R = p_cauer[2*k]
        Lin = p_cauer[2*k + 1]
        L = (1 / Lin) if abs(Lin) > mp.mpf("1e-300") else mp.mpf("nan")
        if abs(R) > mp.mpf("1e-300") and mp.isfinite(L):
            tau_us = (L / R) * mp.mpf("1e6")
        else:
            tau_us = mp.mpf("nan")
        bem_R_2k.append(R)
        bem_L_2k_plus_1.append(L)
        bem_tau_pair_us.append(tau_us)
        print(f"{k:>3} {mp.nstr(R, 8):>16} {mp.nstr(L, 8):>16} "
              f"{mp.nstr(tau_us, 8):>18}")

    print()

    # ---- Breakdown / sign-flip detection ----
    def first_sign_flip(taus):
        for k, t in enumerate(taus):
            if not mp.isfinite(t) or t < 0:
                return k
        return None

    def first_jump(taus, ratio=2.0):
        """First k where |tau[k] - tau[k-1]| / tau[k-1] is large positive jump."""
        for k in range(1, len(taus)):
            prev = taus[k-1]
            curr = taus[k]
            if not (mp.isfinite(prev) and mp.isfinite(curr) and prev > 0):
                return k
            if curr > prev * ratio:
                return k
        return None

    k_sign = first_sign_flip(bem_tau_pair_us)
    k_jump = first_jump(bem_tau_pair_us, 2.0)
    print(f"BEM first sign flip / NaN  : k = {k_sign}")
    print(f"BEM first 2x upward jump   : k = {k_jump}")

    # ---- Save ----
    out = {
        "config": {
            "n_modes": len(eigvals),
            "n_moments": n_moments,
            "mp_dps": mp.mp.dps,
            "sigma": sigma,
        },
        "bem_R_2k":         [mp.nstr(r, 30) for r in bem_R_2k],
        "bem_L_2k_plus_1":  [mp.nstr(l, 30) for l in bem_L_2k_plus_1],
        "tau_pair_us":      [mp.nstr(t, 30) for t in bem_tau_pair_us],
        "tau_pair_us_float": [float(t) if mp.isfinite(t) else None
                              for t in bem_tau_pair_us],
        "alphas":           [mp.nstr(a, 30) for a in alphas],
        "breakdown_k": {"sign_flip": k_sign, "upward_jump": k_jump},
    }
    out_path = HERE / "bem_cauer_high_moments_results.json"
    out_path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"\nSaved: {out_path}")


if __name__ == "__main__":
    main()
