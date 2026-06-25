"""Cuboid 5x2x1 mm BEM HDiv order-3 -> Cauer extraction at high moments.

Companion to `bem_cauer_high_moments.py` (disk axisym). Reads the JSON that
`hex_bem_hdiv_order3_export_json.wls` produces from the overnight .mx file
(81 DOFs, 6 physically-coupled modes) and runs 60-moment Cauer extraction in
50-digit mpmath.

Output: scripts/hex_bem_cauer_high_moments_results.json
"""
from __future__ import annotations

import json
from pathlib import Path

import mpmath as mp

mp.mp.dps = 60

HERE = Path(__file__).resolve().parent
BEM_JSON = Path("W:/30_CauerLadderNetwork/2026_04_01_長方形CLN/"
                "ngsolve_validation/hex_bem_hdiv_order3_export.json")


def compute_alphas_hp(eigvals_lam, gampl, sigma, n_moments):
    """alpha_n = (-1)^n * Sum_k g_k^2 * tau_k^(n+1),  tau_k = sigma * lambda_k."""
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
    print(f"BEM JSON  : {BEM_JSON.name}")
    print(f"  cuboid   : a={bem['params']['a_mm']} mm  b={bem['params']['b_mm']} mm  "
          f"c={bem['params']['c_mm']} mm")
    print(f"  nDOF     : {bem['params']['nDOF']}")
    print(f"  n_modes  : {len(eigvals)}  (sortedPhysical filtered)")
    print(f"  mp.dps   : {mp.mp.dps}")
    print()
    print(f"  Foster spectrum (rank by tau desc):")
    print(f"  {'k':>3} {'tau (us)':>13} {'g':>16} {'g^2 * tau':>16}")
    for k, (tau, g) in enumerate(zip(bem["bem_tau_us"], gampl)):
        g2t = float(g)**2 * float(tau)
        print(f"  {k:>3} {float(tau):>13.4f} {float(g):>+16.4e} {g2t:>16.4e}")
    print()

    n_moments = 60
    print(f"Computing alpha_n at {mp.mp.dps}-digit precision, "
          f"n_moments = {n_moments} ...")
    alphas = compute_alphas_hp(eigvals, gampl, sigma, n_moments)
    tau_lead_estimate = -alphas[1] / alphas[0]
    print(f"  alpha_0  = {mp.nstr(alphas[0], 10)}")
    print(f"  alpha_1  = {mp.nstr(alphas[1], 10)}")
    print(f"  alpha_19 = {mp.nstr(alphas[19], 10)}")
    print(f"  alpha_59 = {mp.nstr(alphas[-1], 10)}")
    print(f"  -alpha_1/alpha_0 = {mp.nstr(tau_lead_estimate * 1e6, 10)} us  "
          f"(weighted-mean tau)")
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

    def first_sign_flip(taus):
        for k, t in enumerate(taus):
            if not mp.isfinite(t) or t < 0:
                return k
        return None

    k_sign = first_sign_flip(bem_tau_pair_us)
    print(f"BEM cuboid first sign flip / NaN  : k = {k_sign}")
    print()

    # -- Cross-check: tau_pair[0] should equal max tau among physical modes
    max_tau_us = max(float(t) for t in bem["bem_tau_us"])
    print(f"Cross-check: max tau across physical modes = {max_tau_us:.4f} us")
    print(f"             tau_pair[0] from Cauer extr.  = "
          f"{float(bem_tau_pair_us[0]):.4f} us")
    print()

    out = {
        "config": {
            "n_dof": bem["params"]["nDOF"],
            "n_modes_physical": len(eigvals),
            "n_moments": n_moments,
            "mp_dps": mp.mp.dps,
            "sigma": sigma,
            "geometry_mm": [bem["params"]["a_mm"], bem["params"]["b_mm"],
                            bem["params"]["c_mm"]],
        },
        "bem_R_2k":         [mp.nstr(r, 30) for r in bem_R_2k],
        "bem_L_2k_plus_1":  [mp.nstr(l, 30) for l in bem_L_2k_plus_1],
        "tau_pair_us":      [mp.nstr(t, 30) for t in bem_tau_pair_us],
        "tau_pair_us_float": [float(t) if mp.isfinite(t) else None
                              for t in bem_tau_pair_us],
        "alphas":           [mp.nstr(a, 30) for a in alphas],
        "physical_modes": {
            "tau_us": [float(t) for t in bem["bem_tau_us"]],
            "g":      [float(g) for g in gampl],
        },
        "breakdown_k": {"sign_flip": k_sign},
    }
    out_path = HERE / "hex_bem_cauer_high_moments_results.json"
    out_path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
