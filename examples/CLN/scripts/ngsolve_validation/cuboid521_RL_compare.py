"""cuboid521_RL_compare.py — Compare R_2n, L_{2n+1}, tau_pair across methods
for cuboid 5x2x1 mm Cu, B0 = 1T applied (z direction).

Methods:
  1. radia-vim Phase F-4 BEM (Spherical Duffy, 81 modes) -> Kameari Cauer-I
  2. NGSolve cln_team28_kelvin (HCurl FEM + snapshot fix) -> Kameari direct

Both expressed in canonical Joule-loss normalization (R_0 = 1/integral(sigma|A_s|^2)).
Investigate where each method breaks down at high stages.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import mpmath as mp
mp.mp.dps = 250   # very high precision to push high-stage extraction

sys.stdout.reconfigure(encoding="utf-8")


# Canonical R_0 from Joule loss formula (B0 = 1T, sigma_Cu, cuboid Lx*Ly*Lz)
# R_0 = 1 / integral_cond(sigma |A_s|^2 dV)
# For uniform B_z, A_s = (B0/2)(-y, x, 0), |A_s|^2 = (B0/2)^2 (x^2+y^2)
# int (x^2+y^2) over cuboid = LxLyLz/12 * (Lx^2 + Ly^2)
def canonical_R0_cuboid(Lx, Ly, Lz, sigma, B0=1.0):
    integral = sigma * (B0/2)**2 * (Lx*Ly*Lz/12) * (Lx**2 + Ly**2)
    return 1.0 / integral


# Hankel-Pade Cauer-I extraction (mpmath, returns p-coefficients)
def hankel_pade_cauer(alphas, max_stages):
    c = [mp.mpf(a) for a in alphas]
    p_list = []
    for _ in range(max_stages):
        if len(c) < 2 or abs(c[0]) < mp.mpf(10) ** (-mp.mp.dps + 10):
            break
        n = len(c)
        e = [mp.mpf(0)] * n
        e[0] = 1 / c[0]
        for k in range(1, n):
            e[k] = -sum(c[j+1] * e[k-j-1] for j in range(k)
                        if k-j-1 >= 0) / c[0]
        p_list.append(e[0])
        c = e[1:]
    return p_list


def bem_cauer_kameari(json_path, n_modes=None, n_moments=120, max_pairs=40):
    with open(json_path) as f:
        d = json.load(f)
    modes = d["top_modes_by_g2"]
    if n_modes is not None:
        modes = modes[:n_modes]
    tau_list = [mp.mpf(e["tau_us"]) * mp.mpf("1e-6") for e in modes]
    g2_list = [mp.mpf(e["g2"]) for e in modes]

    # Kameari moments: alpha_n = (-1)^n sum g_k^2 tau_k^n
    alphas = []
    for n in range(n_moments):
        a = sum(g * t**n for g, t in zip(g2_list, tau_list))
        alphas.append((mp.mpf(-1))**n * a)
    p = hankel_pade_cauer(alphas, 2 * max_pairs)
    rungs = []
    for k in range(len(p) // 2):
        R = p[2*k]
        Linv = p[2*k+1]
        if abs(Linv) < mp.mpf(10)**(-50):
            break
        L = 1 / Linv
        tau_us = float(L/R) * 1e6
        rungs.append({"k": k, "R_2k": float(R), "L_2k_plus_1": float(L),
                      "tau_pair_us": tau_us})
    return rungs, len(modes), d.get("leading_tau_us", None)


def main():
    Lx, Ly, Lz = 5e-3, 2e-3, 1e-3
    sigma = 5.8e7
    R0_can = canonical_R0_cuboid(Lx, Ly, Lz, sigma)
    print("=" * 78)
    print(" Cuboid 5x2x1 mm Cu, B_0 = 1 T applied (z) — R_2n, L_{2n+1} compare")
    print("=" * 78)
    print(f"  Canonical R_0 = 1/integral(sigma|A_s|^2) = {R0_can:.4e} Ohm")
    print(f"    (analytical Joule-loss inverse, sigma = {sigma:.2e} S/m)")
    print()

    # --- BEM Kameari extraction with progressively more modes ---
    base = Path(__file__).parent
    candidates = [
        ("F-4 top 20 modes", base / "cuboid521_radia_vim_F4_foster.json"),
    ]
    full = base / "cuboid521_radia_vim_F4_full81.json"
    if full.exists():
        candidates.append(("F-4 all 81 modes", full))

    bem_results = {}
    for label, path in candidates:
        if not path.exists():
            continue
        rungs, n_modes, leading = bem_cauer_kameari(path)
        bem_results[label] = rungs
        print(f"--- {label} (n_modes={n_modes}, leading Foster pole={leading} us)"
              f" -> Cauer rungs: ---")
        print(f"  {'k':>3} | {'R_2k [Ohm]':>14} | {'L_{2k+1} [H]':>14}"
              f" | {'tau_pair [us]':>14}")
        for r in rungs[:12]:
            print(f"  {r['k']:>3} | {r['R_2k']:>14.6e}"
                  f" | {r['L_2k_plus_1']:>14.6e} | {r['tau_pair_us']:>14.6f}")
        print()

    # --- NGSolve cln_team28 results (Kameari accumulation, snapshot fix) ---
    ngs_files = {
        "p=2, h=1.0mm": "2026-05-10-cln_team28_cuboid521_O2_SNAP_results.json",
        "p=3, h=0.5mm": "2026-05-10-cln_team28_cuboid521_O3_h05_SNAP_results.json",
    }
    ngs_results = {}
    for label, fname in ngs_files.items():
        path = base / fname
        if not path.exists():
            continue
        with open(path) as f:
            d = json.load(f)
        ngs_results[label] = d["stages"]
        print(f"--- NGSolve cln_team28 ({label}) ---")
        print(f"  {'k':>3} | {'R_n [Ohm]':>14} | {'L_n [H]':>14}"
              f" | {'tau_n [us]':>14}")
        for s in d["stages"][:8]:
            print(f"  {s['stage']:>3} | {s['R_n']:>14.6e}"
                  f" | {s['L_n']:>14.6e} | {s['tau_n_us']:>14.6f}")
        print()

    # --- Side-by-side R,L comparison (raw, no normalization) ---
    if bem_results and ngs_results:
        print("=" * 78)
        print(" Side-by-side R, L, tau (raw, no scaling) — leading Cauer rung k=0")
        print("=" * 78)
        for blabel, rungs in bem_results.items():
            print(f"  BEM {blabel}: R_0 = {rungs[0]['R_2k']:.4e},"
                  f" L_1 = {rungs[0]['L_2k_plus_1']:.4e},"
                  f" tau = {rungs[0]['tau_pair_us']:.4f} us")
        for nlabel, stages in ngs_results.items():
            print(f"  NGS {nlabel}: R_0 = {stages[0]['R_n']:.4e},"
                  f" L_0 = {stages[0]['L_n']:.4e},"
                  f" tau = {stages[0]['tau_n_us']:.4f} us")
        print()

        # Use NGSolve p=3,h=0.5 as reference. Scale BEM R, L so R_0 matches.
        ref_label = "p=3, h=0.5mm" if "p=3, h=0.5mm" in ngs_results else list(ngs_results)[0]
        ref_R0 = ngs_results[ref_label][0]["R_n"]
        print(f" Scaling BEM rungs so R_0 matches NGS {ref_label} R_0 = {ref_R0:.4e}")
        print(" (R, L absolute compare in canonical Joule-loss normalization)")
        print()
        for blabel, rungs in bem_results.items():
            scale = ref_R0 / rungs[0]["R_2k"]
            print(f"--- BEM {blabel} scaled (factor {scale:.4e}) ---")
            print(f"  {'k':>3} | {'R_2k [Ohm]':>14} | {'L_{2k+1} [H]':>14}"
                  f" | {'tau_pair [us]':>14}")
            for r in rungs[:8]:
                Rs = r['R_2k'] * scale
                Ls = r['L_2k_plus_1'] * scale
                print(f"  {r['k']:>3} | {Rs:>14.6e} | {Ls:>14.6e}"
                      f" | {r['tau_pair_us']:>14.6f}")
            print()

        # Direct compare (leading rungs k=0..6)
        print("=" * 78)
        print(" tau_pair high-stage breakdown investigation")
        print("=" * 78)
        rl = list(bem_results.items())[0]
        if "F-4 all 81 modes" in bem_results:
            rl = ("F-4 all 81 modes", bem_results["F-4 all 81 modes"])
        bem_label, bem_rungs = rl
        ngs_stages = ngs_results.get("p=3, h=0.5mm", list(ngs_results.values())[0])
        print(f"  (BEM = {bem_label} Kameari conv;"
              f" NGS = cln_team28 p=3 h=0.5mm Kameari)")
        print(f"  {'k':>3} | {'BEM tau [us]':>14} | {'NGS tau [us]':>14} | {'gap':>8} | breakdown?")
        for k in range(min(len(bem_rungs), len(ngs_stages))):
            bt = bem_rungs[k]['tau_pair_us']
            nt = ngs_stages[k]['tau_n_us']
            gap = (nt - bt) / bt * 100 if bt != 0 else float('nan')
            mono_bem = "" if k == 0 else (
                "OK" if bem_rungs[k]['tau_pair_us'] < bem_rungs[k-1]['tau_pair_us']
                else "↑NM"
            )
            mono_ngs = "" if k == 0 else (
                "OK" if ngs_stages[k]['tau_n_us'] < ngs_stages[k-1]['tau_n_us']
                else "↑NM"
            )
            print(f"  {k:>3} | {bt:>14.4f} | {nt:>14.4f} | {gap:>+7.1f}%"
                  f" | BEM:{mono_bem} NGS:{mono_ngs}")


if __name__ == "__main__":
    main()
