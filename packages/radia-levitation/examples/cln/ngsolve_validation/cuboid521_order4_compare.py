"""cuboid521_order4_compare.py — Compare BEM order=3 vs order=4 + NGSolve high-mesh + ELF
for cuboid 5x2x1 mm Cu, B0 = 1T.

Goal: Does BEM order=4 (n_dofs=176) give more reliable Cauer rungs than order=3 (n_dofs=81)?
Are higher Cauer stages more consistent across methods?
"""
from __future__ import annotations

import json, sys
sys.stdout.reconfigure(encoding="utf-8")
from pathlib import Path
import mpmath as mp

mp.mp.dps = 100


def hankel_pade_cauer(alphas, max_stages):
    c = [mp.mpf(a) for a in alphas]
    p = []
    for _ in range(max_stages):
        if len(c) < 2 or abs(c[0]) < mp.mpf(10) ** (-mp.mp.dps + 10):
            break
        n = len(c)
        e = [mp.mpf(0)] * n
        e[0] = 1 / c[0]
        for k in range(1, n):
            e[k] = -sum(c[j+1] * e[k-j-1] for j in range(k)
                        if k-j-1 >= 0) / c[0]
        p.append(e[0])
        c = e[1:]
    return p


def bem_to_cauer(modes, n_moments=80, max_pairs=20):
    tau = [mp.mpf(e["tau_us"]) * mp.mpf("1e-6") for e in modes]
    g2 = [mp.mpf(e["g2"]) for e in modes]
    out = {}
    for label, off in [("Hiruma", 1), ("Kameari", 0)]:
        a = [(mp.mpf(-1))**n * sum(g * t**(n+off) for g,t in zip(g2, tau))
             for n in range(n_moments)]
        p = hankel_pade_cauer(a, 2*max_pairs)
        rungs = []
        for k in range(len(p) // 2):
            R = p[2*k]; Linv = p[2*k+1]
            if abs(Linv) < mp.mpf(10)**(-50): break
            L = 1/Linv
            rungs.append({"k": k, "R_2k": float(R),
                          "L_2k_plus_1": float(L),
                          "tau_us": float(L/R) * 1e6})
        out[label] = rungs
    return out


def main():
    base = Path(__file__).parent
    files = {
        "order=3 (n_dofs=81)":  "cuboid521_radia_vim_F4_full81.json",
        "order=4 (n_dofs=176)": "cuboid521_radia_vim_F4_order4.json",
    }
    results = {}
    for label, f in files.items():
        path = base / f
        if not path.exists():
            print(f"  SKIP (missing): {label}")
            continue
        with open(path) as fp:
            d = json.load(fp)
        modes = d["top_modes_by_g2"]
        results[label] = (d, modes, bem_to_cauer(modes))
        print(f"  {label}: n_modes saved = {len(modes)},"
              f" leading Foster = {d['leading_tau_us']:.4f} us")

    print()
    print("=" * 90)
    print(" BEM Phase F-4 order=3 vs order=4 — Kameari Cauer rungs")
    print("=" * 90)
    print(f"  {'k':>3} | "
          f"{'order=3 tau':>13} | {'order=3 R_2k':>14} | {'order=3 L_2k+1':>14}"
          f" || {'order=4 tau':>13} | {'order=4 R_2k':>14} | {'order=4 L_2k+1':>14}")
    o3 = results.get("order=3 (n_dofs=81)", (None, None, {}))[2].get("Kameari", [])
    o4 = results.get("order=4 (n_dofs=176)", (None, None, {}))[2].get("Kameari", [])
    for k in range(min(10, max(len(o3), len(o4)))):
        r3 = o3[k] if k < len(o3) else {"tau_us": float('nan'),
                                         "R_2k": float('nan'),
                                         "L_2k_plus_1": float('nan')}
        r4 = o4[k] if k < len(o4) else {"tau_us": float('nan'),
                                         "R_2k": float('nan'),
                                         "L_2k_plus_1": float('nan')}
        print(f"  {k:>3} | {r3['tau_us']:>13.4f} | {r3['R_2k']:>14.4e}"
              f" | {r3['L_2k_plus_1']:>14.4e} || {r4['tau_us']:>13.4f}"
              f" | {r4['R_2k']:>14.4e} | {r4['L_2k_plus_1']:>14.4e}")

    # Compare with NGSolve fine + ELF
    ngsf = base / "2026-05-10-cln_team28_cuboid521_O3_h025_air6_results.json"
    if ngsf.exists():
        with open(ngsf) as f:
            ngs = json.load(f)["stages"]
    else:
        ngs = []

    elf_taus = [11.506766, 6.430171, 3.278698, 1.497229]
    elf_amps = [1.28e-7, 1.55e-8, 2.10e-8, 9.24e-10]
    elf_modes = [{"tau_us": t, "g2": a} for t, a in zip(elf_taus, elf_amps)]
    elf_rungs = bem_to_cauer(elf_modes, n_moments=12, max_pairs=4)["Kameari"]

    print()
    print("=" * 90)
    print(" 4-way comparison: BEM order=3, BEM order=4, NGSolve fine, ELF")
    print("=" * 90)
    print(f"  {'k':>3} | {'BEM o3 tau':>13} | {'BEM o4 tau':>13}"
          f" | {'NGS h0.25':>13} | {'ELF':>13}")
    for k in range(8):
        v3 = o3[k]['tau_us'] if k < len(o3) else float('nan')
        v4 = o4[k]['tau_us'] if k < len(o4) else float('nan')
        vn = ngs[k]['tau_n_us'] if k < len(ngs) else float('nan')
        ve = elf_rungs[k]['tau_us'] if k < len(elf_rungs) else float('nan')
        print(f"  {k:>3} | {v3:>13.4f} | {v4:>13.4f} | {vn:>13.4f} | {ve:>13.4f}")


if __name__ == "__main__":
    main()
