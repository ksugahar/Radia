"""hex_vim_kameari_from_foster.py — apply Kameari Hankel-Padé to a Foster
spectrum saved by radia-vim extract_tau_*.py.

Reads {idx, tau_us, g2} list from JSON, computes alpha_n^Kameari moments,
and extracts Cauer-I rungs at 80-digit precision.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

import mpmath as mp
mp.mp.dps = 80


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
            e[k] = -sum(c[j + 1] * e[k - j - 1] for j in range(k)
                        if k - j - 1 >= 0) / c[0]
        p.append(e[0])
        c = e[1:]
    return p


def kameari_from_foster_json(json_path, n_use=None, n_moments=40, n_stages=12):
    data = json.loads(Path(json_path).read_text(encoding="utf-8"))
    modes = data["top_modes_by_g2"]
    if n_use is None:
        n_use = len([m for m in modes if abs(m["g2"]) > 1e-32])
    n_use = min(n_use, len(modes))
    print(f"  Using {n_use} Foster modes (g2 > 1e-32)")

    tau_seconds = [mp.mpf(float(m["tau_us"])) * mp.mpf("1e-6")
                   for m in modes[:n_use]]
    g2_mp = [mp.mpf(float(m["g2"])) for m in modes[:n_use]]

    n_moments = min(n_moments, 2 * n_use - 4)
    alphas = []
    for n in range(n_moments):
        a = mp.mpf(0)
        for gv, tv in zip(g2_mp, tau_seconds):
            a += gv * tv ** n
        alphas.append((mp.mpf(-1)) ** n * a)

    p = hankel_pade_cauer(alphas, n_stages)
    print(f"  Kameari Cauer-I rungs (80 dig):")
    print(f"  k | R_2k       | L_2k+1     | tau_pair us")
    rungs = []
    for k in range(min(n_stages // 2, len(p) // 2)):
        Rv = p[2 * k]
        Linv = p[2 * k + 1]
        if abs(Linv) < mp.mpf(10) ** (-50):
            break
        Lv = 1 / Linv
        tau_p = float(Lv / Rv) * 1e6
        print(f"  {k} | {float(Rv):>10.4e} | {float(Lv):>10.4e} | {tau_p:.4f}")
        rungs.append({"k": k, "R_2k": float(Rv),
                      "L_2k_plus_1": float(Lv), "tau_pair_us": tau_p})

    return rungs


def main():
    if len(sys.argv) < 2:
        print("Usage: python hex_vim_kameari_from_foster.py <foster.json>")
        sys.exit(1)
    json_path = sys.argv[1]
    print(f"Reading: {json_path}")
    rungs = kameari_from_foster_json(json_path)
    out_path = Path(json_path).with_suffix("").parent / (
        Path(json_path).stem + "_cauer_kameari.json")
    out_path.write_text(json.dumps({"rungs": rungs}, indent=2),
                        encoding="utf-8")
    print(f"\nSaved: {out_path}")


if __name__ == "__main__":
    main()
