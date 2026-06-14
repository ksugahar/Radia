"""hex_vim_cauer_interval.py — Verified-interval Cauer-I extraction (Nagamine手法).

Implements interval-arithmetic Hankel-Padé continued fraction from a Foster
spectrum, propagating rigorous error bounds through every operation. Tracks
interval width per stage to identify the truncation point where the
extraction becomes numerically unreliable.

Algorithm:
  Input: Foster {tau_k, g_k} from VIM eigenvalue problem (FP64 from GPU)
  Step 1: convert to interval [tau_k*(1-eps), tau_k*(1+eps)] with eps = 1e-14
          (FP64 rounding noise floor)
  Step 2: compute alpha_n moments in interval arithmetic
  Step 3: Hankel-Padé recurrence in interval arithmetic
  Step 4: report Cauer rungs with rigorous lower/upper bounds + width

Reference: Nagamine, Yamaguchi, Sugahara, Hiruma, Mifune, Matsuo,
"Verified Numerical Computations of the Cauer Network Representation of
a Square Prism Conductor", JJIAM 2026 (submitted).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

import mpmath as mp

# Interval arithmetic precision (mantissa digits)
mp.iv.dps = 80


def hankel_pade_cauer_interval(alphas_iv, max_stages, width_warn=1e-3):
    """Hankel-Padé continued-fraction extraction with interval arithmetic.

    Args:
      alphas_iv: list of mpmath interval (mp.iv.mpf) moments
      max_stages: maximum Cauer rungs to attempt
      width_warn: relative interval width (b-a)/|midpoint| above which to
                  flag unreliable

    Returns: list of dicts with rung index, midpoint, low, high, width.
    """
    c = [a for a in alphas_iv]
    out = []
    for stage in range(max_stages):
        if len(c) < 2:
            break
        if 0 in c[0]:   # interval contains zero -> singular
            print(f"  Stage {stage}: c[0] contains zero, truncating")
            break
        n = len(c)
        e = [mp.iv.mpf(0)] * n
        e[0] = mp.iv.mpf(1) / c[0]
        for k in range(1, n):
            acc = mp.iv.mpf(0)
            for j in range(k):
                if k - j - 1 >= 0:
                    acc = acc + c[j + 1] * e[k - j - 1]
            e[k] = -acc / c[0]
        # Save e[0] as the next Cauer coefficient (use scalar mpf bounds)
        rung = e[0]
        rung_low = mp.mpf(rung.a)
        rung_high = mp.mpf(rung.b)
        midpoint = (rung_low + rung_high) / 2
        width = rung_high - rung_low
        rel_width = (width / abs(midpoint)
                     if abs(midpoint) > 0 else mp.mpf("inf"))
        out.append({
            "stage": stage,
            "midpoint": float(midpoint),
            "low": float(rung.a),
            "high": float(rung.b),
            "width": float(width),
            "rel_width": float(rel_width),
        })
        if rel_width > width_warn:
            print(f"  Stage {stage}: rel_width {float(rel_width):.2e} > "
                  f"{width_warn:.0e} -- truncating (verified limit)")
            break
        c = e[1:]
    return out


def cauer_rungs_interval(rungs_raw):
    """Group continued-fraction coefficients into (R_2k, L_{2k+1}) pairs."""
    pairs = []
    n_pairs = len(rungs_raw) // 2
    for k in range(n_pairs):
        Rv = rungs_raw[2 * k]
        Linv = rungs_raw[2 * k + 1]
        # L = 1 / Linv (interval inverse)
        Linv_mid = (Linv["low"] + Linv["high"]) / 2
        # interval inverse: if Linv = [a, b] both positive,
        # 1/Linv = [1/b, 1/a]
        if Linv["low"] > 0:
            L_low = 1.0 / Linv["high"]
            L_high = 1.0 / Linv["low"]
        elif Linv["high"] < 0:
            L_low = 1.0 / Linv["high"]
            L_high = 1.0 / Linv["low"]
        else:
            L_low = float("-inf"); L_high = float("inf")
        L_mid = 1.0 / Linv_mid if abs(Linv_mid) > 0 else float("inf")
        tau_mid = L_mid / Rv["midpoint"] * 1e6
        pairs.append({
            "k": k,
            "R_2k_mid": Rv["midpoint"],
            "R_2k_low": Rv["low"], "R_2k_high": Rv["high"],
            "R_2k_rel_width": Rv["rel_width"],
            "L_2k1_mid": L_mid,
            "L_2k1_low": L_low, "L_2k1_high": L_high,
            "tau_pair_mid_us": tau_mid,
            "L_inv_rel_width": Linv["rel_width"],
        })
    return pairs


def extract_from_foster_json(json_path, eps_rel=1e-14, n_moments=40,
                              n_stages=12, width_warn=1e-2):
    """Run verified-interval Cauer extraction from radia-vim Foster JSON.

    The JSON should contain "top_modes_by_g2" with tau_us and g2.
    """
    data = json.loads(Path(json_path).read_text(encoding="utf-8"))
    modes = data["top_modes_by_g2"]
    n_use = len([m for m in modes if abs(m["g2"]) > 1e-32])
    print(f"  Using {n_use} Foster modes (g2 > 1e-32)")

    # Build interval Foster spectrum
    tau_iv = []
    g2_iv = []
    eps = mp.iv.mpf((1 - eps_rel, 1 + eps_rel))   # interval [1-eps, 1+eps]
    for m in modes[:n_use]:
        t = mp.iv.mpf(m["tau_us"]) * mp.iv.mpf("1e-6")
        t_iv = t * eps   # widen by eps_rel
        tau_iv.append(t_iv)
        g = mp.iv.mpf(m["g2"]) * eps
        g2_iv.append(g)

    # Compute Kameari moments alpha_n = (-1)^n sum g_k^2 tau_k^n in interval
    n_moments = min(n_moments, 2 * n_use - 4)
    print(f"  Computing {n_moments} moments in interval arithmetic...")
    alphas_iv = []
    for n in range(n_moments):
        acc = mp.iv.mpf(0)
        for gv, tv in zip(g2_iv, tau_iv):
            term = gv * (tv ** n)
            acc = acc + term
        sign = mp.iv.mpf(1) if n % 2 == 0 else mp.iv.mpf(-1)
        alphas_iv.append(sign * acc)
        # Show interval width at this moment (use mp.mpf scalars)
        a_low = mp.mpf(alphas_iv[-1].a)
        a_high = mp.mpf(alphas_iv[-1].b)
        a_mid = (a_low + a_high) / 2
        a_width = a_high - a_low
        a_rel = a_width / abs(a_mid) if abs(a_mid) > 0 else mp.mpf("inf")
        if n < 8:
            print(f"    alpha_{n} midpoint = {float(a_mid):+.4e}, "
                  f"rel_width = {float(a_rel):.2e}")

    # Hankel-Padé in interval
    print(f"\n  Hankel-Padé continued fraction (interval arithmetic):")
    rungs_raw = hankel_pade_cauer_interval(alphas_iv, n_stages,
                                            width_warn=width_warn)
    pairs = cauer_rungs_interval(rungs_raw)

    print(f"\n  Cauer-I verified rungs:")
    print(f"  k | R_2k mid           | R_2k rel_width | L_2k+1 mid    | "
          f"L_2k+1 rel_width | tau_pair mid us | reliable?")
    for p in pairs:
        reliable = "yes" if (p["R_2k_rel_width"] < width_warn
                              and p["L_inv_rel_width"] < width_warn) else "NO"
        print(f"  {p['k']} | {p['R_2k_mid']:>15.4e}   | "
              f"{p['R_2k_rel_width']:>10.2e}     | "
              f"{p['L_2k1_mid']:>10.4e}    | "
              f"{p['L_inv_rel_width']:>10.2e}        | "
              f"{p['tau_pair_mid_us']:>10.4f}      | {reliable}")

    return {"foster_json": str(json_path),
            "n_use": n_use, "eps_rel": eps_rel,
            "verified_rungs": pairs,
            "raw_continued_fraction": rungs_raw}


def main():
    if len(sys.argv) < 2:
        print("Usage: python hex_vim_cauer_interval.py <foster.json> "
              "[eps_rel=1e-14] [width_warn=1e-2]")
        # Default test: run on cuboid 5x2x1 order=5 results
        json_path = (Path("S:/Radia/01_GitHub/packages/radia-vim/scripts")
                     / "cuboid521_order4_results.json")
        if not json_path.exists():
            print(f"Default JSON not found: {json_path}")
            sys.exit(1)
    else:
        json_path = Path(sys.argv[1])

    eps_rel = float(sys.argv[2]) if len(sys.argv) > 2 else 1e-14
    width_warn = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-2

    print(f"\n=== Verified-interval Cauer extraction ===")
    print(f"  Foster JSON: {json_path}")
    print(f"  eps_rel (input noise floor): {eps_rel:.0e}")
    print(f"  width_warn (truncation threshold): {width_warn:.0e}")

    result = extract_from_foster_json(json_path, eps_rel=eps_rel,
                                       width_warn=width_warn)

    out_path = json_path.parent / (json_path.stem + "_verified_cauer.json")
    out_path.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(f"\n  Saved: {out_path}")


if __name__ == "__main__":
    main()
