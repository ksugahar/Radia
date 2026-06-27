"""disk_bem_cauer.py -- Cauer-I CFE from BEM Foster moments, Nagamine pipeline.

Implements the algorithm from
   H. Nagamine, T. Yamaguchi, K. Sugahara, S. Hiruma, T. Mifune, T. Matsuo,
   "Verified Numerical Computations of the Cauer Network Representation of a
    Square Prism Conductor"
   (2026-05-04 manuscript), Fig. 5, eqs (10)-(67).

Pipeline (Nagamine Fig 5):
  (1) Summation:  Foster Form II  ->  Taylor series  c_n = Sum_m g_m^2 tau_m^(n+1)
                  + truncation error analysis (NOT IMPLEMENTED here -- requires
                   MPFR + interval arithmetic; this script uses mpmath 50-digit
                   floats only, so it inherits Nagamine's *structure* but NOT
                   his "verified" interval bounds).
  (2) QD algorithm:  Taylor coefficients -> continued-fraction (CFE)
  (3) Equivalence transform:  CFE -> Cauer ladder (R_0, L_1, R_2, L_3, ...)

Steps (2) + (3) are merged in `cauer_extract` below via the classical Cauer
extraction by repeated Taylor inversion -- mathematically equivalent to
Nagamine eqs (64)-(67), but computed in one pass instead of two.

The output is the Cauer ladder of Fig. 3:
                                                   .
   in --R_0-+--R_2-+--R_4-+-...                    .
            |       |       |                      .
           L_1     L_3     L_5  ...                .
            |       |       |                      .
           gnd     gnd     gnd                     .

with R_{2k} (k=0,1,2,...) the series resistors (even subscripts) and
L_{2k+1} (k=0,1,2,...) the shunt inductors (odd subscripts). Per-stage
time constant: tau_pair[k] = L_{2k+1}/R_{2k}.

Comparison endpoint: tau_pair[k] is normalisation-invariant and directly
comparable with the FE Hiruma 3-term output (radia-core axifem). Absolute R, L
values depend on the Foster-amplitude normalisation chosen on the BEM side.
"""
from __future__ import annotations

import json
from pathlib import Path

import mpmath as mp

mp.mp.dps = 50

REPO = Path(__file__).resolve().parents[0]
BEM_JSON = REPO / "bem_disk_axisym_cauer.json"
Q2_JSON = Path("S:/Radia/01_GitHub/validation_test/axifem/research/verification/test_hiruma_disk_q2_results.json")
Q1_JSON = Path("S:/Radia/01_GitHub/validation_test/axifem/research/verification/test_hiruma_disk_q1_results.json")


def invert_taylor(c):
    """Invert Taylor series 1 / (c_0 + c_1 s + c_2 s² + ...)."""
    n = len(c)
    e = [mp.mpf(0)] * n
    e[0] = 1 / c[0]
    for k in range(1, n):
        sum_term = mp.fsum(c[j] * e[k - j] for j in range(1, k + 1))
        e[k] = -sum_term / c[0]
    return e


def cauer_extract(c_taylor, max_stages):
    """Cauer-I continued-fraction expansion via repeated Taylor inversion.
    Each stage extracts the leading coefficient of the inverted series and
    drops it; remainder is fed back as next Taylor series.
    Returns p_1, p_2, ... (alternating R^-1 / L rung-like values)."""
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
    print("=" * 76)
    print("Phase 3-(3) BEM Cauer-I (Nagamine-style)  vs  radia-core axifem Hiruma 3-term")
    print("=" * 76)
    print()

    bem = json.loads(BEM_JSON.read_text())
    print(f"BEM: {bem['method']}")
    print(f"  ne = {bem['params']['ne']}, n_modes = {bem['params']['n_modes']}, "
          f"n_moments = {bem['params']['n_moments']}")

    # alphas are stored as floats (Mathematica N[..., 16]); load to mpf.
    alphas = [mp.mpf(repr(a)) for a in bem["alphas"]]
    bem_taus = [mp.mpf(repr(t)) for t in bem["bem_tau_us"]]

    print(f"  BEM Foster tau_1..tau_6 (us) = "
          f"{[float(t) for t in bem_taus[:6]]}")
    print(f"  alpha_0 = {mp.nstr(alphas[0], 6)}")
    print(f"  alpha_1 = {mp.nstr(alphas[1], 6)}")
    print(f"  alpha_2 = {mp.nstr(alphas[2], 6)}")
    print()

    # === Sanity: leading mode estimate from -alpha_1/alpha_0 ===
    tau1_moment = -alphas[1] / alphas[0]
    print(f"Leading-mode estimate tau_1 = -alpha_1/alpha_0 = "
          f"{mp.nstr(tau1_moment * 1e6, 8)} us  (BEM Foster: "
          f"{mp.nstr(bem_taus[0], 8)} us)")
    print()

    # === Cauer extraction (Nagamine eqs 64-67 mathematical equivalent) ===
    # c_taylor = [alpha_0, alpha_1, alpha_2, ...] (Foster Taylor moments).
    # The classical Cauer extraction by repeated Taylor inversion produces a
    # sequence p_1, p_2, p_3, ... corresponding to:
    #     p_{2k+1} = R_{2k}        (Cauer series resistor at even subscript)
    #     p_{2k+2} = 1 / L_{2k+1}  (inverse Cauer shunt inductor at odd subscript)
    # These are mathematically the same R, L that Nagamine's QD + equivalence-
    # transform pipeline produces (Fig 5, steps 2+3); the extraction route
    # differs.
    print(f"Cauer extraction on {len(alphas)}-term Taylor series ...")
    p_cauer = cauer_extract(alphas, len(alphas) - 1)
    print(f"  Extracted {len(p_cauer)} p-coefficients")
    print()

    # === Convert p_cauer to Nagamine R_{2k}, L_{2k+1} naming ===
    n_pairs = len(p_cauer) // 2
    bem_R_2k = []          # R_0, R_2, R_4, ...
    bem_L_2k_plus_1 = []   # L_1, L_3, L_5, ...
    bem_tau_pair = []      # tau_pair[k] = L_{2k+1} / R_{2k}
    print(f"{'k':>3} {'R_{2k}':>16} {'L_{2k+1}':>16} {'tau_pair (us)':>16}")
    for k in range(n_pairs):
        R   = p_cauer[2*k]
        Lin = p_cauer[2*k + 1]                    # = 1 / L_{2k+1}
        L   = 1 / Lin if abs(Lin) > mp.mpf("1e-300") else mp.mpf("nan")
        tau = L / R if abs(R) > mp.mpf("1e-300") else mp.mpf("nan")
        bem_R_2k.append(R)
        bem_L_2k_plus_1.append(L)
        bem_tau_pair.append(tau)
        print(f"{k:>3} {mp.nstr(R, 8):>16} {mp.nstr(L, 8):>16} "
              f"{mp.nstr(tau * 1e6, 8):>16}")
    print()
    print("(R, L absolute values use the Foster-amplitude normalisation of the")
    print(" BEM script; tau_pair = L_{2k+1}/R_{2k} is normalisation-invariant)")
    print()

    # === Compare with FE Hiruma 3-term (axihenrotte order=1, order=2) ===
    q2 = json.loads(Q2_JSON.read_text())["results"]["fine"]
    q1 = json.loads(Q1_JSON.read_text())["results"]["very fine"]

    q2_tau_pair = [mp.mpf(repr(s["tau_per_stage_us"])) for s in q2["stages"]]
    q1_tau_pair = [mp.mpf(repr(s["tau_per_stage_us"])) for s in q1["stages"]]

    print("=" * 76)
    print("Cross-validation: tau_pair[k] = L_{2k+1} / R_{2k}  (us)")
    print("=" * 76)
    n_compare = min(6, len(q2_tau_pair), len(bem_tau_pair))
    print(f"{'k':>3} {'BEM Cauer':>15} {'p=2 fine':>12} {'p=1 vfine':>14} "
          f"{'p=2/BEM%':>10} {'p=1/BEM%':>10}")
    for k in range(n_compare):
        bem_v = float(bem_tau_pair[k] * 1e6)
        q2_v = float(q2_tau_pair[k])
        q1_v = float(q1_tau_pair[k])
        gap_q2 = (q2_v - bem_v) / bem_v * 100
        gap_q1 = (q1_v - bem_v) / bem_v * 100
        print(f"{k:>3} {bem_v:>15.6f} {q2_v:>12.4f} {q1_v:>14.4f} "
              f"{gap_q2:>+10.4f} {gap_q1:>+10.4f}")
    print()

    # === Save BEM Cauer output (Nagamine R_{2k}, L_{2k+1} naming) ===
    out = {
        "method": ("BEM Foster -> Cauer ladder (Nagamine pipeline structure: "
                   "summation -> Taylor -> classical Cauer extraction; "
                   "no verified interval arithmetic)"),
        "reference": ("Nagamine, Yamaguchi, Sugahara, Hiruma, Mifune, Matsuo, "
                      "'Verified Numerical Computations of the Cauer Network "
                      "Representation of a Square Prism Conductor', 2026-05-04"),
        "n_moments": len(alphas),
        "n_pairs": n_pairs,
        "naming_convention": ("R_{2k} = series resistor (even subscript), "
                              "L_{2k+1} = shunt inductor (odd subscript), "
                              "tau_pair[k] = L_{2k+1}/R_{2k}"),
        "p_cauer":         [mp.nstr(p, 30) for p in p_cauer],
        "R_2k":            [mp.nstr(r, 30) for r in bem_R_2k],
        "L_2k_plus_1":     [mp.nstr(l, 30) for l in bem_L_2k_plus_1],
        "tau_pair_us":     [mp.nstr(t * 1e6, 30) for t in bem_tau_pair],
        # Back-compat alias for older test scripts (= tau_pair_us):
        "bem_tau_rung_us": [mp.nstr(t * 1e6, 30) for t in bem_tau_pair],
    }
    out_path = REPO / "bem_disk_axisym_cauer_python_results.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(f"Saved BEM Cauer extraction: {out_path}")


if __name__ == "__main__":
    main()
