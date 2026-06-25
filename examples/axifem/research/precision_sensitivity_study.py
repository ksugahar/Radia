"""Precision sensitivity study for BEM-Foster -> Cauer extraction pipeline.

Goal: quantitatively answer "is double precision enough for K matrix
assembly in radia-vim, or do we need quad / mpfr for the cuboid Cauer
ladder paper?"

Method:
  1. Synthetic Foster spectrum with known {(g_k, tau_k)}, k = 1..N.
  2. Compute alpha_n moments at "ground-truth" 200-digit precision.
  3. Run Cauer extraction at high precision -> reference Cauer rungs.
  4. TRUNCATE the input g_k, tau_k to D digits (simulating K assembly at D
     digits) and recompute Cauer rungs.
  5. Count "clean" stages (relative error < 1e-3 vs reference) for each D.

Sweep D in {16 (double), 24, 32 (quad-like), 48, 60, 100} and report the
break-even between assembly precision and number of usable Cauer pairs.
"""
from __future__ import annotations

import json
from pathlib import Path
import mpmath as mp


REF_DPS = 200
mp.mp.dps = REF_DPS


def make_synthetic_spectrum(n_modes, ratio=0.5, g_decay=0.5, seed=42):
    """Foster spectrum mimicking a 3D eddy-current problem.

    tau_k = tau_1 * ratio^(k-1) (geometric decay)
    g_k   = g_1 * g_decay^(k-1) * (-1)^k (alternating signs)
    Default: tau_1 = 1e-4 sec, ratio 0.5, g_decay 0.5.
    """
    mp.mp.dps = REF_DPS
    tau_1 = mp.mpf("1e-4")
    g_1   = mp.mpf("1e-3")
    taus = [tau_1 * mp.mpf(ratio)**k for k in range(n_modes)]
    gs   = [g_1   * mp.mpf(g_decay)**k * ((-1)**k) for k in range(n_modes)]
    return taus, gs


def truncate_to_digits(value, n_digits):
    """Round mpf value to n_digits significant decimal digits."""
    if value == 0:
        return mp.mpf(0)
    s = mp.nstr(value, n_digits)
    return mp.mpf(s)


def compute_alphas(taus, gs, n_moments):
    alphas = []
    for n in range(n_moments):
        s = mp.mpf(0)
        for k in range(len(taus)):
            s += gs[k]**2 * taus[k]**(n + 1)
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


def cauer_to_pairs(p_cauer):
    """Return list of tau_pair[k] = L_{2k+1} / R_{2k} from p_cauer."""
    pairs = []
    n = len(p_cauer) // 2
    for k in range(n):
        R = p_cauer[2*k]
        Lin = p_cauer[2*k + 1]
        if abs(Lin) < mp.mpf("1e-300") or abs(R) < mp.mpf("1e-300"):
            pairs.append(mp.mpf("nan"))
            continue
        L = 1 / Lin
        pairs.append(L / R)
    return pairs


def n_clean_stages(ref_pairs, test_pairs, rtol=mp.mpf("1e-3")):
    n = min(len(ref_pairs), len(test_pairs))
    for k in range(n):
        if not mp.isfinite(test_pairs[k]) or test_pairs[k] <= 0:
            return k
        if abs(test_pairs[k] - ref_pairs[k]) / abs(ref_pairs[k]) > rtol:
            return k
    return n


def study(n_modes, n_moments_factor=4, label=""):
    print("=" * 78)
    print(f"Precision study: {label}  (n_modes={n_modes}, "
          f"n_moments={n_modes * n_moments_factor})")
    print("=" * 78)

    mp.mp.dps = REF_DPS
    taus_ref, gs_ref = make_synthetic_spectrum(n_modes)
    n_moments = n_modes * n_moments_factor

    alphas_ref = compute_alphas(taus_ref, gs_ref, n_moments)
    p_ref = cauer_extract(alphas_ref, n_moments - 1)
    pairs_ref = cauer_to_pairs(p_ref)
    n_pairs_ref = sum(1 for p in pairs_ref if mp.isfinite(p) and p > 0)
    print(f"  Reference (dps={REF_DPS}): "
          f"{len(p_ref)} p-coeffs, {n_pairs_ref} positive pairs")
    print(f"  Reference tau_pair[0..min(5,n)] = "
          f"{[float(p)*1e6 if mp.isfinite(p) else None for p in pairs_ref[:5]]} us")
    print()

    print(f"{'D digits':>10}  {'clean stages':>14}  {'rel.err k=0':>14}  "
          f"{'rel.err k=N/2':>15}")
    results = []
    for D in [10, 14, 16, 20, 24, 32, 48, 60, 100, REF_DPS]:
        # truncate INPUT data to D digits, then recompute everything at REF_DPS
        # (simulating: K assembly precision = D digits, downstream pipeline
        # at full precision)
        taus_d = [truncate_to_digits(t, D) for t in taus_ref]
        gs_d   = [truncate_to_digits(g, D) for g in gs_ref]
        alphas_d = compute_alphas(taus_d, gs_d, n_moments)
        p_d = cauer_extract(alphas_d, n_moments - 1)
        pairs_d = cauer_to_pairs(p_d)
        n_clean = n_clean_stages(pairs_ref, pairs_d)

        rel_k0 = (abs(pairs_d[0] - pairs_ref[0]) / abs(pairs_ref[0])
                  if mp.isfinite(pairs_d[0]) else mp.mpf("inf"))
        mid = min(n_pairs_ref // 2, len(pairs_d) - 1)
        rel_kmid = (abs(pairs_d[mid] - pairs_ref[mid]) / abs(pairs_ref[mid])
                    if mp.isfinite(pairs_d[mid]) and abs(pairs_ref[mid]) > 0
                    else mp.mpf("inf"))
        print(f"{D:>10}  {n_clean:>14}  {mp.nstr(rel_k0, 4):>14}  "
              f"{mp.nstr(rel_kmid, 4):>15}")
        results.append({"D": D, "n_clean": n_clean,
                        "rel_k0": float(rel_k0) if mp.isfinite(rel_k0) else None,
                        "rel_kmid": float(rel_kmid) if mp.isfinite(rel_kmid) else None})
    print()
    return results


def main():
    print()
    print("Precision sensitivity: how many clean Cauer pairs as function of")
    print("input data precision (K matrix entries / eigenvalues / Foster ampl)")
    print()

    out = {}
    out["small_6modes"]   = study(6,  n_moments_factor=10, label="cuboid order 3 analogue (6 modes)")
    out["medium_15modes"] = study(15, n_moments_factor=8,  label="cuboid order 4 analogue (15 modes)")
    out["large_30modes"]  = study(30, n_moments_factor=6,  label="cuboid order 5 analogue (30 modes)")
    out["xlarge_50modes"] = study(50, n_moments_factor=5,  label="cuboid order 6 analogue (50 modes)")

    out_path = Path(__file__).parent / "precision_sensitivity_results.json"
    out_path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"Saved: {out_path}")
    print()
    print("Interpretation:")
    print("  - 'clean stages' = Cauer pairs with rel.err < 1e-3 vs the dps=200")
    print("    reference (input data truncated to D digits, then full-precision")
    print("    moments + Cauer extraction).")
    print("  - D=16 simulates  IEEE double for K assembly")
    print("  - D=32 simulates  __float128 / cpp_bin_float_quad")
    print("  - D=60 simulates  cpp_dec_float<60> or MPFR 60-digit")


if __name__ == "__main__":
    main()
