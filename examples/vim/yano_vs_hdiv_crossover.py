"""Net-cost crossover model: improved yano-MSC (star + cell-graph-Laplacian AMG) vs HDiv-VIM
(charge-Gram + M_mass^-1 GMRES), from MEASURED component scalings.  Per-solve cost ~ iters x matvec,
matvec ~ H-matrix storage (every stored block touched once).  This is a MODEL combining measured
scalings -- NOT a direct head-to-head (that needs the C++ port); it estimates WHERE yano's bounded
iterations overtake HDiv's growing M_mass^-1 iterations despite yano's heavier (worse-compressing)
H-matrix.

Measured inputs:
  HDiv  (examples/vim/hdiv_matvec_scaling.json, unit-cube tet): charge-Gram rank 13->10 (DROPS),
        storage ~N^1.07, M_mass^-1 GMRES iters FLAT ~40 to 10k then GROW 47/80/398 @ 12k/21k/62k.
  yano  (1b C:\temp\yano_1b_hacapk_scaling, unit-cube hex): collocation rank 47->93 (GROWS), storage
        ~N^1.65 (PRE-ASYMPTOTIC -- only to 10k; like HDiv it should flatten at larger N), star+AMG iters
        BOUNDED ~16-26 (1a, real pyamg).

HONEST caveats: (1) cube is the WORST case for H-matrix compression for BOTH; (2) tet (HDiv) vs hex
(yano) DOF differ for the same geometry -- N is the problem-size axis, not matched accuracy; (3) yano's
N^1.65 storage is fit on 4 PRE-ASYMPTOTIC points and is likely pessimistic at large N; (4) this MODEL is
not a benchmark.  It also IGNORES yano's per-element accuracy edge (C-yoke 52 vs 416 hex) which would
shift the crossover further in yano's favor.
"""
import os
import json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

# --- measured HDiv (charge-Gram + M_mass^-1 GMRES), from hdiv_matvec_scaling.json ---
H_ndof = np.array([364, 1331, 2302, 7061, 11964, 21102, 62643], float)
H_iters = np.array([39, 36, 41, 38, 47, 80, 398], float)
H_mem = np.array([0.60, 5.62, 13.36, 44.48, 94.66, 201.07, 456.66], float)   # H-matrix MB ~ matvec cost

# --- measured improved-yano (collocation HACApK + star + AMG), from 1b + 1a ---
Y_ndof = np.array([750, 3072, 6000, 10368], float)
Y_mem = np.array([4.1, 49.1, 139.8, 309.7], float)
Y_iters_const = 20.0     # bounded (1a real-AMG: 16..26 over 18..600 elem)


def _fit_pow(x, y, lo=None):
    """log-log power-law fit y ~ a x^b (optionally only x>=lo, to capture the large-N regime)."""
    m = np.ones_like(x, bool) if lo is None else (x >= lo)
    b, la = np.polyfit(np.log(x[m]), np.log(y[m]), 1)
    return float(np.exp(la)), float(b)


# fit the LARGE-N regimes (HDiv iters grow only past ~10k; storage power laws over the full range)
ai, bi = _fit_pow(H_ndof, H_iters, lo=7061)      # HDiv iter growth (large N)
am, bm = _fit_pow(H_ndof, H_mem)                 # HDiv storage ~N^1.07
aY, bY = _fit_pow(Y_ndof, Y_mem)                 # yano storage ~N^1.65 (pre-asymptotic)

print(f"HDiv  iters ~ {ai:.3g}*N^{bi:.2f} (large N)   storage ~ {am:.3g}*N^{bm:.2f} MB")
print(f"yano  iters ~ {Y_iters_const:.0f} (bounded)     storage ~ {aY:.3g}*N^{bY:.2f} MB (pre-asymptotic)")

# per-solve cost ~ iters x matvec(=storage)
def cost_H(N):
    return max(ai * N ** bi, 36.0) * (am * N ** bm)      # iters floored at the flat ~36-40 regime
def cost_Y(N):
    return Y_iters_const * (aY * N ** bY)

print(f"\n{'N_dof':>8} | {'HDiv iters':>10} {'HDiv cost':>11} | {'yano iters':>10} {'yano cost':>11} | {'yano/HDiv':>9}")
rows = []
cross = None
for N in (1e3, 3e3, 1e4, 3e4, 1e5, 3e5, 1e6):
    iH = max(ai * N ** bi, 36.0)
    cH, cY = cost_H(N), cost_Y(N)
    ratio = cY / cH
    if cross is None and ratio < 1.0:
        cross = N
    print(f"{N:>8.0f} | {iH:>10.0f} {cH:>11.3g} | {Y_iters_const:>10.0f} {cY:>11.3g} | {ratio:>9.2f}")
    rows.append(dict(N=N, hdiv_iters=iH, hdiv_cost=cH, yano_iters=Y_iters_const, yano_cost=cY, ratio=ratio))

# refine crossover by bisection on log N
if cross is not None:
    loN, hiN = cross / 10, cross
    for _ in range(60):
        mid = np.sqrt(loN * hiN)
        (loN, hiN) = (loN, mid) if cost_Y(mid) < cost_H(mid) else (mid, hiN)
    cross = float(np.sqrt(loN * hiN))

print(f"\nCROSSOVER (yano net-cost < HDiv): N_dof ~ {cross:.2e}" if cross else "\nno crossover in range")
print("asymptotics: HDiv cost ~ N^{:.2f} (iters x storage), yano cost ~ N^{:.2f} -> yano wins for large N"
      .format(bi + bm, bY))
with open(os.path.join(HERE, "yano_vs_hdiv_crossover.json"), "w") as f:
    json.dump({"hdiv_iter_exp": bi, "hdiv_storage_exp": bm, "yano_storage_exp": bY,
               "hdiv_cost_exp": bi + bm, "yano_cost_exp": bY, "crossover_Ndof": cross, "rows": rows},
              f, indent=2, default=float)
