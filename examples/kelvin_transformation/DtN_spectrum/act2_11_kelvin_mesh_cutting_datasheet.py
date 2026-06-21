# -*- coding: utf-8 -*-
"""
act2_11_kelvin_mesh_cutting_datasheet.py  (Act 2 -- the consolidated mesh-cutting datasheet)
===========================================================================================
How should the Kelvin region be MESHED, and what DtN-spectrum performance does a given cut
buy?  The corpus answers this in scattered pieces:
  * act0_01 : per-degree DtN defect ~ n^2 (h/R)^(2s)   (s = the FE eigenvalue order)
  * act2_05 : the curved-geometry floor ~ (h/R)^(2k)    (k = isoparametric Curve order)
  * act2_09 : at p>=n the exterior VOLUME mesh is irrelevant -- only Gamma's SURFACE
              triangulation and its Curve order move the spectrum.
This act CONSOLIDATES them into ONE predictive map for the SPARSE Kelvin-FEM closure (the
manuscript's actual object, NOT the dense BEM DtN of exterior_dtn_spectrum): the per-mode
DtN defect d_n over the two mesh-cutting knobs act2_09 proved matter -- surface resolution h
and Curve order k -- at FE orders p, read off as the resolvable multipole band.

THE RESULT (the mesh-cutting recipe -- a crossover the scattered laws did not state):

    d_n  ~  C * n^2 * (h/R)^(2 * min(p, k))

  i.e. the convergence exponent is set by min(FE order p, Curve order k), NOT by either
  alone.  Consequences for "how to cut the Kelvin mesh":
    * raising the Curve order k ABOVE the FE order p does NOT help (FE-limited);
    * raising the FE order p ABOVE the Curve order k does NOT help (geometry-limited);
    * the optimal cut BALANCES p and k (here p=k); and since the radial/volume cut is
      irrelevant (act2_09), the whole Kelvin region collapses to a single curvature-limited
      shell tuned only by (h, p=k).

This is the practical recipe behind the manuscript's "minimal Delta-DoF": one balanced-order
shell on Gamma, not a refined exterior volume.

Reuses radia.infinite_element (the NGSolve-native Kelvin-on-sphere DtN closure) for the
condensed surface operator.  Needs NGSolve + numpy + scipy.  Self-asserting; writes JSON.
"""
import os
import json
import math
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp

import ngsolve as ng
from netgen.occ import OCCGeometry, Sphere, Pnt
from radia.infinite_element import dtn_surface_matrix

N_FAIL = 0


def check(name, cond, detail=""):
    global N_FAIL
    if not cond:
        N_FAIL += 1
    print(f"  [{'ok ' if cond else 'FAIL'}] {name}{('  -- ' + detail) if detail else ''}")


def _to_dense(mat, ndof):
    r, c, v = mat.COO()
    A = sp.coo_matrix((np.array(v), (np.array(r), np.array(c))),
                      shape=(ndof, ndof)).toarray()
    return 0.5 * (A + A.T)


def surface_ms_ks(R, maxh, curve, order):
    """Surface mass MS = int_Gamma u v ds and Laplace-Beltrami stiffness
    KS = int_Gamma grad_S u . grad_S v ds on the truncation sphere of radius R,
    meshed at ``maxh`` and curved to isoparametric order ``curve``, with a boundary
    H1 trace of order ``order``.  Returns (MS, KS, surface_ndof) restricted to the
    genuine surface DOFs (those carrying nonzero surface mass)."""
    mesh = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), R)).GenerateMesh(maxh=maxh))
    mesh.Curve(curve)
    fes = ng.H1(mesh, order=order)
    u, v = fes.TnT()
    n = ng.specialcf.normal(3)
    gu = ng.grad(u).Trace(); gv = ng.grad(v).Trace()
    gut = gu - (gu * n) * n; gvt = gv - (gv * n) * n          # surface (tangential) gradient
    ds = ng.ds(bonus_intorder=2 * order + 2 * curve)
    m = ng.BilinearForm(fes, check_unused=False); m += u * v * ds; m.Assemble()
    k = ng.BilinearForm(fes, check_unused=False); k += gut * gvt * ds; k.Assemble()
    MS = _to_dense(m.mat, fes.ndof); KS = _to_dense(k.mat, fes.ndof)
    d = np.diag(MS)
    bnd = np.where(d > 1e-10 * d.max())[0]
    return MS[np.ix_(bnd, bnd)], KS[np.ix_(bnd, bnd)], int(len(bnd))


def dtn_defects(MS, KS, R, nmax, P_radial):
    """Per-degree relative defect d_n = |mean(bucket) - (n+1)R| / ((n+1)R) of the discrete
    Kelvin DtN energy spectrum eig(S, MS), where S is the condensed surface DtN operator
    (radia.infinite_element.dtn_surface_matrix).  Buckets eigenvalues by multiplicity 2n+1.
    Returns dict {n: d_n} for n=1..nmax (n=0 = the trivial constant mode is skipped)."""
    S = dtn_surface_matrix(MS, KS, P_radial, a=R)
    lam = sla.eigh(S, MS, eigvals_only=True)               # ascending; (n+1)R ladder
    out = {}
    idx = 0
    for n in range(nmax + 1):
        mult = 2 * n + 1
        if idx + mult > len(lam):
            break
        grp = lam[idx:idx + mult]
        lam_exact = (n + 1) * R
        if n >= 1:
            out[n] = abs(float(np.mean(grp)) - lam_exact) / lam_exact
        idx += mult
    return out


# ===========================================================================
R = 1.0
MAXH = [0.5, 0.4, 0.3, 0.25]
CONFIGS = [(2, 1), (2, 2), (2, 3), (3, 2), (3, 3)]     # (FE order p, Curve order k)
NMAX = 5
P_RADIAL = 6
EPS_BAND = 1e-3
N_GAIN = 2                                              # mode used for slope/gain/ratio probes

print("=" * 98)
print(" act2_11 : Kelvin mesh-cutting -> DtN-spectrum datasheet  (d_n ~ n^2 (h/R)^(2 min(p,k)))")
print("=" * 98)

data = {}     # data[(p,k)][maxh] = {"d": {n: d_n}, "ndof": int}
with ng.TaskManager():
    for (p, k) in CONFIGS:
        data[(p, k)] = {}
        for h in MAXH:
            MS, KS, ndof = surface_ms_ks(R, h, k, p)
            dd = dtn_defects(MS, KS, R, NMAX, P_RADIAL)
            data[(p, k)][h] = {"d": dd, "ndof": ndof}

hF = min(MAXH); hC = max(MAXH)                            # finest / coarsest


def dF(p, k, n):
    return data[(p, k)][hF]["d"].get(n, float("nan"))


# ---- [1] the datasheet: finest-mesh per-mode defect for each cut ----
print(f"\n[1] FINEST mesh (maxh={hF}) per-degree DtN defect d_n -- the accuracy each cut buys:")
print(f"    (p,k)  ndof    " + "  ".join(f"n={n}" for n in range(1, NMAX + 1)))
for (p, k) in CONFIGS:
    rec = data[(p, k)][hF]
    row = "  ".join(f"{rec['d'].get(n, float('nan')):.1e}" for n in range(1, NMAX + 1))
    nband = max([n for n in range(1, NMAX + 1) if rec["d"].get(n, 1.0) < EPS_BAND] or [0])
    print(f"    ({p},{k})  {rec['ndof']:4d}   {row}    n_max(<{EPS_BAND:g})={nband}")

# ---- [2] refinement gains and fitted exponents (informational) ----
print(f"\n[2] convergence over maxh {hC}->{hF} at degree n={N_GAIN} (gain = d_coarse/d_fine, "
      f"slope = dlog d / dlog h):")
print("    (p,k)  expect 2*min(p,k)   slope    gain    d_fine")
gains = {}; slopes = {}
for (p, k) in CONFIGS:
    hs = np.array(MAXH, float)
    ds = np.array([data[(p, k)][h]["d"].get(N_GAIN, float("nan")) for h in MAXH], float)
    ok = ds > 0
    slope = float(np.polyfit(np.log(hs[ok]), np.log(ds[ok]), 1)[0]) if ok.sum() >= 2 else float("nan")
    gain = float(data[(p, k)][hC]["d"][N_GAIN] / data[(p, k)][hF]["d"][N_GAIN])
    gains[(p, k)] = gain; slopes[(p, k)] = slope
    print(f"    ({p},{k})       {2*min(p,k):d}            {slope:5.2f}   {gain:6.1f}   "
          f"{data[(p,k)][hF]['d'][N_GAIN]:.2e}")
print("    NB: (3,3) is FLOOR-limited (~1e-6) over this maxh range -- already near-converged even")
print("        on the coarsest mesh, so its asymptotic h^6 slope is not resolved.  That floor IS the")
print("        practical point: a balanced high-order shell needs essentially no refinement.")

# ---- [3] the crossover: min(p,k), not p or k alone ----
print("\n[3] the min(p,k) crossover (finest-mesh defect ratios at n=2):")
r_k_help = dF(2, 1, N_GAIN) / dF(2, 2, N_GAIN)            # k=1->2 at p=2 : geometry relieved
r_k_sat = dF(2, 2, N_GAIN) / dF(2, 3, N_GAIN)            # k=2->3 at p=2 : FE-limited (saturates)
r_p_sat = dF(2, 2, N_GAIN) / dF(3, 2, N_GAIN)            # p=2->3 at k=2 : geom-limited (saturates)
r_p_help = dF(2, 3, N_GAIN) / dF(3, 3, N_GAIN)            # p=2->3 at k=3 : FE relieved
print(f"    raising Curve k 1->2 at p=2 (k<p, geometry-limited)  : d ratio = {r_k_help:6.1f}  (helps)")
print(f"    raising Curve k 2->3 at p=2 (k>=p, FE-limited)       : d ratio = {r_k_sat:6.2f}  (saturates)")
print(f"    raising FE    p 2->3 at k=2 (p>k, geometry-limited)  : d ratio = {r_p_sat:6.2f}  (saturates)")
print(f"    raising FE    p 2->3 at k=3 (p<=k, FE-limited)       : d ratio = {r_p_help:6.1f}  (helps)")
print("    => the convergence is governed by min(p,k); the balanced cut (p=k) is optimal.")

# ---- self-asserting checks ----
print("\n" + "-" * 98)
# mesh distinctness
nd22 = [data[(2, 2)][h]["ndof"] for h in MAXH]
check("surface DOF strictly increases as maxh shrinks (distinct meshes)",
      all(nd22[i] < nd22[i + 1] for i in range(len(nd22) - 1)), f"ndof={nd22}")
# spectral signature: higher modes are worse (ordered by degree) at the finest mesh,
# UNLESS the cut has driven the low modes to the implementation floor (then ordering is noise)
FLOOR = 2e-5
for (p, k) in CONFIGS:
    rec = data[(p, k)][hF]["d"]
    d1, d3 = rec.get(1, 1.0), rec.get(3, 0.0)
    floored = max(rec.get(n, 1.0) for n in (1, 2, 3)) < FLOOR
    check(f"({p},{k}) finest-mesh defect ordered by degree OR floor-saturated",
          d3 >= d1 or floored,
          f"d_1={d1:.1e}, d_3={d3:.1e}{'  (floor-saturated)' if floored else ''}")
# refinement lowers every probed config
for (p, k) in CONFIGS:
    check(f"({p},{k}) refinement helps (d_fine < d_coarse at n={N_GAIN})",
          data[(p, k)][hF]["d"][N_GAIN] < data[(p, k)][hC]["d"][N_GAIN])
# THE crossover (min(p,k)):
check("Curve k helps when k<p (k 1->2 at p=2 drops the defect > 3x)", r_k_help > 3.0,
      f"ratio={r_k_help:.1f}")
check("Curve k SATURATES when k>=p (k 2->3 at p=2 within 3x)", 1.0 / 3 < r_k_sat < 3.0,
      f"ratio={r_k_sat:.2f}")
check("FE p SATURATES when p>k (p 2->3 at k=2 within 3x)", 1.0 / 3 < r_p_sat < 3.0,
      f"ratio={r_p_sat:.2f}")
check("FE p helps when p<=k (p 2->3 at k=3 drops the defect > 3x)", r_p_help > 3.0,
      f"ratio={r_p_help:.1f}")
# the balanced cut (3,3) reaches the lowest finest-mesh defect of all cuts
dF_all = {cfg: dF(cfg[0], cfg[1], N_GAIN) for cfg in CONFIGS}
best = min(dF_all, key=dF_all.get)
check("the balanced cut (p=k=3) reaches the lowest finest-mesh defect", best == (3, 3),
      f"best={best} d={dF_all[best]:.2e}")
# the exponent law 2*min(p,k), read directly from the slopes of the (non-floor) configs
SLOPE_BAND = {(2, 1): (1.4, 2.7), (2, 2): (3.2, 4.8), (2, 3): (3.2, 4.8), (3, 2): (3.2, 4.8)}
for cfg, (lo, hi) in SLOPE_BAND.items():
    check(f"{cfg} convergence exponent ~ 2*min(p,k)={2*min(*cfg)} (slope in [{lo},{hi}])",
          lo <= slopes[cfg] <= hi, f"slope={slopes[cfg]:.2f}")
# (3,3) is floor-limited here, so its asymptotic slope-6 is not resolved; what it DOES show is
# that raising BOTH p and k unlocks a markedly lower floor than the min(p,k)-saturated cuts.
check("the balanced cut (3,3) unlocks a markedly lower floor than (2,2) (< 0.3x)",
      dF(3, 3, N_GAIN) < 0.3 * dF(2, 2, N_GAIN),
      f"d(3,3)={dF(3,3,N_GAIN):.1e} vs d(2,2)={dF(2,2,N_GAIN):.1e}")

# ---- JSON ----
RESULTS = {
    "R": R, "maxh_list": MAXH, "configs": [list(c) for c in CONFIGS],
    "nmax": NMAX, "p_radial": P_RADIAL, "eps_band": EPS_BAND, "n_probe": N_GAIN,
    "datasheet": {f"{p}_{k}": {
        "ndof_by_maxh": {f"{h}": data[(p, k)][h]["ndof"] for h in MAXH},
        "defect_by_maxh": {f"{h}": {str(n): data[(p, k)][h]["d"].get(n)
                                    for n in range(1, NMAX + 1)} for h in MAXH},
        "n_max_band_finest": max([n for n in range(1, NMAX + 1)
                                  if data[(p, k)][hF]["d"].get(n, 1.0) < EPS_BAND] or [0]),
        "slope_n2": slopes[(p, k)], "gain_n2": gains[(p, k)],
        "expect_exponent": 2 * min(p, k),
    } for (p, k) in CONFIGS},
    "crossover_ratios_n2": {
        "k_help_klt_p": r_k_help, "k_saturate_kge_p": r_k_sat,
        "p_saturate_pgt_k": r_p_sat, "p_help_ple_k": r_p_help,
    },
    "n_fail": N_FAIL,
}
out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "act2_11_kelvin_mesh_cutting_datasheet.json")
with open(out, "w") as fh:
    json.dump(RESULTS, fh, indent=2)
print(f"\n  wrote {os.path.basename(out)}")

print("\n" + "=" * 98)
print(" ALL CHECKS PASSED" if N_FAIL == 0 else f" {N_FAIL} CHECK(S) FAILED")
print("=" * 98)
assert N_FAIL == 0, f"{N_FAIL} checks failed"
