# -*- coding: utf-8 -*-
"""
act7_38_ie_spheroidal_allm.py  (Act 7 -- the ALL-m prolate-spheroidal IE: axial + transverse demag)
==================================================================================================
act7_37 built the FULL FE spheroidal IE for the axisymmetric (m=0) demag.  This file establishes the
IE for ALL azimuthal modes m (the genuinely 3-D, non-axisymmetric capability) by validating its radial
closure against the exact associated-Legendre functions Q_n^m AND deriving BOTH the axial (m=0) and
transverse (m=1) Osborn demag factors from the IE's radial Steklov.

ALL-m RADIAL KERNEL.  For azimuthal mode m the exterior energy is (the d/dphi term adds m^2; the
(xi^2-eta^2) factor splits as 1/(1-eta^2)+1/(xi^2-1), so it stays separable):

      E_m = pi f INT INT [ (xi^2-1)(d_xi u)^2 + (1-eta^2)(d_eta u)^2 + m^2( 1/(1-eta^2) + 1/(xi^2-1) ) u^2 ]

-> per angular mode P_n^m(eta) the radial energy is  A_xi + n(n+1) M_xi + m^2 W2_xi  (W2_xi = INT
rho rho /(xi^2-1)), whose decaying eigenfunction is Q_n^m(xi) and whose exact Steklov is
D_{n,m}(xi0) = -(xi0^2-1) Q_n^m'(xi0)/Q_n^m(xi0).  The IE (decay basis in s=xi0/xi) reproduces D_{n,m}.

THE TWO DEMAGS FROM THE IE (matching the permeable spheroid; derived):
      axial      (m=0):  N_a = (xi0^2-1)/(xi0^2-1 + xi0 D_{1,0})
      transverse (m=1):  N_b =  xi0 / (xi0 + D_{1,1})
both == the Osborn (1945) factors (with the sum rule N_a + 2 N_b = 1), so the IE closes the spheroid
for BOTH axial and transverse fields -- the all-m capability, validated modally.

NOTE on the full all-m FEM (vs the modal validation here).  The m=0 axial demag has a clean full FE
solve (act7_37).  For m != 0 the field carries the associated-Legendre factor sqrt((xi^2-1)(1-eta^2))
(vanishing on the focal segment xi=1 and the poles eta=+-1), so a NAIVE nodal FE converges slowly at
those coordinate singularities (a sqrt behaviour); a clean full m!=0 FE solve needs the
sqrt-lifted basis u = sqrt((xi^2-1)(1-eta^2)) v.  That is a DISCRETISATION refinement -- it does NOT
affect the IE's correctness, which is the radial closure proven here (exact Q_n^m, Osborn N_a + N_b).

Pure numpy + scipy.special (associated Legendre Q_n^m).  Self-asserting; writes JSON.
"""
import os
import json
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import numpy as np
from scipy.special import lqmn

N_FAIL = 0


def check(name, cond, detail=""):
    global N_FAIL
    if not cond:
        N_FAIL += 1
    print(f"  [{'ok ' if cond else 'FAIL'}] {name}{('  -- ' + detail) if detail else ''}")


def osborn_Na(e):
    if e < 1e-9:
        return 1.0 / 3.0
    return (1.0 - e * e) / e ** 3 * (np.arctanh(e) - e)


def D_exact(n, m, xi0):
    """Exact spheroidal radial Steklov D_{n,m} = -(xi0^2-1) Q_n^m'(xi0)/Q_n^m(xi0)."""
    Q, Qp = lqmn(m, n, xi0)            # Q[mi,ni], Qp = dQ/dx, x>1
    return -(xi0 ** 2 - 1.0) * Qp[m, n] / Q[m, n]


# ---------------------------------------------------------------------------
# the IE radial operators (nodal decay basis in s = xi0/xi): A, M, W2
# ---------------------------------------------------------------------------
def _legval(j, x):
    c = np.zeros(j + 1); c[j] = 1.0
    return np.polynomial.legendre.legval(x, c)


def ie_ops(P, xi0, nq=800):
    x, w = np.polynomial.legendre.leggauss(nq); s = 0.5 * (x + 1.0); w = 0.5 * w
    R = np.zeros((P, s.size)); Rp = np.zeros((P, s.size))
    R[0] = s; Rp[0] = 1.0
    xi = 2.0 * s - 1.0
    for k in range(2, P + 1):
        R[k - 1] = (_legval(k, xi) - _legval(k - 2, xi)) / (2.0 * k - 1.0)
        Rp[k - 1] = _legval(k - 1, xi) * 2.0
    A = xi0 * (Rp * ((1.0 - s ** 2 / xi0 ** 2) * w)) @ Rp.T      # int (xi^2-1) rho' rho'
    M = xi0 * (R * (w / s ** 2)) @ R.T                          # int rho rho
    W2 = xi0 * (R * (w / (xi0 ** 2 - s ** 2))) @ R.T            # int rho rho /(xi^2-1)
    return A, M, W2


def ie_steklov(n, m, P, xi0):
    A, M, W2 = ie_ops(P, xi0)
    E = A + n * (n + 1) * M + m * m * W2                        # radial energy for mode (n,m)
    g = np.zeros(P); g[0] = 1.0                                 # trace at xi0
    return 1.0 / (g @ np.linalg.solve(E, g))


print("=" * 96)
print(" act7_38 : the ALL-m prolate-spheroidal IE -- radial Steklov Q_n^m + the axial & transverse demag")
print("=" * 96)

CASES = [(3.0, "near-sphere"), (1.5, "moderate"), (1.1, "elongated"), (1.03, "needle")]

# ---- [1] the IE radial Steklov reproduces the exact Q_n^m for m=0,1,2 ----
print("\n[1] IE discrete radial Steklov -> exact D_{n,m}(xi0) = -(xi0^2-1)Q_n^m'/Q_n^m (P=14):")
print("    xi0    (n,m)    D_exact     IE(P=14)    rel.err")
for xi0, _ in CASES[:3]:
    for (n, m) in [(1, 0), (1, 1), (2, 1), (2, 2)]:
        De = D_exact(n, m, xi0); di = ie_steklov(n, m, 14, xi0)
        re = abs(di - De) / abs(De)
        print(f"   {xi0:4.2f}   ({n},{m})   {De:8.4f}   {di:8.4f}    {re:.1e}")
        check(f"xi0={xi0} (n,m)=({n},{m}): IE Steklov == exact Q_n^m", re < 1e-3, f"{re:.1e}")

# ---- [2] BOTH demags from the IE radial Steklov: axial N_a (m=0) AND transverse N_b (m=1) == Osborn ----
print("\n[2] both demags from the IE radial Steklov -> Osborn (axial m=0 + transverse m=1):")
print("    xi0    AR      N_a(IE)     Osborn N_a    N_b(IE)     Osborn N_b    sum N_a+2N_b")
for xi0, _ in CASES:
    AR = xi0 / np.sqrt(xi0 ** 2 - 1.0); e = 1.0 / xi0
    P_dem = 20                                                 # AR-aware: elongated needs higher radial P
    D1 = ie_steklov(1, 0, P_dem, xi0)
    D11 = ie_steklov(1, 1, P_dem, xi0)
    Na_ie = (xi0 ** 2 - 1.0) / (xi0 ** 2 - 1.0 + xi0 * D1)      # axial
    Nb_ie = xi0 / (xi0 + D11)                                   # transverse
    Na_o = osborn_Na(e); Nb_o = (1.0 - Na_o) / 2.0
    sm = Na_ie + 2.0 * Nb_ie
    print(f"   {xi0:4.2f}  {AR:5.3f}   {Na_ie:.6f}    {Na_o:.6f}    {Nb_ie:.6f}    {Nb_o:.6f}    {sm:.6f}")
    check(f"xi0={xi0}: IE axial demag N_a == Osborn", abs(Na_ie - Na_o) < 1e-4, f"{abs(Na_ie-Na_o):.1e}")
    check(f"xi0={xi0}: IE transverse demag N_b == Osborn", abs(Nb_ie - Nb_o) < 1e-4, f"{abs(Nb_ie-Nb_o):.1e}")
    check(f"xi0={xi0}: demag sum rule N_a + 2 N_b == 1", abs(sm - 1.0) < 1e-4, f"{abs(sm-1.0):.1e}")

print("\n" + "-" * 96)
print(" THE ALL-m IE:")
print("   - the spheroidal IE closes the exterior for EVERY azimuthal mode m: its radial Steklov")
print("     reproduces the exact associated-Legendre Q_n^m(xi) (m=0,1,2), so the open boundary is exact")
print("     for axial (m=0), transverse (m=1) and general fields;")
print("   - BOTH Osborn demags follow from the IE radial Steklov: N_a=(xi0^2-1)/(xi0^2-1+xi0 D_{1,0}) and")
print("     N_b=xi0/(xi0+D_{1,1}), == Osborn with the sum rule N_a+2N_b=1 -- the IE is a genuine 3-D")
print("     (all-m) open-boundary closure for the spheroid.")
print("   - (the full m!=0 FE solve needs the sqrt-lifted basis at the focal/pole singularities -- a")
print("     discretisation refinement; the IE's correctness is the radial closure proven here.)")
print("-" * 96)

RESULTS = {
    "cases_xi0": [c[0] for c in CASES],
    "demag": {str(c[0]): dict(
        AR=c[0] / np.sqrt(c[0] ** 2 - 1.0),
        Na_ie=(c[0] ** 2 - 1.0) / (c[0] ** 2 - 1.0 + c[0] * ie_steklov(1, 0, 20, c[0])),
        Nb_ie=c[0] / (c[0] + ie_steklov(1, 1, 20, c[0])),
        Na_osborn=osborn_Na(1.0 / c[0]),
        Nb_osborn=(1.0 - osborn_Na(1.0 / c[0])) / 2.0,
    ) for c in CASES},
    "n_fail": N_FAIL,
}
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "act7_38_ie_spheroidal_allm.json")
with open(out, "w") as fh:
    json.dump(RESULTS, fh, indent=2)
print(f"\n  wrote {os.path.basename(out)}")

print("\n" + "=" * 96)
print(" ALL CHECKS PASSED" if N_FAIL == 0 else f" {N_FAIL} CHECK(S) FAILED")
print("=" * 96)
assert N_FAIL == 0, f"{N_FAIL} checks failed"
