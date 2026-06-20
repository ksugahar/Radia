# -*- coding: utf-8 -*-
"""
act7_36_ie_prolate_spheroidal_radial.py  (Act 7 -- Gate-2 milestone 2: the TIGHT non-spherical radial IE)
========================================================================================================
act7_35 validated the elongated-body physics with the SPHERICAL IE on an enclosing sphere (correct but
NOT tight).  This file builds the radial kernel of the TIGHT (surface-conforming) non-spherical IE for
a prolate spheroid -- the prolate-SPHEROIDAL-coordinate infinite element (the spheroidal analog of the
spherical radial IE act7_25) -- and proves it closes the ellipsoid exactly.

PROLATE SPHEROIDAL COORDINATES (xi >= 1 radial-like; xi=xi0 is the spheroid surface, xi->inf is
infinity).  For an axisymmetric exterior harmonic of mode n (m=0), the exterior Laplace energy
separates into a radial (xi) part times an angular (eta) part; the radial part is

      E_n[R] = int_{xi0}^inf [ (xi^2 - 1) R'(xi)^2 + n(n+1) R(xi)^2 ] dxi ,

whose DECAYING eigenfunction is the Legendre function of the 2nd kind Q_n(xi) (Q_n->0 as xi->inf,
the spheroidal analog of r^-(n+1)).  The exact radial Steklov (DtN) is

      D_n(xi0) = -(xi0^2 - 1) Q_n'(xi0) / Q_n(xi0)  =  -n xi0 + n Q_{n-1}(xi0)/Q_n(xi0).

THE TIGHT IE: discretize E_n on the decay basis rho_k in s = xi0/xi in (0,1] (the spheroidal analog of
(a/r)^k), condense the radial DOFs -> a discrete Steklov s_P, and show s_P -> D_n(xi0).  Unlike the
SPHERE (where (a/r)^(n+1) is IN the basis so the IE is finite-EXACT), Q_n(xi) is NOT a finite
inverse-power sum, so the spheroidal IE CONVERGES (geometrically; more elongated = slower) -- the
honest non-spherical reality (the sphere is the special exact case).

PHYSICS (the demag): matching a permeable spheroid's interior uniform field to the exterior n=1
spheroidal dipole gives  N_a = (xi0^2 - 1) / (xi0^2 - 1 + xi0 D_1(xi0)) ,  which equals the Osborn
(1945) factor at e = 1/xi0 -- so the tight IE's n=1 radial Steklov D_1 IS the demag closure.

Pure numpy.  Self-asserting; writes JSON.
"""
import os
import json
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import numpy as np

N_FAIL = 0


def check(name, cond, detail=""):
    global N_FAIL
    if not cond:
        N_FAIL += 1
    print(f"  [{'ok ' if cond else 'FAIL'}] {name}{('  -- ' + detail) if detail else ''}")


# ---------------------------------------------------------------------------
# Legendre Q_n(xi) (xi>1) + exact spheroidal radial Steklov D_n(xi0)
# ---------------------------------------------------------------------------
def Qn_all(N, xi):
    Q = np.zeros(N + 1)
    Q[0] = 0.5 * np.log((xi + 1.0) / (xi - 1.0))            # arccoth(xi)
    if N >= 1:
        Q[1] = xi * Q[0] - 1.0
    for n in range(1, N):
        Q[n + 1] = ((2 * n + 1) * xi * Q[n] - n * Q[n - 1]) / (n + 1)
    return Q


def D_exact(n, xi0):
    Q = Qn_all(n, xi0)
    return -n * xi0 + n * Q[n - 1] / Q[n]                    # = -(xi0^2-1) Q_n'/Q_n


def osborn_Na(e):
    if e < 1e-9:
        return 1.0 / 3.0
    return (1.0 - e * e) / e ** 3 * (np.arctanh(e) - e)


# ---------------------------------------------------------------------------
# radial decay bases on s = xi0/xi in (0,1]
# ---------------------------------------------------------------------------
def _legval(j, x):
    c = np.zeros(j + 1); c[j] = 1.0
    return np.polynomial.legendre.legval(x, c)


def basis_eval(P, s, kind):
    """rho_k(s), rho_k'(s).  'mono': s^k ; 'nodal': s (vertex) + integrated-Legendre bubbles."""
    s = np.asarray(s, float)
    R = np.zeros((P, s.size)); Rp = np.zeros((P, s.size))
    if kind == "mono":
        for k in range(1, P + 1):
            R[k - 1] = s ** k; Rp[k - 1] = k * s ** (k - 1)
    else:
        R[0] = s; Rp[0] = np.ones_like(s)
        xi = 2.0 * s - 1.0
        for k in range(2, P + 1):
            R[k - 1] = (_legval(k, xi) - _legval(k - 2, xi)) / (2.0 * k - 1.0)
            Rp[k - 1] = _legval(k - 1, xi) * 2.0
    return R, Rp


def E_spheroidal(n, P, xi0, kind, nq=400):
    """E_n,kl = xi0 int_0^1 [ (1 - s^2/xi0^2) rho_k' rho_l' + n(n+1) rho_k rho_l / s^2 ] ds ."""
    x, w = np.polynomial.legendre.leggauss(nq)
    s = 0.5 * (x + 1.0); w = 0.5 * w
    R, Rp = basis_eval(P, s, kind)
    wr = (1.0 - s ** 2 / xi0 ** 2) * w
    wm = n * (n + 1) / s ** 2 * w
    return xi0 * ((Rp * wr) @ Rp.T + (R * wm) @ R.T)


def trace_vec(P, kind):
    R, _ = basis_eval(P, np.array([1.0]), kind)     # rho_k(s=1) = rho_k(xi=xi0)
    return R[:, 0].copy()


def discrete_steklov(n, P, xi0, kind):
    E = E_spheroidal(n, P, xi0, kind)
    g = trace_vec(P, kind)
    return 1.0 / (g @ np.linalg.solve(E, g))


print("=" * 98)
print(" act7_36 : prolate-SPHEROIDAL radial IE -- the tight ellipsoid closure (radial kernel) vs Q_n + Osborn")
print("=" * 98)

# aspect ratios via xi0 (AR = xi0/sqrt(xi0^2-1)); xi0 large -> sphere, small -> needle
CASES = [(3.0, "near-sphere"), (1.5, "moderate"), (1.1, "elongated")]

# ---- [0] exact D_n + the Osborn demag connection ----
print("\n[0] exact spheroidal Steklov D_n(xi0) and the demag connection N_a = (xi0^2-1)/(xi0^2-1 + xi0 D_1):")
print("    xi0    AR      e       D_1        N_a(from D_1)   Osborn N_a    match")
for xi0, _ in CASES:
    AR = xi0 / np.sqrt(xi0 ** 2 - 1); e = 1.0 / xi0
    D1 = D_exact(1, xi0)
    Na_from_D1 = (xi0 ** 2 - 1) / (xi0 ** 2 - 1 + xi0 * D1)
    Na_osb = osborn_Na(e)
    print(f"   {xi0:4.1f}  {AR:5.3f}  {e:5.3f}   {D1:8.4f}    {Na_from_D1:.6f}      {Na_osb:.6f}    "
          f"{abs(Na_from_D1 - Na_osb):.1e}")
    check(f"xi0={xi0}: n=1 spheroidal Steklov D_1 -> Osborn N_a (the tight demag closure)",
          abs(Na_from_D1 - Na_osb) < 1e-9, f"{abs(Na_from_D1 - Na_osb):.1e}")

# ---- [1] the tight IE radial kernel CONVERGES to the exact spheroidal Steklov D_n ----
print("\n[1] spheroidal radial IE discrete Steklov -> exact D_n(xi0) (CONVERGES; more elongated = higher P):")
conv = {}
for xi0, label in CASES:
    print(f"\n    xi0={xi0} ({label}, AR={xi0/np.sqrt(xi0**2-1):.3f}):")
    for n in (1, 2, 3):
        Dn = D_exact(n, xi0)
        errs = [abs(discrete_steklov(n, P, xi0, "nodal") - Dn) / abs(Dn) for P in (2, 4, 6, 8, 12)]
        conv[(xi0, n)] = errs
        print(f"      n={n}: D_n={Dn:8.4f}  relerr(P=2,4,6,8,12) = " + ", ".join(f"{x:.1e}" for x in errs))
        # monomial and nodal span the same space -> same discrete Steklov
        sm = discrete_steklov(n, 8, xi0, "mono"); sn = discrete_steklov(n, 8, xi0, "nodal")
        check(f"xi0={xi0} n={n}: nodal == monomial Steklov (same space)", abs(sm - sn) / abs(Dn) < 1e-6)
    check(f"xi0={xi0}: spheroidal IE CONVERGES to D_n (relerr drops >=100x from P=2 to P=12, n=1)",
          conv[(xi0, 1)][-1] < conv[(xi0, 1)][0] / 100.0,
          f"{conv[(xi0,1)][0]:.1e} -> {conv[(xi0,1)][-1]:.1e}")

# the sphere is the special FINITE-EXACT case; the spheroid converges (Q_n not a finite power sum)
check("near-sphere (xi0=3) converges much faster than elongated (xi0=1.1) -- the honest non-sphere reality",
      conv[(3.0, 1)][-1] < conv[(1.1, 1)][-1], f"xi0=3 {conv[(3.0,1)][-1]:.1e} vs xi0=1.1 {conv[(1.1,1)][-1]:.1e}")

# ---- [2] conditioning: monomial spheroidal basis Hilbert-ill, nodal bounded (the act7_28 lesson) ----
print("\n[2] conditioning of the spheroidal radial energy (xi0=1.5, n=1): monomial vs nodal:")
print("    P     monomial cond     nodal cond")
for P in (2, 4, 6, 8, 12):
    cm = np.linalg.cond(E_spheroidal(1, P, 1.5, "mono"))
    cn = np.linalg.cond(E_spheroidal(1, P, 1.5, "nodal"))
    print(f"   {P:2d}     {cm:.2e}      {cn:.2e}")
check("[2] monomial spheroidal basis is Hilbert-ill at high P (P=12 cond > 1e6)",
      np.linalg.cond(E_spheroidal(1, 12, 1.5, "mono")) > 1e6)
check("[2] nodal (orthogonalized) spheroidal basis well-conditioned (P=12 cond < 1e3)",
      np.linalg.cond(E_spheroidal(1, 12, 1.5, "nodal")) < 1e3,
      f"{np.linalg.cond(E_spheroidal(1, 12, 1.5, 'nodal')):.1e}")

# ---- [3] the demag from the IE: discrete D_1 (P=12) -> N_a -> Osborn ----
print("\n[3] demag from the TIGHT IE radial kernel (discrete D_1, P=12) -> N_a == Osborn:")
print("    xi0    AR      N_a(IE)      Osborn N_a     rel.err")
for xi0, _ in CASES:
    AR = xi0 / np.sqrt(xi0 ** 2 - 1)
    D1_ie = discrete_steklov(1, 12, xi0, "nodal")
    Na_ie = (xi0 ** 2 - 1) / (xi0 ** 2 - 1 + xi0 * D1_ie)
    Na_osb = osborn_Na(1.0 / xi0)
    re = abs(Na_ie - Na_osb) / Na_osb
    print(f"   {xi0:4.1f}  {AR:5.3f}   {Na_ie:.6f}     {Na_osb:.6f}    {re:.1e}")
    check(f"xi0={xi0}: tight-IE demag N_a == Osborn (P=12)", re < 1e-4, f"{re:.1e}")

print("\n" + "-" * 98)
print(" GATE-2 MILESTONE 2 (the tight non-spherical radial kernel):")
print("   - the prolate-SPHEROIDAL-coordinate IE (decay basis in s=xi0/xi) closes the TIGHT ellipsoid:")
print("     its discrete radial Steklov CONVERGES to the exact spheroidal D_n(xi0) = -(xi0^2-1)Q_n'/Q_n,")
print("     and its n=1 value gives N_a == the Osborn demag -- the tight conforming closure is correct.")
print("   - it CONVERGES (not finite-exact like the sphere: Q_n is not a finite inverse-power sum; more")
print("     elongated = higher radial order P), well-conditioned in the orthogonal nodal basis.")
print("   - REMAINING (milestone 3, the big build): the FULL FE spheroidal IE -- the angular (eta,phi)")
print("     surface FE x this radial kernel with the spheroidal metric -- to realise the ~AR^2 DOF edge")
print("     (act7_34) on a real mesh. This modal kernel is its correctness foundation.")
print("-" * 98)

RESULTS = {
    "cases_xi0": [c[0] for c in CASES],
    "D_exact_n1": {str(c[0]): D_exact(1, c[0]) for c in CASES},
    "Na_osborn": {str(c[0]): osborn_Na(1.0 / c[0]) for c in CASES},
    "convergence_relerr": {f"xi0={k[0]},n={k[1]}": v for k, v in conv.items()},
    "n_fail": N_FAIL,
}
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "act7_36_ie_prolate_spheroidal_radial.json")
with open(out, "w") as f:
    json.dump(RESULTS, f, indent=2)
print(f"\n  wrote {os.path.basename(out)}")

print("\n" + "=" * 98)
print(" ALL CHECKS PASSED" if N_FAIL == 0 else f" {N_FAIL} CHECK(S) FAILED")
print("=" * 98)
assert N_FAIL == 0, f"{N_FAIL} checks failed"
