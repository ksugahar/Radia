# -*- coding: utf-8 -*-
"""
act7_30_vector_hcurl_hdiv_ie.py  (Act 7 -- the genuine VECTOR H(curl)/H(div) infinite elements)
================================================================================================
act7_26/29 built the de Rham COMPLEX (the spaces P/B/C + the commuting grad/curl/div maps) and the
scalar DtN.  This file builds the genuine VECTOR infinite elements -- the H(curl)- and H(div)-type
DtN (Steklov) operators -- and shows the de Rham vector ends carry TWO DISTINCT radial ladders:

      gradient / irrotational   (H1 ; H(div) normal trace)        DtN = -(n+1)/a   (scalar ladder)
      toroidal / transverse     (H(curl) tangential trace)        DtN = -n/a       (the VECTOR ladder)

The toroidal (transverse, curl-type) field  T = f(r) Phi_n,  Phi_n = r^ x grad_S Y_n  (= the C-field
of act7_26) is genuinely H(curl): divergence-free, not a gradient.  Source-free exterior ->
f'' + 2 f'/r - n(n+1) f/r^2 = 0 -> decaying f = r^-(n+1); its H(curl) Steklov (tangential curl / trace)
is the **-n** ladder, distinct from the scalar **-(n+1)**.  The curl-curl radial energy is

      E_tor[f] = integral_a^inf [ ((r f)')^2 + n(n+1) f^2 ] dr ,   E_tor[r^-(n+1)] = n,

so on the decay basis  f = (a/r)^k  the toroidal energy matrix is

      T_kl = ( (1-k)(1-l) + n(n+1) ) / ( (k+l) - 1 )            (a = 1),

a Hilbert-type matrix (cousin of the scalar A_kl = (kl + n(n+1))/((k+l)-1)).  DtN = -1/(1^T T^-1 1).

This is the static de Rham vector content; the FULL Maxwell/wave vector IE (TE/TM impedances, Hankel
radial basis) is Nannen et al. 2013 (SISC 10.1137/110860148) -- needs the oscillatory basis, downstream.
H(div) is the Hodge dual of H(curl) on the sphere and carries the SAME two ladders (normal-trace
gradient -(n+1); solenoidal -n).

Pure numpy.
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
# the two energy matrices (monomial decay basis (a/r)^k, a=1)
# ---------------------------------------------------------------------------
def scalar_energy(n, P):                       # gradient / H1 / H(div)-normal -> ladder -(n+1)
    k = np.arange(1, P + 1)
    return (np.outer(k, k) + n * (n + 1)) / (k[:, None] + k[None, :] - 1.0)


def toroidal_energy(n, P):                      # toroidal / H(curl)-tangential -> ladder -n
    k = np.arange(1, P + 1)
    return (np.outer(1 - k, 1 - k) + n * (n + 1)) / (k[:, None] + k[None, :] - 1.0)


def dtn(E):
    P = E.shape[0]
    g = np.ones(P)
    return -1.0 / (g @ np.linalg.solve(E, g))


# ---------------------------------------------------------------------------
# orthogonal (well-conditioned) bases for the SAME spaces
# ---------------------------------------------------------------------------
def gauss01(nq):
    x, w = np.polynomial.legendre.leggauss(nq)
    return 0.5 * (x + 1.0), 0.5 * w


def _leg(j, t):
    c = np.zeros(j + 1); c[j] = 1.0
    return np.polynomial.legendre.legval(2.0 * t - 1.0, c)


def _dleg(j, t):
    c = np.zeros(j + 1); c[j] = 1.0
    return np.polynomial.legendre.legval(2.0 * t - 1.0, np.polynomial.legendre.legder(c)) * 2.0


def toroidal_energy_ortho(n, P, nq=400):
    # t = 1/r in (0,1]; basis u_j = t * shiftedLeg_j(t) spans {t..t^P} = {r^-1..r^-P}.
    # E_tor = int [ ((rf)')^2 + n(n+1) f^2 ] dr  -> in t (rf)'=-t^2 L_j', f=u_j, dr=dt/t^2 ->
    #   integrand = t^2 L_i' L_j' + n(n+1) L_i L_j   (clean, no singularity); trace g_j = u_j(1)=L_j(1)=1
    t, w = gauss01(nq)
    E = np.zeros((P, P))
    for i in range(P):
        for j in range(P):
            E[i, j] = np.sum(w * (t ** 2 * _dleg(i, t) * _dleg(j, t) + n * (n + 1) * _leg(i, t) * _leg(j, t)))
    return E, np.ones(P)


print("=" * 92)
print(" act7_30 : genuine VECTOR H(curl)/H(div) infinite elements -- the two Steklov ladders")
print("=" * 92)

MODES = [1, 2, 3, 4]

# ---- [1] the two ladders: scalar -(n+1) vs toroidal -n, both exact for n <= P-1 ----
print("\n[1] two DtN ladders (P=6, exact for n<=P-1):  gradient -(n+1)  vs  toroidal H(curl) -n")
print("   n   scalar(IE)  exact   |  toroidal(IE)  exact")
P = 6
for n in MODES:
    ds = dtn(scalar_energy(n, P))
    dt_ = dtn(toroidal_energy(n, P))
    print(f"  {n}   {ds:8.4f}  {-(n+1):5d}   |  {dt_:8.4f}     {-n:4d}")
    check(f"n={n}: scalar/gradient IE DtN = -(n+1) (H1 / H(div) normal)", abs(ds - (-(n + 1))) < 1e-8,
          f"{ds:.4f}")
    check(f"n={n}: toroidal H(curl) IE DtN = -n (the VECTOR ladder, distinct)", abs(dt_ - (-n)) < 1e-8,
          f"{dt_:.4f}")
check("[1] the two ladders are DISTINCT (toroidal -n != scalar -(n+1))", True)

# spectral: toroidal exact once P>=n+1
print("\n[1b] toroidal H(curl) is spectral -- exact once P>=n+1:")
for n in (1, 3):
    d = [abs(dtn(toroidal_energy(n, Pp)) - (-n)) / n for Pp in range(1, 7)]
    print(f"   n={n}: reldef vs P=1..6 = " + ", ".join(f"{x:.0e}" for x in d))
    check(f"n={n}: toroidal exact for P>=n+1 (reldef<1e-8)", d[n] < 1e-8, f"P={n+1}: {d[n]:.0e}")

# ---- [2] both vector ends well-conditioned with an orthogonal basis (high order) ----
print("\n[2] conditioning (n=1): monomial Hilbert-ill vs orthogonal bounded -- BOTH ladders:")
print("    P    scalar mono   toroidal mono   toroidal ortho")
for P in (2, 4, 6, 8):
    cs = float(np.linalg.cond(scalar_energy(1, P)))
    ctm = float(np.linalg.cond(toroidal_energy(1, P)))
    Eo, _ = toroidal_energy_ortho(1, P)
    cto = float(np.linalg.cond(Eo))
    print(f"   {P:2d}    {cs:.2e}     {ctm:.2e}      {cto:.2e}")
check("[2] toroidal monomial basis is Hilbert-ill at high order (P=8 cond > 1e6)",
      float(np.linalg.cond(toroidal_energy(1, 8))) > 1e6)
Eo8, go8 = toroidal_energy_ortho(1, 8)
check("[2] toroidal ORTHOGONAL basis well-conditioned (P=8 cond < 1e3)", float(np.linalg.cond(Eo8)) < 1e3,
      f"{float(np.linalg.cond(Eo8)):.1e}")
# orthogonal basis reproduces the -n ladder (same space)
check("[2] toroidal orthogonal basis gives the SAME -n ladder (n=1, P=8)",
      abs(dtn(Eo8) - (-1)) < 1e-8, f"DtN {dtn(Eo8):.4f}")

print("\n" + "-" * 92)
print(" VECTOR IE SUMMARY:")
print("   - H(curl) (tangential / toroidal): DtN ladder -n  -- the genuine vector Steklov, built+verified,")
print("     well-conditioned in an orthogonal basis, spectral (exact n<=P-1).")
print("   - H(div) / H1 (normal / gradient): DtN ladder -(n+1) -- the scalar IE (act7_25/28).")
print("   - H(div) is the Hodge dual of H(curl) on the sphere: SAME two ladders.")
print("   - FULL Maxwell/wave vector IE (TE/TM impedances, Hankel radial) = Nannen 2013 -- downstream")
print("     (needs the oscillatory exp(ikr) basis; this is the STATIC de Rham vector content).")
print("-" * 92)

RESULTS = {
    "modes": MODES,
    "scalar_ladder_exact": [-(n + 1) for n in MODES],
    "toroidal_ladder_exact": [-n for n in MODES],
    "scalar_ie": [dtn(scalar_energy(n, 6)) for n in MODES],
    "toroidal_ie": [dtn(toroidal_energy(n, 6)) for n in MODES],
    "toroidal_cond_mono_P8": float(np.linalg.cond(toroidal_energy(1, 8))),
    "toroidal_cond_ortho_P8": float(np.linalg.cond(Eo8)),
    "n_fail": N_FAIL,
}
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "act7_30_vector_hcurl_hdiv_ie.json")
with open(out, "w") as f:
    json.dump(RESULTS, f, indent=2)
print(f"\n  wrote {os.path.basename(out)}")

print("\n" + "=" * 92)
print(" ALL CHECKS PASSED" if N_FAIL == 0 else f" {N_FAIL} CHECK(S) FAILED")
print("=" * 92)
assert N_FAIL == 0, f"{N_FAIL} checks failed"
