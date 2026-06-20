# -*- coding: utf-8 -*-
"""
act7_26_derham_infinite_element.py  (Act 7 -- the open-boundary comparison; the de Rham / exact-sequence IE)
===========================================================================================================
Extends act7_25 (the SCALAR Bettess infinite element) to a **de Rham / exact-sequence infinite
element**: the  H1 --grad--> H(curl) --curl--> H(div) --div--> L2  complex on the spherical exterior
r >= a, with the RADIAL decay families chosen so that grad / curl / div COMMUTE (the commuting
diagram).  This is the construction of Demkowicz & Pal, "An infinite element for Maxwell's equations"
(CMAME 164, 1998) and the exact-sequence hp infinite element (Demkowicz, "Computing with hp-Adaptive
Finite Elements").

WHAT IS VERIFIED (symbolically with sympy -- EXACT, not sampled):
  Exterior, spherical-harmonic mode n (m=0 representative; the radial structure is m-independent).
  Building blocks (orthonormal spherical components (r^, th^, ph^),  Y = P_n(cos th)):
      0-form  scalar    u    = r^{-k} Y
      1-/2-form vectors  P(k) = (r^{-k} Y,  0,          0)             radial
                         B(k) = (0,         r^{-k} dY,  0)             poloidal-tangential
                         C(k) = (0,         0,          r^{-k} dY)     toroidal
      3-form  scalar         = r^{-k} Y
  RADIAL FAMILIES shift by +1 per form degree (THE design that makes the diagram commute):
      S0 = {n+1, ..., n+P}     S1 = S0+1     S2 = S0+2     S3 = S0+3
  The exact structure constants (each asserted == 0 by sympy):
      grad( r^{-k} Y ) = -k * P(k+1)  +  1 * B(k+1)                    ( -> V1, powers in S1)
      curl( P(k) )     = -1 * C(k+1)                                   ( -> V2)
      curl( B(k) )     = (1-k) * C(k+1)                                ( -> V2)
      curl( C(k) )     = -n(n+1) * P(k+1)  +  (k-1) * B(k+1)           ( -> V2)   [Legendre eig]
      div ( P(k) )     = (2-k) * (r^{-(k+1)} Y)                        ( -> V3)
      div ( B(k) )     = -n(n+1) * (r^{-(k+1)} Y)                      ( -> V3)   [Legendre eig]
      div ( C(k) )     = 0                                             ( -> V3)
  + the complex property  curl(grad u) == 0,  div(curl A) == 0  (operator sanity, exact).

The  -n(n+1)  that keeps the TOROIDAL curl (C) and div(B) closed is the angular-Laplace eigenvalue
  (1/sin th) d/dth( sin th  dY/dth ) = -n(n+1) Y     (Legendre's equation)
-- the heart of why the decay families commute.  So the de Rham infinite element EXISTS and its
radial families are exhibited + verified here.

SCOPE / honesty: this is the SEPARABLE-sphere, m=0 demonstration of the exact-sequence RADIAL
construction (the part the question names: "choose the radial families so grad/curl/div commute").
The full element (all m; a general FE truncation surface instead of analytic harmonics; the eddy/wave
decay basis exp(-q r)) is Demkowicz-Pal -- cited, not re-implemented.  This is a KNOWN construction;
we implement + verify it, no novelty claimed.  The practical SHIPPED de Rham open boundary is instead
the coordinate-MAPPING family (the Kelvin transformation here; the coordinate-scaling infinite-element
domain in commercial FE) -- standard Nedelec / RT elements on the mapped region, de Rham inherited for
free; the decay-BASIS exact-sequence IE below is the academic alternative.

Pure sympy.  Mathematica twin (cross-check via Mathematica's built-in ORTHONORMAL spherical
Grad/Curl/Div, plus a general-m angular-eigenvalue check):
packages/radia-mcp/src/radia_mcp/mathematica/basis_functions/infinite_element_derham.wls
"""
import os
import json
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import sympy as sp

N_FAIL = 0


def check(name, cond, detail=""):
    global N_FAIL
    if not cond:
        N_FAIL += 1
    print(f"  [{'ok ' if cond else 'FAIL'}] {name}{('  -- ' + detail) if detail else ''}")


r, th = sp.symbols("r theta", positive=True)


# ---- spherical vector calculus, orthonormal components (r^, th^, ph^), axisymmetric m=0 (d/dphi = 0) ----
def grad(u):
    """scalar -> vector (1-form)"""
    return (sp.diff(u, r), sp.diff(u, th) / r, sp.Integer(0))


def curl(A):
    """vector -> vector"""
    Ar, Ath, Aph = A
    cr = sp.diff(Aph * sp.sin(th), th) / (r * sp.sin(th))
    cth = -sp.diff(r * Aph, r) / r
    cph = (sp.diff(r * Ath, r) - sp.diff(Ar, th)) / r
    return (cr, cth, cph)


def div(A):
    """vector -> scalar (3-form density)"""
    Ar, Ath, Aph = A
    return sp.diff(r ** 2 * Ar, r) / r ** 2 + sp.diff(Ath * sp.sin(th), th) / (r * sp.sin(th))


def is_zero(e):
    e = sp.simplify(e)
    if e == 0:
        return True
    return sp.simplify(sp.trigsimp(sp.expand(e))) == 0


def vzero(V):
    return all(is_zero(c) for c in V)


def vsub(U, V):
    return tuple(U[i] - V[i] for i in range(3))


def smul(a, V):
    return tuple(a * c for c in V)


def vadd(*Vs):
    return tuple(sum(V[i] for V in Vs) for i in range(3))


# ---- the de Rham building blocks for mode n ----
def Yfun(n):
    return sp.legendre(n, sp.cos(th))


def P_field(k, n):   # radial vector harmonic, decay r^-k
    return (r ** (-k) * Yfun(n), sp.Integer(0), sp.Integer(0))


def B_field(k, n):   # poloidal-tangential vector harmonic
    return (sp.Integer(0), r ** (-k) * sp.diff(Yfun(n), th), sp.Integer(0))


def C_field(k, n):   # toroidal vector harmonic
    return (sp.Integer(0), sp.Integer(0), r ** (-k) * sp.diff(Yfun(n), th))


def scalar0(k, n):   # 0-form / 3-form scalar, decay r^-k
    return r ** (-k) * Yfun(n)


print("=" * 90)
print(" act7_26_derham_infinite_element : exact-sequence (de Rham) infinite element, sympy-verified")
print("=" * 90)
print("  exterior r>=a, mode n, m=0 representative; radial families shift +1 per form degree")

P_ORDER = 3
RESULTS = {
    "n_values": [1, 2, 3],
    "P": P_ORDER,
    "radial_families": {},
    "structure_constants": {
        "grad(r^-k Y)": "-k * P(k+1) + 1 * B(k+1)",
        "curl P(k)": "-1 * C(k+1)",
        "curl B(k)": "(1-k) * C(k+1)",
        "curl C(k)": "-n(n+1) * P(k+1) + (k-1) * B(k+1)   [Legendre eig]",
        "div P(k)": "(2-k) * r^-(k+1) Y",
        "div B(k)": "-n(n+1) * r^-(k+1) Y   [Legendre eig]",
        "div C(k)": "0",
    },
}

for n in (1, 2, 3):
    print(f"\n[n={n}]  Y = P_{n}(cos th)")

    # 0) Legendre angular-Laplace eigenvalue -- the identity that closes the complex
    angL = sp.diff(sp.sin(th) * sp.diff(Yfun(n), th), th) / sp.sin(th) + n * (n + 1) * Yfun(n)
    check(f"n={n}: Legendre angular-Laplace eig  (1/sin)d(sin dY) = -n(n+1) Y", is_zero(angL))

    # 1) operator sanity (the complex property, exact)
    check(f"n={n}: curl(grad u) == 0           (complex: d o d = 0)", vzero(curl(grad(scalar0(n + 1, n)))))
    check(f"n={n}: div(curl A) == 0            (complex: d o d = 0)", is_zero(div(curl(C_field(n + 1, n)))))

    S0 = list(range(n + 1, n + 1 + P_ORDER))
    RESULTS["radial_families"][f"n={n}"] = {
        "S0_0form": S0, "S1_Hcurl": [k + 1 for k in S0],
        "S2_Hdiv": [k + 2 for k in S0], "S3_L2": [k + 3 for k in S0],
    }

    okg = okc = okd = True
    for k in S0:
        # grad : V0 -> V1   (radial power k -> k+1)
        g = grad(scalar0(k, n))
        g_exp = vadd(smul(-k, P_field(k + 1, n)), B_field(k + 1, n))
        okg = okg and vzero(vsub(g, g_exp))

        # curl : V1 -> V2
        cP = vsub(curl(P_field(k, n)), smul(-1, C_field(k + 1, n)))
        cB = vsub(curl(B_field(k, n)), smul(1 - k, C_field(k + 1, n)))
        cC = vsub(curl(C_field(k, n)),
                  vadd(smul(-n * (n + 1), P_field(k + 1, n)), smul(k - 1, B_field(k + 1, n))))
        okc = okc and vzero(cP) and vzero(cB) and vzero(cC)

        # div : V2 -> V3
        dP = div(P_field(k, n)) - (2 - k) * scalar0(k + 1, n)
        dB = div(B_field(k, n)) - (-n * (n + 1)) * scalar0(k + 1, n)
        dC = div(C_field(k, n))
        okd = okd and is_zero(dP) and is_zero(dB) and is_zero(dC)

    check(f"n={n}: grad(V0) subset V1 EXACTLY  (-k*P + B; powers S0->S1)", okg,
          f"S0={S0} -> S1={[k+1 for k in S0]}")
    check(f"n={n}: curl(V1) subset V2 EXACTLY  (toroidal closes via Legendre -n(n+1))", okc,
          f"S1={[k+1 for k in S0]} -> S2={[k+2 for k in S0]}")
    check(f"n={n}: div(V2)  subset V3 EXACTLY  (Legendre -n(n+1) in div B)", okd,
          f"S2={[k+2 for k in S0]} -> S3={[k+3 for k in S0]}")

print("\n[summary] the de Rham (exact-sequence) infinite element:")
print("    H1 --grad--> H(curl) --curl--> H(div) --div--> L2  on the exterior r>=a,")
print("    with radial decay families  S0={n+1..n+P}, S1=S0+1, S2=S0+2, S3=S0+3.")
print("    grad/curl/div map each space EXACTLY into the next (commuting diagram), the")
print("    toroidal/div closure relying on the Legendre eigenvalue -n(n+1).")
print("    => the de Rham IE EXISTS; its radial families are the +1-shift towers above.")
print("    This is Demkowicz-Pal (CMAME 1998); known construction, implemented+verified here.")
print("    (The 0-form tower is act7_25's scalar IE -> static DtN -(n+1)/a.)")

RESULTS["n_fail"] = N_FAIL
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "act7_26_derham_infinite_element.json")
with open(out, "w") as f:
    json.dump(RESULTS, f, indent=2)
print(f"\n  wrote {os.path.basename(out)}")

print("\n" + "=" * 90)
print(" ALL CHECKS PASSED" if N_FAIL == 0 else f" {N_FAIL} CHECK(S) FAILED")
print("=" * 90)
assert N_FAIL == 0, f"{N_FAIL} checks failed"
