# -*- coding: utf-8 -*-
"""
act7_25_infinite_element_dtn.py  (Act 7 -- the open-boundary comparison; the INFINITE-ELEMENT column)
=====================================================================================================
Puts the classic **infinite element** (Zienkiewicz-Bettess / mapped decay-basis) on the SAME
DtN-spectral yardstick as Kelvin / ballooning / Robin / PML / BEM.

WHY this file exists: NGSolve has **no** infinite element (its open-boundary feature is PML only --
verified: `ngsolve.pml` = Radial/Cartesian/BrickRadial/HalfSpace/Custom, zero "infinite" in the
package).  So to rank the infinite element on the yardstick we implement the radial infinite element
directly, in the same closed-form radial style as act7_22's other static closures.  This also
corrects a naming slip: act7_22's "ballooning" entry is actually the crudest member of this family --
a finite truncation WALL (Dirichlet at r=R) -- NOT a true mapped infinite element.

METHOD (radial static Laplace, spherical-harmonic mode n, truncation radius a):
  - the exterior decays as `r^{-(n+1)}`; a mapped / decay-basis infinite element spans it with the
    reciprocal-power shape functions  phi_k(r) = (a/r)^k,  k = 1..P  (P = the element "order");
  - the radial Laplace bilinear form  int_a^inf [ r^2 u' v' + n(n+1) u v ] dr  is CLOSED FORM on this
    basis (the integrals converge because the functions decay):
        A_kl = a * ( k*l + n(n+1) ) / ( (k+l) - 1 )
  - the exterior DtN is the energy-minimising harmonic extension with u(a)=1 (the Dirichlet datum):
        Lambda_n = - 1 / ( a * 1^T A^{-1} 1 )            (since min energy = 1/(1^T A^{-1} 1))
  - exact = -(n+1)/a  (3D static ladder).

KEY RESULT (asserted): the infinite element is **EXACT for every mode n <= P-1** (its decay basis
contains the exact `r^{-(n+1)}`), and degrades gracefully for n >= P -- i.e. it is a p-method in the
decay order P.  This is the **OPPOSITE failure mode** of the truncation wall (act7_22 "ballooning"),
which fails the LOW (slow-decaying) modes `~(a/R)^{2n+1}` and is good at high n.  So the infinite
element joins Kelvin / BEM in the CONVERGENT class (exact-operator, here a modal/decay-basis
realisation valid for the separable exterior), distinct from the finite-reach wall.

CLASSIFICATION (three honest distinctions):
  - TREFFTZ FAMILY.  The decay basis `(a/r)^k` CONTAINS the one Trefftz function `r^-(n+1)` (the exact
    exterior harmonic), so the Galerkin fit is EXACT for `n <= P-1`; a PURE Trefftz element (exact
    multipoles as the basis) would reproduce the exact DtN `-(n+1)/a` trivially per mode.  This is a
    Bettess-Galerkin decay-basis IE, i.e. Trefftz-FAMILY, not pure Trefftz.  (NGSolve core has NO
    infinite element AND no Trefftz space; `ngstrefftz` is a separate add-on = interior Trefftz-DG,
    not an open-boundary element.)
  - SPARSE vs APPROXIMATE.  This MODAL measurement is a tiny DENSE `P x P` block per mode (`P~3-6`).
    The PRACTICAL Bettess element (local decay-shape-function elements on the truncation surface)
    assembles SPARSELY -- a thin FE layer, each element coupling only its own nodes (its appeal vs
    dense BEM).  And it is APPROXIMATE: EXACT for `n <= P-1`, Galerkin best-fit for `n >= P`, converging
    as `P -> inf` (a spectral / p-method).  The exact-and-no-approximation alternative is BEM / pure
    Trefftz -- but those are DENSE.  No free lunch: sparse (Kelvin / IE) => approximate; exact (BEM /
    Trefftz multipole) => dense.
  - de Rham / VECTOR extension is act7_26_derham_infinite_element (the exact-sequence
    H1->H(curl)->H(div)->L2 IE, Demkowicz-Pal CMAME 1998); THIS scalar tower is its 0-form.  NOTE the
    practical SHIPPED de Rham open boundary is instead the coordinate-MAPPING family (the Kelvin
    transformation here; the coordinate-scaling infinite-element domain in commercial FE) -- standard
    Nedelec / RT elements on the mapped region, de Rham inherited for FREE; the decay-basis
    exact-sequence IE (act7_26) is the academic alternative.

Pure numpy.  (The eddy / radiation infinite element would need an `exp(-q r)` decay basis -- out of
scope here; this is the classic static infinite element, which is where the method is used.)
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


A_RAD = 1.0


def infinite_element_dtn(n, P, a=A_RAD):
    """Radial mapped/Bettess infinite element DtN for static mode n, decay order P."""
    k = np.arange(1, P + 1)
    A = a * (np.outer(k, k) + n * (n + 1)) / (k[:, None] + k[None, :] - 1.0)
    g = np.ones(P)
    return -1.0 / (a * (g @ np.linalg.solve(A, g)))


def wall_dtn(n, R, a=A_RAD):
    """Finite truncation wall (Dirichlet at r=R) = act7_22's 'ballooning' (the crudest member)."""
    C = (R / a) ** (2 * n + 1)
    return (n + (n + 1) * C) / (1 - C)


print("=" * 80)
print(" act7_25_infinite_element_dtn : the infinite element on the static DtN yardstick (a=1)")
print("=" * 80)
print("  exact static ladder  lambda_n = -(n+1)")

MODES = list(range(0, 7))
TABLE = {"a": A_RAD, "exact": [-(n + 1) for n in MODES], "infinite_element": {}, "wall": {}}

# ---- infinite element: exact for n <= P-1, degrades for n >= P ----
print("\n[infinite element] reldef d_n vs -(n+1), decay order P (exact for n <= P-1):")
print("   P     n=0       n=1       n=2       n=3       n=4       n=5       n=6")
for P in (2, 4, 6):
    d = [abs(infinite_element_dtn(n, P) - (-(n + 1))) / (n + 1) for n in MODES]
    TABLE["infinite_element"][f"P={P}"] = d
    print(f"  {P}   " + "  ".join(f"{x:.2e}" for x in d))
    check(f"infinite element P={P}: EXACT for every mode n <= P-1 (d_n < 1e-10)",
          max(d[n] for n in MODES if n <= P - 1) < 1e-10,
          f"max(n<=P-1) {max(d[n] for n in MODES if n <= P - 1):.1e}")

# convergence: raising P captures higher modes exactly
check("infinite element CONVERGES: a fixed high mode (n=5) becomes exact when P reaches n+1",
      abs(infinite_element_dtn(5, 4) - (-6)) / 6 > 1e-6 and abs(infinite_element_dtn(5, 6) - (-6)) / 6 < 1e-10,
      f"n=5: P=4 {abs(infinite_element_dtn(5,4)+6)/6:.1e} -> P=6 {abs(infinite_element_dtn(5,6)+6)/6:.1e}")

# ---- truncation wall: fails the LOW modes (the OPPOSITE) ----
print("\n[truncation wall] (act7_22 'ballooning', R=4) reldef d_n -- fails the LOW modes:")
R_wall = 4.0
dw = [abs(wall_dtn(n, R_wall) - (-(n + 1))) / (n + 1) for n in MODES]
TABLE["wall"]["R=4"] = dw
print("        " + "  ".join(f"{x:.2e}" for x in dw))

# ---- the headline: OPPOSITE failure modes ----
dP4 = TABLE["infinite_element"]["P=4"]
check("OPPOSITE failure: infinite element NAILS the low modes (P=4: d[n=0]<1e-10) while the wall FAILS them (d[n=0]>0.1)",
      dP4[0] < 1e-10 and dw[0] > 0.1, f"IE[0] {dP4[0]:.1e}, wall[0] {dw[0]:.2f}")
check("OPPOSITE failure: the wall is good at HIGH modes (d[n=6]<1e-4) while the infinite element degrades there (P=4: d[n=6]>1e-4)",
      dw[6] < 1e-4 and dP4[6] > 1e-4, f"wall[6] {dw[6]:.1e}, IE-P4[6] {dP4[6]:.1e}")

print("\n[summary] the infinite element is a p-method in the DECAY order P:")
print("    - EXACT for n <= P-1 (its (a/r)^k basis contains the exact r^-(n+1) decay);")
print("    - graceful degradation for n >= P; raising P captures more modes exactly.")
print("    => CONVERGENT class (exact-operator, modal/decay-basis realisation for a separable")
print("       exterior) -- alongside Kelvin (conformal-FE) and BEM.  OPPOSITE failure to the")
print("       finite truncation WALL (act7_22 'ballooning'), which fails the LOW modes.")
print("    NGSolve has NO infinite element (PML only); this is implemented directly here.")

out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "act7_25_infinite_element_dtn.json")
with open(out, "w") as f:
    json.dump(TABLE, f, indent=2)
print(f"\n  wrote {os.path.basename(out)}")

print("\n" + "=" * 80)
print(" ALL CHECKS PASSED" if N_FAIL == 0 else f" {N_FAIL} CHECK(S) FAILED")
print("=" * 80)
assert N_FAIL == 0, f"{N_FAIL} checks failed"
