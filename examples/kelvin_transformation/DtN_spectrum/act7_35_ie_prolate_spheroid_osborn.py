# -*- coding: utf-8 -*-
"""
act7_35_ie_prolate_spheroid_osborn.py  (Act 7 -- Gate-2 milestone 1: the elongated-body physics)
================================================================================================
acts 7_32/33 validated the NGSolve-native infinite element (IE) on a SPHERE (== Kelvin).  act7_34
MEASURED the geometry edge a body-conforming closure would buy (Kelvin sphere-lock ~AR^2.7 DOF vs a
conforming ~AR^1.3).  Gate-2 is the tight NON-SPHERICAL (surface-conforming) IE -- the only config
that beats Kelvin.  This file is its FIRST milestone: validate the elongated-body PHYSICS and set the
DOF BASELINE, using the (already-built, spherical) IE on an ENCLOSING sphere.

A permeable PROLATE SPHEROID (semi-axes a=AR*b > b, relative permeability mu_r) in a uniform axial
field has the Osborn (1945) demagnetising factor N_a along its long axis (sphere AR=1 -> 1/3; needle
-> 0); its interior field is uniform, H_in = H0 / (1 + N_a (mu_r-1)).  Here the spheroid (iron) sits
inside an enclosing air sphere of radius R; the reduced scalar potential phi_red (the mu-jump on the
spheroid surface is the source) is closed on the OUTER sphere by the IE.  This is a CORRECT but NOT
tight closure (the enclosing sphere is the expensive, Kelvin-class cost) -- it proves the physics and
gives the DOF the future tight non-spherical IE must beat (~AR^2 fewer, act7_34).

CHECKS (self-asserting):
  - the IE-closed solve reproduces Osborn N_a for AR = 1, 2, 4 and mu_r = 100, 1000
    (d(phi_red)/dz at the centre == H0 N_a (mu_r-1)/(1 + N_a (mu_r-1)), i.e. interior H_in == Osborn);
  - the enclosing-sphere DOF grows with AR (the baseline the tight non-spherical IE must beat).

Needs NGSolve + radia.infinite_element + numpy.  Writes JSON.
"""
import os
import json
import math
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import numpy as np

import ngsolve as ng
from netgen.occ import Sphere, Ellipsoid, Axes, Pnt, X, Z, OCCGeometry, Glue

import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../src"))
from radia import infinite_element as ie

N_FAIL = 0


def check(name, cond, detail=""):
    global N_FAIL
    if not cond:
        N_FAIL += 1
    print(f"  [{'ok ' if cond else 'FAIL'}] {name}{('  -- ' + detail) if detail else ''}")


def osborn_Na(AR):
    """Osborn (1945) demag factor along the LONG axis of a prolate spheroid (a=AR*b > b)."""
    if abs(AR - 1.0) < 1e-9:
        return 1.0 / 3.0
    e = math.sqrt(1.0 - 1.0 / AR ** 2)
    return (1.0 - e * e) / e ** 3 * (math.atanh(e) - e)


def solve_spheroid(AR, mu_r, b=1.0, P=6, order=3, R_frac=1.8, maxh_frac=0.4):
    """RSP solve: permeable prolate spheroid (iron) in an enclosing air sphere, IE-closed on the outer
    sphere.  Returns (d(phi_red)/dz at the centre, total DOF, n_elements)."""
    c = AR * b
    R = R_frac * c
    iron = Ellipsoid(Axes(Pnt(0, 0, 0), Z, X), c, b, b)
    iron.mat("iron"); iron.faces.name = "spheroid"
    outer = Sphere(Pnt(0, 0, 0), R); outer.faces.name = "outer"
    air = outer - iron; air.mat("air")
    mesh = ng.Mesh(OCCGeometry(Glue([iron, air])).GenerateMesh(maxh=maxh_frac * c))
    with ng.TaskManager():
        mesh.Curve(order)
        mu = mesh.MaterialCF({"iron": mu_r, "air": 1.0})
        Xc = ie.ie_compound_space(mesh, P, order=order, definedon=mesh.Boundaries("outer"))
        tr, te = Xc.TrialFunction(), Xc.TestFunction()
        n = ng.specialcf.normal(3)
        a = ng.BilinearForm(Xc, symmetric=True, check_unused=False)
        a += mu * ng.grad(tr[0]) * ng.grad(te[0]) * ng.dx                  # interior (iron + air)
        ie.add_exterior_ie(a, Xc, P, a=R, definedon=mesh.Boundaries("outer"))  # IE open boundary
        a.Assemble()
        f = ng.LinearForm(Xc)
        f += (mu_r - 1.0) * n[2] * te[0].Trace() * ng.ds(definedon=mesh.Boundaries("spheroid"))
        f.Assemble()
        gf = ng.GridFunction(Xc)
        gf.vec.data = a.mat.Inverse(Xc.FreeDofs(), inverse="sparsecholesky") * f.vec
    gz = ng.grad(gf.components[0])(mesh(0, 0, 0))[2]
    return float(gz), Xc.ndof, mesh.ne


print("=" * 96)
print(" act7_35 : permeable PROLATE SPHEROID, IE-closed -- Osborn (1945) demag (Gate-2 milestone 1)")
print("=" * 96)

P = 6
RESULTS = {"P": P, "cases": []}

print("\n  permeable prolate spheroid in an enclosing sphere (IE open boundary), vs Osborn N_a:")
print("   AR   mu_r     d(phi_red)/dz   Osborn target    rel.err    H_in(meas)  H_in(Osborn)   DOF")
for AR in (1.0, 2.0, 4.0):
    Na = osborn_Na(AR)
    for mu_r in (100.0, 1000.0):
        gz, ndof, ne = solve_spheroid(AR, mu_r, P=P)
        gz_target = Na * (mu_r - 1.0) / (1.0 + Na * (mu_r - 1.0))    # = H0 - H_in (H0=1)
        H_in_meas = 1.0 - gz
        H_in_osb = 1.0 / (1.0 + Na * (mu_r - 1.0))
        relerr = abs(gz - gz_target) / gz_target
        print(f"   {AR:<3.0f}  {mu_r:6.0f}    {gz:.6f}      {gz_target:.6f}     {relerr:.2e}   "
              f"{H_in_meas:.5f}    {H_in_osb:.5f}     {ndof}")
        check(f"AR={AR:.0f} mu_r={mu_r:.0f}: interior demag == Osborn N_a={Na:.4f}", relerr < 2e-2,
              f"relerr {relerr:.2e}")
        RESULTS["cases"].append(dict(AR=AR, mu_r=mu_r, Na=Na, gz=gz, gz_target=gz_target,
                                     relerr=relerr, H_in_meas=H_in_meas, H_in_osborn=H_in_osb,
                                     ndof=ndof, ne=ne))

# DOF baseline: the enclosing-sphere closure cost grows with AR (the tight non-spherical IE must beat it)
dofs = {AR: next(c["ndof"] for c in RESULTS["cases"] if c["AR"] == AR and c["mu_r"] == 100.0)
        for AR in (1.0, 2.0, 4.0)}
print(f"\n  enclosing-sphere DOF baseline: AR=1 {dofs[1.0]}, AR=2 {dofs[2.0]}, AR=4 {dofs[4.0]}"
      f"  (grows with AR -- the Kelvin-class cost the tight non-spherical IE must beat, act7_34)")
check("enclosing-sphere DOF grows with AR (the non-tight baseline)", dofs[4.0] > dofs[1.0],
      f"{dofs[1.0]} -> {dofs[4.0]}")

print("\n" + "-" * 96)
print(" GATE-2 MILESTONE 1 (physics + baseline):")
print("   - the NGSolve-native IE, closing an enclosing sphere around a permeable PROLATE SPHEROID,")
print("     reproduces the Osborn (1945) long-axis demag N_a for AR=1,2,4 -- the elongated-body physics;")
print("   - the closure is CORRECT but NOT tight (enclosing-sphere DOF = the Kelvin-class cost).")
print("   - milestone 2 (the hard part): the TIGHT non-spherical (surface-conforming) IE -- per-direction")
print("     radial R0,R1 ~ a_s + a ray-vs-normal Jacobian (prolate-spheroidal-coordinate or Astley mapped")
print("     IE, = Nannen 2013) -- reproducing this Osborn N_a at the ~AR^2-fewer DOF measured in act7_34.")
print("-" * 96)

RESULTS["dof_baseline"] = {str(int(k)): v for k, v in dofs.items()}
RESULTS["n_fail"] = N_FAIL
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "act7_35_ie_prolate_spheroid_osborn.json")
with open(out, "w") as f:
    json.dump(RESULTS, f, indent=2)
print(f"\n  wrote {os.path.basename(out)}")

print("\n" + "=" * 96)
print(" ALL CHECKS PASSED" if N_FAIL == 0 else f" {N_FAIL} CHECK(S) FAILED")
print("=" * 96)
assert N_FAIL == 0, f"{N_FAIL} checks failed"
