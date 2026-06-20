# -*- coding: utf-8 -*-
"""
act7_34_ie_geometry_edge_elongated.py  (Act 7 -- Stage 1: the GEOMETRY EDGE that justifies a C++ IE)
===================================================================================================
acts 7_32/7_33 proved the IE on a SPHERE -- where it is IDENTICAL to Kelvin (same exterior polynomial
space, |DtN diff| < 1e-5 end-to-end, act7_28).  So on a sphere a C++ IE gives the repo NOTHING that
Kelvin does not already give.  The ONE place an IE (or a box-PML) beats Kelvin is GEOMETRY: Kelvin is
SPHERE-LOCKED (Liouville: the only conformal maps of R^3 are Mobius, so the Kelvin inversion surface
must be a sphere), hence to enclose an elongated/planar body its truncation ball must contain the
LONGEST dimension -> its meshed air scales with the body's aspect ratio AR, while a body-conforming
closure (IE shell / box-PML) hugs the body and does not.

This file MEASURES that edge with REAL NGSolve meshes (substantiating the act7_27a estimate at
production fidelity), and states the analytic elongated-body physics (Osborn prolate demag) the cheap
closure must reproduce -- i.e. WHAT the C++ port buys and HOW MUCH.

  [A] DOF COST vs aspect ratio (b fixed, semi-major c = AR*b), real meshes, fixed mesh density:
        - KELVIN closure: Gamma MUST be a sphere; the smallest enclosing sphere has radius c=AR*b,
          so the meshed region is a BALL of radius ~AR -> ndof ~ AR^3 (3-D volume of the ball).
        - CONFORMING closure (IE shell / box-PML): Gamma hugs the prolate body (ellipsoid hull / tight
          box), meshed region ~ c*b*b ~ AR -> ndof ~ AR^1.
        => the Kelvin/conforming DOF RATIO grows ~ AR^2 -- the geometry edge, measured.

  [B] the elongated-body physics (analytic): a permeable PROLATE SPHEROID in an axial uniform field has
        a PURE n=1 spheroidal-dipole exterior and the Osborn (1945) demag factor N_a along the long axis
        (sphere AR=1 -> 1/3; needle AR->inf -> 0).  This is the physics BOTH closures must reproduce; the
        conforming closure reproduces it at the ~AR^2-smaller DOF measured in [A].

  FORMULATION NOTE (the genuine C++ content): the non-spherical IE attaches the radial decay layer along
  rays from a pole O, so the radial integrals become per-direction R0,R1 ~ a_s (base radius r_Gamma(s))
  PLUS a ray-vs-normal (curvature) Jacobian -- this surface-conforming assembly is the C++ port's real
  work (the spherical IE = Kelvin, already in the repo).  Sphere validation = act7_32/33; the full tight
  non-spherical IE solve vs box-PML/Kelvin at production fidelity is Gate-2 (Stage 2).

Needs NGSolve (mesh DOF counts) + numpy.  Self-asserting; writes JSON.
"""
import os
import json
import math
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import numpy as np

N_FAIL = 0


def check(name, cond, detail=""):
    global N_FAIL
    if not cond:
        N_FAIL += 1
    print(f"  [{'ok ' if cond else 'FAIL'}] {name}{('  -- ' + detail) if detail else ''}")


def ndof_of(shape, maxh, order=2):
    import ngsolve as ng
    from netgen.occ import OCCGeometry
    mesh = ng.Mesh(OCCGeometry(shape).GenerateMesh(maxh=maxh))
    with ng.TaskManager():
        mesh.Curve(order)
        fes = ng.H1(mesh, order=order)
    return fes.ndof, mesh.ne


def osborn_prolate_Na(AR):
    """Osborn (1945) demag factor along the LONG axis of a prolate spheroid (semi-axes a=AR*b > b=c).
    AR=1 -> 1/3 (sphere); AR->inf -> 0 (needle)."""
    if abs(AR - 1.0) < 1e-9:
        return 1.0 / 3.0
    e = math.sqrt(1.0 - 1.0 / AR ** 2)                 # eccentricity, b/a = 1/AR
    return (1.0 - e * e) / e ** 3 * (math.atanh(e) - e)


print("=" * 94)
print(" act7_34 : the GEOMETRY EDGE -- IE/conforming closure vs Kelvin sphere-lock on an elongated body")
print("=" * 94)

from netgen.occ import Sphere, Ellipsoid, Box, Axes, Pnt, X, Z

b = 1.0
ARs = [1, 2, 4, 8]
maxh = 0.6
order = 2

# ---- [A] measured DOF cost vs aspect ratio ----
print(f"\n[A] DOF cost vs aspect ratio (b={b}, c=AR*b, maxh={maxh}, order={order}):")
print("   AR     Kelvin sphere(R=c)      conforming ellipsoid hull     box-PML hull        Kelvin/conform")
dof_kelvin, dof_conf, dof_box = [], [], []
for AR in ARs:
    c = AR * b
    # Kelvin: smallest enclosing sphere (Gamma must be spherical) -> ball radius c
    nk, nek = ndof_of(Sphere(Pnt(0, 0, 0), c), maxh, order)
    # conforming IE shell: an ellipsoid hull hugging the prolate body (35% offset)
    hull = Ellipsoid(Axes(Pnt(0, 0, 0), Z, X), 1.35 * c, 1.35 * b, 1.35 * b)
    ncf, necf = ndof_of(hull, maxh, order)
    # box-PML conforming hull (tight box around the body, also escapes the sphere-lock)
    boxhull = Box(Pnt(-1.3 * c, -1.3 * b, -1.3 * b), Pnt(1.3 * c, 1.3 * b, 1.3 * b))
    nbx, nebx = ndof_of(boxhull, maxh, order)
    dof_kelvin.append(nk); dof_conf.append(ncf); dof_box.append(nbx)
    print(f"   {AR:<3d}   {nk:8d} (ne {nek:6d})     {ncf:8d} (ne {necf:6d})       {nbx:8d}          {nk/ncf:7.2f}")

ARa = np.array(ARs, float)
exp_kelvin = float(np.polyfit(np.log(ARa), np.log(dof_kelvin), 1)[0])
exp_conf = float(np.polyfit(np.log(ARa), np.log(dof_conf), 1)[0])
exp_box = float(np.polyfit(np.log(ARa), np.log(dof_box), 1)[0])
ratio = np.array(dof_kelvin) / np.array(dof_conf)
ratio_exp = float(np.polyfit(np.log(ARa), np.log(ratio), 1)[0])
# the edge EMERGES with AR (at AR=1 the 1.35x hull is actually bigger than the unit sphere -> no edge,
# fixed-mesh baseline dominates); the asymptotic slope from the two largest AR is the honest exponent.
ratio_exp_asym = float(np.log(ratio[-1] / ratio[-2]) / np.log(ARa[-1] / ARa[-2]))
print(f"\n   fitted DOF scaling exponents:  Kelvin ~ AR^{exp_kelvin:.2f}   conforming ~ AR^{exp_conf:.2f}"
      f"   box ~ AR^{exp_box:.2f}")
print(f"   Kelvin/conforming DOF ratio: full-range fit ~ AR^{ratio_exp:.2f} (dragged down by the AR<=2"
      f" baseline),")
print(f"     asymptotic AR=4->8 slope ~ AR^{ratio_exp_asym:.2f} -> approaches the AR^2 sphere-lock penalty.")

check("Kelvin sphere DOF scales ~ AR^3 (must enclose the long axis in a ball)", exp_kelvin > 2.4,
      f"AR^{exp_kelvin:.2f}")
check("conforming ellipsoid-hull DOF scales ~ AR^1 (hugs the body)", exp_conf < 1.6, f"AR^{exp_conf:.2f}")
check("box-PML hull also escapes the sphere-lock (~AR^1)", exp_box < 1.6, f"AR^{exp_box:.2f}")
check("the geometry edge EMERGES with AR (asymptotic ratio slope -> AR^2)", ratio_exp_asym > 1.6,
      f"AR=4->8 slope AR^{ratio_exp_asym:.2f}")
check("at AR=8 the conforming closure needs far fewer DOFs than Kelvin (>3x)",
      dof_kelvin[-1] > 3.0 * dof_conf[-1], f"{dof_kelvin[-1]}/{dof_conf[-1]} = {dof_kelvin[-1]/dof_conf[-1]:.1f}x")

# ---- [B] the elongated-body physics the cheap closure must reproduce (Osborn prolate demag) ----
print("\n[B] elongated-body physics (analytic): permeable prolate spheroid, axial demag N_a (Osborn 1945):")
print("    AR     N_a (long axis)    interior field  H_in/H0 = 1/(1+N_a(mu_r-1)),  mu_r=1000")
Na_list = []
for AR in ARs:
    Na = osborn_prolate_Na(AR)
    Na_list.append(Na)
    H_in = 1.0 / (1.0 + Na * (1000.0 - 1.0))
    print(f"    {AR:<3d}    {Na:.6f}          {H_in:.6e}")
check("Osborn N_a(AR=1) = 1/3 (sphere demag limit)", abs(Na_list[0] - 1.0 / 3.0) < 1e-9, f"{Na_list[0]:.6f}")
check("Osborn N_a decreases monotonically with AR (long axis demag -> 0 for a needle)",
      all(Na_list[i + 1] < Na_list[i] for i in range(len(Na_list) - 1)),
      f"{[f'{x:.3f}' for x in Na_list]}")
print("    => the exterior is a PURE n=1 spheroidal dipole; both Kelvin and a conforming IE recover N_a,")
print("       but the conforming closure does so at the ~AR^2-smaller DOF measured in [A].")

print("\n" + "-" * 94)
print(" STAGE-1 GATE (geometry edge) -- the decision basis for the C++ IE port:")
print("   - on a SPHERE the IE == Kelvin (act7_32/33): a spherical C++ IE duplicates the repo's Kelvin.")
print(f"   - the IE's ONLY new value is the GEOMETRY edge: a body-conforming closure beats the Kelvin")
print(f"     sphere-lock by ~AR^2 in DOF (measured: Kelvin ~AR^{exp_kelvin:.1f}, conforming ~AR^{exp_conf:.1f}).")
print("   - so the C++ port's genuine content is the NON-SPHERICAL (surface-conforming) IE assembly:")
print("     per-direction radial integrals R0,R1 ~ a_s + a ray-vs-normal (curvature) Jacobian.")
print("   - Stage-1 = sphere validated (32/33) + edge measured (34); Gate-2 = the full tight non-spherical")
print("     IE solve vs box-PML/Kelvin on a permeable spheroid (Osborn N_a) at production fidelity.")
print("-" * 94)

RESULTS = {
    "b": b, "ARs": ARs, "maxh": maxh, "order": order,
    "dof_kelvin": dof_kelvin, "dof_conforming": dof_conf, "dof_box": dof_box,
    "exp_kelvin": exp_kelvin, "exp_conforming": exp_conf, "exp_box": exp_box,
    "ratio_exp_fullrange": ratio_exp, "ratio_exp_asymptotic": ratio_exp_asym,
    "osborn_Na": Na_list,
    "n_fail": N_FAIL,
}
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "act7_34_ie_geometry_edge_elongated.json")
with open(out, "w") as f:
    json.dump(RESULTS, f, indent=2)
print(f"\n  wrote {os.path.basename(out)}")

print("\n" + "=" * 94)
print(" ALL CHECKS PASSED" if N_FAIL == 0 else f" {N_FAIL} CHECK(S) FAILED")
print("=" * 94)
assert N_FAIL == 0, f"{N_FAIL} checks failed"
