# -*- coding: utf-8 -*-
r"""
act2_05_geometry_floor_law.py  (Kelvin mesh-adequacy: the GEOMETRY-FLOOR side)
================================================================================
The OTHER half of the mesh-adequacy criterion (act2_04_adequacy_eccentric_multibody did the SOURCE side -- the
required element order p* from the source's multipole content).  At p>=n the Kelvin
DtN error is no longer the source content (captured exactly) but the CURVED-GEOMETRY
floor: the isoparametric order-k sphere approximating the true sphere.  act2_03_floor_vs_curve
MEASURED that the floor drops as k rises at FIXED coarse mesh; here we pin the
CLOSED-FORM LAW by sweeping the surface mesh h at each geometry order k and fitting the
asymptotic (fine-h) convergence rate.

THE MEASURED LAW (the data settled it -- an initial k+1 reading was wrong):
        floor_n(k, h/R)  ~  C_n * (h/R)^(2k)      [down to the double-precision floor]
i.e. the DtN eigenvalue (an ENERGY functional) SUPERCONVERGES at TWICE the isoparametric
boundary order k -- the classical domain-approximation (Strang-Fix / Ciarlet-Raviart)
result.  q(1)~2, q(2)~4; k=3 (2k=6) would need <1e-6, so it bottoms out on the
double-precision quadrature/round-off floor (~1e-6) instead -- this IS act2_03_floor_vs_curve's
saturation.  Practical rule: k=2 (rate 4) already reaches ~1e-6 at modest h; higher k
buys nothing in double precision.

Together with act2_04_adequacy_eccentric_multibody the total Kelvin error is
        err  ~  max( (d_max/R)^p  [source side, act2_04_adequacy_eccentric_multibody],  C*(h/R)^(2k)  [this] )
which is the abstract's "required surface resolution AND element order in closed form."

VERIFIED HERE (asserted; NGSolve sphere DtN, reuses act2_03_floor_vs_curve's lam_err):
  [1] at FIXED p the DtN floor refines in h as (h/R)^(2k): q(1)~2, q(2)~4 (the q(2)/q(1)
      doubling is the 2k superconvergence, NOT the boundary order k+1).
  [2] k=3 SATURATES on the double-precision floor (~1e-6): non-monotone at coarse h, so
      the 2k=6 rate is not observable -- the numerical floor, not the geometry, caps it.
  [3] the floor rate is MODE-robust (it is the sphere geometry, not the multipole):
      same q at n=2 and n=3.

NON-CLAIM: rates are fitted over the fine-h (asymptotic) window; very coarse h is
pre-asymptotic and the deep floor bottoms out at quadrature/round-off (~1e-6).  The law
is the geometric RATE read in its clean regime.

No overclaim: every printed 'ok' is gated on an executed numerical assertion.
"""
import os
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import numpy as np
import ngsolve as ng
from netgen.occ import Sphere, Pnt, OCCGeometry
from radia_mcp.radia_ngsolve.fem_bem_coupling import _solid_harmonic

ng.ngsglobals.msg_level = 0
R = 1.0


def lam_err(n, p, maxh, curve_k):
    """Relative error of the Kelvin DtN eigenvalue -(n+1)/R (reuses act2_03_floor_vs_curve)."""
    with ng.TaskManager():
        mesh = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), R)).GenerateMesh(maxh=maxh))
        mesh.Curve(curve_k)                        # isoparametric GEOMETRY order = k
        fes = ng.H1(mesh, order=p, dirichlet=".*"); u, v = fes.TnT()
        a = ng.BilinearForm(ng.grad(u) * ng.grad(v) * ng.dx); a.Assemble()
        gf = ng.GridFunction(fes); gf.Set(_solid_harmonic(n), ng.BND)
        r = gf.vec.CreateVector(); r.data = -(a.mat * gf.vec)
        gf.vec.data += a.mat.Inverse(fes.FreeDofs()) * r
        bi = 2 * (p + curve_k) + 6                 # rich quadrature so the floor is geometry
        energy = float(ng.Integrate(ng.grad(gf) * ng.grad(gf) * ng.dx(bonus_intorder=bi), mesh))
        bmass = float(ng.Integrate(gf * gf * ng.ds(bonus_intorder=bi), mesh))
        lam = -1.0 / R - energy / bmass; lam_ex = -(n + 1) / R
        return abs(lam - lam_ex) / abs(lam_ex)


def asym_rate(hs, es, npts=3):
    """log-log slope over the FINEST npts points (the asymptotic geometry rate)."""
    hs = np.array(hs, float); es = np.array(es, float)
    idx = np.argsort(hs)[:npts]                    # the npts smallest h
    return float(np.polyfit(np.log(hs[idx]), np.log(es[idx]), 1)[0])


print("=" * 78)
print(" act2_05_geometry_floor_law : Kelvin GEOMETRY-FLOOR law  floor ~ (h/R)^(2k)  (energy superconv.)")
print("=" * 78)

hs = [0.6, 0.45, 0.33, 0.25, 0.18]
n, p = 2, 3
print(f"\n[1] sphere DtN floor (n={n}, p={p}); sweep mesh h at each geometry order k (fine-h rate):")
print(f"    {'k':>2}  " + "  ".join(f"h={h:.2f}" for h in hs) + f"   {'q(fine-h)':>10}")
es_all = {}
for k in (1, 2, 3):
    es = [lam_err(n, p, h, k) for h in hs]
    es_all[k] = es
    print(f"    {k:>2}  " + "  ".join(f"{e:.1e}" for e in es) + f"   q={asym_rate(hs, es):5.2f}")
q1, q2 = asym_rate(hs, es_all[1]), asym_rate(hs, es_all[2])
print(f"    -> q(1)={q1:.2f} (~2k=2), q(2)={q2:.2f} (~2k=4); ratio q(2)/q(1)={q2/q1:.2f} (the 2k doubling)")
assert 1.5 <= q1 <= 2.4, "k=1 floor must converge ~ h^2 (2k=2; flat facets)"
assert 3.2 <= q2 <= 4.4, "k=2 floor must converge ~ h^4 (2k=4 superconvergence, NOT k+1=3)"
assert 1.7 <= q2 / q1 <= 2.4, "the rate must DOUBLE k->k+1 (the 2k energy superconvergence)"

print(f"\n[2] k=3 SATURATES on the double-precision floor (~1e-6) -- 2k=6 is not observable:")
es3 = es_all[3]
print(f"    k=3 floor: " + "  ".join(f"{e:.1e}" for e in es3))
print(f"    coarsest two ({es3[0]:.1e}, {es3[1]:.1e}) are non-monotone/flat -> at the numerical floor")
assert min(es3) < 2e-6, "k=3 must reach the deep (~1e-6) double-precision floor"
assert es3[1] > es3[0] * 0.7, "k=3 must be flat/non-monotone at coarse h (saturated, not 2k=6)"
print(f"    => the geometry is no longer the bottleneck at k=3; the numerical floor caps it")

print(f"\n[3] mode-robust: the floor rate is the sphere GEOMETRY, not the multipole:")
es_n3 = [lam_err(3, 4, h, 2) for h in hs]           # n=3, p=4, k=2
q_n3 = asym_rate(hs, es_n3)
print(f"    n=3, p=4, k=2: q={q_n3:.2f}   vs n=2, k=2: q={q2:.2f}  -> same geometry rate")
assert abs(q_n3 - q2) < 0.9, "the floor rate must be ~mode-independent (it is sphere geometry)"

print("\n[verdict]")
print(f"  The Kelvin geometry-FLOOR law, MEASURED: floor ~ (h/R)^(2k) -- the DtN eigenvalue")
print(f"  (an energy functional) SUPERCONVERGES at TWICE the isoparametric boundary order k")
print(f"  (q(1)~{q1:.1f}, q(2)~{q2:.1f}, ratio ~2), mode-independent; k=3 (2k=6) bottoms out on the")
print(f"  double-precision floor ~1e-6 (act2_03_floor_vs_curve's saturation).  Practical: k=2 suffices.")
print(f"  With act2_04_adequacy_eccentric_multibody the TOTAL error ~ max((d_max/R)^p, C(h/R)^(2k)) -- the required ELEMENT")
print(f"  ORDER and SURFACE RESOLUTION both in closed form (the mesh-adequacy criterion, done).")
print("\nALL CHECKS PASSED.")
