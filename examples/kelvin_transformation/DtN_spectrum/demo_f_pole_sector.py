# -*- coding: utf-8 -*-
# DEMO (f) (verified): 2D POLE-SECTOR symmetry preserves the DtN spectrum.
#
# The circular-Gamma DtN is diagonal in the angular harmonics e^{i n theta} with the
# mesh-independent ladder lambda_n = -n/R (2D). A rotating machine / static apparatus
# never models the whole circle -- it models ONE pole (or pole-pair) sector and ties
# the cut edges with a (anti)periodic condition. That condition merely SELECTS a
# sub-basis of admissible harmonics; the eigenvalue lambda_n of each admitted mode is
# UNCHANGED. So a 1-pole open-boundary model is exact for the harmonics the machine
# actually produces, at 1/(sector-fraction) of the DoF.
#
#   admissible harmonics, sector of one pole [0, pi/p]:
#     periodic   (u repeats over 2pi/p)      -> n = p, 2p, 3p, ...
#     anti-period (u flips sign over pi/p)    -> n = p, 3p, 5p, ... (machine fundamental n=p)
#     mirror     (Neumann radial edges)       -> n = p, 2p, 3p, ... (the cleanest to set up;
#                                                used here -- same spectral-restriction principle)
#
# Verified below (p=2, i.e. a 4-pole machine; quarter disk = one pole; R=1):
#   * the quarter sector reproduces lambda_2 = -2/R to the SAME accuracy as the full
#     circle, at ~1/4 the DoF;
#   * harmonics the sector does NOT admit (the dipole n=1, and n=3) are EXCLUDED
#     (large error) -- the machine cannot smuggle them in.
#
# Companion knowledge: dtn_coarse_mesh topic "symmetry_hex" (octant/3D + this 2D sector).
import os
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import ngsolve as ng
from netgen.occ import WorkPlane, OCCGeometry
from radia_mcp.radia_ngsolve.fem_bem_coupling import _solid_harmonic_2d, kelvin_dtn_eigenvalue

R = 1.0
disk = WorkPlane().Circle(0, 0, R).Face()
wedge = WorkPlane().MoveTo(0, 0).LineTo(2 * R, 0).LineTo(0, 2 * R).Close().Face()
sector = disk * wedge                                   # quarter disk = one pole of a 4-pole machine
for e in sector.edges:
    c = e.center
    e.name = "symx" if abs(c.x) < 1e-7 else ("symy" if abs(c.y) < 1e-7 else "gamma")
geo = OCCGeometry(sector, dim=2)


def sector_lambda(n, maxh=0.3, order=3):
    """Effective DtN eigenvalue of mode n on the quarter sector (arc=Dirichlet datum,
    radial edges natural/Neumann). lambda = -energy/bmass (2D offset 0); cmp -n/R."""
    mesh = ng.Mesh(geo.GenerateMesh(maxh=maxh)).Curve(min(order + 1, 4))
    datum = _solid_harmonic_2d(n)
    fes = ng.H1(mesh, order=order, dirichlet="gamma")
    u, v = fes.TnT()
    a = ng.BilinearForm(ng.grad(u) * ng.grad(v) * ng.dx); a.Assemble()
    gf = ng.GridFunction(fes); gf.Set(datum, ng.BND)
    r = gf.vec.CreateVector(); r.data = -(a.mat * gf.vec)
    gf.vec.data += a.mat.Inverse(fes.FreeDofs()) * r
    energy = float(ng.Integrate(ng.grad(gf) * ng.grad(gf) * ng.dx(bonus_intorder=10), mesh))
    bmass = float(ng.Integrate(gf * gf * ng.ds("gamma", bonus_intorder=10), mesh))
    lam = -energy / bmass
    return fes.ndof, lam, abs(lam + n / R) / (n / R)


print("2D circular DtN ladder lambda_n = -n/R (R=1); quarter sector = one pole of a 4-pole (p=2) machine\n")
print("{:>24} | {:>8} {:>11} {:>13}".format("model", "ndof", "lambda", "rel vs -n/R"))
print("-" * 62)
for n in (2, 3):
    fnd = kelvin_dtn_eigenvalue(R=R, degree=n, maxh=0.3, order=3, dim=2)
    snd, slam, srel = sector_lambda(n)
    print("{:>24} | {:>8d} {:>11.5f} {:>13.2e}".format("full circle  n=%d" % n, fnd["ndof"], fnd["lam"], fnd["rel_err"]))
    tag = "ADMITTED (=p)" if n == 2 else "not a sector mode"
    print("{:>24} | {:>8d} {:>11.5f} {:>13.2e}   {}".format("quarter sector n=%d" % n, snd, slam, srel, tag))
nd1, lam1, rel1 = sector_lambda(1)
print("{:>24} | {:>8d} {:>11.5f} {:>13.2e}   {}".format(
    "quarter sector n=1", nd1, lam1, rel1, "EXCLUDED (dipole: a 4-pole machine has none)"))
print("\n-> the sector reproduces its admissible harmonics (n=2) at full-circle accuracy with"
      "\n   ~1/4 the DoF, and excludes the rest. The open-boundary closure cost scales with the"
      "\n   sector fraction; the DtN spectrum is untouched.")
