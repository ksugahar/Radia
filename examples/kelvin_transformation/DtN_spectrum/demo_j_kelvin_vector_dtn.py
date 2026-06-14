# -*- coding: utf-8 -*-
# DEMO (j) (verified): the Kelvin closure works VOLUMETRICALLY for the VECTOR-potential
# (edge-element / HCurl) formulation too -- an independent vector exterior DtN, not the
# A->B.n->scalar-DtN detour of demo3. The continuation of the alpha singular-material
# thread, via transformation optics of the radial inversion x = R^2 x'/|x'|^2 (Jacobian J):
#   scalar Laplace   : weight    W  = J^{-1} J^{-T} |det J| = (R/rho')^2 I
#   vector curl-curl : reluct.   nu' = J^T J / |det J|      = (rho'/R)^2 I
# The HCurl unknown is the 1-form pullback A' = J^T A; on Gamma (rho'=R) the tangential
# trace A'_t = A_t, so the datum is the physical multipole's A_t. The curl-curl gradient
# null space is fixed by a tiny mass regularization (gauge). Effective DtN eigenvalue =
# (energy) / (boundary tangential-trace mass).
#
# Verified:
#   scalar : lambda = (n+1)/R exactly (n=1->2, 2->3, 3->4) -- the Laplace DtN ladder.
#   vector : dipole lambda_vec -> 1/R. Analytic check: E=int|B|^2=1/(6 pi R^3),
#            bm=oint|A_t|^2=1/(6 pi R^2) -> lambda_vec = 1/R. The vector magnetic Steklov
#            ladder is n/R (a DIFFERENT trace pairing than the scalar (n+1)/R -- both are
#            the exterior DtN, both share the p>=n / coarse-mesh / controllability story).
import os
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import numpy as np
import ngsolve as ng
from netgen.occ import Sphere, Pnt, OCCGeometry
from radia_mcp.radia_ngsolve.fem_bem_coupling import _solid_harmonic

R = 1.0
x, y, z = ng.x, ng.y, ng.z
rho2 = x * x + y * y + z * z
rho = ng.sqrt(rho2)


def mk(maxh, order):
    return ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), R)).GenerateMesh(maxh=maxh)).Curve(min(order + 1, 3))


print("scalar Kelvin material-modulation:  lambda = int (R/rho')^2 |grad v|^2 / oint_G v^2")
w = R * R / rho2
for n in (1, 2, 3):
    mesh = mk(0.20, 3)
    fes = ng.H1(mesh, order=3, dirichlet=".*")
    u, v = fes.TnT()
    a = ng.BilinearForm(w * ng.grad(u) * ng.grad(v) * ng.dx(bonus_intorder=14)); a.Assemble()
    gf = ng.GridFunction(fes); gf.Set(_solid_harmonic(n), ng.BND)
    r = gf.vec.CreateVector(); r.data = -(a.mat * gf.vec); gf.vec.data += a.mat.Inverse(fes.FreeDofs()) * r
    E = float(ng.Integrate(w * ng.grad(gf) * ng.grad(gf) * ng.dx(bonus_intorder=14), mesh))
    bm = float(ng.Integrate(gf * gf * ng.ds(bonus_intorder=14), mesh))
    print("  n=%d: lambda=%.5f  exact (n+1)/R=%.1f  rel=%.1e" % (n, E / bm, (n + 1) / R, abs(E / bm - (n + 1) / R) / ((n + 1) / R)))

print("\nHCurl vector Kelvin (nu'=(rho'/R)^2), dipole A=m x r/(4 pi r^3), m=z:")
nup = rho2 / (R * R)
Adip = ng.CoefficientFunction((-y, x, 0.0)) / (4.0 * np.pi * rho2 * rho)
for order in (2, 3):
    mesh = mk(0.20, order)
    fes = ng.HCurl(mesh, order=order, dirichlet=".*")
    u, v = fes.TnT()
    a = ng.BilinearForm(nup * ng.curl(u) * ng.curl(v) * ng.dx(bonus_intorder=12)
                        + 1e-6 * u * v * ng.dx(bonus_intorder=8)); a.Assemble()
    gf = ng.GridFunction(fes); gf.Set(Adip, ng.BND)
    r = gf.vec.CreateVector(); r.data = -(a.mat * gf.vec); gf.vec.data += a.mat.Inverse(fes.FreeDofs()) * r
    E = float(ng.Integrate(nup * ng.curl(gf) * ng.curl(gf) * ng.dx(bonus_intorder=12), mesh))
    bm = float(ng.Integrate(gf.Trace() * gf.Trace() * ng.ds(bonus_intorder=12), mesh))
    print("  order=%d: lambda_vec=%.5f  analytic 1/R=%.1f  rel=%.1e" % (order, E / bm, 1.0 / R, abs(E / bm - 1.0 / R)))
print("\n-> the HCurl Kelvin ball recovers the vector magnetic Steklov eigenvalue n/R (dipole 1/R);"
      "\n   the singular material is integrable (demo_h) and the gradient null space is gauged.")
