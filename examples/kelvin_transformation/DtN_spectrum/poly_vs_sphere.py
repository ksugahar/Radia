# -*- coding: utf-8 -*-
"""Does a FLAT polyhedron (Curve 1) truncation work as well as the curved sphere?
Measure the effective DtN eigenvalue error for the dipole (n=1, linear image) and
the quadrupole (n=2, quadratic image), flat (Curve 1) vs curved geometry, vs mesh.
Hypothesis: the dipole is geometry-robust (faceting only perturbs it at O(h^2),
benign uniform shift); a higher mode is more geometry-sensitive."""
import ngsolve as ng
from ngsolve import (BilinearForm, GridFunction, Integrate, grad, dx, ds, BND)
from netgen.occ import OCCGeometry, Sphere, Pnt
from radia_mcp.radia_ngsolve.fem_bem_coupling import _solid_harmonic


def dtn_eig(degree, maxh, order, curve_order, intorder=12, R=1.0):
    mesh = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), R)).GenerateMesh(maxh=maxh)).Curve(curve_order)
    fes = ng.H1(mesh, order=order, dirichlet=".*")
    u, v = fes.TnT()
    a = BilinearForm(grad(u) * grad(v) * dx); a.Assemble()
    gfu = GridFunction(fes); gfu.Set(_solid_harmonic(degree), BND)
    rhs = gfu.vec.CreateVector(); rhs.data = -(a.mat * gfu.vec)
    gfu.vec.data += a.mat.Inverse(freedofs=fes.FreeDofs()) * rhs
    energy = float(Integrate(grad(gfu) * grad(gfu) * dx(bonus_intorder=intorder), mesh))
    bmass = float(Integrate(gfu * gfu * ds(bonus_intorder=intorder), mesh))
    lam = -1.0 / R - energy / bmass
    return fes.ndof, abs(lam - (-(degree + 1) / R)) / ((degree + 1) / R)


print("=== DIPOLE  n=1, order=1  (image = LINEAR) :  FLAT polyhedron vs curved ===")
print(f"{'maxh':>6} {'ndof(c1)':>9} {'err FLAT(c1)':>14} {'err CURVED(c2)':>15} {'flat/curved':>12}")
for h in (0.6, 0.4, 0.3, 0.2, 0.14):
    nd, ef = dtn_eig(1, h, 1, 1)
    _,  ec = dtn_eig(1, h, 1, 2)
    print(f"{h:>6.2f} {nd:>9d} {ef:>14.3e} {ec:>15.3e} {ef/ec:>12.1f}")

print("\n=== QUADRUPOLE  n=2, order=2  (image = QUADRATIC) :  FLAT vs curved ===")
print(f"{'maxh':>6} {'ndof':>9} {'err FLAT(c1)':>14} {'err CURVED(c3)':>15} {'flat/curved':>12}")
for h in (0.6, 0.4, 0.3, 0.2):
    nd, ef = dtn_eig(2, h, 2, 1)
    _,  ec = dtn_eig(2, h, 2, 3)
    print(f"{h:>6.2f} {nd:>9d} {ef:>14.3e} {ec:>15.3e} {ef/ec:>12.1f}")
print("\nDONE")
