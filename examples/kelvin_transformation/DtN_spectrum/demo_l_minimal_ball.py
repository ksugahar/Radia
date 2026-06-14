# -*- coding: utf-8 -*-
# DEMO (l) (verified): how coarse can the Kelvin ball be? The interior is irrelevant
# (the volume solve is Galerkin-exact at p>=n -- see kelvin_exterior_mesh.py), so the
# ball is meshed as a curvature-limited thin SHELL: its volume-element count is about the
# same as its Gamma surface-triangle count, NOT a refined volume. So the Kelvin closure's
# size tracks the Gamma surface + order p, like a surface method -- but the matrix stays
# SPARSE (cf. demo_k). Below: the minimal unit-sphere ball (curvature-limited; maxh in
# 0.5..1.5 all give the same mesh) is exact once p>=n.
import os
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import ngsolve as ng
from netgen.occ import Sphere, Pnt, OCCGeometry
from radia_mcp.radia_ngsolve.fem_bem_coupling import _solid_harmonic

R = 1.0


def kelvin_eig(n, maxh, order):
    mesh = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), R)).GenerateMesh(maxh=maxh)).Curve(min(order + 1, 3))
    fes = ng.H1(mesh, order=order, dirichlet=".*")
    u, v = fes.TnT()
    a = ng.BilinearForm(ng.grad(u) * ng.grad(v) * ng.dx); a.Assemble()
    gf = ng.GridFunction(fes); gf.Set(_solid_harmonic(n), ng.BND)
    rr = gf.vec.CreateVector(); rr.data = -(a.mat * gf.vec); gf.vec.data += a.mat.Inverse(fes.FreeDofs()) * rr
    energy = float(ng.Integrate(ng.grad(gf) * ng.grad(gf) * ng.dx(bonus_intorder=10), mesh))
    bmass = float(ng.Integrate(gf * gf * ng.ds(bonus_intorder=10), mesh))
    lam = -1.0 / R - energy / bmass
    return mesh.ne, mesh.GetNE(ng.BND), fes.ndof, abs(lam + (n + 1) / R) / ((n + 1) / R)

# minimal (curvature-limited) ball: maxh=1.0 already gives the coarsest sphere mesh
nev, neb, _, _ = kelvin_eig(1, 1.0, 1)
print("minimal unit-sphere Kelvin ball (curvature-limited): %d volume tets, %d Gamma triangles" % (nev, neb))
print("  -> volume-element count ~ surface-triangle count: the ball is a thin shell, no real overhead.\n")
print("DtN rel error on this minimal ball, sweeping order p (exact once p>=n):")
print("  {:>3} | {:>10} {:>10} {:>10}".format("n", "p=1", "p=2", "p=3"))
for n in (1, 2, 3):
    row = "  {:>3} |".format(n)
    for p in (1, 2, 3):
        _, _, nd, rel = kelvin_eig(n, 1.0, p)
        row += " {:>10.2e}".format(rel)
    print(row)
print("\n-> p>=n makes even the minimal ball exact (down to the curved-sphere geometry floor);")
print("   the cost tracks Gamma + order, like a surface method, but with a SPARSE matrix.")
