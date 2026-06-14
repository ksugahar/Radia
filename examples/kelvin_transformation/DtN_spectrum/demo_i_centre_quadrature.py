# -*- coding: utf-8 -*-
# DEMO (i) (verified): do you have to keep nodes / Gauss points off the Kelvin-ball
# centre, or does proper integration of the (R/rho')^2 weight just work?
#
# Answer: it just works. The weight is integrable (see demo_h), so the weak-form
# integrals are finite; the ONLY thing that breaks is a quadrature point landing
# EXACTLY on the centre (w=inf -> NaN), and standard simplex rules put points strictly
# interior to the element (barycentric > 0), so they never hit a vertex. Below, the
# netgen sphere mesh has NO vertex at the centre (nearest vertex ~0.18 R, so the centre
# sits INSIDE an element), and the solve is still stable: the energy is essentially
# independent of the quadrature order and the field is bounded. Regularising the weight
# as R^2/(rho^2+eps^2) with eps->0 gives the identical answer.
import os
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import numpy as np
import ngsolve as ng
from netgen.occ import Sphere, Pnt, OCCGeometry

R = 1.0
x, y, z = ng.x, ng.y, ng.z
rho2 = x * x + y * y + z * z


def solve(maxh, order, bonus, eps=0.0):
    mesh = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), R)).GenerateMesh(maxh=maxh)).Curve(min(order + 1, 3))
    pts = np.array([mesh[v].point for v in mesh.vertices])
    rmin = float(np.min(np.linalg.norm(pts, axis=1)))         # distance of nearest vertex to the centre
    w = R * R / (rho2 + eps * eps)
    fes = ng.H1(mesh, order=order, dirichlet=".*")
    u, v = fes.TnT()
    a = ng.BilinearForm(w * ng.grad(u) * ng.grad(v) * ng.dx(bonus_intorder=bonus)); a.Assemble()
    gf = ng.GridFunction(fes); gf.Set(z / R, ng.BND)
    r = gf.vec.CreateVector(); r.data = -(a.mat * gf.vec)
    gf.vec.data += a.mat.Inverse(fes.FreeDofs()) * r
    E = float(ng.Integrate(w * ng.grad(gf) * ng.grad(gf) * ng.dx(bonus_intorder=max(bonus, 12)), mesh))
    return rmin, E, float(max(abs(gf.vec.FV().NumPy())))


print("quadrature sensitivity (maxh=0.25, order=3); the centre is INSIDE an element:")
print("  {:>14} | {:>12} {:>11} {:>9}".format("bonus_intorder", "rmin_vertex", "energy", "max|v|"))
for bonus in (0, 2, 4, 8, 14):
    rmin, E, vmax = solve(0.25, 3, bonus)
    print("  {:>14d} | {:>12.3e} {:>11.5e} {:>9.4f}".format(bonus, rmin, E, vmax))

print("\nregularisation w=R^2/(rho^2+eps^2), eps->0 (bonus=8):")
print("  {:>8} | {:>11} {:>9}".format("eps", "energy", "max|v|"))
for eps in (1e-1, 1e-2, 1e-3, 1e-6, 0.0):
    _, E, vmax = solve(0.25, 3, 8, eps)
    print("  {:>8.0e} | {:>11.5e} {:>9.4f}".format(eps, E, vmax))
print("\n-> nearest vertex ~0.18R (no centre node) yet stable & quadrature-insensitive; eps->0")
print("   matches eps=0. Proper integration suffices; only a Gauss point EXACTLY at the")
print("   centre (never produced by standard rules) would poison the sum.")
