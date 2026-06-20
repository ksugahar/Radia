# -*- coding: utf-8 -*-
# DEMO (m) (verified): on the MINIMAL (curvature-limited) Kelvin ball, how many multipoles
# can the closure capture? Answer: raising the order p captures degrees up to n ~ p -- the
# ORDER is the multipole reach, on a fixed coarse surface. The Gamma element count sets a
# much higher ceiling (~N_Gamma x p) that p, not the mesh, hits first in this range.
#
# Measured (112 Gamma triangles): p=3 exact to n~3 then blows up at n>=5; p=5 -> n~5;
# p=7 -> n~7; p=9 -> still 9e-6 at n=10. => the source's highest multipole n_src dictates
# p (p>=n_src), and a single coarse ball serves any n_src by raising p.
import os
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import ngsolve as ng
from netgen.occ import Sphere, Pnt, OCCGeometry

R = 1.0
x, y, z = ng.x, ng.y, ng.z
rho2 = x * x + y * y + z * z


def Hn(n):
    # zonal solid harmonic H_n = rho^n P_n(z/rho); n H_n = (2n-1) z H_{n-1} - (n-1) rho^2 H_{n-2}
    H = [ng.CoefficientFunction(1.0), z]
    for k in range(2, n + 1):
        H.append(((2 * k - 1) * z * H[k - 1] - (k - 1) * rho2 * H[k - 2]) / k)
    return H[n]


def kelvin_eig(n, order):
    mesh = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), R)).GenerateMesh(maxh=1.0)).Curve(min(order + 1, 4))
    fes = ng.H1(mesh, order=order, dirichlet=".*")
    u, v = fes.TnT()
    a = ng.BilinearForm(ng.grad(u) * ng.grad(v) * ng.dx); a.Assemble()
    gf = ng.GridFunction(fes); gf.Set(Hn(n), ng.BND)
    rr = gf.vec.CreateVector(); rr.data = -(a.mat * gf.vec); gf.vec.data += a.mat.Inverse(fes.FreeDofs()) * rr
    energy = float(ng.Integrate(ng.grad(gf) * ng.grad(gf) * ng.dx(bonus_intorder=2 * order), mesh))
    bmass = float(ng.Integrate(gf * gf * ng.ds(bonus_intorder=2 * order), mesh))
    return mesh.GetNE(ng.BND), abs(-1.0 / R - energy / bmass + (n + 1) / R) / ((n + 1) / R)


nb, _ = kelvin_eig(1, 1)
print("minimal Kelvin ball: %d Gamma triangles. DtN rel error vs multipole degree n and order p:" % nb)
ns = list(range(1, 11))
print("\n  p \\ n |" + "".join("%8d" % n for n in ns))
print("  " + "-" * (8 + 8 * len(ns)))
for p in (3, 5, 7, 9):
    row = "  %5d |" % p
    for n in ns:
        _, rel = kelvin_eig(n, p)
        row += "%8.1e" % rel
    print(row)
print("\n-> error stays at the floor up to n ~ p, then blows up once the harmonic out-resolves the")
print("   coarse surface. The order IS the multipole reach; one coarse ball serves any n_src via p.")
