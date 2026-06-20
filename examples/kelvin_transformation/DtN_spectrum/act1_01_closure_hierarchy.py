# -*- coding: utf-8 -*-
# DEMO (g) (verified): the truncation-BC hierarchy as PARTIAL DtN spectra, and why
# the value of Kelvin is controllability, not exactness.
#
# Every closure at Gamma approximates the exterior DtN operator, whose sphere
# eigenvalues are the ladder lambda_n = -(n+1)/R. Rank the cheap closures by how
# much of the ladder they reproduce:
#
#   Dirichlet u=0   : reproduces NO eigenvalue (forces the trace to 0)  -> worst.
#   minimal Robin   : d_n u = lambda_1 u, lambda_1 = -2/R  reproduces the n=1
#                     (DIPOLE) eigenvalue EXACTLY -> dipole-exact at ANY R.
#   Kelvin / BEM    : reproduces the whole ladder (p>=n) -> exact up to the floor.
#
# Kameari's canonical problem (magnetic sphere in a uniform field) has a PURE DIPOLE
# exterior, so the minimal Robin is already exact there -- you do NOT need Kelvin for
# ~10% work, but you must NOT use Dirichlet. A non-spherical source adds n>=2 modes,
# which a single Robin cannot reach (it nails ONE mode) and Kelvin captures mode by
# mode as the order p rises. So Kelvin is the controllable closure: coarse mesh + p
# climbs the ladder; sparse/local; per-mode error predictable from the spectrum.
#
# Verified (3D shell, R_inner=a=1, rel-L2 field error vs the analytic field):
#   DIPOLE  : Dirichlet ~0.46->0.22 (R/a 1.5->5)   Robin(lam1) ~1e-3 (R-independent)
#   QUADRUPOLE: Dirichlet ~0.13->0.36  Robin(lam1, mis-tuned) ~3-7%  Robin(lam2) ~2e-3
# Companion knowledge: dtn_coarse_mesh; companion slide (Kameari report, backup).
import os
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import ngsolve as ng
from netgen.occ import Sphere, Pnt, OCCGeometry
from radia_mcp.radia_ngsolve.fem_bem_coupling import _solid_harmonic

ngx, ngy, ngz = ng.x, ng.y, ng.z


def shell(R_inner, R_outer, maxh, order):
    o = Sphere(Pnt(0, 0, 0), R_outer); o.bc("outer")
    i = Sphere(Pnt(0, 0, 0), R_inner); i.bc("inner")
    return ng.Mesh(OCCGeometry(o - i).GenerateMesh(maxh=maxh)).Curve(min(order + 1, 3))


def solve(mesh, R_inner, degree, closure, lam_robin, order=3, io=8):
    inner_datum = _solid_harmonic(degree) / R_inner ** degree
    fes = ng.H1(mesh, order=order, dirichlet=("inner|outer" if closure == "dirichlet" else "inner"))
    u, v = fes.TnT()
    a = ng.BilinearForm(ng.grad(u) * ng.grad(v) * ng.dx(bonus_intorder=io))
    if closure == "robin":
        a += (-lam_robin) * u * v * ng.ds("outer", bonus_intorder=io)
    a.Assemble()
    gf = ng.GridFunction(fes)
    gf.Set(inner_datum, ng.BND, definedon=mesh.Boundaries("inner"))   # outer -> 0 (Dir) or free (Robin)
    r = gf.vec.CreateVector(); r.data = -(a.mat * gf.vec)
    gf.vec.data += a.mat.Inverse(fes.FreeDofs()) * r
    return gf


def relerr(gf, mesh, R_inner, degree):
    r2 = ngx ** 2 + ngy ** 2 + ngz ** 2
    u_ex = R_inner ** (degree + 1) * _solid_harmonic(degree) / r2 ** ((2 * degree + 1) / 2.0)
    num = float(ng.Integrate((gf - u_ex) ** 2 * ng.dx(bonus_intorder=8), mesh))
    den = float(ng.Integrate(u_ex ** 2 * ng.dx(bonus_intorder=8), mesh))
    return (num / den) ** 0.5


a = 1.0
print("DIPOLE source (magnetic sphere in a uniform field) -- exterior is a PURE n=1 mode.")
print("{:>6} | {:>13} | {:>13}".format("R/a", "Dirichlet u=0", "Robin lam1=-2/R"))
print("-" * 40)
for R in (1.5, 2.0, 3.0, 5.0):
    m = shell(a, R, 0.30, 3)
    print("{:>6.1f} | {:>13.2e} | {:>13.2e}".format(
        R / a, relerr(solve(m, a, 1, "dirichlet", 0.0), m, a, 1),
        relerr(solve(m, a, 1, "robin", -2.0 / R), m, a, 1)))

print("\nQUADRUPOLE source (n=2): a single Robin nails only ONE mode.")
print("{:>6} | {:>13} | {:>16} | {:>15}".format("R/a", "Dirichlet u=0", "Robin lam1(-2/R)", "Robin lam2(-3/R)"))
print("-" * 60)
for R in (1.5, 2.0, 3.0):
    m = shell(a, R, 0.30, 3)
    print("{:>6.1f} | {:>13.2e} | {:>16.2e} | {:>15.2e}".format(
        R / a, relerr(solve(m, a, 2, "dirichlet", 0.0), m, a, 2),
        relerr(solve(m, a, 2, "robin", -2.0 / R), m, a, 2),
        relerr(solve(m, a, 2, "robin", -3.0 / R), m, a, 2)))
print("\n-> Dirichlet drops even the dipole (~40%); minimal Robin is dipole-exact; only "
      "Kelvin/BEM\n   climb the whole ladder (p>=n) -- the controllable closure.")
