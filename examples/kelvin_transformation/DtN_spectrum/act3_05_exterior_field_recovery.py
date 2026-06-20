# -*- coding: utf-8 -*-
# DEMO (p) (verified): the Kelvin closure DOES give the EXTERIOR field at any point --
# correcting the misconception that it only yields the interior field (so it is NOT a
# BEM-only advantage). The Kelvin ball IS the compactified exterior: solving on the ball
# gives u' inside it, and the physical field at any exterior point x (|x|>R) is the inverse
# Kelvin transform of that solution,
#       u(x) = (R/rho) u'(x'),   x' = R^2 x / rho^2,   rho=|x|     (3D weight (R/rho)^{d-2}).
# Every exterior point x maps to an interior ball point x'; the field there is recovered by
# a single FE point-evaluation + the conformal weight. No surface integral, no Green function.
#
# Setup: exterior Laplace (decaying) with Dirichlet datum on Gamma = trace of the solid
# harmonic P_n; the exterior solution is u_ext(x) = R^{2n+1} P_n(x)/rho^{2n+1}  (n=1: z/rho^3
# dipole; n=2: (2z^2-x^2-y^2)/rho^5 quadrupole). The Kelvin image of P_n on the ball is P_n
# itself (degree-n polynomial, FE-exact at p>=n).
#
# Measured (R=1, order 3, coarse ball): the inverse-Kelvin field matches the analytic exterior
# field from just outside Gamma (rho=1.2) to far away (rho=20):
#   dipole    n=1: rel err 9e-7 .. 2.4e-5  (FE-exact polynomial)
#   quadrupole n=2: rel err 6e-4 .. 2.3e-2  -- the larger FAR relative error is the ball's
#     ABSOLUTE geometry floor (~1e-5 on the curved sphere) amplified where the mapped value is
#     tiny (far points map near the centre, value ~rho^-(n+1)); the absolute error stays ~1e-5.
# => Kelvin recovers the exterior field everywhere to the ball's FE accuracy. The genuine
#    BEM/H-matrix edge is NOT exterior-field access (Kelvin has it) but arbitrary NON-spherical
#    truncation Gamma (Kelvin's inversion is spherical) + mature oscillatory-kernel FMM.
import os
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import numpy as np
import ngsolve as ng
from netgen.occ import Sphere, Pnt, OCCGeometry
from radia_mcp.radia_ngsolve.fem_bem_coupling import _solid_harmonic

R = 1.0


def solve_ball(n, order=3, maxh=0.5):
    mesh = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), R)).GenerateMesh(maxh=maxh)).Curve(min(order + 1, 3))
    fes = ng.H1(mesh, order=order, dirichlet=".*")
    u, v = fes.TnT()
    a = ng.BilinearForm(ng.grad(u) * ng.grad(v) * ng.dx); a.Assemble()
    gf = ng.GridFunction(fes); gf.Set(_solid_harmonic(n), ng.BND)
    rr = gf.vec.CreateVector(); rr.data = -(a.mat * gf.vec)
    gf.vec.data += a.mat.Inverse(fes.FreeDofs()) * rr
    return mesh, gf


def u_ext_analytic(n, X, Y, Z):
    rho = np.sqrt(X * X + Y * Y + Z * Z)
    Pn = Z if n == 1 else (2 * Z * Z - X * X - Y * Y)
    return R ** (2 * n + 1) * Pn / rho ** (2 * n + 1)


pts = [(0, 0, 1.2), (0, 0, 2.0), (0, 0, 5.0), (0, 0, 20.0),
       (2.0, 0, 2.0), (3.0, 1.0, 1.0), (0, 4.0, 8.0)]

for n in (1, 2):
    mesh, gf = solve_ball(n, order=3, maxh=0.5)
    print("\nn=%d  (%s):  u(x)=(R/rho) u'(R^2 x/rho^2)  (inverse Kelvin)"
          % (n, "dipole z/rho^3" if n == 1 else "quadrupole (2z^2-x^2-y^2)/rho^5"))
    print("   exterior x            |x|     x'=R^2 x/|x|^2 in ball     u_recovered    u_analytic     rel.err")
    for (X, Y, Z) in pts:
        rho = np.sqrt(X * X + Y * Y + Z * Z)
        xp = (R * R * X / rho**2, R * R * Y / rho**2, R * R * Z / rho**2)
        u_rec = (R / rho) * gf(mesh(*xp))
        u_an = u_ext_analytic(n, X, Y, Z)
        rel = abs(u_rec - u_an) / (abs(u_an) + 1e-300)
        print("  (%4.1f,%4.1f,%5.1f)  %6.2f   (%5.3f,%5.3f,%5.3f)      %12.6e  %12.6e  %.2e"
              % (X, Y, Z, rho, xp[0], xp[1], xp[2], u_rec, u_an, rel))

print("\n-> the exterior field at ANY point is recovered from the ball solution by the inverse Kelvin")
print("   map (far points map near the centre = smooth low-multipole tail). Kelvin is NOT 'interior")
print("   field only' -- the ball IS the exterior; the field is available everywhere to FE accuracy.")
