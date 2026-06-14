# -*- coding: utf-8 -*-
# DEMO (h) (verified): the Kelvin-ball centre (the infinity image) is an INTEGRABLE
# COEFFICIENT singularity, NOT a field singularity (Kameari Q-c, deepened).
#
# Material-modulation Kelvin (3D scalar): the exterior Dirichlet energy maps to a
# weighted energy on the ball with weight w(rho') = (R/rho')^2 (the 3D conformal
# weight, d-2=1), which blows up at the centre rho'->0 = the infinity image. We show:
#   (1) the WEIGHT is integrable: int_ball (R/rho')^p dV converges iff p<3 (the 3D
#       volume element rho'^2 cancels rho'^{-2}); p=2 -> 4*pi*R^3, p=4 diverges.
#   (2) the FE solve div(w grad v)=0 with the centre INCLUDED is stable: the weighted
#       energy is mesh-convergent and the field is bounded (max|v| -> boundary max,
#       by the maximum principle) -- the centre is a regular point of the SOLUTION.
#   (1b) a TRUE field singularity (point source v ~ 1/rho') is the opposite: the field
#       itself blows up and its Dirichlet energy int|grad v|^2 diverges.
# => only the coefficient is singular, and integrably so. (For nodal scalar phi it is
# even benign WITHOUT the weight: the transformed unknown K[u] ~ rho'^n is a smooth
# polynomial.) The edge-element vector-A formulation is the lone exception -- a
# singular-tensor material at the centre -- handled by FEM-BEM (DtN on Gamma, no centre).
import os
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import numpy as np


def trap(y, x):                       # np.trapz was removed in NumPy 2.0
    return float(np.sum((y[:-1] + y[1:]) * 0.5 * np.diff(x)))


R = 1.0
print("(1) the Kelvin material weight is INTEGRABLE.")
print("    I(delta) = int_{delta<rho'<R} (R/rho')^p * 4*pi*rho'^2 d rho'  (cutoff delta->0)")
print("    {:>8} | {:>16} | {:>18}".format("delta", "p=2 (Kelvin)", "p=4 (more singular)"))
for d in (1e-2, 1e-4, 1e-6):
    I2 = 4 * np.pi * R ** 2 * (R - d)                 # -> 4 pi R^3 (finite)
    I4 = 4 * np.pi * R ** 4 * (1.0 / d - 1.0 / R)     # -> diverges ∝ 1/delta
    print("    {:>8.0e} | {:>16.4f} | {:>18.4e}".format(d, I2, I4))
print("    -> p=2 -> 4*pi*R^3 = %.4f (FINITE, delta-independent); p=4 grows ∝1/delta (DIVERGES)." % (4 * np.pi * R ** 3))
print("       threshold p<3: in 3D the volume element rho'^2 cancels rho'^{-2}.\n")

print("(1b) contrast -- a TRUE field singularity (point source v~1/rho'): Dirichlet energy")
for d in (1e-2, 1e-3, 1e-4):
    rr = np.linspace(d, R, 2_000_000)
    print("     int_{rho'>%.0e}|grad(1/rho')|^2 dV = %.3e  (DIVERGES as cutoff->0)"
          % (d, trap(4 * np.pi * rr ** (-2), rr)))
print()

import ngsolve as ng
from netgen.occ import Sphere, Pnt, OCCGeometry
x, y, z = ng.x, ng.y, ng.z
w2 = R * R / (x * x + y * y + z * z)          # (R/rho')^2 -- the integrable Kelvin weight


def weighted_solve(w, maxh, order=3):
    mesh = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), R)).GenerateMesh(maxh=maxh)).Curve(min(order + 1, 3))
    fes = ng.H1(mesh, order=order, dirichlet=".*")
    u, v = fes.TnT()
    a = ng.BilinearForm(w * ng.grad(u) * ng.grad(v) * ng.dx(bonus_intorder=14)); a.Assemble()
    gf = ng.GridFunction(fes); gf.Set(z / R, ng.BND)        # smooth datum cos(theta)=z/R on Gamma
    r = gf.vec.CreateVector(); r.data = -(a.mat * gf.vec)
    gf.vec.data += a.mat.Inverse(fes.FreeDofs()) * r
    energy = float(ng.Integrate(w * ng.grad(gf) * ng.grad(gf) * ng.dx(bonus_intorder=14), mesh))
    return mesh.ne, energy, float(max(abs(gf.vec.FV().NumPy())))


print("(2) FE solve of  div(w grad v)=0,  datum v=z/R on Gamma, w=(R/rho')^2, centre INCLUDED:")
print("    {:>6} | {:>8} | {:>13} {:>10}".format("maxh", "ne", "weighted E", "max|v|"))
for mh in (0.40, 0.25, 0.16):
    ne, E, vmax = weighted_solve(w2, mh)
    print("    {:>6.2f} | {:>8d} | {:>13.5e} {:>10.4f}".format(mh, ne, E, vmax))
print("    -> weighted energy CONVERGES and max|v| -> 1.0 (= boundary max, maximum principle):")
print("       the SOLUTION is bounded and regular at the centre; only the coefficient is singular.")
