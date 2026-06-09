# -*- coding: utf-8 -*-
"""Current wire inside an iron shell -- NGSolve 2D A_z FEM vs the closed form.

A straight wire (radius a, current I) sits inside a coaxial iron shell (r1<r<r2,
relative permeability mu_r).  By Ampere's law the field is cylindrically symmetric,
H = I/(2 pi r), so |B|(r) is exact in every region:

    r < a            : mu0 I r / (2 pi a^2)     (inside the wire)
    a < r < r1, r>r2 : mu0 I / (2 pi r)          (air)
    r1 < r < r2      : mu0 * mu_r * I / (2 pi r)  (IRON -- the field is mu_r times larger)

We solve the planar A_z magnetostatic problem  div( (1/mu) grad A_z ) = -J  with
nu = 1/(mu0 mu_r) in the iron and 1/mu0 elsewhere, and read |B| = |grad A_z|.  The
shell carries mu_r-times the flux density of the surrounding air at the same radius --
the basis of flux concentration / magnetic shielding.

Self-contained, headless:  python iron_shell_wire_field.py
"""
import math
from ngsolve import (Mesh, H1, GridFunction, BilinearForm, LinearForm,
                     grad, dx, sqrt)
from netgen.geom2d import SplineGeometry

A_RAD, R1, R2, R_OUT = 0.01, 0.04, 0.08, 0.20
CURRENT, MU_R = 1.0, 100.0
MU0 = 4e-7 * math.pi


def solve_iron_shell(order=4, maxh=0.006):
    Je = CURRENT / (math.pi * A_RAD ** 2)
    geo = SplineGeometry()
    geo.AddCircle((0, 0), R_OUT, leftdomain=1, bc="outer")
    geo.AddCircle((0, 0), R2, leftdomain=2, rightdomain=1)
    geo.AddCircle((0, 0), R1, leftdomain=3, rightdomain=2)
    geo.AddCircle((0, 0), A_RAD, leftdomain=4, rightdomain=3)
    for i, n in enumerate(("airout", "iron", "gap", "cond"), start=1):
        geo.SetMaterial(i, n)
    mesh = Mesh(geo.GenerateMesh(maxh=maxh))
    mesh.Curve(order)                                   # curved circular interfaces
    fes = H1(mesh, order=order, dirichlet="outer")
    u, v = fes.TnT()
    A = GridFunction(fes)
    nu = mesh.MaterialCF({"iron": 1.0 / (MU0 * MU_R)}, default=1.0 / MU0)
    K = BilinearForm(nu * grad(u) * grad(v) * dx).Assemble()
    src = mesh.MaterialCF({"cond": Je}, default=0.0)
    f = LinearForm(src * v * dx).Assemble()
    A.vec.data = K.mat.Inverse(fes.FreeDofs()) * f.vec
    return mesh, sqrt(grad(A) * grad(A))


def analytic(r):
    if r < A_RAD:
        return MU0 * CURRENT * r / (2 * math.pi * A_RAD ** 2)
    if R1 <= r <= R2:
        return MU0 * MU_R * CURRENT / (2 * math.pi * r)
    return MU0 * CURRENT / (2 * math.pi * r)


def main():
    mesh, Bcf = solve_iron_shell()
    print(f"Iron-shell wire field (mu_r={MU_R}):")
    worst = 0.0
    for r, where in ((0.03, "air gap"), (0.06, "IRON"), (0.12, "outer air")):
        val = Bcf(mesh(r, 0.0))
        ref = analytic(r)
        err = abs(val - ref) / ref * 100.0
        worst = max(worst, err)
        print(f"  r={r:.3f} ({where:<9}) NGSolve|B|={val:.6e}  analytic={ref:.6e}  err={err:.3f}%")
    assert worst < 0.5, f"max err {worst:.3f}% should be < 0.5%"
    print(f"PASS: iron shell carries mu_r x the air flux density; FE matches closed form (max {worst:.3f}%)")


if __name__ == "__main__":
    main()
