# -*- coding: utf-8 -*-
"""Coaxial line -- NGSolve 2D A_z FEM vs the closed form, and self-shielding.

An inner wire (radius a, current +I) and a coaxial return shell (r1<r<r2, current -I)
carry equal and opposite currents.  By Ampere's law:

    a < r < r1 (gap)        : |B| = mu0 I / (2 pi r)
    r > r2 (outside)        : |B| = 0     (enclosed current = +I - I = 0)

The return conductor SHIELDS the outside: the field is confined between the conductors.
We solve the planar A_z problem (J = +I/(pi a^2) in the wire, -I/(pi (r2^2-r1^2)) in the
shell) and read |B| = |grad A_z|.

Self-contained, headless:  python coaxial_line_field.py
"""
import math
from ngsolve import (Mesh, H1, GridFunction, BilinearForm, LinearForm,
                     grad, dx, sqrt)
from netgen.geom2d import SplineGeometry

A_RAD, R1, R2, R_OUT = 0.01, 0.05, 0.06, 0.20
CURRENT = 1.0
MU0 = 4e-7 * math.pi


def solve_coax(order=4, maxh=0.006):
    Jw = CURRENT / (math.pi * A_RAD ** 2)
    Js = -CURRENT / (math.pi * (R2 ** 2 - R1 ** 2))
    geo = SplineGeometry()
    geo.AddCircle((0, 0), R_OUT, leftdomain=1, bc="outer")
    geo.AddCircle((0, 0), R2, leftdomain=2, rightdomain=1)
    geo.AddCircle((0, 0), R1, leftdomain=3, rightdomain=2)
    geo.AddCircle((0, 0), A_RAD, leftdomain=4, rightdomain=3)
    for i, n in enumerate(("airout", "shell", "gap", "cond"), start=1):
        geo.SetMaterial(i, n)
    mesh = Mesh(geo.GenerateMesh(maxh=maxh))
    mesh.Curve(order)
    fes = H1(mesh, order=order, dirichlet="outer")
    u, v = fes.TnT()
    A = GridFunction(fes)
    K = BilinearForm((1.0 / MU0) * grad(u) * grad(v) * dx).Assemble()
    src = mesh.MaterialCF({"cond": Jw, "shell": Js}, default=0.0)
    f = LinearForm(src * v * dx).Assemble()
    A.vec.data = K.mat.Inverse(fes.FreeDofs()) * f.vec
    return mesh, sqrt(grad(A) * grad(A))


def main():
    mesh, Bcf = solve_coax()
    B_gap = MU0 * CURRENT / (2 * math.pi * 0.03)
    gap = Bcf(mesh(0.03, 0.0))
    out = Bcf(mesh(0.10, 0.0))
    print("Coaxial line (inner +I, return shell -I):")
    print(f"  r=0.03 (gap)     NGSolve|B|={gap:.6e}  analytic={B_gap:.6e}  err={abs(gap-B_gap)/B_gap*100:.3f}%")
    print(f"  r=0.10 (outside) NGSolve|B|={out:.3e}  analytic=0   (|B_out|/|B_gap|={out/B_gap:.2e})")
    assert abs(gap - B_gap) / B_gap < 0.01, "gap field should match mu0 I/2pi r"
    assert out / B_gap < 1e-4, "the return current should shield the outside"
    print("PASS: field confined between conductors; outside is shielded to ~0")


if __name__ == "__main__":
    main()
