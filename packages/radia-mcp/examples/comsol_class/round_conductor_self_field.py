# -*- coding: utf-8 -*-
"""Round-conductor self-field -- NGSolve 2D A_z FEM vs the closed form, and why the
*curved* element boundary matters.

A long straight conductor of radius a carrying a uniform current I has the exact field

    |B|(r) = mu0 I r / (2 pi a^2)   (r < a, inside),
    |B|(r) = mu0 I / (2 pi r)       (r > a, outside -- Ampere's law).

We solve the planar A_z Poisson problem  -div(grad A_z) = mu0 J  (J = I/(pi a^2) in the
conductor, 0 in air, A_z = 0 on a far circle) and read |B| = |grad A_z|.

LESSON (high-order geometry).  The source term integrates mu0 J over the *meshed*
conductor area.  With a straight-faceted circle that area is < pi a^2, so the enclosed
current -- and hence the whole exterior field -- comes out a few percent LOW *uniformly*
(this is NOT a far-boundary truncation error; it does not grow with r).  Calling
``mesh.Curve(order)`` restores the exact curved boundary and the field matches the closed
form to ~1e-4.  This is the same curved-element capability that lets NGSolve consume
high-order hex meshes (see packages/cubit-mesh-export).

Self-contained, headless:  python round_conductor_self_field.py
"""
import math
from ngsolve import (Mesh, H1, GridFunction, BilinearForm, LinearForm,
                     grad, dx, sqrt)
from netgen.geom2d import SplineGeometry

A_RAD, R_OUT, CURRENT = 0.01, 0.20, 1.0
MU0 = 4e-7 * math.pi


def solve_self_field(order=4, maxh=0.006, curve=True):
    """Return (mesh, |B| CoefficientFunction) for the round-conductor self-field."""
    Je = CURRENT / (math.pi * A_RAD ** 2)
    geo = SplineGeometry()
    geo.AddCircle((0, 0), R_OUT, leftdomain=1, bc="outer")
    geo.AddCircle((0, 0), A_RAD, leftdomain=2, rightdomain=1)
    geo.SetMaterial(1, "air")
    geo.SetMaterial(2, "cond")
    mesh = Mesh(geo.GenerateMesh(maxh=maxh))
    if curve:
        mesh.Curve(order)               # exact curved conductor boundary (the whole point)
    fes = H1(mesh, order=order, dirichlet="outer")
    u, v = fes.TnT()
    A = GridFunction(fes)
    K = BilinearForm(grad(u) * grad(v) * dx).Assemble()
    src = mesh.MaterialCF({"cond": MU0 * Je, "air": 0.0}, default=0.0)
    f = LinearForm(src * v * dx).Assemble()
    A.vec.data = K.mat.Inverse(fes.FreeDofs()) * f.vec
    return mesh, sqrt(grad(A) * grad(A))


def analytic(r):
    if r < A_RAD:
        return MU0 * CURRENT * r / (2 * math.pi * A_RAD ** 2)
    return MU0 * CURRENT / (2 * math.pi * r)


def main():
    radii = (0.005, 0.02, 0.05, 0.10)
    mesh, Bcf = solve_self_field(curve=True)
    print("Round-conductor self-field (curved boundary):")
    worst = 0.0
    for r in radii:
        val = Bcf(mesh(r, 0.0))
        ref = analytic(r)
        err = abs(val - ref) / ref * 100.0
        worst = max(worst, err)
        print(f"  r={r:.3f}  NGSolve|B|={val:.6e}  analytic={ref:.6e}  err={err:.3f}%")
    assert worst < 0.5, f"curved max err {worst:.3f}% should be < 0.5%"

    # contrast: the faceted (un-curved) circle under-integrates the source -> uniform deficit
    mesh0, B0 = solve_self_field(curve=False)
    r = 0.10
    err0 = abs(B0(mesh0(r, 0.0)) - analytic(r)) / analytic(r) * 100.0
    print(f"  [faceted, no mesh.Curve] r={r}: err={err0:.2f}%  -- curving removes this")
    print(f"PASS: curved field matches the closed form (max {worst:.3f}%)")


if __name__ == "__main__":
    main()
