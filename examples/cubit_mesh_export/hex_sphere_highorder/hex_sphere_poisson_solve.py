# -*- coding: utf-8 -*-
"""Solve a PDE on the curved high-order HEX mesh -- not just integrate it.

The companion `hex_sphere_curved_ngsolve.py` shows that the curved hex .vol gives the
right VOLUME.  This one goes further: it SOLVES a Poisson problem on the same curved hex
sphere and shows the FE error fall as the element order rises -- i.e. NGSolve genuinely
assembles and solves on Cubit high-order hexahedra.

Manufactured solution  u_exact = sin(8x) cos(6y)  =>  -div(grad u) = f,
f = 100 sin(8x) cos(6y)  (in 3D; no z-dependence).  Dirichlet u = u_exact on the sphere
surface.  We report the L2 error vs order on the SAME 56-hex mesh (p-refinement):

    order 1 :  ~8e-6
    order 2 :  ~1.3e-6
    order 3 :  ~3.9e-7

Pure NGSolve, no Cubit needed to RUN (the .vol files are committed).  Load the high-order
.vol AS-IS and do NOT call mesh.Curve() (it would reset the curved nodes -- see README).

    python hex_sphere_poisson_solve.py
"""
import os
from ngsolve import (Mesh, H1, GridFunction, BilinearForm, LinearForm, Integrate,
                     grad, dx, x, y, sin, cos, sqrt, VOL, BND)

HERE = os.path.dirname(os.path.abspath(__file__))


def l2_error(order):
    mesh = Mesh(os.path.join(HERE, f"hexsph_o{order}.vol"))     # curved HEX, AS-IS
    u_ex = sin(8 * x) * cos(6 * y)
    f_rhs = 100 * sin(8 * x) * cos(6 * y)
    fes = H1(mesh, order=order, dirichlet=".*")
    gfu = GridFunction(fes)
    gfu.Set(u_ex, BND)
    u, v = fes.TnT()
    a = BilinearForm(grad(u) * grad(v) * dx).Assemble()
    L = LinearForm(f_rhs * v * dx).Assemble()
    r = L.vec - a.mat * gfu.vec
    gfu.vec.data += a.mat.Inverse(fes.FreeDofs()) * r
    nhex = sum(1 for el in mesh.Elements(VOL))
    return nhex, fes.ndof, sqrt(Integrate((gfu - u_ex) ** 2, mesh))


def main():
    print("Poisson on the curved hex sphere (p-refinement, fixed 56-hex mesh):")
    errs = []
    for order in (1, 2, 3):
        nhex, ndof, err = l2_error(order)
        errs.append(err)
        print(f"  order {order}  HEX={nhex}  ndof={ndof}  ||u_h - u_exact||_L2 = {err:.3e}")
    assert errs[1] < errs[0] and errs[2] < errs[1], "L2 error must fall with order"
    print(f"PASS: FE error falls with order on the curved hex mesh "
          f"({errs[0]:.2e} -> {errs[2]:.2e}) -- NGSolve solves on high-order Cubit hexes")


if __name__ == "__main__":
    main()
