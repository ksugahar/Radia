"""
NGSolve FEM ground truth for 3D cube scalar diffusion.

Solve: (-Lap + sMS) u = -sMS in [0,L]^3, u=0 on bnd.
Then Y(s) = Y_DC * (1 + <u>/V).

Uses 1/8 octant [0,L/2]^3 by symmetry, with Neumann at the three coordinate
planes x=0,y=0,z=0 and Dirichlet u=0 on the three outer faces x=L/2,y=L/2,z=L/2.

The Aitken Foster N=799 reference for the FULL cube at f=1e4 is |Y|=3.4319e+00.
The mixed Galerkin (closed K_ss) at rank=20 gives 3.443e+00 (0.33% above).
NGSolve FEM with refinement should tell us which is closer to true.
"""
import cmath
import math
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path("S:/Radia/01_GitHub/examples/hierarchical_cauer_sibc/mixed_galerkin").resolve()))
from _references.cube3d_foster import K_SIBC_cube3d, Y_DC_cube3d

from ngsolve import (
    Mesh, H1, BilinearForm, LinearForm, GridFunction,
    CoefficientFunction, grad, dx, Integrate, TaskManager,
    SetNumThreads
)
from netgen.occ import OCCGeometry, Box, Pnt

SIGMA = 5.8e7
MU = 4 * math.pi * 1e-7
L = 5e-3
MS = MU * SIGMA
Y_DC = Y_DC_cube3d(L, SIGMA)
V_CUBE = L**3

# References at f=1e4 (from earlier Aitken extrapolation + rank-N sweep)
Y_AITKEN_F4 = 3.431919  # |Y| from Foster N=799 Aitken
Y_MIXED_RANK20_F4 = 3.443338  # |Y| from rank-20 closed K_ss


def build_octant_mesh(maxh):
    """1/8 octant [0,L/2]^3.  Dirichlet on outer; Neumann on symmetry planes."""
    cube = Box(Pnt(0, 0, 0), Pnt(L/2, L/2, L/2))
    cube.mat("cube")
    for f in cube.faces:
        c = f.center
        # Outer faces have center at L/2 in their normal direction
        if (abs(c[0] - L/2) < 1e-9 or
            abs(c[1] - L/2) < 1e-9 or
            abs(c[2] - L/2) < 1e-9):
            f.name = "outer"
        else:
            f.name = "sym"
    geo = OCCGeometry(cube)
    ngmesh = geo.GenerateMesh(maxh=maxh)
    return Mesh(ngmesh)


def solve_cube_at_f(mesh, f, order):
    """Solve (-Lap + sMS) u = -sMS, u|outer=0.  Return Y, ndof, V_octant."""
    s = 1j * 2 * math.pi * f
    sMS = s * MS

    fes = H1(mesh, order=order, complex=True, dirichlet="outer")
    u, v = fes.TnT()
    a = BilinearForm(fes, symmetric=True)
    a += grad(u) * grad(v) * dx
    a += sMS * u * v * dx

    F = LinearForm(fes)
    F += -sMS * v * dx

    a.Assemble()
    F.Assemble()

    gfu = GridFunction(fes)
    gfu.vec.data = a.mat.Inverse(
        freedofs=fes.FreeDofs(), inverse="sparsecholesky"
    ) * F.vec

    u_int = Integrate(gfu, mesh)        # complex integral
    V_oct = Integrate(CoefficientFunction(1), mesh)
    u_avg = u_int / V_oct                # = <u>_full by symmetry
    Y = Y_DC * (1 + u_avg)
    return Y, fes.ndof, V_oct


def main():
    SetNumThreads(8)
    print("=== NGSolve cube ground truth (octant, scalar diffusion) ===")
    print(f"L = {L*1e3} mm, sigma = {SIGMA:.2e}")
    print(f"Y_DC = {Y_DC:.6e}, V_cube = {V_CUBE:.4e}")
    print()
    print(f"Reference at f=1e4:")
    print(f"  Aitken Foster N=799:    |Y| = {Y_AITKEN_F4:.6e}")
    print(f"  mixed Galerkin rank=20: |Y| = {Y_MIXED_RANK20_F4:.6e}  ({(Y_MIXED_RANK20_F4/Y_AITKEN_F4-1)*100:+.4f}%)")
    print()

    # Refinement sweep at f=1e4
    print("Refinement sweep at f=1e4 Hz (gL=10.7, skin depth=L/10.7):")
    print(f"  {'maxh(mm)':>10}  {'ord':>3}  {'ne':>6}  {'ndof':>8}  {'t_solve':>8}  {'|Y_FEM|':>14}  {'Y/Y_A':>10}  {'rel(%)':>10}")
    for (maxh, order) in [
        (L/4, 2),
        (L/4, 3),
        (L/6, 3),
        (L/8, 3),
        (L/10, 3),
        (L/12, 3),
        (L/8, 4),
        (L/10, 4),
    ]:
        try:
            mesh = build_octant_mesh(maxh)
            with TaskManager():
                t0 = time.time()
                Y, ndof, V_oct = solve_cube_at_f(mesh, 1e4, order)
                dt = time.time() - t0
            ratio = abs(Y) / Y_AITKEN_F4
            err = (ratio - 1) * 100
            print(f"  {maxh*1e3:10.3f}  {order:3d}  {mesh.ne:6d}  {ndof:8d}  {dt:7.2f}s  {abs(Y):14.6e}  {ratio:10.5f}  {err:+9.4f}%")
        except Exception as e:
            print(f"  {maxh*1e3:10.3f}  {order:3d}  FAILED: {e}")


if __name__ == "__main__":
    main()
