"""
NGSolve cube FEM frequency sweep.  Now that we've established convergence
at order=4, maxh=L/10 (8263 DOFs), use that mesh to settle ALL mid-f
points and compare to (a) Aitken Foster, (b) mixed Galerkin rank-20.
"""
import cmath
import math
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path("S:/Radia/01_GitHub/examples/hierarchical_cauer_sibc/mixed_galerkin").resolve()))
from _references.cube3d_foster import K_SIBC_cube3d, Y_DC_cube3d, Y_mellin_cube3d

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


# References (from earlier sweeps)
Y_AITKEN = {
    1e3: 6.892354,
    1e4: 3.431919,
    1e5: 1.218962,
    1e6: 0.399611,
    1e7: 0.126191,  # underconverged Aitken
}
# Mellin asymptote |Y_M|
def Y_mellin_abs(f):
    return abs(Y_mellin_cube3d(1j * 2 * math.pi * f, L, SIGMA, MU))


def build_octant_mesh(maxh):
    cube = Box(Pnt(0, 0, 0), Pnt(L/2, L/2, L/2))
    cube.mat("cube")
    for f in cube.faces:
        c = f.center
        if (abs(c[0] - L/2) < 1e-9 or
            abs(c[1] - L/2) < 1e-9 or
            abs(c[2] - L/2) < 1e-9):
            f.name = "outer"
        else:
            f.name = "sym"
    geo = OCCGeometry(cube)
    return Mesh(geo.GenerateMesh(maxh=maxh))


def solve_at(mesh, f, order):
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
    gfu.vec.data = a.mat.Inverse(freedofs=fes.FreeDofs(), inverse="sparsecholesky") * F.vec
    u_avg = Integrate(gfu, mesh) / Integrate(CoefficientFunction(1), mesh)
    return Y_DC * (1 + u_avg), fes.ndof


def main():
    SetNumThreads(8)
    print("=== NGSolve cube FEM frequency sweep (high-resolution) ===")
    print(f"L = {L*1e3} mm")
    print()

    # Two resolutions to check convergence at each frequency
    for (maxh, order, label) in [(L/10, 4, "L/10, P4"), (L/12, 4, "L/12, P4")]:
        mesh = build_octant_mesh(maxh)
        print(f"\nMesh: {label}, ne={mesh.ne}")
        print(f"  {'f (Hz)':>10}  {'|gL|':>8}  {'|Y_FEM|':>14}  {'|Y_Aitken|':>14}  {'|Y_Mellin|':>14}  {'FEM/A':>10}  {'FEM/Mellin':>12}")
        for f in [1e3, 1e4, 1e5, 1e6, 1e7]:
            with TaskManager():
                t0 = time.time()
                Y, ndof = solve_at(mesh, f, order)
                dt = time.time() - t0
            gL = abs(cmath.sqrt(1j * 2 * math.pi * f * MS) * L)
            Y_A = Y_AITKEN[f]
            Y_M = Y_mellin_abs(f)
            print(f"  {f:10.2e}  {gL:8.2f}  {abs(Y):14.6e}  {Y_A:14.6e}  {Y_M:14.6e}  {abs(Y)/Y_A:10.5f}  {abs(Y)/Y_M:11.5f}")


if __name__ == "__main__":
    main()
