"""
NGSolve FEM scalar diffusion on a 3D L-shape.

Geometry:
  Horizontal bar:  [0, 2L] x [0, L] x [0,  L]
  Vertical bar:    [0,  L] x [0, L] x [L, 2L]
  Union = 3D L (volume 3 L^3, surface 14 L^2, one concave dihedral edge
                along x=L, z=L, y in [0,L]).

Solve (-Lap + sMS) u = -sMS  in V_L,  u = 0 on partial V_L.
Then Y(s) = sigma * V_L * (1 + <u>/V_L).

Test:
  - Y(0) = sigma * V_L  (DC)
  - |Y(s)| sqrt(s) / K_SIBC -> 1  at high f  (Mellin universality)
    with K_SIBC = S_L * sqrt(sigma / mu) = 14 L^2 sqrt(sigma/mu).
"""
import cmath
import math
import sys
import time
from pathlib import Path

import numpy as np

from ngsolve import (
    Mesh, H1, BilinearForm, LinearForm, GridFunction,
    CoefficientFunction, grad, dx, Integrate, TaskManager,
    SetNumThreads
)
from netgen.occ import OCCGeometry, Box, Pnt, Glue

SIGMA = 5.8e7
MU = 4 * math.pi * 1e-7
L = 5e-3
MS = MU * SIGMA

V_L = 3 * L**3              # L-shape volume
S_L = 14 * L**2             # L-shape surface area
Y_DC_L = SIGMA * V_L        # 21.75 for L=5mm, sigma=5.8e7
K_SIBC_L = S_L * math.sqrt(SIGMA / MU)


def build_lshape(maxh):
    """3D L-shape: union of two boxes meeting at a shared face."""
    horiz = Box(Pnt(0, 0, 0), Pnt(2*L, L, L))
    vert = Box(Pnt(0, 0, L), Pnt(L, L, 2*L))
    shape = Glue([horiz, vert])
    shape.mat("Lshape")
    # All outer faces are Dirichlet (u=0).  Internal shared face (x in [0,L],
    # y in [0,L], z=L) is interior, no BC.
    for f in shape.faces:
        f.name = "outer"
    geo = OCCGeometry(shape)
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
    V_int = Integrate(CoefficientFunction(1), mesh)
    u_avg = Integrate(gfu, mesh) / V_int
    Y = Y_DC_L * (1 + u_avg)
    return Y, fes.ndof, V_int


def main():
    SetNumThreads(8)
    print("=== NGSolve 3D L-shape Y(s) + Mellin universality ===")
    print(f"L = {L*1e3} mm")
    print(f"V_L = 3 L^3 = {V_L:.4e} m^3")
    print(f"S_L = 14 L^2 = {S_L:.4e} m^2")
    print(f"Y_DC = sigma V_L = {Y_DC_L:.4e} S")
    print(f"K_SIBC^L = S_L sqrt(sigma/mu) = {K_SIBC_L:.4e}")
    print()

    # Mesh convergence at f = 1e4 (wall-band, gL = 10.7 on L-side)
    print("Mesh convergence at f = 1e4 (test ne / ndof vs |Y_FEM|):")
    print(f"  {'maxh(mm)':>10}  {'ord':>3}  {'ne':>6}  {'ndof':>8}  {'t':>6}  {'|Y_FEM|':>14}  {'V_int':>14}")
    Y_conv = None
    for (maxh, order) in [(L/4, 3), (L/6, 3), (L/8, 3), (L/10, 3), (L/8, 4)]:
        try:
            mesh = build_lshape(maxh)
            with TaskManager():
                t0 = time.time()
                Y, ndof, V_int = solve_at(mesh, 1e4, order)
                dt = time.time() - t0
            print(f"  {maxh*1e3:10.3f}  {order:3d}  {mesh.ne:6d}  {ndof:8d}  {dt:5.2f}s  {abs(Y):14.6e}  {V_int:14.6e}")
            Y_conv = Y
        except Exception as e:
            print(f"  {maxh*1e3:10.3f}  {order:3d}  FAILED: {e}")

    # Frequency sweep at converged resolution
    print()
    print("Frequency sweep at maxh=L/8, P3 (converged):")
    print(f"  {'f (Hz)':>10}  {'|gL|':>8}  {'|Y_FEM|':>14}  {'arg':>8}  {'|Y_FEM|*sqrt|s|':>16}  {'/K_SIBC':>10}")
    mesh = build_lshape(L/8)
    for f in [1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10]:
        try:
            with TaskManager():
                Y, _, _ = solve_at(mesh, f, 3)
            s = 1j * 2 * math.pi * f
            gL = abs(cmath.sqrt(s * MS) * L)
            ratio = abs(Y) * math.sqrt(abs(s)) / K_SIBC_L
            arg = math.degrees(cmath.phase(Y))
            print(f"  {f:10.2e}  {gL:8.2f}  {abs(Y):14.6e}  {arg:7.2f}d  {abs(Y)*math.sqrt(abs(s)):16.4e}  {ratio:10.5f}")
        except Exception as e:
            print(f"  {f:10.2e}  FAILED: {e}")
    print()
    print("Expected: |Y_FEM|*sqrt|s| / K_SIBC -> 1 at high f (Mellin universality)")


if __name__ == "__main__":
    main()
