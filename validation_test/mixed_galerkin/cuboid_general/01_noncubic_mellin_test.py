"""
Phase 1 experiment: non-cubic cuboid Lx × Ly × Lz, all dihedral angles
still pi/2 but tensor-product structure broken (Lx != Ly != Lz).

Goal: verify the GENERALIZED Mellin asymptote
    c_0 = S_total * sqrt(sigma/mu)              (S = 2(Lx Ly + Ly Lz + Lx Lz))
    c_1 = -(16/(pi mu)) * (Lx + Ly + Lz)        (sum of edge half-lengths)
    c_2 = +48 / (pi mu^1.5 sqrt(sigma))         (vertex term, dim-independent)

against NGSolve FEM Y(s) on a non-cubic cuboid.

If FEM matches the generalized Mellin at high f to within ~1%, the
codim decomposition is validated for any rectangular cuboid -- not
just cubes.  This is the simplest scope-test of the Mixed Galerkin
framework beyond the canonical cube benchmark.

Reference for cube case: cube3d_foster.Y_mellin_cube3d (which fixes
a=b=c=L).  The cuboid generalization just substitutes general (Lx,Ly,Lz).

Test sizes: pick aspect ratios that break tensor symmetry but stay
within meshable bounds:
    case A: 5 mm x 7 mm x 3 mm   (moderate aspect)
    case B: 4 mm x 4 mm x 4 mm   (cube control, must match cube3d_foster)
"""
from __future__ import annotations

import cmath
import math
import sys
import time
from pathlib import Path

import numpy as np

# Reference modules (we reuse the cube one for sanity check at Lx=Ly=Lz)
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from _references.cube3d_foster import Y_mellin_cube3d

from ngsolve import (
    Mesh, H1, BilinearForm, LinearForm, GridFunction,
    CoefficientFunction, grad, dx, Integrate, TaskManager,
    SetNumThreads
)
from netgen.occ import OCCGeometry, Box, Pnt


SIGMA = 5.8e7
MU = 4 * math.pi * 1e-7
MS = MU * SIGMA


# ---------------------------------------------------------------------------
# Generalized Mellin asymptote for cuboid Lx x Ly x Lz
# ---------------------------------------------------------------------------


def Y_DC_cuboid(Lx, Ly, Lz, sigma):
    """DC admittance Y_DC = V * sigma = Lx Ly Lz * sigma (S*m^2)."""
    return Lx * Ly * Lz * sigma


def K_SIBC_cuboid(Lx, Ly, Lz, sigma, mu):
    """c_0 = S * sqrt(sigma/mu), where S = 2(Lx Ly + Ly Lz + Lx Lz)."""
    S = 2.0 * (Lx * Ly + Ly * Lz + Lx * Lz)
    return S * math.sqrt(sigma / mu)


def Y_mellin_cuboid(s, Lx, Ly, Lz, sigma, mu):
    """Generalized Mellin asymptote:
        c_0/sqrt(s) + c_1/s + c_2/s^(3/2)
    with
        c_0 = S_total * sqrt(sigma/mu)
        c_1 = -(16/(pi mu)) * (Lx + Ly + Lz)
        c_2 = +48 / (pi mu^1.5 sqrt(sigma))
    """
    c_0 = K_SIBC_cuboid(Lx, Ly, Lz, sigma, mu)
    c_1 = -(16.0 / (math.pi * mu)) * (Lx + Ly + Lz)
    c_2 = 48.0 / (math.pi * mu**1.5 * math.sqrt(sigma))
    return c_0 / cmath.sqrt(s) + c_1 / s + c_2 / s**1.5


# ---------------------------------------------------------------------------
# NGSolve FEM ground truth
# ---------------------------------------------------------------------------


def build_cuboid_mesh(Lx, Ly, Lz, maxh):
    """Full cuboid mesh (no symmetry reduction)."""
    box = Box(Pnt(0, 0, 0), Pnt(Lx, Ly, Lz))
    box.mat("cuboid")
    for f in box.faces:
        f.name = "outer"
    geo = OCCGeometry(box)
    return Mesh(geo.GenerateMesh(maxh=maxh))


def build_cuboid_octant(Lx, Ly, Lz, maxh):
    """Octant cuboid: Lx/2 x Ly/2 x Lz/2 with sym/outer face labels.

    Outer faces are at x=Lx/2, y=Ly/2, z=Lz/2.
    Symmetry faces are at x=0, y=0, z=0 (Neumann).
    """
    box = Box(Pnt(0, 0, 0), Pnt(Lx / 2, Ly / 2, Lz / 2))
    box.mat("cuboid")
    for f in box.faces:
        c = f.center
        if (abs(c[0] - Lx / 2) < 1e-9 or
            abs(c[1] - Ly / 2) < 1e-9 or
            abs(c[2] - Lz / 2) < 1e-9):
            f.name = "outer"
        else:
            f.name = "sym"
    geo = OCCGeometry(box)
    return Mesh(geo.GenerateMesh(maxh=maxh))


def solve_at(mesh, f, order, Y_DC):
    """Scalar diffusion FEM: (-Lap + s mu sigma) u = -s mu sigma in V,
       u = 0 on outer.  Returns Y = Y_DC (1 + <u>) where <u> = Integral(u)/V."""
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


def run_case(Lx, Ly, Lz, label, maxh_frac=10, order=3):
    """Run one cuboid case: FEM at several f vs generalized Mellin."""
    print(f"\n=== Case {label}: {Lx*1e3:.1f} x {Ly*1e3:.1f} x {Lz*1e3:.1f} mm ===")
    Lmin = min(Lx, Ly, Lz)
    maxh = Lmin / maxh_frac
    Y_DC = Y_DC_cuboid(Lx, Ly, Lz, SIGMA)
    K_SIBC = K_SIBC_cuboid(Lx, Ly, Lz, SIGMA, MU)
    print(f"  V = {Lx*Ly*Lz*1e9:.3f} mm^3, S = {2*(Lx*Ly+Ly*Lz+Lx*Lz)*1e6:.3f} mm^2")
    print(f"  Y_DC = {Y_DC:.4e}, K_SIBC = {K_SIBC:.4e}")
    print(f"  Mesh: octant maxh = Lmin/{maxh_frac} = {maxh*1e3:.4f} mm, order = {order}")

    mesh = build_cuboid_octant(Lx, Ly, Lz, maxh)
    print(f"  ne = {mesh.ne}")
    print(f"  {'f (Hz)':>10}  {'|gL|':>8}  {'|Y_FEM|':>14}  {'|Y_Mellin|':>14}  {'rel err':>10}  {'time (s)':>9}")

    results = []
    for f in [1e3, 1e4, 1e5, 1e6, 1e7, 1e8]:
        s = 1j * 2 * math.pi * f
        with TaskManager():
            t0 = time.time()
            Y_FEM, ndof = solve_at(mesh, f, order, Y_DC)
            dt = time.time() - t0
        Y_M = Y_mellin_cuboid(s, Lx, Ly, Lz, SIGMA, MU)
        gL = abs(cmath.sqrt(s * MS) * Lmin)
        rel = abs(Y_FEM - Y_M) / abs(Y_M) * 100
        print(f"  {f:10.2e}  {gL:8.2f}  {abs(Y_FEM):14.6e}  {abs(Y_M):14.6e}  {rel:9.4f}%  {dt:9.2f}")
        results.append((f, abs(Y_FEM), abs(Y_M), rel))

    return results


def main():
    SetNumThreads(8)
    print("=== Phase 1: non-cubic cuboid generalized Mellin verification ===")
    print(f"sigma = {SIGMA:.2e} S/m, mu = {MU:.4e} H/m")
    print("Octant mesh: 8x speedup vs full cuboid; Dirichlet on outer (3 faces),")
    print("Neumann on sym (3 faces).  scalar diffusion form, Y = Y_DC (1 + <u>).")

    # Case B: cube control (must reduce to existing cube3d Mellin)
    print()
    print("--- Sanity check: cube 4x4x4 mm (Lx=Ly=Lz) ---")
    print("Should reduce to Y_mellin_cube3d exactly.")
    L_cube = 4e-3
    # Compare both Mellin formulas at one f
    s_test = 1j * 2 * math.pi * 1e6
    Y_M_general = Y_mellin_cuboid(s_test, L_cube, L_cube, L_cube, SIGMA, MU)
    Y_M_cube = Y_mellin_cube3d(s_test, L_cube, SIGMA, MU)
    rel = abs(Y_M_general - Y_M_cube) / abs(Y_M_cube) * 100
    print(f"  Y_mellin_cuboid(L,L,L) = {abs(Y_M_general):.6e}")
    print(f"  Y_mellin_cube3d(L)     = {abs(Y_M_cube):.6e}")
    print(f"  rel diff               = {rel:.2e}%  (should be 0)")

    # Case A: non-cubic 5 x 7 x 3 mm
    results_A = run_case(5e-3, 7e-3, 3e-3, label="A", maxh_frac=8, order=3)

    # Case B FEM: cube 4x4x4 (gives the cube reference value direct from FEM)
    results_B = run_case(4e-3, 4e-3, 4e-3, label="B (cube control)", maxh_frac=8, order=3)

    # Case C: another non-cubic shape, more extreme aspect
    results_C = run_case(3e-3, 8e-3, 2e-3, label="C (3:8:2 aspect)", maxh_frac=6, order=3)

    print()
    print("=== Summary at high f (where Mellin should dominate) ===")
    print(f"{'Case':>20}  {'f=1e7':>10}  {'f=1e8':>10}")
    for label, results in [("A (5x7x3)", results_A), ("B (cube 4x4x4)", results_B), ("C (3x8x2)", results_C)]:
        err_1e7 = next((r[3] for r in results if abs(r[0] - 1e7) < 1), None)
        err_1e8 = next((r[3] for r in results if abs(r[0] - 1e8) < 1), None)
        print(f"  {label:>20}  {err_1e7:9.4f}%  {err_1e8:9.4f}%")


if __name__ == "__main__":
    main()
