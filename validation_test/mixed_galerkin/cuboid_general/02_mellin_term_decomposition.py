"""
Phase 1 follow-up: term-by-term Mellin decomposition.

Verify that each Mellin term (c_0, c_1, c_2) contributes the right amount
on non-cubic cuboids.  Specifically:

    Y_c0_only  = c_0 / sqrt(s)
    Y_c0_c1    = c_0 / sqrt(s) + c_1 / s
    Y_full     = c_0 / sqrt(s) + c_1 / s + c_2 / s^(3/2)

If the generalization (Lx + Ly + Lz) in c_1 is correct, then:
    - cube case (Lx=Ly=Lz=L): Y_c0_c1 should match the published cube Mellin
    - non-cubic case (Lx != Ly != Lz): Y_c0_c1 with new (Lx+Ly+Lz) should
      beat Y_c0_only by a measurable margin (proving c_1 has the right shape)

If using the WRONG c_1 (e.g., 3*L_avg) for a non-cubic case, the formula
would deviate from FEM by a few percent at mid-f.  We will measure this.
"""
from __future__ import annotations

import cmath
import math
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from ngsolve import (
    Mesh, H1, BilinearForm, LinearForm, GridFunction,
    CoefficientFunction, grad, dx, Integrate, TaskManager,
    SetNumThreads
)
from netgen.occ import OCCGeometry, Box, Pnt


SIGMA = 5.8e7
MU = 4 * math.pi * 1e-7
MS = MU * SIGMA


def K_SIBC_cuboid(Lx, Ly, Lz):
    S = 2.0 * (Lx * Ly + Ly * Lz + Lx * Lz)
    return S * math.sqrt(SIGMA / MU)


def c1_correct(Lx, Ly, Lz):
    """Correct: c_1 = -16 (Lx+Ly+Lz) / (pi mu)."""
    return -16.0 * (Lx + Ly + Lz) / (math.pi * MU)


def c1_wrong_avg(Lx, Ly, Lz):
    """Wrong: pretend cube shape, c_1 = -48 * L_avg / (pi mu) where L_avg = (Lx+Ly+Lz)/3."""
    L_avg = (Lx + Ly + Lz) / 3.0
    return -48.0 * L_avg / (math.pi * MU)
    # Note: this is actually the SAME as the correct formula!  c_1 only depends on Lx+Ly+Lz.
    # Let me try a TRULY wrong c_1 below to distinguish.


def c1_wrong_geom_mean(Lx, Ly, Lz):
    """Wrong: c_1 = -48 * (Lx*Ly*Lz)^(1/3) / (pi mu).
    Geometric mean instead of arithmetic sum -- breaks for non-cubic shapes.
    """
    L_gm = (Lx * Ly * Lz) ** (1.0 / 3.0)
    return -48.0 * L_gm / (math.pi * MU)


def c2_value():
    """c_2 = +48 / (pi mu^1.5 sqrt(sigma)).  Dim-independent (8 vertices, 90 deg solid angle)."""
    return 48.0 / (math.pi * MU**1.5 * math.sqrt(SIGMA))


def Y_mellin_terms(s, Lx, Ly, Lz, c1_func):
    c0 = K_SIBC_cuboid(Lx, Ly, Lz)
    c1 = c1_func(Lx, Ly, Lz)
    c2 = c2_value()
    Y_c0 = c0 / cmath.sqrt(s)
    Y_c0_c1 = Y_c0 + c1 / s
    Y_full = Y_c0_c1 + c2 / s**1.5
    return Y_c0, Y_c0_c1, Y_full


def build_cuboid_octant(Lx, Ly, Lz, maxh):
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
    return Y_DC * (1 + u_avg)


def run_case(Lx, Ly, Lz, label, maxh_frac=8, order=3):
    print(f"\n=== Case {label}: {Lx*1e3:.1f} x {Ly*1e3:.1f} x {Lz*1e3:.1f} mm ===")
    Lmin = min(Lx, Ly, Lz)
    maxh = Lmin / maxh_frac
    Y_DC = Lx * Ly * Lz * SIGMA

    mesh = build_cuboid_octant(Lx, Ly, Lz, maxh)

    print(f"  Y_DC = {Y_DC:.4e}, V = {Lx*Ly*Lz*1e9:.3f} mm^3")
    print(f"  c_0 (K_SIBC)        = {K_SIBC_cuboid(Lx, Ly, Lz):.4e}")
    print(f"  c_1 (correct)       = {c1_correct(Lx, Ly, Lz):.4e}")
    print(f"  c_1 (geom mean WRONG) = {c1_wrong_geom_mean(Lx, Ly, Lz):.4e}")
    print(f"  c_2                 = {c2_value():.4e}")

    print(f"\n  Term-by-term comparison at f = 1e5 (where c_0 dominates):")
    print(f"  {'f (Hz)':>10}  {'delta_skin/L':>13}  {'|Y_FEM|':>13}  {'err_c0only':>11}  {'err_c0+c1_corr':>15}  {'err_c0+c1_wrong':>16}  {'err_full':>10}")

    for f in [3e4, 1e5, 3e5]:
        s = 1j * 2 * math.pi * f
        delta = 1.0 / math.sqrt(math.pi * f * MS)  # skin depth
        delta_over_L = delta / Lmin

        with TaskManager():
            Y_FEM = solve_at(mesh, f, order, Y_DC)

        Y_c0_corr, Y_c0c1_corr, Y_full_corr = Y_mellin_terms(s, Lx, Ly, Lz, c1_correct)
        _, Y_c0c1_wrong, _ = Y_mellin_terms(s, Lx, Ly, Lz, c1_wrong_geom_mean)

        err_c0 = abs(Y_FEM - Y_c0_corr) / abs(Y_FEM) * 100
        err_c0c1_corr = abs(Y_FEM - Y_c0c1_corr) / abs(Y_FEM) * 100
        err_c0c1_wrong = abs(Y_FEM - Y_c0c1_wrong) / abs(Y_FEM) * 100
        err_full = abs(Y_FEM - Y_full_corr) / abs(Y_FEM) * 100

        print(f"  {f:10.2e}  {delta_over_L:13.4f}  {abs(Y_FEM):13.6e}  {err_c0:10.4f}%  {err_c0c1_corr:14.4f}%  {err_c0c1_wrong:15.4f}%  {err_full:9.4f}%")


def main():
    SetNumThreads(8)
    print("=== Phase 1 follow-up: Mellin term decomposition ===")
    print("Test whether each term (c_0, c_1, c_2) of the generalized Mellin")
    print("contributes correctly for non-cubic cuboids.")
    print()
    print("WRONG variant tested: c_1 with geometric mean of (Lx,Ly,Lz)")
    print("(instead of arithmetic sum).  Non-cubic shapes should reveal the bug.")

    run_case(4e-3, 4e-3, 4e-3, label="B (cube control)", maxh_frac=8, order=3)
    run_case(5e-3, 7e-3, 3e-3, label="A (5x7x3 moderate)", maxh_frac=8, order=3)
    run_case(3e-3, 8e-3, 2e-3, label="C (3x8x2 extreme)", maxh_frac=6, order=3)


if __name__ == "__main__":
    main()
