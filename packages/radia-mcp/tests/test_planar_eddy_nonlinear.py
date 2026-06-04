"""Nonlinear time-harmonic eddy (solve.solve_planar_eddy_nonlinear, FEMM
nonlinear AC via amplitude-based effective permeability). Two checks:

  1. REDUCTION: with a constant nu_of_B it reproduces the linear solve_planar_eddy
     (Jz mode) to ~machine precision -- validates the Picard machinery.
  2. SATURATION: with a saturating nu(|B|), the iron's mean |B| grows SUB-LINEARLY
     with drive (5x current -> <5x flux) and the effective permeability drops --
     the physical hallmark of nonlinear AC.
"""
import math
import os
import sys

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from ngsolve import (Mesh, CoefficientFunction, grad, Conj, InnerProduct, sqrt,
                     IfPos, Integrate, TaskManager)
from netgen.occ import OCCGeometry, MoveTo, Glue, X, Y
from radia_mcp.radia_ngsolve.solve import (solve_planar_eddy,
                                           solve_planar_eddy_nonlinear,
                                           reluctivity, NU0)

MU0 = 4e-7 * math.pi
SIGMA_FE = 2.0e6
OMEGA = 2 * math.pi * 50.0
MUR0, BSAT = 2000.0, 1.8


def build_mesh():
    box = MoveTo(-1, -1).Rectangle(2, 2).Face()
    iron = MoveTo(-0.3, -0.3).Rectangle(0.6, 0.6).Face()
    iron.faces.name = "iron"; iron.maxh = 0.05
    coil = MoveTo(0.5, -0.3).Rectangle(0.2, 0.6).Face()
    coil.faces.name = "coil"; coil.maxh = 0.05
    air = box - iron - coil; air.faces.name = "air"
    for s in (air.edges.Max(X), air.edges.Min(X), air.edges.Max(Y), air.edges.Min(Y)):
        s.name = "outer"
    return Mesh(OCCGeometry(Glue([air, iron, coil]), dim=2).GenerateMesh(maxh=0.18))


def iron_mean_B(Az, mesh):
    B = CoefficientFunction((grad(Az)[1], -grad(Az)[0]))
    Bmag = sqrt((B[0] * Conj(B[0]) + B[1] * Conj(B[1])).real)
    area = Integrate(CoefficientFunction(1.0), mesh, definedon=mesh.Materials("iron"))
    return Integrate(Bmag, mesh, definedon=mesh.Materials("iron")) / area


def main():
    mesh = build_mesh()
    sigma = mesh.MaterialCF({"iron": SIGMA_FE}, default=0.0)
    iron_ind = mesh.MaterialCF({"iron": 1.0}, default=0.0)
    nu_lin = reluctivity(mesh, {"iron": MUR0})

    with TaskManager():
        # 1. reduction to linear: constant nu_of_B == nu_lin
        J0 = 1.0e5
        Jz = mesh.MaterialCF({"coil": J0}, default=0.0)
        A_lin = solve_planar_eddy(mesh, nu_lin, sigma, OMEGA, Jz=Jz)
        A_red = solve_planar_eddy_nonlinear(mesh, lambda Bmag: nu_lin, sigma,
                                            OMEGA, Jz=Jz, tol=1e-10, max_iter=200)
        dA = A_lin - A_red
        num = Integrate((dA * Conj(dA)).real, mesh)
        den = Integrate((A_lin * Conj(A_lin)).real, mesh)
        rel = math.sqrt(num / den)
        print(f"  [reduction] ||A_nl - A_lin|| / ||A_lin|| = {rel:.2e}")
        assert rel < 1e-3, f"nonlinear w/ const curve != linear ({rel:.2e})"

        # 2. saturation: physical saturating nu(|B|)
        def nu_of_B(Bmag):
            s = (Bmag / BSAT) ** 4 / (1.0 + (Bmag / BSAT) ** 4)
            nu_fe = NU0 * (1.0 / MUR0 + (1.0 - 1.0 / MUR0) * s)
            return iron_ind * nu_fe + (1.0 - iron_ind) * NU0

        Jlo = mesh.MaterialCF({"coil": 2.0e6}, default=0.0)
        Jhi = mesh.MaterialCF({"coil": 1.0e7}, default=0.0)   # 5x drive
        A_lo = solve_planar_eddy_nonlinear(mesh, nu_of_B, sigma, OMEGA, Jz=Jlo,
                                           relax=0.4, max_iter=150)
        A_hi = solve_planar_eddy_nonlinear(mesh, nu_of_B, sigma, OMEGA, Jz=Jhi,
                                           relax=0.4, max_iter=150)
        Blo, Bhi = iron_mean_B(A_lo, mesh), iron_mean_B(A_hi, mesh)
        ratio = Bhi / Blo
        print(f"  [saturation] <|B|>_iron: drive x1 -> {Blo:.4f} T, x5 -> {Bhi:.4f} T; "
              f"flux ratio = {ratio:.2f} (linear would be 5.0)")
        assert 1.0 < ratio < 4.5, f"expected sub-linear saturation, got {ratio:.2f}"

    print(f"\n[OK] nonlinear AC validated: reduces to linear ({rel:.1e}); "
          f"saturation makes flux sub-linear in drive (x5 -> x{ratio:.1f}).")


if __name__ == "__main__":
    main()
