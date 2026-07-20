"""FEMM "stranded" / litz conductor model (solve.stranded_source): an imposed
UNIFORM current density (no eddy redistribution). Validated against the DC
uniform-current field and contrasted with a SOLID conductor's skin effect.

For a round wire of radius a carrying total current I with UNIFORM density, the
interior field is the DC Ampere profile  B_theta(r) = mu0 I r / (2 pi a^2)
at ANY frequency (a stranded bundle stays uniform). A solid conductor at the
same high frequency pushes current to the surface, so its interior B is lower.
"""
import math
import os
import sys

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from ngsolve import (Mesh, CoefficientFunction, grad, Conj, InnerProduct, sqrt,
                     Integrate, TaskManager)
from netgen.occ import OCCGeometry, MoveTo, WorkPlane, Glue, X, Y
from radia_mcp.radia_ngsolve.solve import solve_planar_eddy, stranded_source

MU0 = 4e-7 * math.pi
A, SIGMA, I_TOT, Q = 1e-3, 5.8e7, 1.0, 4.0
OMEGA = Q * Q / (A * A * MU0 * SIGMA)
R_FAR = 12 * A


def build_mesh():
    delta = math.sqrt(2.0 / (OMEGA * MU0 * SIGMA))
    wire = WorkPlane().Circle(0, 0, A).Face(); wire.faces.name = "wire"
    wire.maxh = delta / 4.0
    box = MoveTo(-R_FAR, -R_FAR).Rectangle(2 * R_FAR, 2 * R_FAR).Face()
    air = box - wire; air.faces.name = "air"
    for s in (air.edges.Max(X), air.edges.Min(X), air.edges.Max(Y), air.edges.Min(Y)):
        s.name = "outer"
    return Mesh(OCCGeometry(Glue([air, wire]), dim=2).GenerateMesh(maxh=R_FAR / 8))


def bmag_at(Az, mesh, px, py):
    B = CoefficientFunction((grad(Az)[1], -grad(Az)[0]))
    return sqrt((B[0] * Conj(B[0]) + B[1] * Conj(B[1])).real)(mesh(px, py))


def main():
    mesh = build_mesh()
    nu = CoefficientFunction(1.0 / MU0)
    B_dc_half = MU0 * I_TOT * (A / 2) / (2 * math.pi * A * A)  # uniform-J profile

    with TaskManager():
        # stranded: uniform imposed current, sigma=0 (no eddy)
        Jz = stranded_source(mesh, {"wire": I_TOT})
        I_check = Integrate(Jz, mesh, definedon=mesh.Materials("wire"))
        gfu_s = solve_planar_eddy(mesh, nu, CoefficientFunction(0.0), OMEGA, Jz=Jz)
        b_strand = bmag_at(gfu_s, mesh, A / 2, 0.0)

        # solid: current-driven, full eddy (skin effect) at the same omega
        sigma = mesh.MaterialCF({"wire": SIGMA}, default=0.0)
        gfu_solid = solve_planar_eddy(mesh, nu, sigma, OMEGA,
                                      driven_region="wire", total_current=I_TOT)
        b_solid = bmag_at(gfu_solid.components[0], mesh, A / 2, 0.0)

    e_strand = (b_strand - B_dc_half) / B_dc_half
    print(f"  total current check: int Jz = {I_check:.6f} (want {I_TOT})")
    print(f"  |B|(a/2): stranded={b_strand:.6e} (DC exact {B_dc_half:.6e}, "
          f"{100*e_strand:+.2f}%)  solid(skin)={b_solid:.6e} "
          f"(ratio solid/DC={b_solid/B_dc_half:.3f})")
    assert abs(I_check - I_TOT) < 1e-6, "stranded total current not preserved"
    assert abs(e_strand) < 0.02, f"stranded B off DC profile {100*e_strand:.2f}%"
    assert b_solid < 0.9 * B_dc_half, "solid should show skin-effect field drop"
    print(f"\n[OK] stranded conductor = uniform-current DC field ({100*abs(e_strand):.2f}%), "
          f"distinct from solid skin effect (solid/DC={b_solid/B_dc_half:.2f}).")


if __name__ == "__main__":
    main()
