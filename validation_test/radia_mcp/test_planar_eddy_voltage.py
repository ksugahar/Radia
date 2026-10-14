"""VOLTAGE-DRIVEN 2D planar eddy conductor (FEMM voltage-driven circuit):
prescribe the axial field Vc = applied_Ez and let the net current follow. The
per-length impedance Z = Vc / I has Re(Z) = Rac, validated against the round-wire
Kelvin (ber/bei) AC resistance -- the same wire as the current-driven test,
confirming both drive modes agree.
"""
import math
import os
import sys

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from scipy.special import kelvin
from ngsolve import Mesh, CoefficientFunction, Integrate, dx, TaskManager
from netgen.occ import OCCGeometry, MoveTo, WorkPlane, Glue, X, Y
from radia_mcp.radia_ngsolve.solve import solve_planar_eddy

MU0 = 4e-7 * math.pi
A, SIGMA, Q = 1e-3, 5.8e7, 4.0
OMEGA = Q * Q / (A * A * MU0 * SIGMA)
R_FAR = 12 * A
VC = 1.0   # applied axial E-field [V/m]


def rac_over_rdc_exact(q):
    be, ke, bep, kep = kelvin(q)
    return (q / 2.0) * (be.real * bep.imag - be.imag * bep.real) / (bep.real ** 2 + bep.imag ** 2)


def build_mesh():
    delta = math.sqrt(2.0 / (OMEGA * MU0 * SIGMA))
    wire = WorkPlane().Circle(0, 0, A).Face(); wire.faces.name = "wire"; wire.maxh = delta / 4.0
    box = MoveTo(-R_FAR, -R_FAR).Rectangle(2 * R_FAR, 2 * R_FAR).Face()
    air = box - wire; air.faces.name = "air"
    for sel in (air.edges.Max(X), air.edges.Min(X), air.edges.Max(Y), air.edges.Min(Y)):
        sel.name = "outer"
    return Mesh(OCCGeometry(Glue([air, wire]), dim=2).GenerateMesh(maxh=R_FAR / 8))


def main():
    mesh = build_mesh()
    sigma = mesh.MaterialCF({"wire": SIGMA}, default=0.0)
    nu = CoefficientFunction(1.0 / MU0)
    with TaskManager():
        Az = solve_planar_eddy(mesh, nu, sigma, OMEGA, applied_Ez=VC, order=3)
        Ez = -1j * OMEGA * Az + VC          # in the wire (sigma=0 elsewhere)
        I = Integrate(sigma * Ez * dx, mesh)   # net current [A]
    Z = VC / I                                  # per-length impedance [ohm/m]
    Rac = Z.real
    Rdc = 1.0 / (SIGMA * math.pi * A * A)
    ratio_fem = Rac / Rdc
    ratio_exact = rac_over_rdc_exact(Q)
    err = (ratio_fem - ratio_exact) / ratio_exact
    print(f"  voltage-driven: I={I:.4e} A   Z={Z:.5e} ohm/m")
    print(f"  Rac/Rdc  FEM={ratio_fem:.5f}  exact(Kelvin)={ratio_exact:.5f}  err={100*err:+.2f}%")
    assert abs(err) < 0.01, f"voltage-driven Rac off by {100*err:.2f}% (>1%)"
    print(f"\n[OK] voltage-driven eddy conductor validated vs Kelvin Rac ({100*abs(err):.2f}%).")


if __name__ == "__main__":
    main()
