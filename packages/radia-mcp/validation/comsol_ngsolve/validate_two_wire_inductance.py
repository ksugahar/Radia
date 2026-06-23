r"""Two-wire line external inductance -- regression test (transmission-line parameter).

Two parallel wires (radius a, separation D) carrying +-I; the line-dipole field energy
converges, so L_ext = 2W/I^2 = (mu0/pi) acosh(D/2a) (two_wire_external_inductance). With the
two-wire C = pi eps0/acosh(D/2a): v = 1/sqrt(LC) = c. radia gets L from force.magnetic_energy_2d.
Exact analytic reference -> tool-independent. The magnetic twin of the two-wire capacitance.
"""
import math
import os
import sys

import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen")

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from ngsolve import Mesh, CoefficientFunction, grad
from netgen.occ import OCCGeometry, WorkPlane, Glue
from radia_mcp.radia_ngsolve.solve import (solve_planar_magnetostatic, reluctivity,
                                           two_wire_external_inductance)
from radia_mcp.radia_ngsolve.force import magnetic_energy_2d, MU0

EPS0 = 8.8541878128e-12
a, D, Rout, I = 0.001, 0.010, 0.12, 1.0


def validate_two_wire_external_inductance():
    w1 = WorkPlane().Circle(-D/2, 0, a).Face(); w1.faces.name = "w1"
    w2 = WorkPlane().Circle(D/2, 0, a).Face(); w2.faces.name = "w2"
    box = WorkPlane().Circle(0, 0, Rout).Face(); box.edges.name = "outer"
    air = box - w1 - w2; air.faces.name = "air"
    w1.faces.maxh = a/4; w2.faces.maxh = a/4; air.faces.maxh = Rout/8
    for w in (w1, w2):
        for e in w.edges:
            e.maxh = a/4
    mesh = Mesh(OCCGeometry(Glue([w1, w2, air]), dim=2).GenerateMesh(maxh=Rout/8))
    Jz = CoefficientFunction([{"w1": I/(math.pi*a*a), "w2": -I/(math.pi*a*a)}.get(m, 0.0)
                             for m in mesh.GetMaterials()])
    A = solve_planar_magnetostatic(mesh, reluctivity(mesh, {}), Jz=Jz, order=3, dirichlet="outer")
    B = CoefficientFunction((grad(A)[1], -grad(A)[0]))
    L_fe = 2.0 * magnetic_energy_2d(B, mesh, "air") / (I*I)
    L_cf = two_wire_external_inductance(D, a)
    C_cf = math.pi*EPS0/math.acosh(D/(2*a))
    v_fe = 1.0/math.sqrt(L_fe*C_cf)
    relL = abs(L_fe - L_cf)/L_cf
    print(f"[two-wire L] FE={L_fe:.4e} exact={L_cf:.4e} ({100*relL:.2f}%)  v/c={v_fe/(1/math.sqrt(MU0*EPS0)):.4f}")
    assert relL < 0.03, f"L_ext FE {L_fe:.3e} vs exact {L_cf:.3e} ({100*relL:.1f}%)"
    assert abs(v_fe/(1/math.sqrt(MU0*EPS0)) - 1.0) < 0.03, "v = 1/sqrt(LC) should be ~c"


if __name__ == "__main__":
    validate_two_wire_external_inductance()
    print("[OK] two-wire line external inductance validated.")
