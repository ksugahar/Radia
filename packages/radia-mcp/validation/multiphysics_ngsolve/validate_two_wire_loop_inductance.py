r"""Two-wire line TOTAL loop inductance (external + internal) -- regression test (#48).

L_loop = (mu0/pi) acosh(D/2a) + 2*(mu0/8pi) = external (#30) + both wires' internal (#36).
radia partitions the FE field energy (air -> external, wire interiors -> internal); the closed
form is the tool-independent gate. The +-I pair's energy converges, so the internal part is
exact (mu0/8pi each, radius-independent) and the external is the acosh term to ~1.8 %."""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.solve import (two_wire_external_inductance,
                                           two_wire_loop_inductance,
                                           internal_inductance_round_wire)

MU0 = 4e-7 * math.pi
A, D, ROUT, I = 0.001, 0.010, 0.12, 1.0


def validate_loop_inductance_formula():
    Le = two_wire_external_inductance(D, A)
    Ll = two_wire_loop_inductance(D, A)
    # loop = external + 2 * internal (mu0/8pi each)
    assert math.isclose(Ll, Le + MU0 / (4 * math.pi), rel_tol=1e-12)
    assert math.isclose(Ll, Le + 2 * internal_inductance_round_wire(), rel_tol=1e-12)
    assert Ll > Le                                   # loop exceeds external-only
    # the internal term is radius-independent (mu0/4pi)
    assert math.isclose(two_wire_loop_inductance(D, A) - two_wire_external_inductance(D, A),
                        two_wire_loop_inductance(2 * D, 2 * A) - two_wire_external_inductance(2 * D, 2 * A),
                        rel_tol=1e-12)


def validate_two_wire_loop_fe():
    pytest.importorskip("ngsolve")
    pytest.importorskip("netgen")
    from ngsolve import Mesh, CoefficientFunction, grad
    from netgen.occ import OCCGeometry, WorkPlane, Glue
    from radia_mcp.radia_ngsolve.solve import solve_planar_magnetostatic, reluctivity
    from radia_mcp.radia_ngsolve.force import magnetic_energy_2d

    w1 = WorkPlane().Circle(-D/2, 0, A).Face(); w1.faces.name = "w1"
    w2 = WorkPlane().Circle(D/2, 0, A).Face(); w2.faces.name = "w2"
    box = WorkPlane().Circle(0, 0, ROUT).Face(); box.edges.name = "outer"
    air = box - w1 - w2; air.faces.name = "air"
    w1.faces.maxh = A/4; w2.faces.maxh = A/4; air.faces.maxh = ROUT/8
    for w in (w1, w2):
        for e in w.edges:
            e.maxh = A/4
    mesh = Mesh(OCCGeometry(Glue([w1, w2, air]), dim=2).GenerateMesh(maxh=ROUT/8))
    Jz = CoefficientFunction([{"w1": I/(math.pi*A*A), "w2": -I/(math.pi*A*A)}.get(m, 0.0)
                             for m in mesh.GetMaterials()])
    Az = solve_planar_magnetostatic(mesh, reluctivity(mesh, {}), Jz=Jz, order=3, dirichlet="outer")
    B = CoefficientFunction((grad(Az)[1], -grad(Az)[0]))
    L_int = 2 * (magnetic_energy_2d(B, mesh, "w1") + magnetic_energy_2d(B, mesh, "w2")) / (I*I)
    L_loop = 2 * (magnetic_energy_2d(B, mesh, "air")
                  + magnetic_energy_2d(B, mesh, "w1") + magnetic_energy_2d(B, mesh, "w2")) / (I*I)
    # internal part exact (boundary-independent), total within the #30 energy-truncation margin
    assert abs(L_int - MU0 / (4 * math.pi)) / (MU0 / (4 * math.pi)) < 0.01
    assert abs(L_loop - two_wire_loop_inductance(D, A)) / two_wire_loop_inductance(D, A) < 0.03
    assert L_loop > two_wire_external_inductance(D, A)   # includes the internal


if __name__ == "__main__":
    validate_loop_inductance_formula()
    validate_two_wire_loop_fe()
    print("[OK] two-wire loop inductance = external acosh + internal mu0/4pi validated.")
