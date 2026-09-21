r"""Internal inductance of a round wire = mu0/(8 pi) per unit length -- regression test.

A wire with a UNIFORM (DC) current has L_int = 2 W_in/I^2 = mu0/(8 pi) [H/m] from the
field energy INSIDE the conductor (internal_inductance_round_wire) -- INDEPENDENT of radius.
radia uses force.magnetic_energy_2d on the wire region (boundary-independent: Ampere fixes the
interior B). Checked at two radii (radius-independence). Tool-independent (exact analytic).
The internal companion to the external coax/two-wire inductances.
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
                                           internal_inductance_round_wire)
from radia_mcp.radia_ngsolve.force import magnetic_energy_2d

I = 1.0


def _L_internal(R, Rout):
    wire = WorkPlane().Circle(0, 0, R).Face(); wire.faces.name = "wire"
    box = WorkPlane().Circle(0, 0, Rout).Face(); box.edges.name = "outer"
    air = box - wire; air.faces.name = "air"
    wire.faces.maxh = R/12; air.faces.maxh = Rout/8
    mesh = Mesh(OCCGeometry(Glue([wire, air]), dim=2).GenerateMesh(maxh=Rout/8))
    Jz = CoefficientFunction([{"wire": I/(math.pi*R*R)}.get(m, 0.0) for m in mesh.GetMaterials()])
    A = solve_planar_magnetostatic(mesh, reluctivity(mesh, {}), Jz=Jz, order=3, dirichlet="outer")
    B = CoefficientFunction((grad(A)[1], -grad(A)[0]))
    return 2.0 * magnetic_energy_2d(B, mesh, "wire") / (I*I)


def validate_internal_inductance_round_wire():
    L_cf = internal_inductance_round_wire()
    Ls = {R: _L_internal(R, 20*R) for R in (0.001, 0.002)}
    for R, L in Ls.items():
        rel = abs(L - L_cf) / L_cf
        print(f"[internal L] R={R*1e3:.1f}mm  FE={L:.4e}  exact mu0/8pi={L_cf:.4e}  ({100*rel:.2f}%)")
        assert rel < 0.02, f"R={R}: L_int {L:.3e} vs mu0/8pi {L_cf:.3e} ({100*rel:.1f}%)"
    # radius-independence
    vals = list(Ls.values())
    assert abs(vals[0] - vals[1]) / L_cf < 0.01, "L_int must be independent of wire radius"


if __name__ == "__main__":
    validate_internal_inductance_round_wire()
    print("[OK] internal inductance mu0/(8 pi) validated (radius-independent).")
