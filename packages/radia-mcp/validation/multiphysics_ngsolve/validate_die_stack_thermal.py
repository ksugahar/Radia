r"""Multilayer thermal resistance / die-stack junction temperature -- regression test (#50).

Heat-generating top layer -> conduction down through passive layers (series thermal R=L/k) ->
cooled base. Drop across passive layer i = flux*L_i/k_i; junction rise = q L1^2/(2k1) +
Q*(L2/k2+L3/k3), Q=q L1. radia's solve_heat_steady (region-wise k) reproduces the interface
temperatures; the closed-form thermal network is the tool-independent gate. The
electronics-cooling junction-to-ambient block."""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.multiphysics import thermal_resistance_series

WD = 0.02
L1, L2, L3 = 0.002, 0.002, 0.002
K1, K2, K3 = 200.0, 10.0, 400.0
Q_VOL = 1.0e7
QF = Q_VOL * L1


def validate_thermal_resistance_series():
    R = thermal_resistance_series([L2, L3], [K2, K3])
    assert math.isclose(R, L2 / K2 + L3 / K3, rel_tol=1e-12)
    # the low-k layer (solder, k2=10) dominates over the high-k spreader (k3=400)
    assert (L2 / K2) > 0.9 * R
    # series adds; a thicker low-k layer raises R
    assert thermal_resistance_series([2 * L2, L3], [K2, K3]) > R


def validate_die_stack_fe():
    pytest.importorskip("ngsolve")
    pytest.importorskip("netgen")
    from ngsolve import Mesh, CoefficientFunction
    from netgen.occ import OCCGeometry, MoveTo, Glue
    from radia_mcp.radia_ngsolve.multiphysics import solve_heat_steady

    r3 = MoveTo(-WD/2, 0).Rectangle(WD, L3).Face(); r3.faces.name = "lay3"
    r2 = MoveTo(-WD/2, L3).Rectangle(WD, L2).Face(); r2.faces.name = "lay2"
    r1 = MoveTo(-WD/2, L3+L2).Rectangle(WD, L1).Face(); r1.faces.name = "lay1"
    stack = Glue([r3, r2, r1])
    for e in stack.edges:
        if abs(e.center[1]) < 1e-9:
            e.name = "base"
    mesh = Mesh(OCCGeometry(stack, dim=2).GenerateMesh(maxh=L3/6))
    kcf = mesh.MaterialCF({"lay1": K1, "lay2": K2, "lay3": K3}, default=K1)
    qcf = mesh.MaterialCF({"lay1": Q_VOL}, default=0.0)
    gT = solve_heat_steady(mesh, qcf, kcf, conductor="lay1|lay2|lay3", dirichlet="base", order=3)

    # interface rises = flux * cumulative series resistance below
    dT_32 = QF * (L3 / K3)
    dT_21 = QF * (L3 / K3 + L2 / K2)
    dT_j = dT_21 + Q_VOL * L1 * L1 / (2 * K1)
    assert abs(gT(mesh(0.0, L3)) - dT_32) / dT_32 < 5e-3
    assert abs(gT(mesh(0.0, L3 + L2)) - dT_21) / dT_21 < 5e-3
    assert abs(gT(mesh(0.0, (L3+L2+L1) * 0.999)) - dT_j) / dT_j < 5e-3
    # monotone increasing from base to junction
    assert gT(mesh(0.0, L3)) < gT(mesh(0.0, L3 + L2)) < gT(mesh(0.0, (L3+L2+L1)*0.999))


if __name__ == "__main__":
    validate_thermal_resistance_series()
    validate_die_stack_fe()
    print("[OK] die-stack thermal resistance network validated (series R, junction temperature).")
