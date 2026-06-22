r"""Slot-leakage inductance of a machine slot -- regression test (#54).

The leakage flux crossing a slot from tooth to tooth: specific permeance
lambda_s = h_c/(3 w_s) + h_o/w_o, L = mu0 N^2 l_stk lambda_s. The 1/3 conductor term (linear
cross-slot MMF buildup) and a same-width current-free extension are EXACT; a narrowed tooth-tip
opening adds shoulder fringing the simple series form omits. Closed-form helper (tool-independent)
+ an FE check on the well-posed 1-D-equivalent BVP (top Dirichlet yoke / Neumann tooth walls)."""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.solve import slot_leakage_permeance, slot_leakage_inductance, MU0


def test_slot_leakage_closed_form():
    # closed slot: lambda = h_c/(3 w_s), the exact 1/3 factor
    assert math.isclose(slot_leakage_permeance(0.020, 0.004), 0.020 / (3 * 0.004), rel_tol=1e-12)
    # opening adds a series term h_o/w_o
    lam = slot_leakage_permeance(0.024, 0.006, 0.004, 0.002)
    assert math.isclose(lam, 0.024 / (3 * 0.006) + 0.004 / 0.002, rel_tol=1e-12)
    # default opening width = slot width
    assert math.isclose(slot_leakage_permeance(0.02, 0.004, 0.01),
                        slot_leakage_permeance(0.02, 0.004, 0.01, 0.004), rel_tol=1e-12)
    # scaling: lambda ~ h_c, ~ 1/w_s
    assert math.isclose(slot_leakage_permeance(0.04, 0.004), 2 * slot_leakage_permeance(0.02, 0.004))
    assert math.isclose(slot_leakage_permeance(0.02, 0.008), 0.5 * slot_leakage_permeance(0.02, 0.004))
    # inductance: mu0 N^2 l_stk lambda, with N^2 and l_stk scaling
    L = slot_leakage_inductance(0.020, 0.004, axial_length=0.10, turns=12)
    assert math.isclose(L, MU0 * 12**2 * 0.10 * (0.020 / (3 * 0.004)), rel_tol=1e-12)
    assert math.isclose(slot_leakage_inductance(0.02, 0.004, 0.1, 24),
                        4 * slot_leakage_inductance(0.02, 0.004, 0.1, 12), rel_tol=1e-12)


def _slot_L_fe(w_s, h_c, h_o=0.0, w_o=None, order=3):
    from ngsolve import Mesh, CoefficientFunction, grad
    from netgen.occ import OCCGeometry, MoveTo, Glue, Y
    from radia_mcp.radia_ngsolve.solve import solve_planar_magnetostatic, NU0
    from radia_mcp.radia_ngsolve.force import magnetic_energy_2d
    cond = MoveTo(0, 0).Rectangle(w_s, h_c).Face(); cond.faces.name = "cond"
    if h_o > 0:
        w_o = w_o or w_s
        op = MoveTo((w_s - w_o) / 2, h_c).Rectangle(w_o, h_o).Face(); op.faces.name = "air"
        shape = Glue([cond, op])
    else:
        shape = cond
    shape.edges.Max(Y).name = "top"
    mesh = Mesh(OCCGeometry(shape, dim=2).GenerateMesh(maxh=min(w_s, w_o or w_s) / 6))
    Jz = CoefficientFunction([1.0 / (w_s * h_c) if m == "cond" else 0.0 for m in mesh.GetMaterials()])
    A = solve_planar_magnetostatic(mesh, CoefficientFunction(NU0), Jz=Jz, order=order, dirichlet="top")
    B = CoefficientFunction((grad(A)[1], -grad(A)[0]))
    return 2 * sum(magnetic_energy_2d(B, mesh, m) for m in set(mesh.GetMaterials()))
