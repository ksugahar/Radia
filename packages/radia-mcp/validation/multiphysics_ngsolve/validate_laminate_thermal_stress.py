r"""Composite (laminate) residual thermal stress / rule-of-mixtures CTE -- #38.

A FREE bonded sandwich B|A|B under uniform dT develops internal axial stress from CTE
mismatch alone: sigma_i = E_i (alpha_eff - alpha_i) dT, alpha_eff the stiffness-weighted
rule of mixtures, self-equilibrated (sum A_i sigma_i = 0). The radia elasticity solver
(region-wise E, alpha thermal eigenstrain) reproduces sigma_xx per layer and the common
central strain alpha_eff*dT. Tool-independent (closed form is the reference). The FREE-bar
counterpart of the externally blocked bar (#20) / clamped plate (#35); symmetric -> no bend.
"""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.elasticity import (rule_of_mixtures_cte,
                                                laminate_residual_stress)

EA, alA = 70e9, 23e-6
EB, alB = 200e9, 12e-6
NU, DT, L, A, B = 0.3, 100.0, 40e-3, 2e-3, 1e-3   # core width A, each face B (2B=A)


def validate_rule_of_mixtures_self_equilibrium():
    """Pure closed-form identities: alpha_eff bounded by the layer CTEs, and the
    free composite is self-equilibrated sum(A_i sigma_i) = 0 (no net axial force)."""
    a_eff = rule_of_mixtures_cte([EA, EB], [A, 2 * B], [alA, alB])
    assert alB < a_eff < alA                      # stiffness-weighted, between the two
    sA = laminate_residual_stress(EA, alA, a_eff, DT)
    sB = laminate_residual_stress(EB, alB, a_eff, DT)
    assert sA < 0.0 < sB                          # high-alpha core compressed, faces in tension
    assert abs(A * sA + 2 * B * sB) < 1e-6 * abs(A * sA)   # force balance


def validate_laminate_fe_matches_rule_of_mixtures():
    pytest.importorskip("ngsolve")
    pytest.importorskip("netgen")
    from ngsolve import Mesh, Integrate, IfPos, grad, x
    from netgen.occ import OCCGeometry, WorkPlane, Glue
    from radia_mcp.radia_ngsolve.elasticity import solve_linear_elasticity, stress_2d

    a_eff = rule_of_mixtures_cte([EA, EB], [A, 2 * B], [alA, alB])
    sA = laminate_residual_stress(EA, alA, a_eff, DT)
    sB = laminate_residual_stress(EB, alB, a_eff, DT)

    core = WorkPlane().MoveTo(0, -A / 2).Rectangle(L, A).Face(); core.faces.name = "core"
    ftop = WorkPlane().MoveTo(0, A / 2).Rectangle(L, B).Face(); ftop.faces.name = "face"
    fbot = WorkPlane().MoveTo(0, -A / 2 - B).Rectangle(L, B).Face(); fbot.faces.name = "face"
    strip = Glue([fbot, core, ftop])
    for e in strip.edges:
        if abs(e.center[0]) < 1e-9:
            e.name = "left"
    mesh = Mesh(OCCGeometry(strip, dim=2).GenerateMesh(maxh=B / 2))
    assert set(mesh.GetMaterials()) == {"core", "face"}

    E_cf = mesh.MaterialCF({"core": EA, "face": EB}, default=EA)
    al_cf = mesh.MaterialCF({"core": alA, "face": alB}, default=alA)
    gu = solve_linear_elasticity(mesh, E_cf, NU, dirichletx="left", dirichlety="left",
                                 thermal=(al_cf, DT), plane="stress", order=2)
    sxx = stress_2d(gu, E_cf, NU, plane="stress", thermal=(al_cf, DT))[0, 0]

    win = IfPos(x - 0.4 * L, IfPos(0.6 * L - x, 1.0, 0.0), 0.0)   # uniform central window
    for name, sc in (("core", sA), ("face", sB)):
        m = mesh.MaterialCF({name: 1.0}, default=0.0) * win
        smean = Integrate(sxx * m, mesh) / Integrate(m, mesh)
        assert abs(smean - sc) / abs(sc) < 5e-3, f"{name}: FE {smean:.4e} vs closed {sc:.4e}"

    # central common strain == alpha_eff*dT (rule-of-mixtures CTE)
    exx = Integrate(grad(gu)[0, 0] * win, mesh) / Integrate(win, mesh)
    assert abs(exx - a_eff * DT) / (a_eff * DT) < 1e-2


if __name__ == "__main__":
    validate_rule_of_mixtures_self_equilibrium()
    validate_laminate_fe_matches_rule_of_mixtures()
    print("[OK] laminate residual thermal stress = rule of mixtures (0.00%).")
