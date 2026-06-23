r"""2D biaxial thermal stress of a fully-constrained plate (plane stress) -- regression test.

A plate fully constrained in-plane (u=0 all edges) under a uniform rise dT -> biaxial stress
sigma_x=sigma_y=-E alpha dT/(1-nu) (plane stress, biaxial_thermal_stress); u=0 is exact so the
FE is exact. The 2D constraint factor 1/(1-nu) vs the 1D uniaxial -E alpha dT. Tool-independent.
(Plane STRAIN -> -E alpha dT/(1-2 nu); radia's plane='strain' thermal load is currently (1+nu)
too small -- flagged separately -- so this test covers plane STRESS only.)
"""
import os
import sys

import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen")

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from ngsolve import Mesh, TaskManager
from netgen.occ import OCCGeometry, MoveTo
from radia_mcp.radia_ngsolve.elasticity import (solve_linear_elasticity, stress_2d,
                                                biaxial_thermal_stress)

L, E, nu, alpha, dT = 0.05, 200e9, 0.30, 12e-6, 50.0


def validate_biaxial_thermal_stress_plane_stress():
    plate = MoveTo(0, 0).Rectangle(L, L).Face(); plate.faces.name = "plate"
    for e in plate.edges:
        e.name = "fix"
    mesh = Mesh(OCCGeometry(plate, dim=2).GenerateMesh(maxh=L/12))
    with TaskManager():
        u = solve_linear_elasticity(mesh, E, nu, dirichlet="fix",
                                    thermal=(alpha, dT), plane="stress", order=2)
    sig = stress_2d(u, E, nu, "stress", thermal=(alpha, dT))
    sxx = sig[0](mesh(L/2, L/2)); syy = sig[3](mesh(L/2, L/2))
    cf = biaxial_thermal_stress(E, alpha, dT, nu)
    print(f"[biaxial] sxx={sxx/1e6:.3f} MPa  syy={syy/1e6:.3f} MPa  analytic={cf/1e6:.3f} MPa")
    assert abs(sxx - cf) / abs(cf) < 0.01, f"sxx {sxx:.3e} vs -E a dT/(1-nu) {cf:.3e}"
    assert abs(syy - cf) / abs(cf) < 0.01, "biaxial: syy must equal sxx"
    # the 2D constraint amplifies the 1D uniaxial by 1/(1-nu)
    assert abs(abs(cf) - abs(-E*alpha*dT)/(1-nu)) < 1.0


if __name__ == "__main__":
    validate_biaxial_thermal_stress_plane_stress()
    print("[OK] biaxial (plane-stress) thermal stress validated.")
