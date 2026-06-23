r"""Convective electro-thermal (Joule heat + Newton cooling) -- regression test (#44).

A slab with a uniform Joule source q, convectively cooled (film h) on both faces, has the
centre rise dT = qL^2/(8k) + qL/(2h) = conduction parabola + convective film. radia's Robin
heat solver reproduces it; the closed form is the tool-independent gate. A reference
commercial FE code matches the same closed form to 0.000% live (recorded internally). As
h -> inf the film vanishes and the fixed-T parabola (#19/#41) is recovered."""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.multiphysics import convective_slab_peak_dT

L, WD = 0.02, 0.10
Q, K = 5.0e5, 2.0


def validate_peak_dt_formula():
    # dT = conduction qL^2/8k + film qL/2h
    assert math.isclose(convective_slab_peak_dT(Q, L, K, 500.0),
                        Q*L*L/(8*K) + Q*L/(2*500.0), rel_tol=1e-12)
    # for these values: conduction 12.5 K + film 10 K = 22.5 K
    assert math.isclose(convective_slab_peak_dT(Q, L, K, 500.0), 22.5, rel_tol=1e-12)
    # h -> inf : film vanishes -> the fixed-T parabola qL^2/8k = 12.5 K
    assert math.isclose(convective_slab_peak_dT(Q, L, K, 1e12), Q*L*L/(8*K), rel_tol=1e-6)
    # smaller h (worse cooling) -> larger film -> hotter
    assert convective_slab_peak_dT(Q, L, K, 250.0) > convective_slab_peak_dT(Q, L, K, 500.0)


def validate_convective_slab_fe():
    pytest.importorskip("ngsolve")
    pytest.importorskip("netgen")
    from ngsolve import Mesh, CoefficientFunction
    from netgen.occ import OCCGeometry, WorkPlane, X, Y
    from radia_mcp.radia_ngsolve.multiphysics import solve_heat_steady_robin

    slab = WorkPlane().MoveTo(-WD/2, -L/2).Rectangle(WD, L).Face(); slab.faces.name = "slab"
    slab.edges.Max(Y).name = "conv"; slab.edges.Min(Y).name = "conv"
    slab.edges.Max(X).name = "side"; slab.edges.Min(X).name = "side"
    mesh = Mesh(OCCGeometry(slab, dim=2).GenerateMesh(maxh=L/20))

    for h in (250.0, 500.0, 2000.0):
        gT = solve_heat_steady_robin(mesh, CoefficientFunction(Q), K, "slab", "conv", h, order=3)
        dT_an = convective_slab_peak_dT(Q, L, K, h)
        assert abs(gT(mesh(0.0, 0.0)) - dT_an) / dT_an < 5e-3, f"h={h}: centre"
        # surface floats q L/(2h) above coolant
        dT_surf = Q * L / (2.0 * h)
        assert abs(gT(mesh(0.0, L/2 * 0.999)) - dT_surf) / dT_surf < 1e-2, f"h={h}: surface"


if __name__ == "__main__":
    validate_peak_dt_formula()
    validate_convective_slab_fe()
    print("[OK] convective electro-thermal dT = qL^2/8k + qL/2h validated.")
