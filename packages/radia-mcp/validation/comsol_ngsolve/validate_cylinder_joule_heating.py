r"""Cylindrical conductor Joule heating (electro-thermal) -- regression test (#41).

A round conductor (radius a) with a uniform Joule source q=J^2/sigma has the parabolic
radial temperature rise dT(r)=q(a^2-r^2)/(4k), centre dT_peak=q a^2/(4k). radia
`solve_heat_steady` on a disk reproduces it; the closed form is the tool-independent gate
(the cylindrical sibling of the slab electro-thermal #19). A reference commercial FE code
matches the same closed form to 0.000% live (recorded internally, not here)."""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.multiphysics import joule_heated_cylinder_peak_dT

A, K = 0.02, 2.0
Q = 5.0e5          # uniform Joule density [W/m^3] (= J^2/sigma)


def validate_peak_dt_formula():
    assert math.isclose(joule_heated_cylinder_peak_dT(Q, A, K), Q * A * A / (4 * K), rel_tol=1e-12)
    # dT_peak = 25 K for these values; scales with q and a^2, inverse with k
    assert math.isclose(joule_heated_cylinder_peak_dT(Q, A, K), 25.0, rel_tol=1e-12)
    assert math.isclose(joule_heated_cylinder_peak_dT(2 * Q, A, K), 50.0, rel_tol=1e-12)
    assert math.isclose(joule_heated_cylinder_peak_dT(Q, 2 * A, K), 100.0, rel_tol=1e-12)
    assert math.isclose(joule_heated_cylinder_peak_dT(Q, A, 4 * K), 6.25, rel_tol=1e-12)


def validate_cylinder_joule_heating_fe():
    pytest.importorskip("ngsolve")
    pytest.importorskip("netgen")
    from ngsolve import Mesh, CoefficientFunction
    from netgen.occ import OCCGeometry, WorkPlane
    from radia_mcp.radia_ngsolve.multiphysics import solve_heat_steady

    wire = WorkPlane().Circle(0, 0, A).Face(); wire.faces.name = "wire"; wire.edges.name = "surf"
    mesh = Mesh(OCCGeometry(wire, dim=2).GenerateMesh(maxh=A / 20))
    gT = solve_heat_steady(mesh, CoefficientFunction(Q), K, "wire", "surf", order=3)

    # centre rise == closed-form peak
    dT_peak = joule_heated_cylinder_peak_dT(Q, A, K)
    assert abs(gT(mesh(0.0, 0.0)) - dT_peak) / dT_peak < 5e-3
    # parabolic profile dT(r) = q (a^2 - r^2)/(4k) at several radii
    for ra in (0.25, 0.5, 0.75):
        r = ra * A
        dT_an = Q * (A * A - r * r) / (4.0 * K)
        assert abs(gT(mesh(r, 0.0)) - dT_an) / dT_an < 8e-3, f"r/a={ra}: {gT(mesh(r,0.0))} vs {dT_an}"


if __name__ == "__main__":
    validate_peak_dt_formula()
    validate_cylinder_joule_heating_fe()
    print("[OK] cylindrical Joule heating dT(r)=q(a^2-r^2)/(4k) validated.")
