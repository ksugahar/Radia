r"""phi-theta (voltage-temperature) relation / constriction supertemperature -- test (#56).

A current-carrying conductor with both terminals at T0 and voltage U has a GEOMETRY-INDEPENDENT
hot spot: T(V)^2 = T0^2 + V(U-V)/L, T_max^2 - T0^2 = U^2/(4L) (L = Lorenz number). The
Wiedemann-Franz-exact generalisation of the constant-property Joule bar (#19, sigma V^2/8k),
linearised by the Kirchhoff transform. Closed-form helpers (tool-independent) + an electro-thermal
FE on two different shapes (the geometry independence)."""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.multiphysics import (phi_theta_local_temperature,
                                                  phi_theta_max_temperature, melting_voltage, LORENZ)

T0 = 300.0


def validate_phi_theta_closed_form():
    U = 0.1
    # endpoints sit at T0; the relation is symmetric about V = U/2
    assert math.isclose(phi_theta_local_temperature(0.0, U, T0), T0, rel_tol=1e-12)
    assert math.isclose(phi_theta_local_temperature(U, U, T0), T0, rel_tol=1e-12)
    assert math.isclose(phi_theta_local_temperature(0.3 * U, U, T0),
                        phi_theta_local_temperature(0.7 * U, U, T0), rel_tol=1e-12)
    # max at V = U/2: T_max^2 - T0^2 = U^2/(4L)
    Tmax = phi_theta_max_temperature(U, T0)
    assert math.isclose(Tmax**2 - T0**2, U**2 / (4 * LORENZ), rel_tol=1e-12)
    assert phi_theta_local_temperature(0.5 * U, U, T0) >= phi_theta_local_temperature(0.4 * U, U, T0)
    # inverse round-trip: melting_voltage(T_max) recovers U
    assert math.isclose(melting_voltage(Tmax, T0), U, rel_tol=1e-12)
    # classic copper contact voltages (geometry-independent): soften ~0.11 V, melt ~0.43 V
    assert math.isclose(melting_voltage(463.0, T0), 0.110, abs_tol=3e-3)
    assert math.isclose(melting_voltage(1356.0, T0), 0.43, abs_tol=2e-2)


def _fe_Tmax(mesh, sigma, U, order=3):
    from ngsolve import CoefficientFunction, grad, InnerProduct, sqrt
    from radia_mcp.radia_ngsolve.scalar_fem2d import solve_poisson_2d
    psi0 = LORENZ * sigma * T0 * T0 / 2.0
    V = solve_poisson_2d(mesh, CoefficientFunction(sigma), {"hi": U, "lo": 0.0}, order=order)
    q = sigma * InnerProduct(grad(V), grad(V))
    psi = solve_poisson_2d(mesh, CoefficientFunction(1.0), {"hi": psi0, "lo": psi0}, source=q, order=order)
    Tcf = sqrt(2.0 * psi / (LORENZ * sigma))
    Lx, W = 0.02, 0.006
    Tmax = T0
    for i in range(61):
        for j in range(31):
            try:
                Tmax = max(Tmax, Tcf(mesh(i / 60.0 * Lx, (-0.5 + j / 30.0) * W)))
            except Exception:
                pass
    return Tmax


def validate_phi_theta_fe_geometry_independent():
    pytest.importorskip("ngsolve")
    pytest.importorskip("netgen")
    from ngsolve import Mesh, TaskManager
    from netgen.occ import OCCGeometry, MoveTo, Glue, X

    sigma, U = 5.8e7, 0.1
    Lx, W = 0.02, 0.006

    def bar():
        f = MoveTo(0, -W / 2).Rectangle(Lx, W).Face()
        f.edges.Min(X).name = "hi"; f.edges.Max(X).name = "lo"
        return Mesh(OCCGeometry(f, dim=2).GenerateMesh(maxh=W / 6))

    def dogbone():
        nw, nh = 0.006, W / 3
        left = MoveTo(0, -W / 2).Rectangle((Lx - nw) / 2, W).Face()
        right = MoveTo((Lx + nw) / 2, -W / 2).Rectangle((Lx - nw) / 2, W).Face()
        neck = MoveTo((Lx - nw) / 2, -nh / 2).Rectangle(nw, nh).Face()
        s = Glue([left, neck, right])
        s.edges.Min(X).name = "hi"; s.edges.Max(X).name = "lo"
        return Mesh(OCCGeometry(s, dim=2).GenerateMesh(maxh=W / 8))

    Tcf = phi_theta_max_temperature(U, T0)
    with TaskManager():
        Tbar = _fe_Tmax(bar(), sigma, U)
        Tdog = _fe_Tmax(dogbone(), sigma, U)
    # each shape matches the closed form ...
    assert math.isclose(Tbar, Tcf, rel_tol=2e-3)
    assert math.isclose(Tdog, Tcf, rel_tol=2e-3)
    # ... and the two shapes agree (the geometry independence)
    assert math.isclose(Tbar, Tdog, rel_tol=2e-3)


if __name__ == "__main__":
    validate_phi_theta_closed_form()
    validate_phi_theta_fe_geometry_independent()
    print("[OK] phi-theta: electro-thermal FE == U^2/4L supertemperature, geometry-independent.")
