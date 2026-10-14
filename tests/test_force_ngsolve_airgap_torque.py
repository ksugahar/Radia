"""Air-gap Maxwell torque: exact analytic value, and line vs Arkkio agreement.

For a planar field prescribed on an annulus as

    B_r   = B0 cos(phi)
    B_phi = C  cos(phi - delta)

the single-contour torque at radius R is exactly

    T = (L / mu) * R^2 * B0 * C * pi * cos(delta)

and Arkkio's thickness-average over r_in < r < r_out is the same integrand
weighted by r^2, so the two agree as the annulus thins about R.  Arkkio is NOT
the contour formula -- this file is what keeps the two from being confused
again, after `_compute_torque_arkkio` shipped for years as a contour integral
under Arkkio's name.
"""

import math

import numpy as np
import pytest

ngsolve = pytest.importorskip("ngsolve")
occ = pytest.importorskip("netgen.occ")

from ngsolve import CoefficientFunction, Mesh, TaskManager, sqrt, x, y  # noqa: E402

from radia.force_ngsolve import (  # noqa: E402
    MU0,
    air_gap_maxwell_torque_2d,
    air_gap_maxwell_torque_arkkio_2d,
    air_gap_maxwell_torque_line_2d,
)

B0 = 0.8
C_AMP = 0.35
DELTA = 0.4
STACK = 0.12


def _annulus_mesh(inner, outer, mid, maxh=0.004):
    disc_outer = occ.Circle(occ.Pnt(0, 0), outer).Face()
    disc_inner = occ.Circle(occ.Pnt(0, 0), inner).Face()
    annulus = disc_outer - disc_inner
    annulus.faces.name = "airgap"
    # a meshed contour at the mid radius so the line route has a boundary
    lower = (occ.Circle(occ.Pnt(0, 0), mid).Face() - disc_inner)
    upper = (disc_outer - occ.Circle(occ.Pnt(0, 0), mid).Face())
    lower.faces.name = "airgap"
    upper.faces.name = "airgap"
    glued = occ.Glue([lower, upper])
    for edge in glued.edges:
        edge.name = "other"
    glued.edges.Nearest(occ.Pnt(mid, 0, 0)).name = "airgap_mid"
    return Mesh(occ.OCCGeometry(glued, dim=2).GenerateMesh(maxh=maxh))


def _prescribed_field():
    """Cartesian components of B_r = B0 cos(phi), B_phi = C cos(phi - delta)."""
    radius = sqrt(x * x + y * y)
    cos_phi = x / radius
    sin_phi = y / radius
    radial = B0 * cos_phi
    tangential = C_AMP * (cos_phi * math.cos(DELTA) + sin_phi * math.sin(DELTA))
    return CoefficientFunction(
        (radial * cos_phi - tangential * sin_phi,
         radial * sin_phi + tangential * cos_phi)
    )


def _exact_line_torque(radius):
    return (STACK / MU0) * radius**2 * B0 * C_AMP * math.pi * math.cos(DELTA)


def test_line_torque_matches_the_closed_form():
    mid = 0.05
    mesh = _annulus_mesh(0.048, 0.052, mid)
    field = _prescribed_field()
    with TaskManager():
        torque = air_gap_maxwell_torque_line_2d(
            field, mesh, mesh.Boundaries("airgap_mid"), stack_length_m=STACK
        )
    assert torque == pytest.approx(_exact_line_torque(mid), rel=2e-3)


def test_arkkio_average_matches_its_own_closed_form():
    inner, outer, mid = 0.045, 0.055, 0.050
    mesh = _annulus_mesh(inner, outer, mid, maxh=0.003)
    field = _prescribed_field()
    with TaskManager():
        torque = air_gap_maxwell_torque_arkkio_2d(
            field, mesh, mesh.Materials("airgap"),
            inner_radius_m=inner, outer_radius_m=outer, stack_length_m=STACK,
        )
    # int r * (B_r B_phi) * r dr dphi / (r_out - r_in)
    radial_factor = (outer**3 - inner**3) / (3.0 * (outer - inner))
    exact = (STACK / MU0) * radial_factor * B0 * C_AMP * math.pi * math.cos(DELTA)
    assert torque == pytest.approx(exact, rel=2e-3)


def test_arkkio_converges_to_the_contour_value_as_the_gap_thins():
    mid = 0.05
    field = _prescribed_field()
    errors = []
    for half_span in (0.006, 0.003, 0.0015):
        inner, outer = mid - half_span, mid + half_span
        mesh = _annulus_mesh(inner, outer, mid, maxh=max(half_span / 2.0, 0.0015))
        with TaskManager():
            arkkio = air_gap_maxwell_torque_arkkio_2d(
                field, mesh, mesh.Materials("airgap"),
                inner_radius_m=inner, outer_radius_m=outer, stack_length_m=STACK,
            )
        errors.append(abs(arkkio / _exact_line_torque(mid) - 1.0))
    assert errors[-1] < errors[0], f"no convergence: {errors}"
    assert errors[-1] < 5e-3, errors


def test_selector_dispatches_and_rejects_mismatched_arguments():
    inner, outer, mid = 0.045, 0.055, 0.050
    mesh = _annulus_mesh(inner, outer, mid, maxh=0.003)
    field = _prescribed_field()
    with TaskManager():
        by_selector = air_gap_maxwell_torque_2d(
            field, mesh, mesh.Boundaries("airgap_mid"),
            method="line", stack_length_m=STACK,
        )
        direct = air_gap_maxwell_torque_line_2d(
            field, mesh, mesh.Boundaries("airgap_mid"), stack_length_m=STACK
        )
    assert by_selector == pytest.approx(direct, rel=1e-12)

    with pytest.raises(ValueError):
        air_gap_maxwell_torque_2d(
            field, mesh, mesh.Materials("airgap"), method="arkkio"
        )
    with pytest.raises(ValueError):
        air_gap_maxwell_torque_2d(
            field, mesh, mesh.Boundaries("airgap_mid"),
            method="line", inner_radius_m=inner, outer_radius_m=outer,
        )
    with pytest.raises(ValueError):
        air_gap_maxwell_torque_2d(field, mesh, mesh.Materials("airgap"), method="nope")


def test_nonpositive_inputs_fail_loudly():
    mesh = _annulus_mesh(0.048, 0.052, 0.05)
    field = _prescribed_field()
    with pytest.raises(ValueError):
        air_gap_maxwell_torque_line_2d(
            field, mesh, mesh.Boundaries("airgap_mid"), stack_length_m=0.0
        )
    with pytest.raises(ValueError):
        air_gap_maxwell_torque_arkkio_2d(
            field, mesh, mesh.Materials("airgap"),
            inner_radius_m=0.052, outer_radius_m=0.048,
        )
