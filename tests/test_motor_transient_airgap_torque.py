"""calc_motor_transient's air-gap torque, against a closed form.

The script used to inline the contour integral under the name
`_compute_torque_arkkio` and fed `grad(A_op)` straight into a `BND` integral on
the internal "airgap_mid" edge.  That evaluates a volume CoefficientFunction in
the wrong space and returned a near-zero, mesh-dependent torque -- measured
0.0116 -> 0.0027 -> 0.0007 against an exact 28.125 N m as the gap thinned.  The
unused `BoundaryFromVolumeCF` import the function carried was the missing piece.

The route now goes through radia.force_ngsolve, which takes the volume trace, so
`T_em` changes: it used to be roughly three orders of magnitude too small.

Closed form: for A_z = f1(r) cos(phi) + f2(r) sin(phi) the cycle mean of
B_r B_phi is (f2 f1' - f1 f2')/2.  With f1 = r, f2 = r^2 that is r/2, so

    T = (L / mu) * R * (R/2) * 2 pi R = (L / mu) * pi * R^3 .
"""

import math

import pytest

ngsolve = pytest.importorskip("ngsolve")
occ = pytest.importorskip("netgen.occ")

from ngsolve import (BND, CoefficientFunction, GridFunction, H1,  # noqa: E402
                     Integrate, Mesh, TaskManager, grad, sqrt, x, y)

from radia.panels.calc_motor_transient import (MU_0,  # noqa: E402
                                               _compute_airgap_torque)

STACK = 0.09
R_MID = 0.05
EXACT = (STACK / MU_0) * math.pi * R_MID**3          # 28.125 N m


def _annulus_mesh(inner=0.0475, outer=0.0525, mid=R_MID, maxh=0.0012):
    disc_outer = occ.Circle(occ.Pnt(0, 0), outer).Face()
    disc_inner = occ.Circle(occ.Pnt(0, 0), inner).Face()
    disc_mid = occ.Circle(occ.Pnt(0, 0), mid).Face()
    lower = disc_mid - disc_inner
    upper = disc_outer - disc_mid
    lower.faces.name = "airgap"
    upper.faces.name = "airgap"
    glued = occ.Glue([lower, upper])
    for edge in glued.edges:
        edge.name = "other"
    glued.edges.Nearest(occ.Pnt(mid, 0, 0)).name = "airgap_mid"
    return Mesh(occ.OCCGeometry(glued, dim=2).GenerateMesh(maxh=maxh))


def _vector_potential(mesh, order=4):
    """A_z = r cos(phi) + r^2 sin(phi) = x + y*r, whose torque is (L/mu) pi R^3."""
    gf = GridFunction(H1(mesh, order=order))
    gf.Set(x + y * sqrt(x * x + y * y))
    return gf


def test_contour_route_matches_the_closed_form():
    mesh = _annulus_mesh()
    with TaskManager():
        a_op = _vector_potential(mesh)
        torque = _compute_airgap_torque(mesh, a_op, R_MID, STACK, method="line")
    assert torque == pytest.approx(EXACT, rel=2e-3), (torque, EXACT)


def test_raw_gradient_on_the_internal_contour_is_the_defect_being_guarded():
    """Lock the reason the volume trace is required.

    This recomputes what calc_motor_transient did before delegating.  If a
    future edit drops BoundaryFromVolumeCF, the production route silently
    returns this number instead -- three orders of magnitude too small.
    """
    mesh = _annulus_mesh()
    with TaskManager():
        a_op = _vector_potential(mesh)
        gA = grad(a_op)
        bx, by = gA[1], -gA[0]
        radius = sqrt(x * x + y * y)
        cos_phi, sin_phi = x / radius, y / radius
        b_r = bx * cos_phi + by * sin_phi
        b_phi = -bx * sin_phi + by * cos_phi
        value = Integrate(radius * b_r * b_phi, mesh, BND,
                          definedon=mesh.Boundaries("airgap_mid"))
        legacy = (STACK / MU_0) * float(getattr(value, "real", value))
    assert abs(legacy) < 0.01 * EXACT, (
        "the raw-gradient contour integral is expected to collapse; if it no "
        f"longer does, NGSolve changed and this guard needs revisiting: {legacy}"
    )


def test_arkkio_route_agrees_with_the_contour():
    inner, outer = 0.0475, 0.0525
    mesh = _annulus_mesh(inner, outer, R_MID)
    with TaskManager():
        a_op = _vector_potential(mesh)
        line = _compute_airgap_torque(mesh, a_op, R_MID, STACK, method="line")
        arkkio = _compute_airgap_torque(
            mesh, a_op, R_MID, STACK, method="arkkio",
            r_inner=inner, r_outer=outer)
    assert line == pytest.approx(EXACT, rel=2e-3)
    # Arkkio weights by r^2 across the gap, so it sits slightly above the
    # mid-radius contour by the known (r_out^3 - r_in^3)/(3 dr R^2) factor.
    weighted = (outer**3 - inner**3) / (3.0 * (outer - inner) * R_MID**2)
    assert arkkio == pytest.approx(EXACT * weighted, rel=3e-3), (arkkio, line)


def test_missing_contour_still_returns_zero():
    """Historical behaviour: no 'airgap_mid' sideset means no torque handle."""
    disc = occ.Circle(occ.Pnt(0, 0), 0.05).Face()
    disc.faces.name = "airgap"
    for edge in disc.edges:
        edge.name = "other"
    mesh = Mesh(occ.OCCGeometry(disc, dim=2).GenerateMesh(maxh=0.01))
    with TaskManager():
        a_op = _vector_potential(mesh)
        assert _compute_airgap_torque(mesh, a_op, R_MID, STACK) == 0.0


def test_arkkio_route_reports_what_it_needs():
    mesh = _annulus_mesh()
    with TaskManager():
        a_op = _vector_potential(mesh)
        with pytest.raises(ValueError, match="r-airgap-inner"):
            _compute_airgap_torque(mesh, a_op, R_MID, STACK, method="arkkio")
        with pytest.raises(ValueError, match="airgap-torque arkkio needs"):
            _compute_airgap_torque(
                mesh, a_op, R_MID, STACK, method="arkkio",
                airgap_region="no_such_region",
                r_inner=0.0475, r_outer=0.0525)
        with pytest.raises(ValueError, match="line.*arkkio"):
            _compute_airgap_torque(mesh, a_op, R_MID, STACK, method="nope")
