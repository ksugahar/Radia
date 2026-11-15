"""Discrete differential-identity checks for the axifem production space.

axifem is not a full H1/HCurl/HDiv/L2 de Rham complex.  It is a Henrotte
axisymmetric lane: the unknown u is A_phi with the V-DOF basis
N_i = r_i psi_i / r, and ``grad(u)`` is the physical meridian gradient entering

    B_r = -d u / d z,      B_z = d u / d r + u / r.

Guards: the fields whose flux function psi = r A_phi lies in every element's
scalar space are reproduced exactly by interpolation, value, gradient and B:
A = 1/r (psi = 1), A = z/r (psi = z) and, except for P1, A = r (psi = s, a
uniform B_z = 2).  P1 is linear in the reference coordinates of the quadratic
(s, z) map, so psi = s is not in its space (stated, not hidden).
Assumes: straight off-axis elements; curved elements are refused.
"""

import math

import pytest

ng = pytest.importorskip("ngsolve")

from _vol_mesh import reload_via_vol, structured_rect_vol_mesh
from netgen.geom2d import SplineGeometry
from netgen.occ import MoveTo, OCCGeometry, X, Y
from ngsolve import GridFunction, Integrate, Mesh, grad, x, y


@pytest.fixture(scope="module", autouse=True)
def _taskmanager():
    with ng.TaskManager():
        yield


def _identity_cases(with_uniform_b=True):
    zero = 0 * x
    cases = [
        ("inv_r", 1 / x, -1 / (x * x), zero, zero, zero),
        ("z_over_r", y / x, -y / (x * x), 1 / x, -1 / x, zero),
    ]
    if with_uniform_b:
        cases.append(("r", x, 1 + zero, zero, zero, 2 + zero))
    return cases


def _assert_axisym_field_identities(mesh, order, *, curvedquad=False, tol=1e-18,
                                    with_uniform_b=True):
    axifem = pytest.importorskip("radia.axifem")
    fes = axifem.H1Henrotte(mesh, order=order, curvedquad=curvedquad)
    gfu = GridFunction(fes)

    for label, expr, du_dr, du_dz, br_exact, bz_exact in _identity_cases(with_uniform_b):
        gfu.Set(expr)
        interp_err = Integrate((gfu - expr) * (gfu - expr), mesh)
        grad_err = Integrate(
            (grad(gfu)[0] - du_dr) * (grad(gfu)[0] - du_dr)
            + (grad(gfu)[1] - du_dz) * (grad(gfu)[1] - du_dz),
            mesh,
        )
        field_err = Integrate(
            (-grad(gfu)[1] - br_exact) * (-grad(gfu)[1] - br_exact)
            + (grad(gfu)[0] + gfu / x - bz_exact) * (grad(gfu)[0] + gfu / x - bz_exact),
            mesh,
        )
        assert interp_err < tol, f"{label}: interpolation identity drifted ({interp_err:.3e})"
        assert grad_err < tol, f"{label}: grad identity drifted ({grad_err:.3e})"
        assert field_err < tol, f"{label}: axisymmetric B identity drifted ({field_err:.3e})"


def test_q1_q2_axis_aligned_field_identities():
    mesh = structured_rect_vol_mesh(
        0.6,
        1.7,
        -0.4,
        0.4,
        quads=True,
        nx=2,
        ny=2,
        mapping=lambda xi, eta: (0.6 + 1.1 * xi, -0.4 + 0.8 * eta),
        stem="axifem_identity_quad",
    )
    _assert_axisym_field_identities(mesh, order=1)
    _assert_axisym_field_identities(mesh, order=2)


def test_p1_p2_triangle_field_identities():
    mesh = structured_rect_vol_mesh(
        0.6,
        1.7,
        -0.4,
        0.4,
        quads=False,
        nx=2,
        ny=2,
        mapping=lambda xi, eta: (0.6 + 1.1 * xi, -0.4 + 0.8 * eta),
        stem="axifem_identity_tri",
    )
    _assert_axisym_field_identities(mesh, order=1, with_uniform_b=False)
    _assert_axisym_field_identities(mesh, order=2)


@pytest.mark.parametrize("order", [1, 2])
def test_curved_triangle_mesh_is_refused(order):
    axifem = pytest.importorskip("radia.axifem")
    geo = SplineGeometry()
    geo.AddCircle((1.2, 0.0), 0.3, leftdomain=1, rightdomain=0, bc="outer")
    geo.SetMaterial(1, "conductor")
    mesh = Mesh(geo.GenerateMesh(maxh=0.15))
    mesh.Curve(2)
    gfu = GridFunction(axifem.H1Henrotte(mesh, order=order))
    with pytest.raises(Exception, match="curved element"):
        gfu.Set(1 / x)


def test_q2_general_quad_field_identities():
    def annulus(xi, eta):
        rho = 0.8 + 0.5 * xi
        phi = -0.3 + 0.6 * eta
        return (rho * math.cos(phi), rho * math.sin(phi))

    mesh = structured_rect_vol_mesh(
        0.8,
        1.3,
        -0.3,
        0.3,
        quads=True,
        nx=2,
        ny=4,
        mapping=annulus,
        stem="axifem_identity_q2_curved",
    )
    _assert_axisym_field_identities(mesh, order=2, curvedquad=True, tol=1e-16)


def test_axis_aligned_boundary_value_trace_matches_analytic_edge_integrals():
    axifem = pytest.importorskip("radia.axifem")
    box = MoveTo(1.0, 0.0).Rectangle(1.0, 1.0).Face()
    box.edges.Min(X).name = "left"
    box.edges.Max(X).name = "right"
    box.edges.Min(Y).name = "bottom"
    box.edges.Max(Y).name = "top"
    mesh = Mesh(OCCGeometry(box, dim=2).GenerateMesh(maxh=5.0, quad_dominated=True))
    mesh = reload_via_vol(mesh, "axifem_boundary_trace")

    fes = axifem.H1Henrotte(mesh, order=1)
    gfu = GridFunction(fes)
    gfu.Set(x + 2 * y / x)   # psi = s + 2 z is in the Q1 space

    top = Integrate(gfu, mesh, definedon=mesh.Boundaries("top"), order=24)
    right = Integrate(gfu, mesh, definedon=mesh.Boundaries("right"), order=24)
    top_exact = 1.5 + 2.0 * math.log(2.0)   # INT_1^2 (r + 2 / r) dr at z = 1
    right_exact = 2.5                        # INT_0^1 (2 + z) dz at r = 2
    assert abs(top - top_exact) < 1e-12, (top, top_exact)
    assert abs(right - right_exact) < 1e-12, (right, right_exact)
