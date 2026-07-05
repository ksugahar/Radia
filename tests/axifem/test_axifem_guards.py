"""Fail-fast guards for axifem geometry contracts.

The closed-form Q1/Q2 quad elements are valid only on axis-aligned
rectangles in the (r, z) meridian plane.  General quads must use the Q2
curved path, and degenerate P1 triangles must fail before they can produce
silent garbage matrices.
"""

import pytest

pytest.importorskip("ngsolve")

from ngsolve import BilinearForm, CoefficientFunction as CF, TaskManager

from _vol_mesh import structured_rect_vol_mesh


def _assemble_stiffness(mesh, order, **flags):
    axifem = pytest.importorskip("radia.axifem")
    fes = axifem.H1Henrotte(mesh, order=order, **flags)
    a = BilinearForm(fes, symmetric=True)
    a += axifem.AxiHenrotteStiffnessBFI(CF(1.0))
    with TaskManager():
        a.Assemble()


def _skew_quad_mesh():
    return structured_rect_vol_mesh(
        0.8,
        1.8,
        -0.2,
        0.8,
        quads=True,
        nx=1,
        ny=1,
        mapping=lambda x, y: (0.8 + x + 0.25 * y, -0.2 + y),
        stem="axifem_skew_quad",
    )


@pytest.mark.parametrize("order", [1, 2])
def test_default_closed_form_quads_reject_skewed_geometry(order):
    with pytest.raises(Exception, match="non-axis-aligned quad"):
        _assemble_stiffness(_skew_quad_mesh(), order=order)


def test_q2_curved_quad_accepts_skewed_geometry():
    _assemble_stiffness(_skew_quad_mesh(), order=2, curvedquad=True)


def test_p1_triangle_rejects_singular_r2z_vandermonde():
    axifem = pytest.importorskip("radia.axifem")
    with pytest.raises(Exception, match="singular Vandermonde"):
        axifem.AxiHenrotteFE_P1_Triangle([1.0, 1.0, 1.0], [0.0, 1.0, 2.0])
