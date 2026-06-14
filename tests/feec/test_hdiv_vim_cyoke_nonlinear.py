"""Golden test: NONLINEAR HDiv-VIM on a C-yoke (non-convex, reentrant corners, NON-uniform M),
cross-validated vs the shipped Radia solver -- closes the "C-yoke accuracy verification" gate.

(1) The C-yoke nonlinear (with analytic_gram, the required volume Gram for div M != 0) converges in a
    handful of damped-Newton iters AND its volume-avg M_z matches the shipped Radia MMM solver on the
    SAME flat mesh / M-H law / applied field to < 3% (a C-yoke is flat, so Radia is a valid
    cross-reference -- cross-method agreement, not vs analytic truth).
(2) The surface-only wilton_surface Gram does NOT converge for a non-uniform-M nonlinear problem -> the
    solver now RAISES (fail-loud), instead of silently returning a wrong M -- locked on the cube (also
    non-uniform nonlinear, fast).
"""
import os
import sys

import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.occ")
pytest.importorskip("radia")

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "examples", "vim"))
import hdiv_cyoke_nonlinear as cy  # noqa: E402
from radia.vim import _nonlinear as nl  # noqa: E402
import ngsolve as ng  # noqa: E402


def test_cyoke_nonlinear_converges_and_matches_radia():
    """The gate: C-yoke non-uniform nonlinear converges fast (analytic_gram) and matches shipped Radia."""
    mesh = cy.cyoke_mesh(0.02)
    Mh, nit = cy.hdiv_cyoke_Mz(mesh, 1000.0, 1.0e6, 2.0e5)
    assert nit < 15, f"C-yoke Newton not fast with analytic_gram: {nit} iters"
    Mr = cy.radia_cyoke_Mz(mesh, 1000.0, 1.0e6, 2.0e5)
    assert abs(Mh / Mr - 1.0) < 0.03, f"HDiv M_z {Mh:.1f} vs Radia {Mr:.1f} disagree by {100*(Mh/Mr-1):.2f}%"


def test_nonuniform_nonlinear_wrong_gram_fails_loud():
    """Fail-loud: the surface-only wilton_surface Gram cannot converge a NON-uniform-M nonlinear problem
    (div M != 0); the solver must RAISE pointing at analytic_gram, not silently return a wrong M.  Locked
    on the convex cube (non-uniform nonlinear, fast)."""
    from netgen.csg import CSGeometry, OrthoBrick, Pnt
    g = CSGeometry(); g.Add(OrthoBrick(Pnt(-0.5, -0.5, -0.5), Pnt(0.5, 0.5, 0.5)))
    mesh = ng.Mesh(g.GenerateMesh(maxh=0.4))
    with ng.TaskManager():
        with pytest.raises(RuntimeError, match="analytic_gram"):
            nl.solve_nonlinear_newton(mesh, 1000.0, 1.0e6, 2.0e5, wilton_surface=True, maxit=20)
