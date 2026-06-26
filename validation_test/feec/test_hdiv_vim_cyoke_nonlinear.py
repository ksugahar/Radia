"""Golden test: NONLINEAR HDiv-VIM on a C-yoke (non-convex, reentrant corners, NON-uniform M),
cross-validated vs the shipped Radia solver -- closes the "C-yoke accuracy verification" gate.

The C-yoke nonlinear (with the production C++ analytic charge Gram, the required volume Gram for
div M != 0) converges in a handful of damped-Newton iters AND its volume-avg M_z matches the shipped
Radia MMM solver on the SAME flat mesh / M-H law / applied field to < 3% (a C-yoke is flat, so Radia is
a valid cross-reference -- cross-method agreement, not vs analytic truth).
"""
import os
import sys

import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.occ")
pytest.importorskip("radia")

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "examples", "vim"))
import hdiv_cyoke_nonlinear as cy  # noqa: E402


def test_cyoke_nonlinear_converges_and_matches_radia():
    """The gate: C-yoke non-uniform nonlinear converges (C++ analytic Gram) and matches shipped Radia."""
    mesh = cy.cyoke_mesh(0.02)
    Mh, nit = cy.hdiv_cyoke_Mz(mesh, 1000.0, 1.0e6, 2.0e5)
    assert nit < 60, f"C-yoke Newton not converging with the C++ analytic Gram: {nit} iters"
    Mr = cy.radia_cyoke_Mz(mesh, 1000.0, 1.0e6, 2.0e5)
    assert abs(Mh / Mr - 1.0) < 0.03, f"HDiv M_z {Mh:.1f} vs Radia {Mr:.1f} disagree by {100*(Mh/Mr-1):.2f}%"
