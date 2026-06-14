"""Golden test: 2D bidirectional potential-coordinate map (Tampere, Stage A).

Locks the forward (geometry -> A,phi), the verified Jacobian identity
det[grad A; grad phi] = |B||H| (mu=1, conjugate-harmonic orthogonality), and
the inverse (A,phi -> geometry, x/y harmonic in (A,phi)) + round-trip, on the
analytic complex-potential pair w = zeta^2.

See examples/vim/bidirectional_map_2d.py and the MCP topic
streamfunction("clebsch_3d").
"""
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__),
                                "../../examples/vim"))


@pytest.fixture(scope="module")
def solved():
    from ngsolve import TaskManager
    import bidirectional_map_2d as bm
    with TaskManager():
        mesh_f, gfA, gfP = bm.solve_forward()
        fchk = bm.forward_checks(mesh_f, gfA, gfP)
        mesh_i, gfx, gfy, ichk = bm.solve_inverse()
        rt = bm.roundtrip_error(mesh_i, gfx, gfy)
    return fchk, ichk, rt


def test_forward_recovers_analytic(solved):
    fchk, _, _ = solved
    assert fchk["errA"] < 1e-6, f"A L2 err {fchk['errA']:.2e}"
    assert fchk["errP"] < 1e-6, f"phi L2 err {fchk['errP']:.2e}"


def test_conjugate_harmonic_orthogonality(solved):
    """grad(A) . grad(phi) ~ 0 -- A,phi are conjugate harmonic."""
    fchk, _, _ = solved
    assert fchk["ortho"] < 1e-6, f"orthogonality {fchk['ortho']:.2e}"


def test_jacobian_equals_B_times_H(solved):
    """det[grad A; grad phi] = |grad A||grad phi| = |B||H| (verified identity)."""
    fchk, _, _ = solved
    assert fchk["jac_rel"] < 1e-6, f"Jacobian rel-err {fchk['jac_rel']:.2e}"


def test_inverse_recovers_geometry(solved):
    """x(A,phi), y(A,phi) (harmonic in (A,phi)) match sqrt(A + i phi)."""
    _, ichk, _ = solved
    assert ichk["errx"] < 1e-2, f"x L2 err {ichk['errx']:.2e}"
    assert ichk["erry"] < 1e-2, f"y L2 err {ichk['erry']:.2e}"


def test_roundtrip_identity(solved):
    """p -> (A,phi) -> p recovers the physical point."""
    _, _, rt = solved
    assert rt < 1e-2, f"round-trip error {rt:.2e}"
