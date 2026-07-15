"""RT1/Q2 and RT2/Q3 contracts for the production planar HDiv-VIM path."""

import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
pytest.importorskip("netgen.occ")
from netgen.occ import OCCGeometry, WorkPlane  # noqa: E402

from radia import vim  # noqa: E402


def _ellipse_mesh(curve_order):
    mesh = ng.Mesh(
        OCCGeometry(WorkPlane().Ellipse(0.2, 0.1).Face(), dim=2)
        .GenerateMesh(maxh=0.1))
    mesh.Curve(curve_order)
    return mesh


def test_planar_rt2_q3_matches_ellipse_demag_and_native_field():
    mesh = _ellipse_mesh(3)
    with ng.TaskManager():
        result = vim.Solve(
            mesh, order=2, mu_r=1000.0,
            H_ext=ng.CoefficientFunction((1.0e5, 0.0)), tol=1e-11)
        field = result["body"].H_at(
            np.array([[0.3, 0.0], [0.0, 0.2], [-0.3, 0.1]]), result["m"])

    assert result["order"] == 2
    assert result["geometry_order"] == 3
    np.testing.assert_allclose(result["demag_factors"], (1.0/3.0, 2.0/3.0), atol=1e-5)
    assert np.isfinite(field).all()


def test_planar_rt1_and_rt2_agree_on_the_same_q2_mesh():
    mesh = _ellipse_mesh(2)
    applied = ng.CoefficientFunction((1.0e5/np.sqrt(2.0), 1.0e5/np.sqrt(2.0)))
    with ng.TaskManager():
        rt1 = vim.Solve(mesh, order=1, mu_r=1000.0, H_ext=applied, tol=1e-11)
        rt2 = vim.Solve(mesh, order=2, mu_r=1000.0, H_ext=applied, tol=1e-11)

    np.testing.assert_allclose(rt1["M_avg"], rt2["M_avg"], rtol=5e-5)
    assert rt2["ndof"] > rt1["ndof"]


def test_planar_rt1_rejects_q3_geometry():
    with pytest.raises(ValueError, match="does not support geometry order 3 for 2D tri RT1"):
        vim.PlanarDemagBody(_ellipse_mesh(3), order=1)
