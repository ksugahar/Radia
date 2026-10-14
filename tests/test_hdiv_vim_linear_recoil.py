"""Fast contracts for the linear-recoil permanent-magnet HDiv model."""

import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
from ngsolve.meshes import MakeStructured2DMesh, MakeStructured3DMesh

from radia import vim


MU0 = 4.0e-7 * np.pi


def _hex_mesh():
    return MakeStructured3DMesh(
        hexes=True, nx=1, ny=1, nz=1,
        mapping=lambda x, y, z: (0.4*x-0.2, 0.3*y-0.15, 0.2*z-0.1),
    )


@pytest.mark.parametrize("order", [1, 2])
def test_linear_recoil_is_exactly_the_shifted_symmetric_hdiv_system(order):
    mesh = _hex_mesh()
    mu_rec = 1.05
    B_r = np.array([0.1, -0.2, 1.15])
    H_applied = ng.CoefficientFunction((2.0e4, -1.0e4, 3.0e4))
    equivalent_H = H_applied + ng.CoefficientFunction(tuple(B_r)) / (
        MU0 * (mu_rec - 1.0)
    )

    with ng.TaskManager():
        permanent_magnet = vim.Solve(
            mesh, mu_r=mu_rec, B_r=B_r, H_ext=H_applied, order=order, tol=1.0e-10
        )
        shifted_linear = vim.Solve(
            mesh, mu_r=mu_rec, H_ext=equivalent_H, order=order, tol=1.0e-10
        )

    np.testing.assert_allclose(
        permanent_magnet["_m_coefficients"], shifted_linear["_m_coefficients"],
        rtol=2.0e-13, atol=1.0e-8
    )
    np.testing.assert_allclose(
        permanent_magnet["M_avg"], shifted_linear["M_avg"],
        rtol=2.0e-13, atol=1.0e-8,
    )
    assert permanent_magnet["permanent_magnet_model"] == "linear-recoil"
    assert permanent_magnet["permanent_magnet_level"] == 2
    assert permanent_magnet["recoil_mu_r"] == pytest.approx(mu_rec)
    assert permanent_magnet["B_r_supplied"] is True
    assert permanent_magnet["order"] == order
    assert permanent_magnet["nonlinear"] is False
    field = vim.FieldFromSolution(
        permanent_magnet, np.array([[0.0, 0.0, 0.5], [0.3, -0.2, 0.4]]),
        algorithm="direct",
    )
    assert field.shape == (2, 3)
    assert np.isfinite(field).all()


def test_linear_recoil_accepts_spatial_remanence_and_zero_applied_field():
    mesh = _hex_mesh()
    distributed_B_r = ng.CoefficientFunction((
        0.05 * ng.x, 0.1 * ng.y, 1.1 + 0.2 * ng.z
    ))

    with ng.TaskManager():
        result = vim.Solve(mesh, mu_r=1.08, B_r=distributed_B_r, tol=1.0e-9)

    assert result["_B_r"] is distributed_B_r
    assert np.isfinite(result["_m_coefficients"]).all()
    assert np.linalg.norm(result["M_avg"]) > 1.0e5


def test_planar_linear_recoil_uses_the_same_material_law():
    mesh = MakeStructured2DMesh(quads=True, nx=2, ny=2)
    mu_rec = 1.04
    B_r = np.array([0.15, 0.9])
    equivalent_H = ng.CoefficientFunction(tuple(B_r)) / (MU0 * (mu_rec - 1.0))

    with ng.TaskManager():
        permanent_magnet = vim.Solve(mesh, mu_r=mu_rec, B_r=B_r)
        shifted_linear = vim.Solve(mesh, mu_r=mu_rec, H_ext=equivalent_H)

    np.testing.assert_allclose(
        permanent_magnet["m"], shifted_linear["m"], rtol=2.0e-13, atol=1.0e-8
    )
    assert permanent_magnet["permanent_magnet_model"] == "linear-recoil"
    assert permanent_magnet["permanent_magnet_level"] == 2


@pytest.mark.parametrize(
    "kwargs, message",
    [
        ({"mu_r": 1.0, "B_r": (0.0, 0.0, 1.0)}, "mu_r > 1"),
        ({"mu_r": np.nan, "B_r": (0.0, 0.0, 1.0)}, "mu_r > 1"),
        ({"mu_r": 1.05, "B_r": (0.0, np.inf, 1.0)}, "finite values"),
        ({"mu_r": 1.05, "B_r": (0.0, 1.0)}, "length-3 vector"),
        ({"mu_r": {"body": 1.05}, "B_r": (0.0, 0.0, 1.0)}, "scalar recoil"),
        ({
            "mu_r": 1.05,
            "B_r": (0.0, 0.0, 1.0),
            "bh_table": [[0.0, 0.0], [1.0e5, 1.0]],
        }, "cannot be combined"),
    ],
)
def test_linear_recoil_rejects_ambiguous_material_contracts(kwargs, message):
    with pytest.raises((ValueError, NotImplementedError), match=message):
        vim.Solve(_hex_mesh(), **kwargs)
