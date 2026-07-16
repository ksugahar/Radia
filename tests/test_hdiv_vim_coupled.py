"""Fast contracts for cached and mutually coupled HDiv body solves."""

import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
from ngsolve.meshes import MakeStructured3DMesh  # noqa: E402

from radia import vim  # noqa: E402


def _hex(x0, x1):
    return MakeStructured3DMesh(
        hexes=True, nx=1, ny=1, nz=1,
        mapping=lambda x, y, z: (
            x0 + (x1-x0)*x, 0.2*(y-0.5), 0.2*(z-0.5)))


def test_prepared_operator_reuses_geometry_only_charge_gram():
    mesh = _hex(-0.1, 0.1)
    with ng.TaskManager():
        first = vim.Solve(
            mesh, mu_r=20.0, H_ext=ng.CF((0.0, 0.0, 1.0e4)),
            gram_eps=1.0e-8, tol=1.0e-10)
        second = vim.Solve(
            mesh, mu_r=20.0, H_ext=ng.CF((0.0, 0.0, 2.0e4)),
            gram_eps=1.0e-8, tol=1.0e-10,
            _prepared_operator=first["_prepared_operator"])

    assert first["prepared_operator_reused"] is False
    assert second["prepared_operator_reused"] is True
    assert second["_charge_gram"] is first["_charge_gram"]
    assert second["charge_gram_wall_s"] == 0.0
    np.testing.assert_allclose(
        second["_m_coefficients"], 2.0*first["_m_coefficients"],
        rtol=2.0e-11, atol=1.0e-8)


def test_linear_recoil_pm_and_nonlinear_iron_reach_a_block_fixed_point():
    pm_mesh = _hex(-0.45, -0.25)
    iron_mesh = _hex(0.25, 0.45)
    bh_table = np.array([
        [0.0, 0.0],
        [1.0e2, 0.8],
        [1.0e3, 1.35],
        [1.0e5, 1.8],
    ])
    bodies = [
        vim.CoupledBody(
            pm_mesh, "pm", mu_r=1.05, B_r=(0.0, 0.0, 1.2),
            solve_options={"gram_eps": 1.0e-8, "tol": 1.0e-9}),
        vim.CoupledBody(
            iron_mesh, "iron", bh_table=bh_table,
            solve_options={
                "gram_eps": 1.0e-8, "tol": 1.0e-8,
                "nl_tol": 1.0e-6, "nl_maxit": 40}),
    ]

    with ng.TaskManager():
        result = vim.SolveCoupled(bodies, tol=2.0e-6, maxit=12)

    assert result["converged"] is True
    assert result["permanent_magnet_body_count"] == 1
    assert result["nonlinear_iron_body_count"] == 1
    assert 2 <= result["iterations"] <= 12
    assert result["relative_step"] < 2.0e-6
    assert all(body["prepared_operator_reused"] for body in result["bodies"])
    assert all(body["charge_gram_wall_s"] == 0.0 for body in result["bodies"])
    assert np.linalg.norm(result["bodies"][0]["M_avg"]) > 1.0e5
    assert np.linalg.norm(result["bodies"][1]["M_avg"]) > 1.0e2
    field = vim.FieldFromCoupledSolution(
        result, [[0.0, 0.0, 0.0], [0.0, 0.0, 0.5]], algorithm="direct")
    assert field.shape == (2, 3)
    assert np.isfinite(field).all()


def test_segmented_recoil_magnets_keep_separate_normal_trace_spaces():
    left = _hex(-0.2, 0.0)
    right = _hex(0.0, 0.2)
    bodies = [
        vim.CoupledBody(
            left, "left", mu_r=1.04, B_r=(1.1, 0.0, 0.0),
            solve_options={"gram_eps": 1.0e-8, "tol": 1.0e-9}),
        vim.CoupledBody(
            right, "right", mu_r=1.04, B_r=(-1.1, 0.0, 0.0),
            solve_options={"gram_eps": 1.0e-8, "tol": 1.0e-9}),
    ]

    with ng.TaskManager():
        result = vim.SolveCoupled(bodies, tol=2.0e-7, maxit=12)

    assert result["permanent_magnet_body_count"] == 2
    assert result["nonlinear_iron_body_count"] == 0
    assert result["relative_step"] < 2.0e-7
    left_m, right_m = [body["M_avg"] for body in result["bodies"]]
    assert left_m[0] > 1.0e5
    assert right_m[0] < -1.0e5
    assert all(body["prepared_operator_reused"] for body in result["bodies"])


def test_coupled_bodies_reject_a_shared_conforming_space():
    mesh = _hex(-0.1, 0.1)
    body_a = vim.CoupledBody(mesh, "a", mu_r=1.05, B_r=(1.0, 0.0, 0.0))
    body_b = vim.CoupledBody(mesh, "b", mu_r=1.05, B_r=(-1.0, 0.0, 0.0))
    with pytest.raises(ValueError, match="distinct mesh object"):
        vim.SolveCoupled([body_a, body_b])
