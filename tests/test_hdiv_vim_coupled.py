"""Fast contracts for persistent and mutually coupled HDiv body solves."""

import inspect
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


def test_hdiv_solver_reuses_geometry_without_exposing_private_cache():
    mesh = _hex(-0.1, 0.1)
    with ng.TaskManager():
        solver = vim.HDivSolver(mesh, gram_eps=1.0e-8)
        first = solver.Solve(
            mu_r=20.0, H_ext=ng.CF((0.0, 0.0, 1.0e4)), tol=1.0e-10)
        second = solver.Solve(
            mu_r=20.0, H_ext=ng.CF((0.0, 0.0, 2.0e4)), tol=1.0e-10)

    assert "_operator_cache" not in inspect.signature(vim.Solve).parameters
    assert "_operator_cache" not in inspect.signature(vim.SolveHysteresis).parameters
    assert "_operator_cache" not in first
    assert "_operator_cache" not in second
    with pytest.raises(TypeError, match="HDivSolver"):
        vim.Solve(mesh, mu_r=20.0, H_ext=ng.CF((0.0, 0.0, 1.0)),
                  _operator_cache={})
    assert solver.operator_build_count == 1
    assert solver.charge_gram is first["_charge_gram"]
    assert first["operator_reused"] is False
    assert second["operator_reused"] is True
    assert second["_charge_gram"] is first["_charge_gram"]
    assert second["charge_gram_wall_s"] == 0.0
    np.testing.assert_allclose(
        second["_m_coefficients"], 2.0*first["_m_coefficients"],
        rtol=2.0e-11, atol=1.0e-8)


def test_hdiv_solver_restores_geometry_mass_after_per_region_case():
    mesh = _hex(-0.1, 0.1)
    material = mesh.GetMaterials()[0]
    field = ng.CF((0.0, 0.0, 1.0e4))
    with ng.TaskManager():
        solver = vim.HDivSolver(mesh, gram_eps=1.0e-8)
        reference = solver.Solve(mu_r=20.0, H_ext=field, tol=1.0e-10)
        solver.Solve(mu_r={material: 20.0}, H_ext=field, tol=1.0e-10)
        restored = solver.Solve(mu_r=20.0, H_ext=field, tol=1.0e-10)

    assert solver.operator_build_count == 1
    assert restored["operator_reused"] is True
    np.testing.assert_allclose(
        restored["_m_coefficients"], reference["_m_coefficients"],
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
    assert all(body["operator_reused"] for body in result["bodies"])
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
    assert all(body["operator_reused"] for body in result["bodies"])


def test_energy_stop_pm_and_nonlinear_iron_commit_one_physical_step():
    eta = np.array([0.15, 0.35, 0.65])
    peaks = [2.0e4, 4.0e4, 7.0e4]
    tables = [
        (np.linspace(0.0, radius, 17), np.linspace(0.0, peak, 17))
        for radius, peak in zip(eta, peaks)
    ]
    material = vim.EnergyStopMaterial(
        eta, tables, alpha=5.0, gamma=0.0, b_max=1.5)
    initial_b_path = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.9]])
    pm = vim.CoupledHistoryBody(
        _hex(-0.45, -0.25), "history-pm", material,
        initial_b_path=initial_b_path,
        solve_options={"gram_eps": 1.0e-8, "nl_tol": 3.0e-5})
    iron = vim.CoupledBody(
        _hex(0.25, 0.45), "nonlinear-iron",
        bh_table=np.array([
            [0.0, 0.0], [1.0e2, 0.8], [1.0e3, 1.35], [1.0e5, 1.8]]),
        solve_options={
            "gram_eps": 1.0e-8, "tol": 1.0e-8,
            "nl_tol": 1.0e-6, "nl_maxit": 40})

    with ng.TaskManager():
        result = vim.SolveCoupledHysteresis(
            [pm], [iron], [[0.0, 0.0, 0.0]], tol=2.0e-6, maxit=12)
        replay = vim.SolveHysteresis(
            pm.mesh,
            [vim.FieldCoefficientFromSolution(result["bodies"][0])],
            material=material, initial_b_path=initial_b_path,
            gram_eps=1.0e-8, nl_tol=3.0e-5)

    step = result["steps"][0]
    history_result = result["history_bodies"][0]
    assert result["history_material_models"] == ("b-input-energy-stop",)
    assert result["history_material_levels"] == (4,)
    assert result["operator_build_counts"] == {
        "history-pm": 1, "nonlinear-iron": 1}
    assert result["nonlinear_iron_body_count"] == 1
    assert 2 <= step["iterations"] <= 12
    assert step["relative_step"] < 2.0e-6
    assert history_result["operator_reused"] is True
    assert result["bodies"][0]["operator_reused"] is True
    assert history_result["charge_gram_wall_s"] == 0.0
    assert result["bodies"][0]["charge_gram_wall_s"] == 0.0
    np.testing.assert_allclose(
        history_result["_m_coefficients"],
        replay["_m_coefficients"], rtol=2.0e-11, atol=1.0e-8)
    np.testing.assert_allclose(
        result["states"][0]["material_states"],
        replay["state"]["material_states"], rtol=0.0, atol=2.0e-12)
    field = vim.FieldFromCoupledHysteresis(
        result, [[0.0, 0.0, 0.0], [0.0, 0.0, 0.5]])
    assert np.isfinite(field).all()


class _TrackingLinearMaterial:
    def state0(self):
        return np.zeros(1)

    def nu_bound(self):
        return 1.0/(4.0e-7*np.pi*50.0)

    def forward(self, flux_density, states):
        return np.asarray(flux_density, dtype=float)/(4.0e-7*np.pi*50.0)

    def commit(self, flux_density, states):
        result = np.asarray(states, dtype=float).copy()
        result[:, 0] = np.asarray(flux_density, dtype=float)[:, 2]
        return result


def test_multiple_history_bodies_commit_independent_states_after_global_convergence():
    history_bodies = [
        vim.CoupledHistoryBody(
            _hex(-0.45, -0.25), "pm-left", _TrackingLinearMaterial(),
            solve_options={"gram_eps": 1.0e-8, "nl_tol": 1.0e-8}),
        vim.CoupledHistoryBody(
            _hex(0.25, 0.45), "pm-right", _TrackingLinearMaterial(),
            applied_field=(0.0, 0.0, -3.0e3),
            solve_options={"gram_eps": 1.0e-8, "nl_tol": 1.0e-8}),
    ]

    with ng.TaskManager():
        result = vim.SolveCoupledHysteresis(
            history_bodies, [], [[0.0, 0.0, 1.0e4]],
            tol=1.0e-7, maxit=12)

    assert result["history_body_count"] == 2
    assert result["body_names"] == ("pm-left", "pm-right")
    assert result["operator_build_counts"] == {"pm-left": 1, "pm-right": 1}
    assert len(result["states"]) == 2
    assert all("_operator_cache" not in body for body in result["history_bodies"])
    left_state, right_state = result["states"]
    assert np.mean(left_state["material_states"][:, 0]) > np.mean(
        right_state["material_states"][:, 0])
    field = vim.FieldFromCoupledHysteresis(result, [[0.0, 0.0, 0.0]])
    assert field.shape == (1, 3)
    assert np.isfinite(field).all()


def test_coupled_bodies_reject_a_shared_conforming_space():
    mesh = _hex(-0.1, 0.1)
    body_a = vim.CoupledBody(mesh, "a", mu_r=1.05, B_r=(1.0, 0.0, 0.0))
    body_b = vim.CoupledBody(mesh, "b", mu_r=1.05, B_r=(-1.0, 0.0, 0.0))
    with pytest.raises(ValueError, match="distinct mesh object"):
        vim.SolveCoupled([body_a, body_b])
