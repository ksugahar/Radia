"""BDM2 B-input history state lives on NGSolve integration-rule points."""

import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
from ngsolve.meshes import MakeStructured3DMesh  # noqa: E402

from radia import vim  # noqa: E402


MU0 = 4.0e-7*np.pi


class _TrackingLinearMaterial:
    """Linear B-input law whose committed scalar records local Bz."""

    def state0(self):
        return np.zeros(1)

    def nu_bound(self):
        return 1.0/(MU0*50.0)

    def forward(self, flux_density, states):
        return np.asarray(flux_density, dtype=float)/(MU0*50.0)

    def commit(self, flux_density, states):
        committed = np.asarray(states, dtype=float).copy()
        committed[:, 0] = np.asarray(flux_density, dtype=float)[:, 2]
        return committed


def _curved_hex():
    mesh = MakeStructured3DMesh(
        hexes=True, nx=1, ny=1, nz=1,
        mapping=lambda x, y, z: (x-0.5, y-0.5, z-0.5))
    mesh.Curve(2)
    return mesh


def _play_material():
    radius = np.linspace(0.0, 1.5, 31)
    slope = 1.0/(MU0*80.0)
    eta = np.asarray([0.0, 0.25])
    return vim.PlayHysteresisMaterial(
        2, eta,
        [(radius, slope*radius), (radius, -0.12*slope*radius)])


def _energy_stop_material():
    eta = np.asarray([0.2, 0.5])
    tables = [
        (np.linspace(0.0, radius, 17),
         np.linspace(0.0, peak, 17))
        for radius, peak in zip(eta, (2.0e4, 4.0e4))
    ]
    return vim.EnergyStopMaterial(
        eta, tables, alpha=5.0, gamma=0.0, b_max=1.2)


def test_rt2_history_uses_quadrature_state_and_restarts():
    mesh = _curved_hex()
    material = _TrackingLinearMaterial()
    varying_field = ng.CF((0.0, 0.0, 2.0e4*(1.0 + 0.4*ng.x)))
    with ng.TaskManager():
        solver = vim.HDivSolver(mesh, order=2)
        first = solver.SolveHysteresis(
            [varying_field], material=material,
            tol=1.0e-10, maxit=500, nl_tol=1.0e-8)

    assert first["order"] == 2
    assert first["state_layout"] == "quadrature"
    assert first["state_quadrature_order"] == 3
    assert first["state_points"] > first["n_el"]
    states = first["state"]["material_states"]
    assert states.shape == (first["state_points"], 1)
    assert np.std(states[:, 0]) > 1.0e-3
    assert first["steps"][0]["M"].shape == (first["n_el"], 3)

    with ng.TaskManager():
        continued = solver.SolveHysteresis(
            [[0.0, 0.0, 1.0e4]], material=material,
            initial_state=first["state"], tol=1.0e-10, maxit=500,
            nl_tol=1.0e-8)

    assert continued["state"]["material_states"].shape == states.shape
    assert continued["state"]["state_layout"] == "quadrature"
    assert continued["operator_reused"] is True
    assert solver.operator_build_count == 1
    assert continued["_charge_gram"] is first["_charge_gram"]
    assert continued["charge_gram_wall_s"] == 0.0


@pytest.mark.parametrize(
    ("material_factory", "expected_model", "expected_level"),
    [
        (_play_material, "simplified-play", 3),
        (_energy_stop_material, "b-input-energy-stop", 4),
    ],
)
def test_curved_rt2_production_hysteresis_models_reuse_quadrature_state(
        material_factory, expected_model, expected_level):
    """Actual Play/EnergyStop models run through the curved BDM2 solver."""
    mesh = _curved_hex()
    material = material_factory()
    with ng.TaskManager():
        solver = vim.HDivSolver(mesh, order=2, curve_order=2)
        first = solver.SolveHysteresis(
            [[0.0, 0.0, 1.0e4]], material=material,
            tol=1.0e-10, maxit=500, nl_tol=1.0e-7)
        continued = solver.SolveHysteresis(
            [[0.0, 0.0, -5.0e3]], material=material,
            initial_state=first["state"],
            tol=1.0e-10, maxit=500, nl_tol=1.0e-7)

    assert first["order"] == continued["order"] == 2
    assert first["state_layout"] == continued["state_layout"] == "quadrature"
    assert first["state_points"] > first["n_el"]
    assert first["permanent_magnet_model"] == expected_model
    assert first["permanent_magnet_level"] == expected_level
    assert continued["operator_reused"] is True
    assert solver.operator_build_count == 1
    assert continued["charge_gram_wall_s"] == 0.0
    assert np.isfinite(continued["state"]["material_states"]).all()
