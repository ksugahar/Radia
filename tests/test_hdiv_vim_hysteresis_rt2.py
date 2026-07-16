"""RT2 B-input history state lives on NGSolve integration-rule points."""

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


def test_rt2_history_uses_quadrature_state_and_restarts():
    mesh = _curved_hex()
    material = _TrackingLinearMaterial()
    varying_field = ng.CF((0.0, 0.0, 2.0e4*(1.0 + 0.4*ng.x)))
    with ng.TaskManager():
        first = vim.SolveHysteresis(
            mesh, [varying_field], material=material, order=2,
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
        continued = vim.SolveHysteresis(
            mesh, [[0.0, 0.0, 1.0e4]], material=material, order=2,
            initial_state=first["state"], tol=1.0e-10, maxit=500,
            nl_tol=1.0e-8, _prepared_operator=first["_prepared_operator"])

    assert continued["state"]["material_states"].shape == states.shape
    assert continued["state"]["state_layout"] == "quadrature"
    assert continued["prepared_operator_reused"] is True
    assert continued["_charge_gram"] is first["_charge_gram"]
    assert continued["charge_gram_wall_s"] == 0.0
