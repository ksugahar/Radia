"""Small coupled gates for irreversible hard-magnet state in HDiv-VIM.

The tables below are deliberately synthetic.  They validate the production
state/solver contract; they are not a calibrated commercial magnet grade.
"""

import numpy as np
import ngsolve as ng
from ngsolve.meshes import MakeStructured3DMesh

from radia import vim


def _mesh():
    return MakeStructured3DMesh(
        hexes=True, nx=1, ny=1, nz=1,
        mapping=lambda x, y, z: (x - 0.5, y - 0.5, z - 0.5),
    )


def _material():
    eta = np.array([0.15, 0.35, 0.65])
    peaks = [2.0e4, 4.0e4, 7.0e4]
    tables = [
        (np.linspace(0.0, radius, 17), np.linspace(0.0, peak, 17))
        for radius, peak in zip(eta, peaks)
    ]
    return vim.EnergyStopMaterial(
        eta, tables, alpha=5.0, gamma=0.0, b_max=1.5
    )


_INITIAL_B_PATH = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.9]])
_SOLVE_OPTIONS = dict(gram_eps=1.0e-8, nl_tol=3.0e-5)


def test_reverse_field_causes_irreversible_loss_and_restart_is_reproducible():
    drive = [[0.0, 0.0, 0.0], [0.0, 0.0, -1.2e5], [0.0, 0.0, 0.0]]
    with ng.TaskManager():
        full = vim.SolveHysteresis(
            _mesh(), drive, material=_material(),
            initial_b_path=_INITIAL_B_PATH, **_SOLVE_OPTIONS,
        )
        first = vim.SolveHysteresis(
            _mesh(), drive[:2], material=_material(),
            initial_b_path=_INITIAL_B_PATH, **_SOLVE_OPTIONS,
        )
        restarted = vim.SolveHysteresis(
            _mesh(), drive[2:], material=_material(),
            initial_state=first["state"], **_SOLVE_OPTIONS,
        )

    magnetized = float(full["steps"][0]["M_avg"][2])
    unloaded = float(full["steps"][-1]["M_avg"][2])
    assert full["permanent_magnet_model"] == "b-input-energy-stop"
    assert full["permanent_magnet_level"] == 4
    assert magnetized > 2.0e5
    assert 0.0 < unloaded < 0.5 * magnetized

    np.testing.assert_allclose(
        restarted["steps"][-1]["B"], full["steps"][-1]["B"],
        rtol=0.0, atol=2.0e-13,
    )
    np.testing.assert_allclose(
        restarted["steps"][-1]["M"], full["steps"][-1]["M"],
        rtol=0.0, atol=2.0e-8,
    )
    np.testing.assert_allclose(
        restarted["state"]["material_states"],
        full["state"]["material_states"], rtol=0.0, atol=2.0e-13,
    )
    probes = np.array([[0.0, 0.0, 1.0], [0.2, -0.1, 1.5]])
    field_full = vim.FieldFromSolution(full, probes, algorithm="direct")
    field_restarted = vim.FieldFromSolution(restarted, probes, algorithm="direct")
    assert np.linalg.norm(field_full) > 1.0
    np.testing.assert_allclose(field_restarted, field_full, rtol=0.0, atol=2.0e-8)
    assert full["field_evaluator_stats"]["source_count"] > 0


def test_spatial_coefficient_function_is_an_applied_field_step():
    with ng.TaskManager():
        field = ng.CoefficientFunction((0.0, 0.0, 2.0e4 * (ng.x + 0.5)))
        result = vim.SolveHysteresis(
            _mesh(), [field], material=_material(),
            initial_b_path=_INITIAL_B_PATH, **_SOLVE_OPTIONS,
        )

    step = result["steps"][0]
    assert step["h_applied_is_uniform"] is False
    np.testing.assert_allclose(step["h_applied_avg"], [0.0, 0.0, 1.0e4], atol=1e-10)
