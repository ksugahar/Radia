import copy
import json

from radia_mcp.radia_ngsolve.server import two_terminal_dc_conduction_power_gate
from radia_mcp.radia_ngsolve.two_terminal_dc_conduction_gate import (
    two_terminal_dc_conduction_power_gate as build_gate,
)


def _run(run_id: str) -> dict:
    return {
        "run_id": run_id,
        "solver_complete": True,
        "port_currents_a": [-0.9999703, 0.9999706],
        "port_potentials_v": [0.0, 1013.9277],
        "integrated_loss_w": 1013.8976,
        "loss_parts_w": [1013.8976],
        "adaptive_rows": [
            {"mesh_cells": 502212, "degrees_of_freedom": 1007045, "power_w": 1013.1410},
            {"mesh_cells": 503434, "degrees_of_freedom": 1008873, "power_w": 1013.8976},
        ],
        "final_adaptive_relative_error": 7.46221e-4,
    }


def summary() -> dict:
    return {"runs": [_run("a"), _run("b")]}


def test_accepts_two_fresh_stationary_current_runs():
    result = build_gate(summary())
    assert result["status"] == "ok"
    assert result["solver_ready"] is True
    assert json.loads(two_terminal_dc_conduction_power_gate(json.dumps(summary())))["status"] == "ok"


def test_rejects_current_and_power_imbalance():
    bad = summary()
    bad["runs"][1]["port_currents_a"][1] = 0.8
    bad["runs"][1]["integrated_loss_w"] = 700.0
    result = build_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["opposed_terminal_currents_close"] is False
    assert result["checks"]["terminal_power_matches_integrated_joule_loss"] is False


def test_rejects_nonconverged_or_nonrepeatable_adaptive_result():
    bad = copy.deepcopy(summary())
    bad["runs"][1]["adaptive_rows"][-1]["mesh_cells"] = 500000
    bad["runs"][1]["final_adaptive_relative_error"] = 0.02
    bad["runs"][1]["port_potentials_v"][1] *= 1.05
    result = build_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["adaptive_mesh_and_dof_strictly_increase"] is False
    assert result["checks"]["final_adaptive_error_is_bounded"] is False
    assert result["checks"]["fresh_runs_repeat_scalar_observables"] is False
