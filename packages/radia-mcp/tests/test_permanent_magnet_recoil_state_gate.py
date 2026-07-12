from __future__ import annotations

import copy
import json

from radia_mcp.radia_ngsolve.server import permanent_magnet_recoil_state_gate


def _summary() -> dict:
    fresh = [1.2, 0.006, 0.075, 0.065, 0.29, 0.016]
    return {
        "model_contract": {
            "initial": "nonlinear_in_circuit",
            "out_of_circuit": "nonlinear_open_circuit",
            "recoil_return": "linear_recoil_in_circuit",
            "same_geometry": True,
            "same_mesh": True,
            "same_observation_points": True,
        },
        "units": {"magnetic_flux_density": "T"},
        "states": {
            "initial": {"on_axis": fresh[0], "off_axis": fresh[1]},
            "out_of_circuit": {"on_axis": fresh[2], "off_axis": fresh[3]},
            "recoil_return": {"on_axis": fresh[4], "off_axis": fresh[5]},
        },
        "stored_reference": [value * 0.9999 for value in fresh],
        "fresh_replay": fresh,
    }


def _call(summary: dict) -> dict:
    return json.loads(permanent_magnet_recoil_state_gate(json.dumps(summary)))


def test_accepts_partial_recoil_and_cross_run_replay():
    result = _call(_summary())
    assert result["status"] == "ok"
    assert result["metrics"]["recovery_ratios"]["recoil_to_initial_on_axis"] < 1.0


def test_rejects_full_recovery_claim_and_changed_mesh():
    summary = copy.deepcopy(_summary())
    summary["states"]["recoil_return"]["on_axis"] = 1.2
    summary["model_contract"]["same_mesh"] = False
    result = _call(summary)
    assert result["status"] == "needs_attention"
    assert set(result["issues"]) >= {
        "geometry_mesh_and_points_shared",
        "recoil_is_partial_not_full_recovery",
    }


def test_rejects_replay_drift_and_open_field_concentration():
    summary = copy.deepcopy(_summary())
    summary["fresh_replay"][5] *= 1.1
    summary["states"]["out_of_circuit"]["off_axis"] = 0.005
    result = _call(summary)
    assert result["status"] == "needs_attention"
    assert set(result["issues"]) >= {
        "saved_and_fresh_replay_close",
        "open_field_is_spatially_spread",
    }
