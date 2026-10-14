from __future__ import annotations

import copy
import json
import math

from radia_mcp.radia_ngsolve.server import periodic_unwrapped_pm_machine_replay_gate


def _summary() -> dict:
    circumference = math.pi * 110.0
    run = {
        "element_type": "TL3",
        "node_count": 10600,
        "element_count": 20800,
        "energy": 21.78,
        "coenergy": 22.36,
        "normal_flux_rms": 0.83,
        "normal_flux_mean_relative": 0.043,
        "half_turn_antiperiodicity_relative_error": 0.045,
        "finite_profile_coverage": 0.96,
        "antiperiodic_pair_coverage": 0.92,
    }
    replay = copy.deepcopy(run)
    replay["node_count"] = 10700
    replay["element_count"] = 21000
    replay["energy"] *= 1.0002
    replay["normal_flux_rms"] *= 1.0002
    return {
        "machine": {
            "slot_count": 12,
            "pole_count": 10,
            "circumference": circumference,
            "slot_pitch": circumference / 12.0,
            "pole_pitch": circumference / 10.0,
            "field_symmetry_shift": "half_circumference",
        },
        "periodic_boundary": {"master_count": 33, "slave_count": 33},
        "runs": [run, replay],
    }


def _call(summary: dict) -> dict:
    return json.loads(periodic_unwrapped_pm_machine_replay_gate(json.dumps(summary)))


def test_periodic_pm_machine_gate_accepts_topology_aware_replay():
    result = _call(_summary())
    assert result["status"] == "ok"


def test_periodic_pm_machine_gate_rejects_one_pole_symmetry_claim():
    summary = _summary()
    summary["machine"]["field_symmetry_shift"] = "one_pole_pitch"
    result = _call(summary)
    assert result["status"] == "needs_attention"
    assert "symmetry_shift_is_topology_aware" in result["issues"]


def test_periodic_pm_machine_gate_rejects_observable_drift():
    summary = _summary()
    summary["runs"][1]["energy"] *= 1.1
    result = _call(summary)
    assert result["status"] == "needs_attention"
    assert "energy_observables_replay" in result["issues"]


def test_periodic_pm_machine_gate_rejects_even_half_turn_pole_count():
    summary = _summary()
    summary["machine"]["pole_count"] = 12
    summary["machine"]["pole_pitch"] = summary["machine"]["circumference"] / 12.0
    result = _call(summary)
    assert result["status"] == "needs_attention"
    assert "half_turn_magnetization_reverses" in result["issues"]
