from __future__ import annotations

import copy
import json

from radia_mcp.radia_ngsolve.server import transient_conductor_replay_identity_gate as mcp_gate
from radia_mcp.radia_ngsolve.transient_conductor_replay_gate import (
    transient_conductor_replay_identity_gate,
)


def _variant(label: str) -> dict:
    current = [0.0, 1.0, 0.9]
    resistance = [0.0, 2.0, 2.1]
    inductance = [0.0, 3.0, 3.1]
    return {
        "label": label,
        "times_s": [0.0, 1.0e-6, 2.0e-6],
        "current_a": current,
        "joule_loss_w": [0.0, 2.0, 0.9 * 0.9 * 2.1],
        "circuit_power_w": [0.0, 2.5, 2.0],
        "flux_linkage_wb": [0.0, 3.0, 0.9 * 3.1],
        "resistance_ohm": resistance,
        "inductance_h": inductance,
        "mesh_elements": 100,
        "mesh_vertices": 60,
    }


def _summary() -> dict:
    return {
        "original_step_count": 10000,
        "bounded_step_count": 3,
        "time_step_s": 1.0e-6,
        "replays": [
            {"label": "a", "variants": [_variant("definition"), _variant("alias")]},
            {"label": "b", "variants": [_variant("definition"), _variant("alias")]},
        ],
        "timing_breakdown_s": {"stage": 4.0, "solve": 2.0, "extract": 1.0, "verify": 0.5},
    }


def test_accepts_full_waveform_replay_and_dispatches():
    result = transient_conductor_replay_identity_gate(_summary())
    assert result["status"] == "ok"
    assert result["metrics"]["maximum_joule_i2r_relative_error"] < 1.0e-15
    assert json.loads(mcp_gate(json.dumps(_summary())))["status"] == "ok"


def test_rejects_identity_replay_and_timing_failures():
    row = copy.deepcopy(_summary())
    row["replays"][0]["variants"][1]["joule_loss_w"][2] *= 1.2
    row["replays"][1]["variants"][0]["flux_linkage_wb"][1] *= 0.8
    row["timing_breakdown_s"].pop("verify")
    result = transient_conductor_replay_identity_gate(row)
    assert result["status"] == "needs_attention"
    assert result["checks"]["joule_loss_closes_i_squared_r"] is False
    assert result["checks"]["flux_linkage_closes_l_times_i"] is False
    assert result["checks"]["exactly_four_timing_stages"] is False


def test_rejects_non_monotonic_time_and_mesh_drift():
    row = copy.deepcopy(_summary())
    row["replays"][1]["variants"][1]["times_s"][2] = 0.5e-6
    row["replays"][1]["variants"][1]["mesh_elements"] += 1
    result = transient_conductor_replay_identity_gate(row)
    assert result["status"] == "needs_attention"
    assert result["checks"]["time_axes_are_strictly_increasing"] is False
    assert result["checks"]["mesh_inventory_is_replay_invariant"] is False
