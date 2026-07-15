from __future__ import annotations

import copy
import json

import pytest

from radia_mcp.radia_ngsolve.hysteresis_minor_loop_gate import (
    hysteresis_minor_loop_replay_gate as evaluate_gate,
)
from radia_mcp.radia_ngsolve.server import hysteresis_minor_loop_replay_gate


def _summary() -> dict:
    time_s = [float(index) for index in range(49)]
    drive = [0.0, 1.0, 0.0, -1.0] * 12 + [0.0]
    response = []
    previous = drive[0]
    for index, value in enumerate(drive):
        direction = 1.0 if value > previous else -1.0 if value < previous else (1.0 if index % 2 else -1.0)
        response.append(0.2 * value + 0.05 * direction)
        previous = value
    hysteresis_power = [-1.0 if index % 2 == 0 else 2.0 for index in range(49)]
    return {
        "fresh": {"time_s": time_s, "drive": drive, "response": response},
        "repeat": {"time_s": time_s, "drive": drive, "response": response},
        "saved_reference": {
            "time_s": time_s[1:],
            "drive": drive[1:],
            "response": response[1:],
        },
        "historical_reference": {
            "time_s": time_s[1:],
            "response_magnitude": [abs(value) for value in response[1:]],
        },
        "losses": {
            "time_s": time_s,
            "joule_power_W": [0.0] * len(time_s),
            "hysteresis_power_W": hysteresis_power,
            "iron_power_W": hysteresis_power,
        },
    }


def test_accepts_history_knot_normalization_signed_loss_and_replay() -> None:
    result = evaluate_gate(_summary())
    assert result["status"] == "ok"
    assert result["metrics"]["saved_alignment"]["dropped_fresh_row"] == 0
    assert result["metrics"]["hysteresis_energy_J"] > 0.0


def test_rejects_single_valued_response_control() -> None:
    bad = _summary()
    bad["fresh"]["response"] = [0.2 * value for value in bad["fresh"]["drive"]]
    bad["repeat"]["response"] = list(bad["fresh"]["response"])
    bad["saved_reference"]["response"] = bad["fresh"]["response"][1:]
    bad["historical_reference"]["response_magnitude"] = [abs(value) for value in bad["fresh"]["response"][1:]]
    result = evaluate_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["response_is_multivalued_at_repeated_drive"] is False


def test_rejects_two_missing_serialized_knots() -> None:
    bad = _summary()
    for key in ("time_s", "drive", "response"):
        bad["saved_reference"][key] = bad["saved_reference"][key][1:]
    result = evaluate_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["saved_result_matches_after_one_knot_normalization"] is False


def test_rejects_pointwise_positive_but_negative_integrated_loss() -> None:
    bad = copy.deepcopy(_summary())
    bad["losses"]["hysteresis_power_W"] = [-2.0 if index % 2 == 0 else 1.0 for index in range(49)]
    bad["losses"]["iron_power_W"] = list(bad["losses"]["hysteresis_power_W"])
    result = evaluate_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["signed_hysteresis_power_integrates_to_positive_loss"] is False


def test_stdio_wrapper_returns_json_and_invalid_input() -> None:
    positive = json.loads(hysteresis_minor_loop_replay_gate(json.dumps(_summary())))
    assert positive["status"] == "ok"
    invalid = json.loads(hysteresis_minor_loop_replay_gate("[]"))
    assert invalid["status"] == "invalid_input"
    assert "summary must be an object" in invalid["error"]


def test_rejects_single_knot_iron_hysteresis_power_disagreement() -> None:
    bad = copy.deepcopy(_summary())
    bad["losses"]["iron_power_W"] = list(bad["losses"]["iron_power_W"])
    bad["losses"]["iron_power_W"][20] += 0.01
    result = evaluate_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["iron_power_equals_hysteresis_power"] is False


@pytest.mark.parametrize(
    "case_id",
    ["fresh_repeat", "serialized_knots", "drive_polarity", "joule_loss", "signed_loss"],
)
def test_counterfactual_curriculum90_public(case_id: str) -> None:
    bad = copy.deepcopy(_summary())
    if case_id == "fresh_repeat":
        bad["repeat"]["response"] = list(bad["repeat"]["response"])
        bad["repeat"]["response"][20] *= 1.2
    elif case_id == "serialized_knots":
        bad["saved_reference"]["time_s"].pop(0)
    elif case_id == "drive_polarity":
        bad["fresh"]["drive"] = list(bad["fresh"]["drive"])
        bad["fresh"]["drive"][3] = 1.0
    elif case_id == "joule_loss":
        bad["losses"]["joule_power_W"][20] = 0.1
    else:
        bad["losses"]["hysteresis_power_W"] = [-1.0] * 49
    result = json.loads(hysteresis_minor_loop_replay_gate(json.dumps(bad)))
    assert result["status"] in {"needs_attention", "invalid_input"}


def test_generalization_v3s_rejects_historical_response_drift() -> None:
    bad = copy.deepcopy(_summary())
    bad["historical_reference"]["response_magnitude"][20] *= 1.5
    result = evaluate_gate(bad)
    assert result["status"] == "needs_attention"
