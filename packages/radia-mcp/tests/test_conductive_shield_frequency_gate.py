from __future__ import annotations

import copy
import json

from radia_mcp.radia_ngsolve.conductive_shield_frequency_gate import (
    magnetic_conductive_shield_frequency_gate,
)
from radia_mcp.radia_ngsolve.server import magnetic_conductive_shield_frequency_gate as mcp_gate


def _rows(ratios: list[float]) -> list[dict]:
    rows = []
    for frequency, ratio in zip([10.0, 100.0, 1000.0, 10000.0], ratios, strict=True):
        rows.append(
            {
                "frequency_hz": frequency,
                "primary_response": [1.0, 0.1],
                "secondary_response": [ratio, 0.1 * ratio],
                "secondary_flux_linkage": [ratio * 1.0e-6, ratio * 0.1e-6],
                "maximum_faraday_relative_error": 2.0e-5,
            }
        )
    return rows


def _summary() -> dict:
    baseline = _rows([1.0, 1.0, 1.0, 1.0])
    shielded = _rows([1.10, 1.05, 0.80, 0.50])
    return {
        "models": [
            {"label": "baseline", "replays": [{"rows": copy.deepcopy(baseline)}, {"rows": copy.deepcopy(baseline)}]},
            {"label": "shielded", "replays": [{"rows": copy.deepcopy(shielded)}, {"rows": copy.deepcopy(shielded)}]},
        ],
        "timing_breakdown_s": {"copy": 1.0, "mesh": 2.0, "solve": 3.0, "verify": 1.0},
    }


def test_accepts_dual_regime_shield_and_dispatches():
    result = magnetic_conductive_shield_frequency_gate(_summary())
    assert result["status"] == "ok"
    assert result["metrics"]["crossover_frequency_bracket_hz"] == [100.0, 1000.0]
    assert json.loads(mcp_gate(json.dumps(_summary())))["status"] == "ok"


def test_rejects_missing_attenuation_and_replay_drift():
    row = copy.deepcopy(_summary())
    for sample in row["models"][1]["replays"][0]["rows"]:
        sample["secondary_response"] = [1.1, 0.11]
        sample["secondary_flux_linkage"] = [1.1e-6, 0.11e-6]
    row["models"][1]["replays"][1]["rows"][2]["secondary_response"][0] *= 1.2
    result = magnetic_conductive_shield_frequency_gate(row)
    assert result["status"] == "needs_attention"
    assert result["checks"]["two_replays_are_deterministic_per_model"] is False
    assert result["checks"]["high_frequency_conductive_shield_attenuates_coupling"] is False
    assert result["checks"]["single_gain_to_attenuation_crossover"] is False


def test_rejects_faraday_primary_and_timing_failures():
    row = copy.deepcopy(_summary())
    row["models"][1]["replays"][0]["rows"][0]["maximum_faraday_relative_error"] = 0.1
    row["models"][1]["replays"][0]["rows"][0]["primary_response"] = [2.0, 0.1]
    row["timing_breakdown_s"].pop("verify")
    result = magnetic_conductive_shield_frequency_gate(row)
    assert result["status"] == "needs_attention"
    assert result["checks"]["faraday_identity_is_closed"] is False
    assert result["checks"]["primary_response_is_nearly_invariant"] is False
    assert result["checks"]["exactly_four_timing_stages"] is False
