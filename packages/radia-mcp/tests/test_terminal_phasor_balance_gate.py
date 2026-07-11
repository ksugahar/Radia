import copy
import json
import math

from radia_mcp.radia_ngsolve.server import cyclic_terminal_phasor_balance_gate


def _triplet(magnitude: float, phase_deg: float) -> list[list[float]]:
    return [
        [
            magnitude * math.cos(math.radians(phase_deg - 120.0 * index)),
            magnitude * math.sin(math.radians(phase_deg - 120.0 * index)),
        ]
        for index in range(3)
    ]


def _summary() -> dict:
    return {
        "expected_phase_step_deg": -120.0,
        "voltage_unit": "V",
        "current_unit": "A",
        "groups": [
            {
                "label": "inner",
                "voltage_phasors": _triplet(100.0, 0.0),
                "current_phasors": _triplet(2.0, 90.0),
                "reference_current_magnitude": 2.01,
                "reference_current_unit": "A",
            },
            {
                "label": "return",
                "voltage_phasors": _triplet(10.0, 0.0),
                "current_phasors": _triplet(2.0, -90.0),
            },
        ],
    }


def _gate(summary: dict) -> dict:
    return json.loads(cyclic_terminal_phasor_balance_gate(json.dumps(summary)))


def test_accepts_balanced_cyclic_triplets_and_all_terminal_kcl():
    result = _gate(_summary())
    assert result["status"] == "ok"
    assert result["metrics"]["terminal_count"] == 6
    assert result["metrics"]["all_terminal_kcl_residual"] < 1.0e-15
    assert result["checks"]["group_1_reference_current_matches"] is True


def test_rejects_wrong_phase_sequence_and_open_return_current():
    summary = copy.deepcopy(_summary())
    summary["groups"][0]["voltage_phasors"][1] = _triplet(100.0, 120.0)[0]
    summary["groups"][1]["current_phasors"][2][0] += 0.2
    result = _gate(summary)
    assert result["status"] == "needs_attention"
    assert set(result["issues"]) >= {
        "group_1_voltage_phase_sequence",
        "group_1_voltage_zero_sequence_small",
        "group_2_current_magnitudes_balanced",
        "group_2_current_zero_sequence_small",
        "all_terminal_kcl_closes",
    }


def test_rejects_missing_triplet_instead_of_vacuous_success():
    summary = _summary()
    summary["groups"][1]["current_phasors"].pop()
    result = _gate(summary)
    assert result["status"] == "invalid_input"
    assert "three voltages and currents" in result["error"]


def test_rejects_dimensionally_incompatible_reference_current():
    summary = _summary()
    summary["groups"][0]["reference_current_unit"] = "A/m"
    result = _gate(summary)
    assert result["status"] == "invalid_input"
    assert "unit must match" in result["error"]
