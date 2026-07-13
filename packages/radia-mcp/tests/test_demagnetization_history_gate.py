import json

from radia_mcp.motor.demagnetization_history_gate import permanent_magnet_demagnetization_history_gate
from radia_mcp.motor.server import motor_permanent_magnet_demagnetization_history_gate


def _summary():
    return {
        "state_unit": "fraction_of_reference_remanence",
        "field_unit": "T",
        "stress_step_index": 1,
        "recovery_step_index": 3,
        "steps": [
            {"step": 0, "peak_flux_density_T": 0.58, "magnet_state_fraction": [1.0, 0.96, 1.0]},
            {"step": 1, "peak_flux_density_T": 0.46, "magnet_state_fraction": [0.83, 0.74, 0.99]},
            {"step": 2, "peak_flux_density_T": 0.65, "magnet_state_fraction": [0.83, 0.74, 0.99]},
            {"step": 3, "peak_flux_density_T": 0.57, "magnet_state_fraction": [0.83, 0.74, 0.99]},
        ],
    }


def _family_summary():
    states = [
        {
            "state_id": "initial",
            "step": 0,
            "element_ids": [1, 2],
            "material_ids": [7, 7],
            "remanence_ratio": [1.0, 0.98],
            "flux_density_T": [0.52, 0.49],
        },
        {
            "state_id": "reverse_field",
            "step": 1,
            "element_ids": [1, 2],
            "material_ids": [7, 7],
            "remanence_ratio": [0.81, 0.76],
            "flux_density_T": [0.31, 0.28],
        },
        {
            "state_id": "first_unloaded",
            "step": 2,
            "element_ids": [1, 2],
            "material_ids": [7, 7],
            "remanence_ratio": [0.81, 0.76],
            "flux_density_T": [0.44, 0.40],
        },
        {
            "state_id": "alternate_curve",
            "step": 3,
            "element_ids": [1, 2],
            "material_ids": [7, 7],
            "remanence_ratio": [0.81, 0.76],
            "flux_density_T": [0.50, 0.46],
        },
        {
            "state_id": "final_unloaded",
            "step": 4,
            "element_ids": [1, 2],
            "material_ids": [7, 7],
            "remanence_ratio": [0.81, 0.76],
            "flux_density_T": [0.44, 0.40],
        },
    ]
    return {
        "schema": "permanent-magnet-demagnetization-history/v1",
        "result_authority": ".mao",
        "state_variable": "elementwise_remanence_ratio",
        "instantaneous_field_observable": "elementwise_flux_density_magnitude_T",
        "cases": [
            {
                "case_id": "precondition_then_temperature_cycle",
                "states": states,
                "history_blocks": [
                    {
                        "pre_state": "initial",
                        "stress_state": "reverse_field",
                        "unloaded_state": "first_unloaded",
                        "expect_additional_demagnetization": True,
                    },
                    {
                        "pre_state": "first_unloaded",
                        "stress_state": "alternate_curve",
                        "unloaded_state": "final_unloaded",
                        "expect_additional_demagnetization": False,
                    },
                ],
                "replay_count": 2,
                "replay_max_abs": 0.0,
            }
        ],
    }


def test_demagnetization_history_accepts_irreversible_load_memory():
    result = permanent_magnet_demagnetization_history_gate(json.dumps(_summary()))
    assert result["status"] == "ok"
    assert result["metrics"]["irreversibly_damaged_partition_count"] == 3
    assert json.loads(motor_permanent_magnet_demagnetization_history_gate(json.dumps(_summary())))["status"] == "ok"


def test_demagnetization_history_rejects_spontaneous_recovery():
    summary = _summary()
    summary["steps"][3]["magnet_state_fraction"][0] = 1.0
    result = permanent_magnet_demagnetization_history_gate(json.dumps(summary))
    assert result["status"] == "needs_attention"
    assert result["checks"]["state_never_spontaneously_recovers"] is False
    assert result["checks"]["recovery_retains_stressed_state"] is False


def test_demagnetization_history_rejects_final_field_only_without_damage():
    summary = _summary()
    for step in summary["steps"]:
        step["magnet_state_fraction"] = [1.0, 1.0, 1.0]
    result = permanent_magnet_demagnetization_history_gate(json.dumps(summary))
    assert result["status"] == "needs_attention"
    assert result["checks"]["stress_causes_resolved_damage"] is False


def test_demagnetization_family_accepts_damage_then_no_additional_damage():
    summary = _family_summary()
    result = permanent_magnet_demagnetization_history_gate(json.dumps(summary))
    assert result["status"] == "ok"
    assert result["policy"] == "permanent_magnet_demagnetization_history_gate_v2"
    assert result["metrics"]["damage_block_count"] == 1
    assert result["metrics"]["no_additional_damage_block_count"] == 1
    assert json.loads(
        motor_permanent_magnet_demagnetization_history_gate(json.dumps(summary))
    )["status"] == "ok"


def test_demagnetization_family_rejects_state_healing_after_unload():
    summary = _family_summary()
    summary["cases"][0]["states"][2]["remanence_ratio"][0] = 0.90
    result = permanent_magnet_demagnetization_history_gate(json.dumps(summary))
    checks = result["case_checks"]["precondition_then_temperature_cycle"]
    assert result["status"] == "needs_attention"
    assert checks["state_never_spontaneously_recovers"] is False
    assert checks["unload_preserves_stressed_remanence_state"] is False


def test_demagnetization_family_rejects_unexpected_second_block_damage():
    summary = _family_summary()
    summary["cases"][0]["states"][3]["remanence_ratio"][0] = 0.70
    summary["cases"][0]["states"][4]["remanence_ratio"][0] = 0.70
    result = permanent_magnet_demagnetization_history_gate(json.dumps(summary))
    checks = result["case_checks"]["precondition_then_temperature_cycle"]
    assert result["status"] == "needs_attention"
    assert checks["damage_expectation_matches_each_block"] is False


def test_demagnetization_family_rejects_stale_replay():
    summary = _family_summary()
    summary["cases"][0]["replay_max_abs"] = 1.0e-4
    result = permanent_magnet_demagnetization_history_gate(json.dumps(summary))
    checks = result["case_checks"]["precondition_then_temperature_cycle"]
    assert result["status"] == "needs_attention"
    assert checks["replays_match_within_state_tolerance"] is False
