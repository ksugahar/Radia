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
