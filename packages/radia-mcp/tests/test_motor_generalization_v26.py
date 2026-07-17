from __future__ import annotations

from radia_mcp.radia_ngsolve.pwm_controlled_motor_loss_gate import pwm_controlled_motor_loss_gate
from test_motor_generalization_v25 import _payload_v25


def _payload_v26():
    payload = _payload_v25()
    identity = payload["artifact_identity"]
    identity["iron_loss_hysteresis_eddy_excess_harmonic_frequency_material_volume_generation_identity"] = {
        "loss_generation": "loss-131", "component_loss_generation": "loss-131",
        "harmonic_loss_generation": "loss-131", "frequency_loss_generation": "loss-131",
        "material_loss_generation": "loss-131", "volume_loss_generation": "loss-131",
        "result_loss_generation": "loss-131", "hysteresis_loss_w": 12.0,
        "result_hysteresis_loss_w": 12.0, "eddy_loss_w": 5.0, "result_eddy_loss_w": 5.0,
        "excess_loss_w": 1.5, "result_excess_loss_w": 1.5,
        "total_iron_loss_w": 18.5, "result_total_iron_loss_w": 18.5,
        "harmonic_orders": [1, 3, 5], "result_harmonic_orders": [1, 3, 5],
        "harmonic_frequencies_hz": [50.0, 150.0, 250.0],
        "result_harmonic_frequencies_hz": [50.0, 150.0, 250.0],
        "material_law_sha256": "1" * 64, "result_material_law_sha256": "1" * 64,
        "integration_volume_m3": 0.002, "result_integration_volume_m3": 0.002,
        "mesh_sha256": "2" * 64, "result_mesh_sha256": "2" * 64,
        "result_sha256": "3" * 64, "accepted_result_sha256": "3" * 64,
    }
    identity["skew_slice_torque_phase_angle_weight_periodicity_mesh_generation_identity"] = {
        "skew_generation": "skew-131", "phase_skew_generation": "skew-131",
        "angle_skew_generation": "skew-131", "weight_skew_generation": "skew-131",
        "periodicity_skew_generation": "skew-131", "mesh_skew_generation": "skew-131",
        "result_skew_generation": "skew-131", "slice_phase_deg": [-10.0, 0.0, 10.0],
        "result_slice_phase_deg": [-10.0, 0.0, 10.0], "mechanical_angle_deg": [0.0, 1.0, 2.0],
        "result_mechanical_angle_deg": [0.0, 1.0, 2.0], "slice_weights": [0.25, 0.5, 0.25],
        "result_slice_weights": [0.25, 0.5, 0.25], "periodicity": 8, "result_periodicity": 8,
        "slice_mesh_sha256": ["4" * 64, "5" * 64, "6" * 64],
        "result_slice_mesh_sha256": ["4" * 64, "5" * 64, "6" * 64],
        "slice_torque_nm": [1.0, 1.2, 1.1], "result_slice_torque_nm": [1.0, 1.2, 1.1],
        "skew_averaged_torque_nm": 1.125, "result_skew_averaged_torque_nm": 1.125,
        "result_sha256": "7" * 64, "accepted_result_sha256": "7" * 64,
    }
    return payload


def test_v26_public_positive_iron_loss_and_skew_identity():
    assert pwm_controlled_motor_loss_gate(_payload_v26())["status"] == "ok"


def test_v26_public_iron_loss_hysteresis_eddy_excess_harmonic_frequency_material_volume_mismatch():
    payload = _payload_v26()
    payload["artifact_identity"]["iron_loss_hysteresis_eddy_excess_harmonic_frequency_material_volume_generation_identity"].update({
        "component_loss_generation": "loss-130", "harmonic_loss_generation": "loss-129",
        "result_hysteresis_loss_w": 8.0, "result_eddy_loss_w": 9.0,
        "result_total_iron_loss_w": 17.0, "result_harmonic_orders": [1, 2, 3],
        "result_harmonic_frequencies_hz": [50.0, 100.0, 150.0],
        "result_material_law_sha256": "a" * 64, "result_integration_volume_m3": 0.001,
    })
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["iron_loss_uses_current_components_harmonics_frequency_material_volume_and_result"]


def test_v26_public_skew_slice_torque_phase_angle_weight_periodicity_mesh_generation_mismatch():
    payload = _payload_v26()
    payload["artifact_identity"]["skew_slice_torque_phase_angle_weight_periodicity_mesh_generation_identity"].update({
        "phase_skew_generation": "skew-130", "angle_skew_generation": "skew-129",
        "result_slice_phase_deg": [10.0, 0.0, -10.0],
        "result_mechanical_angle_deg": [0.0, 2.0, 4.0], "result_slice_weights": [0.5, 0.5, 0.5],
        "result_periodicity": 4, "result_slice_mesh_sha256": ["4" * 64, "a" * 64, "6" * 64],
        "result_slice_torque_nm": [0.8, 1.2, 1.4], "result_skew_averaged_torque_nm": 1.7,
    })
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["skew_torque_uses_current_slice_phases_angles_weights_periodicity_meshes_and_result"]
