from __future__ import annotations

from radia_mcp.radia_ngsolve.pwm_controlled_motor_loss_gate import (
    pwm_controlled_motor_loss_gate,
)
from test_motor_generalization_v26 import _payload_v26


def _payload_v27():
    payload = _payload_v26()
    identity = payload["artifact_identity"]
    identity["rotating_sector_pole_pair_periodic_phase_skew_slice_torque_frame_generation_identity"] = {
        "sector_generation": "sector-141",
        "pole_pair_sector_generation": "sector-141",
        "periodic_sector_generation": "sector-141",
        "skew_sector_generation": "sector-141",
        "rotor_frame_sector_generation": "sector-141",
        "torque_sector_generation": "sector-141",
        "mesh_sector_generation": "sector-141",
        "result_sector_generation": "sector-141",
        "pole_pairs": 4,
        "result_pole_pairs": 4,
        "sector_angle_deg": 45.0,
        "result_sector_angle_deg": 45.0,
        "periodic_phase_deg": 180.0,
        "result_periodic_phase_deg": 180.0,
        "periodic_pair_ids": [[101, 201], [102, 202]],
        "result_periodic_pair_ids": [[101, 201], [102, 202]],
        "periodic_pair_orientation": [1, -1],
        "result_periodic_pair_orientation": [1, -1],
        "skew_slice_deg": [-2.0, 0.0, 2.0],
        "result_skew_slice_deg": [-2.0, 0.0, 2.0],
        "skew_slice_weights": [0.25, 0.5, 0.25],
        "result_skew_slice_weights": [0.25, 0.5, 0.25],
        "rotor_mechanical_angle_deg": [0.0, 1.0, 2.0],
        "result_rotor_mechanical_angle_deg": [0.0, 1.0, 2.0],
        "torque_frame": "rotor-mechanical-ccw",
        "result_torque_frame": "rotor-mechanical-ccw",
        "slice_torque_nm": [1.0, 1.2, 1.1],
        "result_slice_torque_nm": [1.0, 1.2, 1.1],
        "torque_average_nm": 1.125,
        "result_torque_average_nm": 1.125,
        "sector_mesh_sha256": "a" * 64,
        "result_sector_mesh_sha256": "a" * 64,
        "torque_result_sha256": "b" * 64,
        "accepted_torque_result_sha256": "b" * 64,
    }
    components = {
        "hysteresis": [3.0, 1.0, 0.5],
        "eddy": [1.0, 0.5, 0.25],
        "excess": [0.3, 0.1, 0.05],
    }
    identity["iron_loss_harmonic_decomposition_model_temperature_frequency_element_volume_result_generation_identity"] = {
        "decomposition_generation": "iron-decomposition-141",
        "model_decomposition_generation": "iron-decomposition-141",
        "temperature_decomposition_generation": "iron-decomposition-141",
        "frequency_decomposition_generation": "iron-decomposition-141",
        "volume_decomposition_generation": "iron-decomposition-141",
        "material_decomposition_generation": "iron-decomposition-141",
        "result_decomposition_generation": "iron-decomposition-141",
        "loss_model": "bertotti-three-term",
        "result_loss_model": "bertotti-three-term",
        "material_temperature_c": 120.0,
        "result_material_temperature_c": 120.0,
        "harmonic_orders": [1, 3, 5],
        "result_harmonic_orders": [1, 3, 5],
        "frequency_hz": [50.0, 150.0, 250.0],
        "result_frequency_hz": [50.0, 150.0, 250.0],
        "harmonic_loss_w": components,
        "result_harmonic_loss_w": components,
        "element_ids": [11, 12, 13],
        "result_element_ids": [11, 12, 13],
        "element_volume_m3": [0.0005, 0.0007, 0.0008],
        "result_element_volume_m3": [0.0005, 0.0007, 0.0008],
        "integration_volume_m3": 0.002,
        "result_integration_volume_m3": 0.002,
        "material_state_sha256": "c" * 64,
        "result_material_state_sha256": "c" * 64,
        "mesh_sha256": "d" * 64,
        "result_mesh_sha256": "d" * 64,
        "loss_result_sha256": "e" * 64,
        "accepted_loss_result_sha256": "e" * 64,
    }
    return payload


def test_v27_public_positive_rotating_sector_and_iron_loss_decomposition_identity():
    assert pwm_controlled_motor_loss_gate(_payload_v27())["status"] == "ok"


def test_v27_public_rotating_sector_pole_pair_periodic_phase_skew_slice_torque_frame_mismatch():
    payload = _payload_v27()
    payload["artifact_identity"][
        "rotating_sector_pole_pair_periodic_phase_skew_slice_torque_frame_generation_identity"
    ].update({
        "pole_pair_sector_generation": "sector-140",
        "periodic_sector_generation": "sector-139",
        "mesh_sector_generation": "sector-138",
        "result_pole_pairs": 2,
        "result_sector_angle_deg": 90.0,
        "result_periodic_phase_deg": 0.0,
        "result_periodic_pair_ids": [[101, 202], [102, 201]],
        "result_periodic_pair_orientation": [1, 1],
        "result_skew_slice_deg": [2.0, 0.0, -2.0],
        "result_rotor_mechanical_angle_deg": [0.0, 2.0, 4.0],
        "result_torque_frame": "stator-clockwise",
        "result_torque_average_nm": 1.7,
        "result_sector_mesh_sha256": "4" * 64,
        "accepted_torque_result_sha256": "5" * 64,
    })
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "rotating_sector_torque_uses_current_pole_pairs_periodic_phase_skew_frame_mesh_and_result"
    ]


def test_v27_public_iron_loss_harmonic_decomposition_material_temperature_frequency_volume_mismatch():
    payload = _payload_v27()
    payload["artifact_identity"][
        "iron_loss_harmonic_decomposition_model_temperature_frequency_element_volume_result_generation_identity"
    ].update({
        "model_decomposition_generation": "iron-decomposition-140",
        "temperature_decomposition_generation": "iron-decomposition-139",
        "volume_decomposition_generation": "iron-decomposition-138",
        "result_loss_model": "two-term",
        "result_material_temperature_c": 20.0,
        "result_harmonic_orders": [1, 2, 3],
        "result_frequency_hz": [50.0, 100.0, 150.0],
        "result_harmonic_loss_w": {"hysteresis": [3.0], "eddy": [2.0], "excess": []},
        "result_element_ids": [13, 12, 11],
        "result_element_volume_m3": [0.0008, 0.0007, 0.0004],
        "result_integration_volume_m3": 0.0019,
        "result_material_state_sha256": "6" * 64,
        "result_mesh_sha256": "7" * 64,
        "accepted_loss_result_sha256": "8" * 64,
    })
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "iron_loss_decomposition_uses_current_model_temperature_frequency_elements_volume_and_result"
    ]
