from __future__ import annotations

import math

from radia_mcp.radia_ngsolve.pwm_controlled_motor_loss_gate import (
    pwm_controlled_motor_loss_gate,
)
from test_motor_generalization_v28 import _payload_v28


_PROMOTED_CASE_IDS = (
    "v29_public_iron_loss_hysteresis_eddy_anomalous_frequency_harmonic_volume_mismatch",
    "v29_public_induction_motor_slip_frequency_rotor_current_torque_power_balance_frame_mismatch",
)


def _payload_v29():
    payload = _payload_v28()
    identity = payload["artifact_identity"]
    identity["iron_loss_component_harmonic_frequency_volume_generation_identity"] = {
        "iron_loss_generation": "iron-loss-161",
        "hysteresis_iron_loss_generation": "iron-loss-161",
        "eddy_iron_loss_generation": "iron-loss-161",
        "anomalous_iron_loss_generation": "iron-loss-161",
        "frequency_iron_loss_generation": "iron-loss-161",
        "harmonic_iron_loss_generation": "iron-loss-161",
        "material_iron_loss_generation": "iron-loss-161",
        "volume_iron_loss_generation": "iron-loss-161",
        "mesh_iron_loss_generation": "iron-loss-161",
        "result_iron_loss_generation": "iron-loss-161",
        "frequency_hz": 50.0,
        "result_frequency_hz": 50.0,
        "harmonic_orders": [1, 3, 5],
        "result_harmonic_orders": [1, 3, 5],
        "flux_density_harmonic_t": [1.2, 0.12, 0.06],
        "result_flux_density_harmonic_t": [1.2, 0.12, 0.06],
        "hysteresis_component_w": [80.0, 12.0, 5.0],
        "result_hysteresis_component_w": [80.0, 12.0, 5.0],
        "eddy_component_w": [30.0, 9.0, 5.0],
        "result_eddy_component_w": [30.0, 9.0, 5.0],
        "anomalous_component_w": [10.0, 3.0, 1.0],
        "result_anomalous_component_w": [10.0, 3.0, 1.0],
        "total_iron_loss_w": 155.0,
        "result_total_iron_loss_w": 155.0,
        "element_ids": [101, 102],
        "result_element_ids": [101, 102],
        "element_volumes_m3": [0.001, 0.002],
        "result_element_volumes_m3": [0.001, 0.002],
        "material_coefficients_sha256": "1" * 64,
        "result_material_coefficients_sha256": "1" * 64,
        "mesh_sha256": "2" * 64,
        "result_mesh_sha256": "2" * 64,
        "result_sha256": "3" * 64,
        "accepted_result_sha256": "3" * 64,
    }
    slip = 0.04
    stator_frequency = 50.0
    pole_pairs = 2
    torque = 48.0
    speed = (1.0 - slip) * 2.0 * math.pi * stator_frequency / pole_pairs
    output = torque * speed
    electrical_input = output + 500.0 + 300.0 + 200.0 + 100.0
    identity["induction_slip_rotor_current_torque_power_frame_generation_identity"] = {
        "induction_generation": "induction-balance-161",
        "stator_frequency_induction_generation": "induction-balance-161",
        "slip_induction_generation": "induction-balance-161",
        "rotor_frequency_induction_generation": "induction-balance-161",
        "rotor_current_induction_generation": "induction-balance-161",
        "frame_induction_generation": "induction-balance-161",
        "torque_induction_generation": "induction-balance-161",
        "power_induction_generation": "induction-balance-161",
        "loss_induction_generation": "induction-balance-161",
        "result_induction_generation": "induction-balance-161",
        "stator_frequency_hz": stator_frequency,
        "result_stator_frequency_hz": stator_frequency,
        "slip": slip,
        "result_slip": slip,
        "rotor_frequency_hz": slip * stator_frequency,
        "result_rotor_frequency_hz": slip * stator_frequency,
        "pole_pairs": pole_pairs,
        "result_pole_pairs": pole_pairs,
        "rotor_current_rms_a": [10.0, 10.0, 10.0],
        "result_rotor_current_rms_a": [10.0, 10.0, 10.0],
        "reference_frame": "stator-mechanical-ccw",
        "result_reference_frame": "stator-mechanical-ccw",
        "torque_nm": torque,
        "result_torque_nm": torque,
        "mechanical_speed_rad_s": speed,
        "result_mechanical_speed_rad_s": speed,
        "mechanical_output_w": output,
        "result_mechanical_output_w": output,
        "stator_copper_loss_w": 500.0,
        "result_stator_copper_loss_w": 500.0,
        "rotor_copper_loss_w": 300.0,
        "result_rotor_copper_loss_w": 300.0,
        "iron_loss_w": 200.0,
        "result_iron_loss_w": 200.0,
        "mechanical_loss_w": 100.0,
        "result_mechanical_loss_w": 100.0,
        "electrical_input_w": electrical_input,
        "result_electrical_input_w": electrical_input,
        "result_sha256": "4" * 64,
        "accepted_result_sha256": "4" * 64,
    }
    return payload


def test_v29_public_positive_iron_loss_and_induction_balance_identities():
    assert pwm_controlled_motor_loss_gate(_payload_v29())["status"] == "ok"


def test_v29_public_iron_loss_hysteresis_eddy_anomalous_frequency_harmonic_volume_mismatch():
    payload = _payload_v29()
    identity = payload["artifact_identity"][
        "iron_loss_component_harmonic_frequency_volume_generation_identity"
    ]
    identity.update(
        {
            "eddy_iron_loss_generation": "iron-loss-160",
            "volume_iron_loss_generation": "iron-loss-159",
            "result_frequency_hz": 100.0,
            "result_harmonic_orders": [1, 5, 7],
            "result_flux_density_harmonic_t": [1.1, 0.2, 0.1],
            "result_hysteresis_component_w": [70.0, 10.0, 4.0],
            "result_eddy_component_w": [60.0, 18.0, 10.0],
            "result_anomalous_component_w": [14.0, 5.0, 2.0],
            "result_total_iron_loss_w": 193.0,
            "result_element_ids": [102, 101],
            "result_element_volumes_m3": [0.002, 0.001],
            "result_material_coefficients_sha256": "8" * 64,
            "result_mesh_sha256": "9" * 64,
            "accepted_result_sha256": "a" * 64,
        }
    )
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "iron_loss_components_use_current_frequency_harmonics_material_volumes_mesh_and_result"
    ]


def test_v29_public_induction_motor_slip_frequency_rotor_current_torque_power_balance_frame_mismatch():
    payload = _payload_v29()
    identity = payload["artifact_identity"][
        "induction_slip_rotor_current_torque_power_frame_generation_identity"
    ]
    identity.update(
        {
            "slip_induction_generation": "induction-balance-160",
            "frame_induction_generation": "induction-balance-159",
            "result_stator_frequency_hz": 60.0,
            "result_slip": -0.04,
            "result_rotor_frequency_hz": 6.0,
            "result_rotor_current_rms_a": [8.0, 9.0, 10.0],
            "result_reference_frame": "rotor-electrical-clockwise",
            "result_torque_nm": 41.0,
            "result_mechanical_speed_rad_s": 170.0,
            "result_mechanical_output_w": 6000.0,
            "result_rotor_copper_loss_w": 900.0,
            "result_electrical_input_w": 7000.0,
            "accepted_result_sha256": "b" * 64,
        }
    )
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "induction_motor_uses_current_slip_rotor_frequency_current_frame_torque_and_power_balance"
    ]
