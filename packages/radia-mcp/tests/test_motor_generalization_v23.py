from __future__ import annotations

from radia_mcp.radia_ngsolve.pwm_controlled_motor_loss_gate import (
    pwm_controlled_motor_loss_gate,
)
from test_motor_generalization_v22 import _payload_v22


def _payload_v23():
    payload = _payload_v22()
    identity = payload["artifact_identity"]
    identity["winding_current_phase_circuit_sequence_torque_generation_identity"] = {
        "motor_sweep_generation": "motor-sweep-51",
        "winding_order_motor_sweep_generation": "motor-sweep-51",
        "phase_convention_motor_sweep_generation": "motor-sweep-51",
        "circuit_sequence_motor_sweep_generation": "motor-sweep-51",
        "rotor_angle_motor_sweep_generation": "motor-sweep-51",
        "torque_result_motor_sweep_generation": "motor-sweep-51",
        "winding_order": ["u", "v", "w"],
        "torque_winding_order": ["u", "v", "w"],
        "current_phase_convention": "abc_positive_sequence",
        "torque_current_phase_convention": "abc_positive_sequence",
        "circuit_sequence_ids": [101, 102, 103],
        "torque_circuit_sequence_ids": [101, 102, 103],
        "rotor_angles_deg": [0.0, 5.0, 10.0, 15.0],
        "torque_rotor_angles_deg": [0.0, 5.0, 10.0, 15.0],
        "phase_current_table_sha256": "1" * 64,
        "torque_phase_current_table_sha256": "1" * 64,
        "torque_nm": [1.0, 1.2, 0.9, 1.1],
        "reported_torque_nm": [1.0, 1.2, 0.9, 1.1],
        "torque_table_sha256": "2" * 64,
        "reported_torque_table_sha256": "2" * 64,
    }
    identity[
        "demagnetization_knee_temperature_recoil_operating_generation_identity"
    ] = {
        "demag_generation": "demag-51",
        "knee_curve_demag_generation": "demag-51",
        "temperature_demag_generation": "demag-51",
        "recoil_line_demag_generation": "demag-51",
        "operating_state_demag_generation": "demag-51",
        "margin_result_demag_generation": "demag-51",
        "knee_curve_sha256": "3" * 64,
        "margin_knee_curve_sha256": "3" * 64,
        "temperature_c": 140.0,
        "margin_temperature_c": 140.0,
        "recoil_line_sha256": "4" * 64,
        "margin_recoil_line_sha256": "4" * 64,
        "operating_state_id": "id=-220A;iq=260A;theta=20deg",
        "margin_operating_state_id": "id=-220A;iq=260A;theta=20deg",
        "demag_margin_a_per_m": [120000.0, 65000.0],
        "reported_demag_margin_a_per_m": [120000.0, 65000.0],
        "demag_state_sha256": "5" * 64,
        "reported_demag_state_sha256": "5" * 64,
    }
    return payload


def test_v23_public_positive_winding_torque_and_demag_knee_identity():
    assert pwm_controlled_motor_loss_gate(_payload_v23())["status"] == "ok"


def test_v23_public_winding_current_phase_convention_circuit_sequence_torque_generation_mismatch():
    payload = _payload_v23()
    payload["artifact_identity"][
        "winding_current_phase_circuit_sequence_torque_generation_identity"
    ].update(
        {
            "winding_order_motor_sweep_generation": "motor-sweep-50",
            "phase_convention_motor_sweep_generation": "motor-sweep-49",
            "circuit_sequence_motor_sweep_generation": "motor-sweep-48",
            "rotor_angle_motor_sweep_generation": "motor-sweep-47",
            "torque_result_motor_sweep_generation": "motor-sweep-46",
            "torque_winding_order": ["u", "w", "v"],
            "torque_current_phase_convention": "acb_negative_sequence",
            "torque_circuit_sequence_ids": [101, 103, 102],
            "torque_rotor_angles_deg": [0.0, 10.0, 5.0, 15.0],
            "torque_phase_current_table_sha256": "c" * 64,
            "reported_torque_nm": [1.0, 0.9, 1.2, 1.1],
            "reported_torque_table_sha256": "d" * 64,
        }
    )
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "motor_torque_uses_current_winding_phase_circuit_sequence_and_angles"
    ]


def test_v23_public_demagnetization_knee_temperature_operating_state_generation_mismatch():
    payload = _payload_v23()
    payload["artifact_identity"][
        "demagnetization_knee_temperature_recoil_operating_generation_identity"
    ].update(
        {
            "knee_curve_demag_generation": "demag-50",
            "temperature_demag_generation": "demag-49",
            "recoil_line_demag_generation": "demag-48",
            "operating_state_demag_generation": "demag-47",
            "margin_result_demag_generation": "demag-46",
            "margin_knee_curve_sha256": "e" * 64,
            "margin_temperature_c": 80.0,
            "margin_recoil_line_sha256": "f" * 64,
            "margin_operating_state_id": "id=-100A;iq=120A;theta=5deg",
            "reported_demag_margin_a_per_m": [250000.0, 190000.0],
            "reported_demag_state_sha256": "0" * 64,
        }
    )
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "demag_margin_uses_current_knee_temperature_recoil_and_operating_state"
    ]
