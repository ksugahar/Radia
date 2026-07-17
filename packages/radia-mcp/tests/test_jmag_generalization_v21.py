from __future__ import annotations

from radia_mcp.radia_ngsolve.pwm_controlled_motor_loss_gate import (
    pwm_controlled_motor_loss_gate,
)
from test_jmag_generalization_v20 import _payload_v20


def _payload_v21():
    payload = _payload_v20()
    identity = payload["artifact_identity"]
    identity["dq_transform_rotor_angle_phase_order_generation_identity"] = {
        "operating_point_generation": "operating-point-31",
        "rotor_angle_operating_point_generation": "operating-point-31",
        "electrical_offset_operating_point_generation": "operating-point-31",
        "phase_order_operating_point_generation": "operating-point-31",
        "dq_result_operating_point_generation": "operating-point-31",
        "rotor_mechanical_angle_deg": 15.0,
        "dq_rotor_mechanical_angle_deg": 15.0,
        "pole_pairs": 4,
        "dq_pole_pairs": 4,
        "electrical_offset_deg": 30.0,
        "dq_electrical_offset_deg": 30.0,
        "phase_order": ["u", "v", "w"],
        "dq_phase_order": ["u", "v", "w"],
        "phase_values": [10.0, -5.0, -5.0],
        "dq_source_phase_values": [10.0, -5.0, -5.0],
        "dq_transform_table_sha256": "1" * 64,
        "resolved_dq_transform_table_sha256": "1" * 64,
    }
    identity["iron_loss_frequency_harmonic_material_curve_generation_identity"] = {
        "loss_study_generation": "iron-loss-31",
        "frequency_loss_study_generation": "iron-loss-31",
        "harmonic_spectrum_loss_study_generation": "iron-loss-31",
        "material_curve_loss_study_generation": "iron-loss-31",
        "loss_result_study_generation": "iron-loss-31",
        "fundamental_frequency_hz": 400.0,
        "loss_frequency_hz": 400.0,
        "harmonic_orders": [1, 3, 5, 7],
        "loss_harmonic_orders": [1, 3, 5, 7],
        "harmonic_amplitudes_t": [1.0, 0.12, 0.05, 0.02],
        "loss_harmonic_amplitudes_t": [1.0, 0.12, 0.05, 0.02],
        "material_curve_ids": ["stator-r3", "rotor-r2"],
        "loss_material_curve_ids": ["stator-r3", "rotor-r2"],
        "loss_input_table_sha256": "2" * 64,
        "resolved_loss_input_table_sha256": "2" * 64,
    }
    return payload


def test_v21_public_positive_dq_and_iron_loss_identity():
    result = pwm_controlled_motor_loss_gate(_payload_v21())
    assert result["status"] == "ok"
    assert result["checks"][
        "dq_transform_uses_current_rotor_angle_offset_and_phase_order"
    ]
    assert result["checks"][
        "iron_loss_uses_current_frequency_harmonics_and_material_curves"
    ]


def test_v21_public_dq_transform_rotor_angle_phase_order_generation_mismatch():
    payload = _payload_v21()
    payload["artifact_identity"][
        "dq_transform_rotor_angle_phase_order_generation_identity"
    ].update(
        {
            "rotor_angle_operating_point_generation": "operating-point-30",
            "electrical_offset_operating_point_generation": "operating-point-29",
            "phase_order_operating_point_generation": "operating-point-28",
            "dq_rotor_mechanical_angle_deg": 10.0,
            "dq_electrical_offset_deg": -30.0,
            "dq_phase_order": ["u", "w", "v"],
            "dq_source_phase_values": [10.0, -5.0, -4.5],
            "resolved_dq_transform_table_sha256": "a" * 64,
        }
    )
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert result["checks"][
        "dq_transform_uses_current_rotor_angle_offset_and_phase_order"
    ] is False


def test_v21_public_iron_loss_frequency_harmonic_material_curve_generation_mismatch():
    payload = _payload_v21()
    payload["artifact_identity"][
        "iron_loss_frequency_harmonic_material_curve_generation_identity"
    ].update(
        {
            "frequency_loss_study_generation": "iron-loss-30",
            "harmonic_spectrum_loss_study_generation": "iron-loss-29",
            "material_curve_loss_study_generation": "iron-loss-28",
            "loss_frequency_hz": 60.0,
            "loss_harmonic_orders": [1, 5, 7],
            "loss_harmonic_amplitudes_t": [1.0, 0.05, 0.02],
            "loss_material_curve_ids": ["stator-r2", "rotor-r1"],
            "resolved_loss_input_table_sha256": "b" * 64,
        }
    )
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert result["checks"][
        "iron_loss_uses_current_frequency_harmonics_and_material_curves"
    ] is False
