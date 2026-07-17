from __future__ import annotations

from radia_mcp.radia_ngsolve.pwm_controlled_motor_loss_gate import (
    pwm_controlled_motor_loss_gate,
)
from test_motor_generalization_v24 import _payload_v24


def _payload_v25():
    payload = _payload_v24()
    identity = payload["artifact_identity"]
    identity[
        "torque_map_current_angle_temperature_speed_interpolation_generation_identity"
    ] = {
        "map_generation": "torque-map-111",
        "current_map_generation": "torque-map-111",
        "angle_map_generation": "torque-map-111",
        "temperature_map_generation": "torque-map-111",
        "speed_map_generation": "torque-map-111",
        "interpolation_map_generation": "torque-map-111",
        "query_map_generation": "torque-map-111",
        "result_map_generation": "torque-map-111",
        "current_axis_a": [0.0, 5.0, 10.0],
        "result_current_axis_a": [0.0, 5.0, 10.0],
        "electrical_angle_axis_deg": [0.0, 30.0, 60.0],
        "result_electrical_angle_axis_deg": [0.0, 30.0, 60.0],
        "temperature_axis_c": [20.0, 80.0],
        "result_temperature_axis_c": [20.0, 80.0],
        "speed_axis_rpm": [1000.0, 3000.0],
        "result_speed_axis_rpm": [1000.0, 3000.0],
        "angle_period_deg": 360.0,
        "result_angle_period_deg": 360.0,
        "interpolation_method": "multilinear_periodic_angle",
        "result_interpolation_method": "multilinear_periodic_angle",
        "torque_tensor_sha256": "1" * 64,
        "result_torque_tensor_sha256": "1" * 64,
        "query_point": [7.5, 45.0, 50.0, 2000.0],
        "result_query_point": [7.5, 45.0, 50.0, 2000.0],
        "interpolated_torque_nm": 1.35,
        "result_interpolated_torque_nm": 1.35,
        "result_sha256": "2" * 64,
        "accepted_result_sha256": "2" * 64,
    }
    identity[
        "demagnetization_margin_operating_point_temperature_recoil_generation_identity"
    ] = {
        "demag_generation": "demag-111",
        "material_demag_generation": "demag-111",
        "temperature_demag_generation": "demag-111",
        "recoil_demag_generation": "demag-111",
        "operating_point_demag_generation": "demag-111",
        "margin_demag_generation": "demag-111",
        "result_demag_generation": "demag-111",
        "magnet_ids": ["pm-1", "pm-2"],
        "result_magnet_ids": ["pm-1", "pm-2"],
        "temperature_c": 120.0,
        "result_temperature_c": 120.0,
        "coercivity_a_m": 720000.0,
        "result_coercivity_a_m": 720000.0,
        "recoil_relative_permeability": 1.05,
        "result_recoil_relative_permeability": 1.05,
        "minimum_operating_h_a_m": -510000.0,
        "result_minimum_operating_h_a_m": -510000.0,
        "demagnetization_margin_a_m": 210000.0,
        "result_demagnetization_margin_a_m": 210000.0,
        "material_curve_sha256": "3" * 64,
        "result_material_curve_sha256": "3" * 64,
        "operating_point_field_sha256": "4" * 64,
        "result_operating_point_field_sha256": "4" * 64,
        "result_sha256": "5" * 64,
        "accepted_result_sha256": "5" * 64,
    }
    return payload


def test_v25_public_positive_torque_map_and_demagnetization_margin_identity():
    assert pwm_controlled_motor_loss_gate(_payload_v25())["status"] == "ok"


def test_v25_public_nonlinear_torque_map_current_angle_temperature_speed_interpolation_mismatch():
    payload = _payload_v25()
    payload["artifact_identity"][
        "torque_map_current_angle_temperature_speed_interpolation_generation_identity"
    ].update(
        {
            "current_map_generation": "torque-map-110",
            "angle_map_generation": "torque-map-109",
            "temperature_map_generation": "torque-map-108",
            "speed_map_generation": "torque-map-107",
            "interpolation_map_generation": "torque-map-106",
            "query_map_generation": "torque-map-105",
            "result_current_axis_a": [0.0, 10.0, 5.0],
            "result_electrical_angle_axis_deg": [60.0, 30.0, 0.0],
            "result_temperature_axis_c": [20.0, 100.0],
            "result_speed_axis_rpm": [1000.0, 6000.0],
            "result_angle_period_deg": 180.0,
            "result_interpolation_method": "nearest",
            "result_torque_tensor_sha256": "a" * 64,
            "result_query_point": [7.5, 45.0, 90.0, 4000.0],
            "result_interpolated_torque_nm": 0.85,
            "accepted_result_sha256": "b" * 64,
        }
    )
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "torque_map_uses_current_axes_interpolation_query_and_result_generation"
    ]


def test_v25_public_demagnetization_margin_operating_point_temperature_recoil_generation_mismatch():
    payload = _payload_v25()
    payload["artifact_identity"][
        "demagnetization_margin_operating_point_temperature_recoil_generation_identity"
    ].update(
        {
            "material_demag_generation": "demag-110",
            "temperature_demag_generation": "demag-109",
            "recoil_demag_generation": "demag-108",
            "operating_point_demag_generation": "demag-107",
            "margin_demag_generation": "demag-106",
            "result_magnet_ids": ["pm-1", "pm-old"],
            "result_temperature_c": 20.0,
            "result_coercivity_a_m": 900000.0,
            "result_recoil_relative_permeability": 1.2,
            "result_minimum_operating_h_a_m": -850000.0,
            "result_demagnetization_margin_a_m": -130000.0,
            "result_material_curve_sha256": "c" * 64,
            "result_operating_point_field_sha256": "d" * 64,
            "accepted_result_sha256": "e" * 64,
        }
    )
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "demagnetization_margin_uses_current_temperature_recoil_material_and_operating_point"
    ]
