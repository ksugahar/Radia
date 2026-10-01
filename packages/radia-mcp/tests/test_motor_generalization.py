"""Motor generalization gate contracts (formerly test_motor_generalization_v19..v56).

Positive closure is asserted once per gate: the v43 positive covers the whole
v19..v43 ``pwm_controlled_motor_loss_gate`` chain; each v45..v56 identity gate keeps
its own positive.  Every negative failure signal is kept verbatim.
"""

from __future__ import annotations

from copy import deepcopy
import math

from radia_mcp.radia_ngsolve.motor_artifact_identity_v49 import (
    DEMAG as DEMAG_V49,
    IRON as IRON_V49,
    validate_public_identity as validate_public_identity_v49,
)
from radia_mcp.radia_ngsolve.motor_artifact_identity_v50 import (
    DQ as DQ_V50,
    THERMAL,
    validate_public_identity as validate_public_identity_v50,
)
from radia_mcp.radia_ngsolve.motor_artifact_identity_v51 import (
    TORQUE as TORQUE_V51,
    WINDING,
    validate_public_identity as validate_public_identity_v51,
)
from radia_mcp.radia_ngsolve.motor_artifact_identity_v52 import (
    COGGING,
    IRON_LOSS,
    validate_public_identity as validate_public_identity_v52,
)
from radia_mcp.radia_ngsolve.motor_artifact_identity_v53 import (
    DEMAG as DEMAG_V53,
    SKEW as SKEW_V53,
    validate_public_identity as validate_public_identity_v53,
)
from radia_mcp.radia_ngsolve.motor_artifact_identity_v54 import (
    DEMAG as DEMAG_V54,
    TORQUE as TORQUE_V54,
    validate_public_identity as validate_public_identity_v54,
)
from radia_mcp.radia_ngsolve.motor_artifact_identity_v55 import (
    DQ as DQ_V55,
    IRON as IRON_V55,
    validate_public_identity as validate_public_identity_v55,
)
from radia_mcp.radia_ngsolve.motor_artifact_lineage_v47 import (
    DQ as DQ_V47,
    WINDOW,
    validate_public_identity as validate_public_identity_v47,
)
from radia_mcp.radia_ngsolve.motor_map_induction_identity_v56 import (
    INDUCTION,
    MAP,
    validate_public_identity as validate_public_identity_v56,
)
from radia_mcp.radia_ngsolve.motor_semantic_identity_v48 import (
    PWM,
    SKEW as SKEW_V48,
    validate_public_identity as validate_public_identity_v48,
)
from radia_mcp.radia_ngsolve.motor_v44_identity import (
    validate_public_identity as validate_public_identity_v45,
)
from radia_mcp.radia_ngsolve.pwm_controlled_motor_loss_gate import pwm_controlled_motor_loss_gate

from _motor_generalization_payloads import (
    _AXIAL,
    _DEMAG,
    _FIELDWEAKENING,
    _INDUCTION,
    _INDUCTION_V40,
    _IPM,
    _PMSM_FORCE,
    _SKEW,
    _SRM,
    _WOUNDFIELD,
    _identity_v45,
    _identity_v47,
    _identity_v48,
    _identity_v49,
    _identity_v50,
    _identity_v51,
    _identity_v52,
    _identity_v53,
    _identity_v56,
    _payload_v19,
    _payload_v20,
    _payload_v21,
    _payload_v22,
    _payload_v23,
    _payload_v24,
    _payload_v25,
    _payload_v26,
    _payload_v27,
    _payload_v28,
    _payload_v29,
    _payload_v30,
    _payload_v31,
    _payload_v32,
    _payload_v33,
    _payload_v34,
    _payload_v35,
    _payload_v36,
    _payload_v37,
    _payload_v38,
    _payload_v39,
    _payload_v40,
    _payload_v41,
    _payload_v42,
    _payload_v43,
    _payload_v54,
    _payload_v55,
)


def test_v19_public_iron_loss_harmonic_frequency_basis_coefficient_unit_mismatch():
    payload = _payload_v19()
    payload["artifact_identity"][
        "iron_loss_harmonic_frequency_coefficient_unit_basis_identity"
    ].update(
        {
            "coefficient_loss_generation": "iron-loss-20",
            "coefficient_frequency_unit": "kHz",
            "coefficient_flux_density_unit": "mT",
            "evaluated_harmonic_frequencies_hz": [0.05, 0.15, 0.25],
            "evaluated_loss_coefficients": [1000.0, 20.0, 1.0],
            "evaluated_loss_basis_sha256": "e" * 64,
        }
    )
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert result["checks"][
        "iron_loss_harmonics_and_coefficients_share_frequency_flux_units"
    ] is False


def test_v19_public_demagnetization_operating_point_temperature_phase_generation_mismatch():
    payload = _payload_v19()
    payload["artifact_identity"][
        "demagnetization_temperature_current_phase_operating_point_identity"
    ].update(
        {
            "temperature_operating_point_generation": "operating-point-20",
            "demag_margin_operating_point_generation": "operating-point-20",
            "demag_margin_temperature_c": 80.0,
            "demag_margin_current_phase_deg": 60.0,
            "demag_margin_operating_point_sha256": "e" * 64,
        }
    )
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert result["checks"][
        "demag_margin_uses_current_temperature_and_current_phase_state"
    ] is False


def test_v20_public_skew_slice_torque_angle_weight_periodicity_generation_mismatch():
    payload = _payload_v20()
    payload["artifact_identity"][
        "skew_slice_torque_angle_weight_periodicity_generation_identity"
    ].update(
        {
            "angle_skew_generation": "skew-21",
            "weight_skew_generation": "skew-21",
            "torque_slice_ids": [3, 2, 1],
            "torque_slice_angles_deg": [5.0, 0.0, -5.0],
            "torque_quadrature_weights": [0.5, 0.25, 0.25],
            "torque_periodic_wrap_deg": 180.0,
            "torque_skew_average_table_sha256": "f" * 64,
        }
    )
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert result["checks"][
        "skew_torque_uses_current_slice_angles_weights_and_periodicity"
    ] is False


def test_v20_public_incremental_inductance_current_perturbation_phase_state_generation_mismatch():
    payload = _payload_v20()
    payload["artifact_identity"][
        "incremental_inductance_current_perturbation_phase_state_generation_identity"
    ].update(
        {
            "perturbation_operating_point_generation": "operating-point-21",
            "phase_state_operating_point_generation": "operating-point-21",
            "matrix_base_solve_generation": "solve-21-base",
            "matrix_perturbation_solve_generations": [
                "solve-22-b",
                "solve-22-c",
                "solve-22-a",
            ],
            "matrix_phase_names": ["b", "c", "a"],
            "matrix_perturbation_currents_a": [
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
                [1.0, 0.0, 0.0],
            ],
            "resolved_incremental_inductance_table_sha256": "f" * 64,
        }
    )
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert result["checks"][
        "incremental_inductance_uses_current_perturbation_phase_and_state"
    ] is False


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


def test_v22_public_motion_skew_force_harmonic_time_angle_phase_generation_mismatch():
    payload = _payload_v22()
    identity = payload["artifact_identity"][
        "motion_skew_force_harmonic_time_angle_phase_generation_identity"
    ]
    identity.update(
        {
            "time_motion_study_generation": "motion-40",
            "angle_motion_study_generation": "motion-39",
            "skew_motion_study_generation": "motion-38",
            "phase_motion_study_generation": "motion-37",
            "force_time_s": [0.0, 0.002, 0.004],
            "force_mechanical_angle_deg": [0.0, 6.0, 12.0],
            "force_slice_weights": [0.5, 0.25, 0.25],
            "force_phase_reference_deg": -30.0,
            "force_harmonic_orders": [1, 5, 7],
            "reported_force_harmonics_n": [120.0, 1.2, 0.7],
            "resolved_force_harmonic_table_sha256": "a" * 64,
        }
    )
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "force_harmonics_use_current_motion_skew_time_angle_and_phase"
    ]


def test_v22_public_ipm_irreversible_demag_recoil_temperature_operating_generation_mismatch():
    payload = _payload_v22()
    identity = payload["artifact_identity"][
        "ipm_irreversible_demag_recoil_temperature_operating_generation_identity"
    ]
    identity.update(
        {
            "recoil_curve_demag_study_generation": "demag-40",
            "temperature_demag_study_generation": "demag-39",
            "operating_point_demag_study_generation": "demag-38",
            "magnet_orientation_demag_study_generation": "demag-37",
            "result_temperature_c": 80.0,
            "result_operating_point_id": "id=-100A;iq=120A;theta=5deg",
            "result_magnet_orientation_vectors": [[-1.0, 0.0], [0.0, 1.0]],
            "result_recoil_curve_sha256": "b" * 64,
            "result_magnet_state_sha256": "c" * 64,
            "reported_demag_margin_a_per_m": [310000.0, 205000.0],
        }
    )
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "irreversible_demag_uses_current_recoil_temperature_operating_state"
    ]


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


def test_v24_public_loss_torque_speed_power_balance_harmonic_window_generation_mismatch():
    payload = _payload_v24()
    payload["artifact_identity"][
        "loss_torque_speed_power_balance_harmonic_window_generation_identity"
    ].update(
        {
            "torque_speed_power_balance_generation": "power-100",
            "harmonic_window_power_balance_generation": "power-99",
            "time_average_power_balance_generation": "power-98",
            "iron_loss_power_balance_generation": "power-97",
            "copper_loss_power_balance_generation": "power-96",
            "mechanical_loss_power_balance_generation": "power-95",
            "power_balance_torque_nm": [1.2, 1.0],
            "power_balance_speed_rad_s": [50.0, 200.0],
            "power_balance_mechanical_output_w": [60.0, 200.0],
            "loss_harmonic_window_samples": [0, 40],
            "loss_time_average_window_s": [0.0, 0.04],
            "power_balance_iron_loss_w": [7.0, 8.0],
            "power_balance_copper_loss_w": [4.0, 5.0],
            "power_balance_mechanical_loss_w": [0.0, 0.0],
            "power_balance_electrical_input_w": [71.0, 213.0],
            "reported_power_balance_sha256": "c" * 64,
        }
    )
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "loss_power_balance_uses_current_torque_speed_windows_and_loss_components"
    ]


def test_v24_public_skew_slice_weight_rotor_angle_phase_periodicity_generation_mismatch():
    payload = _payload_v24()
    payload["artifact_identity"][
        "skew_slice_weight_rotor_angle_phase_periodicity_generation_identity"
    ].update(
        {
            "weight_skew_generation": "skew-100",
            "angle_skew_generation": "skew-99",
            "phase_skew_generation": "skew-98",
            "periodicity_skew_generation": "skew-97",
            "solve_skew_generation": "skew-96",
            "result_quadrature_weights": [0.5, 0.5, 0.5],
            "result_rotor_angles_deg": [5.0, 0.0, -5.0],
            "result_current_phase_ids": ["acb@-5", "abc@0", "abc@5"],
            "result_periodic_map_ids": ["p3", "p2", "p1"],
            "result_slice_solve_sha256": ["d" * 64, "3" * 64, "4" * 64],
            "result_slice_torque_nm": [1.2, 1.0, 0.8],
            "reported_weighted_torque_nm": 1.5,
            "reported_skew_result_sha256": "e" * 64,
        }
    )
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "skew_slice_result_uses_current_weights_angles_phases_and_periodicity"
    ]


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


def test_v28_public_pwm_current_harmonic_time_alignment_electrical_angle_torque_loss_average_mismatch():
    payload = _payload_v28()
    identity = payload["artifact_identity"][
        "pwm_current_harmonic_time_electrical_angle_torque_loss_mesh_result_generation_identity"
    ]
    identity.update(
        {
            "current_pwm_generation": "pwm-observables-150",
            "time_pwm_generation": "pwm-observables-149",
            "result_harmonic_orders": [1, 3, 5],
            "result_current_harmonic_a": [100.0, 20.0, 8.0],
            "result_current_phase_deg": [0.0, 30.0, -20.0],
            "result_time_s": [0.0, 0.001, 0.002],
            "result_electrical_angle_deg": [0.0, 60.0, 120.0],
            "result_pole_pairs": 2,
            "result_torque_window_s": [0.0, 0.002],
            "result_torque_average_nm": 38.0,
            "result_loss_average_w": 500.0,
            "result_mesh_sha256": "8" * 64,
            "accepted_result_sha256": "9" * 64,
        }
    )
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "pwm_observables_use_current_harmonics_time_angle_torque_loss_mesh_and_result"
    ]


def test_v28_public_skew_slice_angular_offset_periodic_weight_rotor_frame_torque_ripple_mismatch():
    payload = _payload_v28()
    identity = payload["artifact_identity"][
        "skew_slice_angle_weight_frame_interpolation_torque_ripple_mesh_generation_identity"
    ]
    identity.update(
        {
            "angle_skew_generation": "skew-average-150",
            "frame_skew_generation": "skew-average-149",
            "result_slice_angles_deg": [3.0, 0.0, -3.0],
            "result_slice_weights": [0.5, 0.5, 0.5],
            "result_rotor_frame": "electrical-clockwise",
            "result_interpolation_rule": "linear",
            "result_slice_torque_nm": [38.0, 45.0, 41.0],
            "result_torque_average_nm": 41.0,
            "result_torque_ripple_nm": 7.0,
            "result_mesh_sha256": "a" * 64,
            "accepted_result_sha256": "b" * 64,
        }
    )
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "skew_average_uses_current_angles_weights_frame_interpolation_torque_mesh_and_result"
    ]


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


def test_v30_public_ipm_dq_inductance_current_angle_frame_saturation_reciprocity_result_mismatch():
    payload = _payload_v30()
    payload["artifact_identity"]["ipm_dq_inductance_current_angle_park_saturation_flux_derivative_reciprocity_mesh_result_identity"].update({
        "current_dq_generation": "ipm-dq-170", "mesh_dq_generation": "ipm-dq-169",
        "result_current_magnitude_a": 80.0, "result_current_angle_electrical_deg": -30.0,
        "result_park_frame": "stator_q_aligned_clockwise",
        "result_saturation_operating_point_a": [50.0, 86.6],
        "result_flux_linkage_derivative_h": [[0.004, 0.001], [-0.0005, 0.005]],
        "result_reciprocity_tolerance_h": 0.1,
        "result_mesh_sha256": "8" * 64, "accepted_result_sha256": "9" * 64,
    })
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["ipm_dq_inductance_uses_current_angle_park_frame_saturation_derivatives_reciprocity_mesh_and_result"]


def test_v30_public_srm_torque_current_position_coenergy_periodicity_phase_sequence_mismatch():
    payload = _payload_v30()
    payload["artifact_identity"]["srm_torque_current_position_coenergy_periodicity_phase_sequence_mesh_result_identity"].update({
        "current_srm_generation": "srm-coenergy-170", "periodicity_srm_generation": "srm-coenergy-169",
        "result_current_a": [0.0, 20.0, 40.0],
        "result_rotor_position_mechanical_deg": [0.0, 1.0, 2.0],
        "result_coenergy_j_at_50a": [1.9, 2.2, 2.1], "result_torque_nm_at_50a": -20.0,
        "result_sector_period_mechanical_deg": 45.0, "result_phase_sequence": ["A", "C", "B"],
        "result_mesh_sha256": "a" * 64, "accepted_result_sha256": "b" * 64,
    })
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["srm_torque_uses_current_positions_coenergy_periodicity_phase_sequence_mesh_and_result"]


def test_v31_public_pwm_iron_loss_sampling_carrier_sideband_angle_alias_energy_balance_mismatch():
    payload = _payload_v31(); record = payload["artifact_identity"]["pwm_iron_loss_sampling_sideband_angle_alias_volume_energy_result_identity"]
    record.update({"sampling_loss_generation": "pwm-loss-180", "result_sample_period_s": 2.5e-5, "result_carrier_sidebands_hz": [9800.0, 10200.0], "result_pole_pairs": 3, "result_alias_filter": "none", "result_active_volume_m3": 0.001, "result_cycle_energy_j": 0.2, "accepted_result_sha256": "b" * 64})
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["pwm_iron_loss_uses_current_sampling_sidebands_angles_alias_volume_energy_mesh_and_result"]


def test_v31_public_skew_slice_torque_weight_phase_offset_periodicity_ripple_mismatch():
    payload = _payload_v31(); record = payload["artifact_identity"]["skew_slice_torque_weight_axial_phase_periodicity_ripple_mesh_result_identity"]
    record.update({"weight_skew_generation": "skew-torque-180", "result_slice_weights": [1.0, 1.0, 1.0], "result_electrical_phase_offsets_deg": [2.0, 0.0, -2.0], "result_periodic_wrap_electrical_deg": 180.0, "result_skew_mean_torque_nm": 30.2, "result_skew_ripple_harmonics_nm": {"6": 0.8}, "accepted_result_sha256": "c" * 64})
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["skew_torque_uses_current_slice_weights_axial_phase_periodicity_ripple_mesh_and_result"]


def test_v32_public_ipm_demagnetization_knee_temperature_current_angle_region_fraction_mismatch():
    payload = _payload_v32()
    record = payload["artifact_identity"][
        "ipm_demagnetization_knee_temperature_current_angle_region_fraction_mesh_result_identity"
    ]
    record.update(
        {
            "knee_demag_generation": "ipm-demag-190",
            "temperature_demag_generation": "ipm-demag-189",
            "result_demag_generation": "ipm-demag-188",
            "result_knee_criterion": "b_magnitude_below_room_temperature_knee",
            "result_magnet_temperature_c": 20.0,
            "result_phase_current_rms_a": 120.0,
            "result_current_angle_electrical_deg": 90.0,
            "result_irreversible_region_labels": ["magnet_3/center"],
            "result_demagnetized_fraction": 0.15,
            "result_operating_point_owner": "ipm/old-case",
            "result_mesh_sha256": "8" * 64,
            "accepted_result_sha256": "9" * 64,
        }
    )
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "ipm_demagnetization_uses_current_knee_temperature_current_angle_regions_fraction_mesh_owner_and_result"
    ]


def test_v32_public_synrm_dq_map_angle_saturation_cross_coupling_mtpa_torque_derivative_mismatch():
    payload = _payload_v32()
    record = payload["artifact_identity"][
        "synrm_dq_map_angle_saturation_cross_coupling_mtpa_torque_mesh_result_identity"
    ]
    record.update(
        {
            "angle_map_generation": "synrm-dq-190",
            "cross_coupling_map_generation": "synrm-dq-189",
            "result_map_generation": "synrm-dq-188",
            "result_electrical_angle_deg": [0.0, 15.0, 30.0, 45.0],
            "result_saturation_branch": "linearized",
            "result_flux_map_rows": [],
            "result_dpsi_d_diq_h": [0.0, 0.0, 0.0],
            "result_dpsi_q_did_h": [-1.0e-3, -1.0e-3, -1.0e-3],
            "result_mtpa_row_indices": [2, 1, 0],
            "result_pole_pairs": 3,
            "result_torque_reconstruction": "1.5*(psi_d*iq-psi_q*id)",
            "result_torque_nm": [10.0, 20.0, 30.0],
            "result_mesh_sha256": "a" * 64,
            "accepted_map_sha256": "b" * 64,
        }
    )
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "synrm_dq_map_uses_current_angles_saturation_cross_coupling_mtpa_torque_mesh_and_result"
    ]


def test_v33_public_srm_commutation_current_chop_dwell_overlap_coenergy_torque_loss_mismatch():
    payload = _payload_v33()
    record = payload["artifact_identity"][
        "srm_commutation_phase_dwell_chop_overlap_coenergy_torque_loss_angle_mesh_result_identity"
    ]
    record.update(
        {
            "phase_generation": "srm-commutation-200",
            "coenergy_generation": "srm-commutation-199",
            "result_generation": "srm-commutation-198",
            "result_phase_sequence": ["C", "B", "A"],
            "result_turn_on_deg": [5.0, 35.0, 65.0],
            "result_turn_off_deg": [10.0, 40.0, 70.0],
            "result_current_chop_a": 50.0,
            "result_overlap_deg": -5.0,
            "result_angle_grid_rad": [0.0, 0.2, 0.1],
            "result_coenergy_j": [0.0, 1.0, 0.5],
            "result_torque_nm": [-5.0, 0.0, 5.0],
            "result_copper_loss_w": 10.0,
            "result_iron_loss_w": 300.0,
            "result_total_loss_w": 100.0,
            "result_mesh_sha256": "8" * 64,
            "accepted_result_sha256": "9" * 64,
        }
    )
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "srm_commutation_uses_current_phases_dwell_chop_overlap_coenergy_torque_loss_mesh_and_result"
    ]


def test_v33_public_axial_flux_pm_sector_airgap_end_effect_torque_axial_force_surface_mismatch():
    payload = _payload_v33()
    record = payload["artifact_identity"][
        "axial_flux_pm_sector_airgap_end_effect_torque_force_surface_direction_frame_mesh_result_identity"
    ]
    record.update(
        {
            "sector_generation": "axial-flux-pm-200",
            "force_generation": "axial-flux-pm-199",
            "result_generation": "axial-flux-pm-198",
            "result_sector_multiplier": 6,
            "result_air_gaps_m": [0.001, 0.002],
            "result_end_effect_factor": 1.2,
            "result_surface_coordinates": [[0.0, 0.0]],
            "result_torque_surface_nm": [5.0],
            "result_axial_force_surface_n": [-18.0],
            "result_force_direction": "-x",
            "result_axial_frame": "stale_local_y",
            "result_mesh_sha256": "a" * 64,
            "accepted_result_lineage_sha256": "b" * 64,
        }
    )
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "axial_flux_pm_uses_current_sector_airgaps_end_effect_torque_force_surface_frame_mesh_and_result"
    ]


def test_v34_public_demagnetization_temperature_recoil_loadline_operating_point_margin_mismatch():
    payload = _payload_v34()
    record = payload["artifact_identity"][
        "pm_demagnetization_temperature_recoil_loadline_operating_point_knee_margin_angle_mesh_owner_result_identity"
    ]
    record.update(
        {
            "temperature_generation": "pm-demag-operating-point-210",
            "margin_generation": "pm-demag-operating-point-209",
            "result_generation": "pm-demag-operating-point-208",
            "result_operating_temperature_c": 20.0,
            "result_remanence_operating_t": 1.2,
            "result_coercivity_operating_a_m": 900000.0,
            "result_recoil_permeability_relative": -1.05,
            "result_loadline_slope_t_per_a_m": -1.0e-6,
            "result_operating_field_a_m": -800000.0,
            "result_operating_flux_density_t": -0.2,
            "result_knee_field_a_m": -600000.0,
            "result_irreversible_margin_a_m": -200000.0,
            "result_rotor_angle_rad": -0.3,
            "result_demag_mesh_sha256": "9" * 64,
            "accepted_demag_result_owner": "stale/demag",
            "accepted_demag_result_sha256": "a" * 64,
        }
    )
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "pm_demagnetization_uses_current_temperature_recoil_loadline_knee_margin_angle_mesh_owner_and_result"
    ]


def test_v34_public_eccentricity_unbalanced_magnetic_pull_harmonic_frame_force_torque_mismatch():
    payload = _payload_v34()
    record = payload["artifact_identity"][
        "eccentricity_static_dynamic_frame_radial_force_harmonic_ump_torque_pole_periodicity_angle_owner_result_identity"
    ]
    record.update(
        {
            "static_generation": "eccentricity-ump-210",
            "harmonic_generation": "eccentricity-ump-209",
            "result_generation": "eccentricity-ump-208",
            "result_static_eccentricity_m": [-1.0e-4, 0.0],
            "result_dynamic_eccentricity_amplitude_m": -5.0e-5,
            "result_mechanical_frame": "rotor_local_yz",
            "result_radial_force_harmonics_n": [[0, 0.0, 0.0], [1, -100.0, 50.0]],
            "result_unbalanced_magnetic_pull_n": [0.0, -100.0],
            "result_torque_nm": -20.0,
            "result_pole_pairs": 2,
            "result_periodicity_angle_rad": math.pi,
            "result_angle_grid_rad": [0.0, 1.0, 0.5],
            "result_eccentricity_mesh_sha256": "b" * 64,
            "accepted_eccentricity_result_owner": "stale/eccentricity",
            "accepted_eccentricity_result_sha256": "c" * 64,
        }
    )
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "eccentricity_ump_uses_current_static_dynamic_frame_harmonics_force_torque_periodicity_angles_owner_and_result"
    ]


def test_v34_public_rejects_self_consistent_but_temperature_wrong_remanence():
    payload = _payload_v34()
    record = payload["artifact_identity"][
        "pm_demagnetization_temperature_recoil_loadline_operating_point_knee_margin_angle_mesh_owner_result_identity"
    ]
    record["remanence_operating_t"] = record["result_remanence_operating_t"] = 1.1
    record["operating_flux_density_t"] = record["result_operating_flux_density_t"] = 0.5
    assert pwm_controlled_motor_loss_gate(payload)["status"] == "needs_attention"


def test_v34_public_rejects_self_consistent_ump_opposite_static_eccentricity():
    payload = _payload_v34()
    record = payload["artifact_identity"][
        "eccentricity_static_dynamic_frame_radial_force_harmonic_ump_torque_pole_periodicity_angle_owner_result_identity"
    ]
    harmonics = [[0, 0.0, 0.0], [1, -100.0, 0.0], [2, 0.0, 0.0]]
    record["radial_force_harmonics_n"] = harmonics
    record["result_radial_force_harmonics_n"] = harmonics
    record["unbalanced_magnetic_pull_n"] = [-100.0, 0.0]
    record["result_unbalanced_magnetic_pull_n"] = [-100.0, 0.0]
    assert pwm_controlled_motor_loss_gate(payload)["status"] == "needs_attention"


def test_v35_public_skew_slice_torque_harmonic_axial_weight_phase_periodicity_mismatch():
    payload = _payload_v35()
    row = payload["artifact_identity"]["skew_slice_torque_angle_axial_weight_harmonic_phase_pole_periodicity_mean_ripple_mesh_owner_result_identity"]
    row.update({"slice_generation": "skew-slice-torque-234", "phase_generation": "skew-slice-torque-233",
                "result_generation": "skew-slice-torque-232", "result_slice_angles_rad": [0.1, 0.0, -0.1],
                "result_axial_weights": [0.8, 0.8, -0.6],
                "result_harmonic_phase_shifts_rad": [[1, [0.0, 0.0, 0.0]], [6, [1.0, 1.0, 1.0]]],
                "result_pole_pairs": 2, "result_pole_periodicity_angle_rad": math.pi,
                "result_slice_mean_torque_nm": [10.0, -20.0, 100.0], "result_weighted_mean_torque_nm": -49.0,
                "result_torque_ripple_spectrum_nm": [[0, -49.0], [5, 20.0]],
                "result_skew_mesh_sha256": "9" * 64, "accepted_skew_result_owner": "stale/skew",
                "accepted_skew_result_sha256": "a" * 64})
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["skew_slice_torque_closes_angles_axial_weights_harmonic_phases_pole_periodicity_mean_ripple_mesh_owner_and_result"]


def test_v35_public_ironloss_hysteresis_eddy_excess_waveform_frequency_volume_temperature_mismatch():
    payload = _payload_v35()
    row = payload["artifact_identity"]["ironloss_hysteresis_eddy_excess_waveform_frequency_coeff_volume_temperature_total_owner_result_identity"]
    row.update({"component_generation": "iron-loss-separation-234", "temperature_generation": "iron-loss-separation-233",
                "result_generation": "iron-loss-separation-232", "result_b_waveform_peak_t": -1.2,
                "result_waveform_factor": 0.0, "result_frequency_hz": -100.0,
                "result_material_coefficients": [2.0, -0.1, 0.0], "result_active_volume_m3": -0.001,
                "result_temperature_c": 20.0, "result_temperature_factor": -1.0,
                "result_loss_components_w_m3": [100.0, -200.0, 0.0], "result_total_iron_loss_w": 99.0,
                "result_waveform_sha256": "b" * 64, "accepted_ironloss_owner": "stale/loss",
                "accepted_ironloss_result_sha256": "c" * 64})
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["iron_loss_closes_hysteresis_eddy_excess_waveform_frequency_coefficients_volume_temperature_total_owner_and_result"]


def test_v35_rejects_self_consistent_skew_weight_sum_error():
    payload = _payload_v35()
    row = payload["artifact_identity"]["skew_slice_torque_angle_axial_weight_harmonic_phase_pole_periodicity_mean_ripple_mesh_owner_result_identity"]
    row["axial_weights"] = row["result_axial_weights"] = [0.5, 0.5, 0.5]
    assert pwm_controlled_motor_loss_gate(payload)["status"] == "needs_attention"


def test_v35_rejects_self_consistent_ironloss_total_error():
    payload = _payload_v35()
    row = payload["artifact_identity"]["ironloss_hysteresis_eddy_excess_waveform_frequency_coeff_volume_temperature_total_owner_result_identity"]
    row["total_iron_loss_w"] = row["result_total_iron_loss_w"] = 99.0
    assert pwm_controlled_motor_loss_gate(payload)["status"] == "needs_attention"


def test_v36_public_dq_flux_inductance_torque_mtpa_current_angle_speed_convention_mismatch():
    payload = _payload_v36()
    row = payload["artifact_identity"]["dq_flux_inductance_torque_mtpa_current_angle_speed_convention_owner_result_identity"]
    row.update({"park_generation": "dq-mtpa-235", "speed_generation": "dq-mtpa-234",
                "result_generation": "dq-mtpa-233", "result_park_convention": "amplitude_invariant_d_leads_q",
                "result_pole_pairs": 2, "result_current_dq_a": [100.0, -50.0],
                "result_current_magnitude_a": -1.0, "result_current_angle_rad": -2.0,
                "result_pm_flux_linkage_wb_turn": -0.075, "result_flux_linkage_dq_wb_turn": [0.15, 0.025],
                "result_differential_inductance_dq_h": [-0.001, 0.0], "result_torque_nm": -60.0,
                "result_mechanical_speed_rad_s": -100.0, "result_electrical_speed_rad_s": 100.0,
                "accepted_dq_owner": "stale/dq", "accepted_dq_result_sha256": "a" * 64})
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["dq_map_closes_park_flux_inductance_torque_mtpa_current_angle_speed_owner_and_result"]


def test_v36_public_iron_loss_hysteresis_eddy_excess_frequency_flux_energy_balance_mismatch():
    payload = _payload_v36()
    row = payload["artifact_identity"]["iron_loss_component_frequency_flux_region_thermal_energy_balance_owner_result_identity"]
    row.update({"component_generation": "iron-loss-energy-235", "thermal_generation": "iron-loss-energy-234",
                "result_generation": "iron-loss-energy-233", "result_frequency_hz": -100.0,
                "result_flux_peak_t": -1.2, "result_loss_coefficients": [2.0, -0.1, 0.0],
                "result_loss_components_w_m3": [1.0, -2.0, 3.0],
                "result_regional_volumes_m3": [["stator", -0.0007], ["old", 0.0]],
                "result_temperature_c": 20.0, "result_temperature_factor": -1.0,
                "result_total_iron_loss_w": -5.0, "result_integration_duration_s": -0.2,
                "result_loss_energy_j": 99.0, "accepted_iron_loss_owner": "stale/loss",
                "accepted_iron_loss_result_sha256": "b" * 64})
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["iron_loss_closes_components_frequency_flux_regions_thermal_power_energy_owner_and_result"]


def test_v36_rejects_self_consistent_dq_torque_error():
    payload = _payload_v36()
    row = payload["artifact_identity"]["dq_flux_inductance_torque_mtpa_current_angle_speed_convention_owner_result_identity"]
    row["torque_nm"] = row["result_torque_nm"] = -60.0
    assert pwm_controlled_motor_loss_gate(payload)["status"] == "needs_attention"


def test_v36_rejects_self_consistent_iron_loss_energy_error():
    payload = _payload_v36()
    row = payload["artifact_identity"]["iron_loss_component_frequency_flux_region_thermal_energy_balance_owner_result_identity"]
    row["loss_energy_j"] = row["result_loss_energy_j"] = 99.0
    assert pwm_controlled_motor_loss_gate(payload)["status"] == "needs_attention"


def test_v37_public_induction_motor_slip_airgap_power_torque_rotor_loss_efficiency_owner_mismatch():
    payload = _payload_v37()
    row = payload["artifact_identity"][
        "induction_motor_slip_synchronous_speed_airgap_power_torque_rotor_loss_mechanical_output_efficiency_owner_result_identity"
    ]
    row.update({
        "slip_generation": "induction-power-245", "torque_generation": "induction-power-244",
        "result_generation": "induction-power-243", "result_slip": 0.2,
        "result_synchronous_speed_rad_s": -1.0, "result_mechanical_speed_rad_s": 200.0,
        "result_airgap_power_w": -1000.0, "result_electromagnetic_torque_nm": -5.0,
        "result_rotor_copper_loss_w": 500.0, "result_mechanical_output_w": 100.0,
        "result_input_power_w": 100.0, "result_efficiency": 2.0,
        "accepted_motor_owner": "stale/motor", "accepted_motor_result_sha256": "a" * 64})
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "induction_motor_closes_slip_speed_airgap_power_torque_rotor_loss_output_efficiency_owner_and_result"
    ]


def test_v37_public_axial_flux_motor_periodicity_airgap_flux_torque_ripple_backemf_owner_mismatch():
    payload = _payload_v37()
    row = payload["artifact_identity"][
        "axial_flux_motor_sector_periodicity_dual_airgap_axial_flux_torque_ripple_backemf_frame_mesh_owner_result_identity"
    ]
    row.update({
        "sector_generation": "axial-flux-245", "backemf_generation": "axial-flux-244",
        "result_generation": "axial-flux-243", "result_sector_factor": 6,
        "result_sector_angle_rad": math.pi, "result_dual_airgap_m": [0.001, -0.001],
        "result_sector_axial_flux_per_gap_wb": [0.002, -0.002],
        "result_sector_torque_samples_nm": [2.0, -2.1, 1.9],
        "result_full_machine_torque_samples_nm": [1.0, 2.0, 3.0],
        "result_average_torque_nm": -1.0, "result_torque_ripple_ratio": -0.5,
        "result_backemf_phase_angles_rad": [0.0, 0.0, 0.0],
        "result_coordinate_frame": "cartesian_x_axial",
        "accepted_mesh_owner": "stale/mesh",
        "accepted_axial_flux_result_sha256": "b" * 64})
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "axial_flux_motor_closes_sector_dual_airgap_flux_torque_ripple_backemf_frame_mesh_owner_and_result"
    ]


def test_v37_public_rejects_self_consistent_induction_power_imbalance():
    payload = _payload_v37()
    row = payload["artifact_identity"][
        "induction_motor_slip_synchronous_speed_airgap_power_torque_rotor_loss_mechanical_output_efficiency_owner_result_identity"
    ]
    row["rotor_copper_loss_w"] = row["result_rotor_copper_loss_w"] = 500.0
    assert pwm_controlled_motor_loss_gate(payload)["status"] == "needs_attention"


def test_v37_public_rejects_self_consistent_axial_flux_phase_collapse():
    payload = _payload_v37()
    row = payload["artifact_identity"][
        "axial_flux_motor_sector_periodicity_dual_airgap_axial_flux_torque_ripple_backemf_frame_mesh_owner_result_identity"
    ]
    row["backemf_phase_angles_rad"] = row["result_backemf_phase_angles_rad"] = [0.0, 0.0, 0.0]
    assert pwm_controlled_motor_loss_gate(payload)["status"] == "needs_attention"


def test_v38_public_wound_field_synchronous_excitation_flux_torque_powerfactor_fieldloss_energy_mismatch():
    payload = _payload_v38()
    row = payload["artifact_identity"][
        "wound_field_synchronous_excitation_flux_torque_angle_powerfactor_field_stator_loss_mechanical_energy_mesh_owner_result_identity"
    ]
    row.update(
        {
            "flux_generation": "wound-field-257",
            "result_field_current_a": -5.0,
            "result_excitation_flux_linkage_wb_turn": -0.2,
            "result_torque_angle_rad": math.pi,
            "result_electromagnetic_torque_nm": -1.0,
            "result_power_factor": 1.5,
            "result_field_copper_loss_w": -25.0,
            "result_energy_balance_residual_w": 100.0,
            "accepted_mesh_owner": "stale:wound-field",
        }
    )
    assert pwm_controlled_motor_loss_gate(payload)["status"] == "needs_attention"


def test_v38_public_flux_switching_pm_slot_pole_polarity_harmonic_backemf_torque_ripple_mismatch():
    payload = _payload_v38()
    row = payload["artifact_identity"][
        "flux_switching_pm_slot_pole_polarity_phase_harmonic_backemf_torque_ripple_periodicity_mesh_owner_result_identity"
    ]
    row.update(
        {
            "polarity_generation": "flux-switching-257",
            "result_slot_count": 10,
            "result_pole_count": 12,
            "result_magnet_polarity_sequence": [1] * 10,
            "result_phase_sequence": "ACB",
            "result_working_harmonic_order": 3,
            "result_backemf_phase_angles_rad": [0.0, 0.0, 0.0],
            "result_torque_ripple_ratio": -1.0,
            "result_periodic_multiplier": 1,
            "accepted_mesh_owner": "stale:flux-switching",
        }
    )
    assert pwm_controlled_motor_loss_gate(payload)["status"] == "needs_attention"


def test_v38_public_rejects_self_consistent_wrong_wound_field_torque():
    payload = _payload_v38()
    row = payload["artifact_identity"][
        "wound_field_synchronous_excitation_flux_torque_angle_powerfactor_field_stator_loss_mechanical_energy_mesh_owner_result_identity"
    ]
    row["electromagnetic_torque_nm"] *= 2.0
    row["result_electromagnetic_torque_nm"] = row["electromagnetic_torque_nm"]
    assert pwm_controlled_motor_loss_gate(payload)["status"] == "needs_attention"


def test_v38_public_rejects_self_consistent_wrong_flux_switching_periodicity():
    payload = _payload_v38()
    row = payload["artifact_identity"][
        "flux_switching_pm_slot_pole_polarity_phase_harmonic_backemf_torque_ripple_periodicity_mesh_owner_result_identity"
    ]
    row["periodic_multiplier"] = 1
    row["result_periodic_multiplier"] = 1
    row["sector_slot_count"] = 12
    row["result_sector_slot_count"] = 12
    row["sector_pole_count"] = 10
    row["result_sector_pole_count"] = 10
    assert pwm_controlled_motor_loss_gate(payload)["status"] == "needs_attention"


def test_v39_public_skewed_rotor_slice_angle_phase_weight_torque_ripple_power_mismatch():
    payload = _payload_v39()
    payload["artifact_identity"][_SKEW].update({"phase_generation": "skewed-rotor-270", "power_generation": "skewed-rotor-269", "result_generation": "skewed-rotor-268", "result_slice_phase_offsets_electrical_deg": [5.0, 0.0, -5.0], "result_axial_weights": [0.5, 0.5, 0.5], "result_slice_mean_torque_nm": [9.0, 9.0, 9.0], "result_weighted_mean_torque_nm": -1.0, "result_weighted_ripple_residual": -1.0, "result_mechanical_power_w": -100.0, "accepted_model_owner": "stale:motor", "accepted_skew_result_sha256": "a" * 64})
    assert pwm_controlled_motor_loss_gate(payload)["status"] == "needs_attention"


def test_v39_public_pm_irreversible_demag_recoil_temperature_operatingpoint_flux_torque_mismatch():
    payload = _payload_v39()
    payload["artifact_identity"][_DEMAG].update({"temperature_generation": "pm-demag-270", "operating_generation": "pm-demag-269", "result_generation": "pm-demag-268", "result_temperature_adjusted_remanence_t": -1.0, "result_operating_h_a_per_m": 8.0e5, "result_operating_b_t": -1.0, "result_irreversible_region": False, "result_remanence_loss_fraction": -0.1, "result_airgap_flux_after_wb": 2.0e-2, "result_torque_after_nm": 20.0, "accepted_mesh_owner": "stale:mesh", "accepted_demag_result_sha256": "b" * 64})
    assert pwm_controlled_motor_loss_gate(payload)["status"] == "needs_attention"


def test_v39_public_rejects_self_consistent_wrong_skew_phase_conversion():
    payload = _payload_v39()
    row = payload["artifact_identity"][_SKEW]
    phases = row["slice_angles_mechanical_deg"]
    row["slice_phase_offsets_electrical_deg"] = phases
    row["result_slice_phase_offsets_electrical_deg"] = phases
    assert pwm_controlled_motor_loss_gate(payload)["status"] == "needs_attention"


def test_v39_public_rejects_self_consistent_reversible_demag_claim():
    payload = _payload_v39()
    row = payload["artifact_identity"][_DEMAG]
    row["irreversible_region"] = False
    row["result_irreversible_region"] = False
    assert pwm_controlled_motor_loss_gate(payload)["status"] == "needs_attention"


def test_v40_public_induction_cage_rotorbar_slip_current_loss_torque_endring_power_mismatch():
    payload = _payload_v40()
    payload["artifact_identity"][_INDUCTION_V40].update({"slip_generation": "induction-cage-279", "result_slip": 0.2, "result_airgap_power_w": -1.0, "accepted_mesh_owner": "stale:mesh"})
    assert pwm_controlled_motor_loss_gate(payload)["status"] == "needs_attention"


def test_v40_public_ipm_fieldweakening_dq_flux_voltage_limit_current_angle_torque_power_mismatch():
    payload = _payload_v40()
    payload["artifact_identity"][_FIELDWEAKENING].update({"flux_generation": "ipm-fieldweakening-279", "result_flux_d_wb": -1.0, "result_voltage_d_v": 200.0, "result_electrical_power_w": -1.0, "accepted_model_owner": "stale:motor"})
    assert pwm_controlled_motor_loss_gate(payload)["status"] == "needs_attention"


def test_v40_public_rejects_self_consistent_wrong_induction_slip():
    payload = _payload_v40()
    row = payload["artifact_identity"][_INDUCTION_V40]
    row["slip"] = row["result_slip"] = 0.2
    assert pwm_controlled_motor_loss_gate(payload)["status"] == "needs_attention"


def test_v40_public_rejects_self_consistent_fieldweakening_power_gap():
    payload = _payload_v40()
    row = payload["artifact_identity"][_FIELDWEAKENING]
    row["electrical_power_w"] = row["result_electrical_power_w"] = row["mechanical_power_w"]
    assert pwm_controlled_motor_loss_gate(payload)["status"] == "needs_attention"


def test_v41_public_srm_inductance_position_current_coenergy_torque_ripple_power_mismatch():
    payload = _payload_v41()
    payload["artifact_identity"][_SRM].update(
        {
            "inductance_generation": "srm-map-723",
            "result_rotor_position_rad": [0.4, 0.3, 0.2, 0.1, 0.0],
            "result_torque_nm": [-1.0],
            "accepted_model_owner": "stale:model",
        }
    )
    assert pwm_controlled_motor_loss_gate(payload)["status"] == "needs_attention"


def test_v41_public_axialflux_sector_periodicity_skew_endeffect_flux_torque_loss_power_mismatch():
    payload = _payload_v41()
    payload["artifact_identity"][_AXIAL].update(
        {
            "periodicity_generation": "axial-flux-723",
            "result_periodicity": "antiperiodic",
            "result_corrected_airgap_flux_wb": -1.0,
            "accepted_mesh_owner": "stale:mesh",
        }
    )
    assert pwm_controlled_motor_loss_gate(payload)["status"] == "needs_attention"


def test_v41_public_rejects_self_consistent_wrong_srm_coenergy():
    payload = _payload_v41()
    row = payload["artifact_identity"][_SRM]
    row["coenergy_j"] = row["result_coenergy_j"] = [1.0] * 5
    assert pwm_controlled_motor_loss_gate(payload)["status"] == "needs_attention"


def test_v41_public_rejects_self_consistent_wrong_axial_flux_power():
    payload = _payload_v41()
    row = payload["artifact_identity"][_AXIAL]
    row["electrical_power_w"] = row["result_electrical_power_w"] = row["mechanical_power_w"]
    assert pwm_controlled_motor_loss_gate(payload)["status"] == "needs_attention"


def test_v42_public_ipm_dq_mismatch():
    payload = _payload_v42()
    payload["artifact_identity"][_IPM].update(
        {
            "flux_generation": "ipm-dq-841",
            "energy_generation": "ipm-dq-840",
            "result_generation": "ipm-dq-839",
            "result_current_d_a": 20.0,
            "result_flux_d_wb": -0.08,
            "result_torque_nm": -36.0,
            "result_voltage_d_v": 101.0,
            "result_voltage_magnitude_v": 200.0,
            "result_power_factor": -0.5,
            "result_mtpv_branch": "stale_branch",
            "result_mtpv_voltage_margin_v": -50.0,
            "result_field_energy_j": -2.7,
            "result_coenergy_j": -0.7,
            "accepted_mesh_owner": "stale:mesh",
            "accepted_ipm_result_sha256": "9" * 64,
        }
    )
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "ipm_dq_maps_close_currents_flux_torque_voltage_powerfactor_mtpv_energy_mesh_and_result"
    ]


def test_v42_public_induction_motor_mismatch():
    payload = _payload_v42()
    payload["artifact_identity"][_INDUCTION].update(
        {
            "slip_generation": "induction-power-841",
            "loss_generation": "induction-power-840",
            "result_generation": "induction-power-839",
            "result_slip": -0.1,
            "result_rotor_electrical_frequency_hz": -5.0,
            "result_torque_nm": -20.0,
            "result_airgap_power_w": -1.0,
            "result_rotor_copper_loss_w": -1.0,
            "result_converted_power_w": -1.0,
            "result_mechanical_power_w": -1.0,
            "result_efficiency": 1.5,
            "accepted_motor_result_sha256": "a" * 64,
        }
    )
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "induction_motors_close_slip_rotor_frequency_losses_torque_airgap_mechanical_power_efficiency_and_result"
    ]


def test_v42_public_rejects_self_consistent_wrong_ipm_torque():
    payload = _payload_v42()
    record = payload["artifact_identity"][_IPM]
    record["torque_nm"] = 30.0
    record["result_torque_nm"] = 30.0
    assert pwm_controlled_motor_loss_gate(payload)["status"] == "needs_attention"


def test_v42_public_rejects_self_consistent_wrong_rotor_loss():
    payload = _payload_v42()
    record = payload["artifact_identity"][_INDUCTION]
    record["rotor_copper_loss_w"] *= 2.0
    record["result_rotor_copper_loss_w"] = record["rotor_copper_loss_w"]
    assert pwm_controlled_motor_loss_gate(payload)["status"] == "needs_attention"


def test_v43_public_positive_pmsm_force_nvh_and_woundfield():
    result = pwm_controlled_motor_loss_gate(_payload_v43())
    assert result["status"] == "ok", result


def test_v43_public_rejects_pmsm_force_nvh_mismatch():
    payload = _payload_v43()
    payload["artifact_identity"][_PMSM_FORCE]["result_radial_force_space_orders"] = [6]
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["pmsm_force_nvh_closes_space_time_orders_phase_torque_modal_power_energy_mesh_and_result"]


def test_v43_public_rejects_woundfield_energy_mismatch():
    payload = _payload_v43()
    payload["artifact_identity"][_WOUNDFIELD]["result_efficiency"] = 1.2
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["woundfield_motor_closes_excitation_flux_torque_copper_losses_efficiency_power_mesh_and_result"]


def test_v45_public_identity_accepts_closed_artifacts():
    checks = validate_public_identity_v45(_identity_v45())
    assert checks and all(checks.values())


def test_v45_public_identity_rejects_torque_mutation():
    identity = _identity_v45()
    identity["v45_public_ipmsm_torque_ripple_radial_force_modal_power_efficiency_energy_mesh_owner_mismatch"]["result_efficiency"] = 1.1
    checks = validate_public_identity_v45(identity)
    assert checks and not all(checks.values())


def test_v47_positive_motor_artifacts_are_accepted() -> None:
    assert all(validate_public_identity_v47(_identity_v47()).values())


def test_v47_dq_mapping_mutation_is_rejected() -> None:
    identity = _identity_v47()
    identity[DQ_V47]["result_phase_order"] = ["A", "C", "B"]
    identity[DQ_V47]["result_electrical_angle_origin_deg"] = 30.0
    identity[DQ_V47]["result_pole_pairs"] = 3
    identity[DQ_V47]["result_transform_identity"] = "amplitude_invariant_park"
    assert not all(validate_public_identity_v47(identity).values())


def test_v47_window_parameter_row_mutation_is_rejected() -> None:
    identity = _identity_v47()
    identity[WINDOW]["result_integration_window_s"] = [0.02, 0.03]
    identity[WINDOW]["result_parameter_row_key"] = "speed=6000rpm,current=50A"
    assert not all(validate_public_identity_v47(identity).values())


def test_v48_positive_skew_and_pwm_artifacts_are_accepted() -> None:
    assert all(validate_public_identity_v48(_identity_v48()).values())


def test_v48_skew_permutation_and_phase_mutations_are_rejected() -> None:
    identity = deepcopy(_identity_v48())
    identity[SKEW_V48]["result_slice_weights"] = [0.50, 0.25, 0.25]
    identity[SKEW_V48]["result_phase_origins_deg"] = [30.0, 0.0, 0.0]
    assert validate_public_identity_v48(identity)["motor_v48_skew_slice_angle_harmonic_phase_owner"] is False


def test_v48_pwm_timeline_and_owner_mutations_are_rejected() -> None:
    identity = deepcopy(_identity_v48())
    identity[PWM]["result_sample_times_s"] = [0.0, 10.0e-6, 5.0e-6, 15.0e-6]
    identity[PWM]["result_timeline_owner"] = "timeline:pwm-v48-old"
    assert validate_public_identity_v48(identity)["motor_v48_pwm_timeline_switch_electrical_loss_owner"] is False


def test_v49_positive_demag_and_iron_loss_artifacts_are_accepted() -> None:
    assert all(validate_public_identity_v49(_identity_v49()).values())


def test_v49_demag_operating_point_and_state_mutations_are_rejected() -> None:
    identity = deepcopy(_identity_v49())
    identity[DEMAG_V49]["result_temperature_c"] = 20.0
    identity[DEMAG_V49]["result_irreversible_magnet_state"] = {"magnet:north": 1.0, "magnet:south": 1.0}
    identity[DEMAG_V49]["result_operating_point_owner"] = "operating-point:old"
    assert validate_public_identity_v49(identity)["motor_v49_demag_temperature_current_angle_state_owner"] is False


def test_v49_iron_loss_window_harmonic_coefficient_and_owner_mutations_are_rejected() -> None:
    identity = deepcopy(_identity_v49())
    identity[IRON_V49]["result_time_window_s"] = [0.0, 0.005]
    identity[IRON_V49]["result_harmonic_rows"] = list(reversed(identity[IRON_V49]["harmonic_rows"]))
    identity[IRON_V49]["result_loss_coefficients"] = {"hysteresis": 0.5, "eddy": 0.04, "excess": 0.0}
    identity[IRON_V49]["result_loss_owner"] = "loss-table:old"
    assert validate_public_identity_v49(identity)["motor_v49_iron_loss_harmonic_time_frequency_coefficient_owner"] is False


def test_v50_positive_dq_and_electrothermal_artifacts_are_accepted() -> None:
    assert all(validate_public_identity_v50(_identity_v50()).values())


def test_v50_dq_angle_current_saliency_operating_point_and_owner_drift_is_rejected() -> None:
    identity = deepcopy(_identity_v50())
    identity[DQ_V50]["result_park_angle_electrical_deg"] = 7.5
    identity[DQ_V50]["result_dq_currents"] = {"id_a": 82.0, "iq_a": -35.0, "phase_order": "uwv"}
    identity[DQ_V50]["result_operating_point_id"] = "operating-point:old"
    identity[DQ_V50]["accepted_result_owner"] = "dq-result:foreign"
    assert validate_public_identity_v50(identity)["motor_v50_dq_park_current_saliency_operating_point_owner"] is False


def test_v50_thermal_loss_boundary_temperature_material_and_owner_drift_is_rejected() -> None:
    identity = deepcopy(_identity_v50())
    identity[THERMAL]["replayed_loss_map"] = {"copper_w": 34.0}
    identity[THERMAL]["replayed_convection_boundaries"] = [{"boundary": "shaft", "h_w_m2k": 5.0, "ambient_c": 40.0}]
    identity[THERMAL]["replayed_material_revisions"] = {"copper": "cu:old"}
    identity[THERMAL]["replayed_thermal_owner"] = "thermal-result:foreign"
    assert validate_public_identity_v50(identity)["motor_v50_thermal_loss_convection_temperature_material_owner"] is False


def test_v51_positive_public_artifacts_are_accepted() -> None:
    assert all(validate_public_identity_v51(_identity_v51()).values())


def test_v51_frozen_counterfactuals_are_rejected() -> None:
    identity = deepcopy(_identity_v51())
    identity[TORQUE_V51].update({"result_electrical_period_deg": 90.0, "result_torque_owner": "torque:stale"})
    identity[WINDING].update({"result_resistance_at_temperature_ohm": 0.08, "result_winding_owner": "winding:stale"})
    assert not all(validate_public_identity_v51(identity).values())


def test_v51_self_consistent_wrong_physics_are_rejected() -> None:
    identity = deepcopy(_identity_v51())
    identity[TORQUE_V51]["electrical_period_deg"] = identity[TORQUE_V51]["result_electrical_period_deg"] = 90.0
    identity[WINDING]["resistance_at_temperature_ohm"] = identity[WINDING]["result_resistance_at_temperature_ohm"] = 0.08
    assert not all(validate_public_identity_v51(identity).values())


def test_v52_positive_public_artifacts_are_accepted() -> None:
    assert all(validate_public_identity_v52(_identity_v52()).values())


def test_v52_frozen_counterfactuals_are_rejected() -> None:
    identity = deepcopy(_identity_v52())
    identity[IRON_LOSS]["result_fft_window"] = "rectangular"
    identity[COGGING]["result_cogging_period_mechanical_deg"] = 30.0
    assert not all(validate_public_identity_v52(identity).values())


def test_v52_self_consistent_wrong_physics_are_rejected() -> None:
    identity = deepcopy(_identity_v52())
    identity[IRON_LOSS]["fft_window"] = identity[IRON_LOSS]["result_fft_window"] = "rectangular"
    identity[COGGING]["cogging_period_mechanical_deg"] = identity[COGGING]["result_cogging_period_mechanical_deg"] = 30.0
    assert not all(validate_public_identity_v52(identity).values())


def test_v53_positive_public_artifacts_are_accepted():
    assert all(validate_public_identity_v53(_identity_v53()).values())


def test_v53_frozen_counterfactuals_are_rejected():
    identity = deepcopy(_identity_v53())
    identity[SKEW_V53]["result_slice_weights"] = [1.0, 0.0, 0.0]
    identity[DEMAG_V53]["result_temperature_c"] = 20.0
    assert not all(validate_public_identity_v53(identity).values())


def test_v53_self_consistent_wrong_physics_is_rejected():
    identity = deepcopy(_identity_v53())
    identity[SKEW_V53]["slice_weights"] = identity[SKEW_V53]["result_slice_weights"] = [0.5, 0.5, 0.5]
    identity[DEMAG_V53]["irreversible_demag_fraction"] = identity[DEMAG_V53]["result_irreversible_demag_fraction"] = 0.0
    assert not all(validate_public_identity_v53(identity).values())


def test_v54_positive_identities_are_accepted():
    assert all(validate_public_identity_v54(_payload_v54()).values())


def test_v54_frozen_mutations_are_rejected():
    payload = deepcopy(_payload_v54())
    payload[TORQUE_V54]["result_pole_pairs"] = 3
    payload[DEMAG_V54]["result_temperature_c"] = 20.0
    assert not all(validate_public_identity_v54(payload).values())


def test_v54_self_consistent_nonphysical_records_are_rejected():
    payload = deepcopy(_payload_v54())
    payload[TORQUE_V54]["electrical_angles_deg"] = payload[TORQUE_V54]["result_electrical_angles_deg"] = [0.0, 15.0, 30.0]
    payload[DEMAG_V54]["current_vector_abc_a"] = payload[DEMAG_V54]["result_current_vector_abc_a"] = [120.0, -50.0, -50.0]
    assert not all(validate_public_identity_v54(payload).values())


def test_v54_malformed_values_reject_without_raising():
    payload = deepcopy(_payload_v54())
    payload[TORQUE_V54]["torque_harmonics"] = [{"order": [6], "amplitude_nm": 0.35, "phase_electrical_deg": 25.0}]
    payload[DEMAG_V54]["current_vector_abc_a"] = [[120.0], -60.0, -60.0]
    assert not all(validate_public_identity_v54(payload).values())


def test_v55_positive_identities_are_accepted():
    assert all(validate_public_identity_v55(_payload_v55()).values())


def test_v55_frozen_mutations_are_rejected():
    payload = deepcopy(_payload_v55()); payload[DQ_V55]["accepted_result_owner"] = "result:stale"; payload[IRON_V55]["result_material_owner"] = "material:stale"
    assert not all(validate_public_identity_v55(payload).values())


def test_v55_self_consistent_nonphysical_records_are_rejected():
    payload = deepcopy(_payload_v55())
    payload[DQ_V55]["flux_linkage_dq_wb"] = payload[DQ_V55]["result_flux_linkage_dq_wb"] = {"d": 0.4, "q": 0.4}
    payload[IRON_V55]["loss_components"] = payload[IRON_V55]["result_loss_components"] = {"hysteresis_w": -1.0, "eddy_w": 8.0, "excess_w": 2.0}
    payload[IRON_V55]["total_iron_loss_w"] = payload[IRON_V55]["result_total_iron_loss_w"] = 9.0
    assert not all(validate_public_identity_v55(payload).values())


def test_v55_malformed_values_reject_without_raising():
    payload = deepcopy(_payload_v55()); payload[DQ_V55]["current_dq_a"] = {"d": [1.0], "q": 2.0}; payload[IRON_V55]["flux_density_waveform_t"] = [0.0, [1.0], 0.0]
    assert not all(validate_public_identity_v55(payload).values())


def test_v55_numeric_sha256_values_are_rejected():
    payload = _payload_v55()
    numeric_digest = int("9" * 64)
    for row in payload.values():
        row["result_sha256"] = numeric_digest
        row["accepted_result_sha256"] = numeric_digest
    assert not all(validate_public_identity_v55(payload).values())


def test_v56_positive_identity_is_accepted() -> None:
    assert all(validate_public_identity_v56(_identity_v56()).values())


def test_v56_frozen_result_mutations_are_rejected() -> None:
    identity = deepcopy(_identity_v56()); identity[MAP]["result_efficiency"] = 2.0; identity[INDUCTION]["result_slip"] = -0.5
    assert not all(validate_public_identity_v56(identity).values())


def test_v56_self_consistent_power_and_slip_errors_are_rejected() -> None:
    identity = deepcopy(_identity_v56()); identity[MAP]["output_power_w"] = identity[MAP]["result_output_power_w"] = 9999.0; identity[INDUCTION]["slip"] = identity[INDUCTION]["result_slip"] = -0.5
    assert not all(validate_public_identity_v56(identity).values())


def test_v56_malformed_losses_reject_without_raising() -> None:
    identity = deepcopy(_identity_v56()); identity[MAP]["loss_components_w"] = []
    assert not all(validate_public_identity_v56(identity).values())


def test_v56_numeric_digests_are_rejected() -> None:
    identity = deepcopy(_identity_v56())
    numeric_digest = int("1" * 64)
    for contract_name in (MAP, INDUCTION):
        identity[contract_name]["result_sha256"] = numeric_digest
        identity[contract_name]["accepted_result_sha256"] = numeric_digest
    assert not all(validate_public_identity_v56(identity).values())
