"""Payload builders and constants for ``test_motor_generalization.py``.

The v19..v43 builders form one cumulative chain on ``pwm_controlled_motor_loss_gate``:
each ``_payload_vN`` extends ``_payload_v(N-1)``.  The v45..v56 builders are
independent fixtures for separate ``validate_public_identity`` gates.  Names that
collided across the former per-version modules carry a version suffix.
"""

from __future__ import annotations

import math

from radia_mcp.radia_ngsolve.motor_artifact_identity_v49 import (
    DEMAG as DEMAG_V49,
    IRON as IRON_V49,
)
from radia_mcp.radia_ngsolve.motor_artifact_identity_v50 import (
    DQ as DQ_V50,
    THERMAL,
)
from radia_mcp.radia_ngsolve.motor_artifact_identity_v51 import (
    TORQUE as TORQUE_V51,
    WINDING,
)
from radia_mcp.radia_ngsolve.motor_artifact_identity_v52 import (
    COGGING,
    IRON_LOSS,
)
from radia_mcp.radia_ngsolve.motor_artifact_identity_v53 import (
    DEMAG as DEMAG_V53,
    SKEW as SKEW_V53,
)
from radia_mcp.radia_ngsolve.motor_artifact_identity_v54 import (
    DEMAG as DEMAG_V54,
    TORQUE as TORQUE_V54,
)
from radia_mcp.radia_ngsolve.motor_artifact_identity_v55 import (
    DQ as DQ_V55,
    IRON as IRON_V55,
)
from radia_mcp.radia_ngsolve.motor_artifact_lineage_v47 import (
    DQ as DQ_V47,
    WINDOW,
)
from radia_mcp.radia_ngsolve.motor_map_induction_identity_v56 import (
    INDUCTION,
    MAP,
)
from radia_mcp.radia_ngsolve.motor_semantic_identity_v48 import (
    PWM,
    SKEW as SKEW_V48,
)
from test_pwm_controlled_motor_loss_gate import _payload, _with_artifact_identity


def _payload_v19():
    payload = _with_artifact_identity(_payload())
    identity = payload["artifact_identity"]
    identity["iron_loss_harmonic_frequency_coefficient_unit_basis_identity"] = {
        "loss_generation": "iron-loss-21",
        "harmonic_loss_generation": "iron-loss-21",
        "coefficient_loss_generation": "iron-loss-21",
        "frequency_unit": "Hz",
        "coefficient_frequency_unit": "Hz",
        "flux_density_unit": "T",
        "coefficient_flux_density_unit": "T",
        "harmonic_frequencies_hz": [50.0, 150.0, 250.0],
        "evaluated_harmonic_frequencies_hz": [50.0, 150.0, 250.0],
        "loss_coefficients": [1.0, 0.02, 0.001],
        "evaluated_loss_coefficients": [1.0, 0.02, 0.001],
        "loss_basis_sha256": "a" * 64,
        "evaluated_loss_basis_sha256": "a" * 64,
    }
    identity["demagnetization_temperature_current_phase_operating_point_identity"] = {
        "operating_point_generation": "operating-point-21",
        "temperature_operating_point_generation": "operating-point-21",
        "current_phase_operating_point_generation": "operating-point-21",
        "demag_margin_operating_point_generation": "operating-point-21",
        "magnet_temperature_c": 120.0,
        "demag_margin_temperature_c": 120.0,
        "current_phase_deg": 90.0,
        "demag_margin_current_phase_deg": 90.0,
        "operating_point_sha256": "b" * 64,
        "demag_margin_operating_point_sha256": "b" * 64,
    }
    return payload


def _payload_v20():
    payload = _payload_v19()
    identity = payload["artifact_identity"]
    identity["skew_slice_torque_angle_weight_periodicity_generation_identity"] = {
        "skew_generation": "skew-22",
        "torque_skew_generation": "skew-22",
        "angle_skew_generation": "skew-22",
        "weight_skew_generation": "skew-22",
        "periodicity_skew_generation": "skew-22",
        "slice_ids": [1, 2, 3],
        "torque_slice_ids": [1, 2, 3],
        "slice_angles_deg": [-5.0, 0.0, 5.0],
        "torque_slice_angles_deg": [-5.0, 0.0, 5.0],
        "quadrature_weights": [0.25, 0.5, 0.25],
        "torque_quadrature_weights": [0.25, 0.5, 0.25],
        "periodic_wrap_deg": 360.0,
        "torque_periodic_wrap_deg": 360.0,
        "skew_average_table_sha256": "1" * 64,
        "torque_skew_average_table_sha256": "1" * 64,
    }
    identity[
        "incremental_inductance_current_perturbation_phase_state_generation_identity"
    ] = {
        "operating_point_generation": "operating-point-22",
        "matrix_operating_point_generation": "operating-point-22",
        "perturbation_operating_point_generation": "operating-point-22",
        "phase_state_operating_point_generation": "operating-point-22",
        "base_solve_generation": "solve-22-base",
        "matrix_base_solve_generation": "solve-22-base",
        "perturbation_solve_generations": ["solve-22-a", "solve-22-b", "solve-22-c"],
        "matrix_perturbation_solve_generations": [
            "solve-22-a",
            "solve-22-b",
            "solve-22-c",
        ],
        "phase_names": ["a", "b", "c"],
        "matrix_phase_names": ["a", "b", "c"],
        "perturbation_currents_a": [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ],
        "matrix_perturbation_currents_a": [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ],
        "incremental_inductance_table_sha256": "2" * 64,
        "resolved_incremental_inductance_table_sha256": "2" * 64,
    }
    return payload


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


def _payload_v22():
    payload = _payload_v21()
    identity = payload["artifact_identity"]
    identity["motion_skew_force_harmonic_time_angle_phase_generation_identity"] = {
        "motion_study_generation": "motion-41",
        "time_motion_study_generation": "motion-41",
        "angle_motion_study_generation": "motion-41",
        "skew_motion_study_generation": "motion-41",
        "phase_motion_study_generation": "motion-41",
        "force_result_motion_study_generation": "motion-41",
        "time_s": [0.0, 0.001, 0.002],
        "force_time_s": [0.0, 0.001, 0.002],
        "mechanical_angle_deg": [0.0, 5.0, 10.0],
        "force_mechanical_angle_deg": [0.0, 5.0, 10.0],
        "skew_slice_angles_deg": [-5.0, 0.0, 5.0],
        "force_skew_slice_angles_deg": [-5.0, 0.0, 5.0],
        "slice_weights": [0.25, 0.5, 0.25],
        "force_slice_weights": [0.25, 0.5, 0.25],
        "phase_reference_deg": 30.0,
        "force_phase_reference_deg": 30.0,
        "harmonic_orders": [1, 3, 5],
        "force_harmonic_orders": [1, 3, 5],
        "force_harmonics_n": [120.0, 4.5, 1.2],
        "reported_force_harmonics_n": [120.0, 4.5, 1.2],
        "force_harmonic_table_sha256": "1" * 64,
        "resolved_force_harmonic_table_sha256": "1" * 64,
    }
    identity[
        "ipm_irreversible_demag_recoil_temperature_operating_generation_identity"
    ] = {
        "demag_study_generation": "demag-41",
        "recoil_curve_demag_study_generation": "demag-41",
        "temperature_demag_study_generation": "demag-41",
        "operating_point_demag_study_generation": "demag-41",
        "magnet_orientation_demag_study_generation": "demag-41",
        "result_demag_study_generation": "demag-41",
        "temperature_c": 120.0,
        "result_temperature_c": 120.0,
        "operating_point_id": "id=-180A;iq=240A;theta=17.5deg",
        "result_operating_point_id": "id=-180A;iq=240A;theta=17.5deg",
        "magnet_orientation_vectors": [[1.0, 0.0], [0.0, 1.0]],
        "result_magnet_orientation_vectors": [[1.0, 0.0], [0.0, 1.0]],
        "recoil_curve_sha256": "2" * 64,
        "result_recoil_curve_sha256": "2" * 64,
        "magnet_state_sha256": "3" * 64,
        "result_magnet_state_sha256": "3" * 64,
        "demag_margin_a_per_m": [175000.0, 82000.0],
        "reported_demag_margin_a_per_m": [175000.0, 82000.0],
    }
    return payload


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


def _payload_v24():
    payload = _payload_v23()
    identity = payload["artifact_identity"]
    identity[
        "loss_torque_speed_power_balance_harmonic_window_generation_identity"
    ] = {
        "power_balance_generation": "power-101",
        "torque_speed_power_balance_generation": "power-101",
        "harmonic_window_power_balance_generation": "power-101",
        "time_average_power_balance_generation": "power-101",
        "iron_loss_power_balance_generation": "power-101",
        "copper_loss_power_balance_generation": "power-101",
        "mechanical_loss_power_balance_generation": "power-101",
        "result_power_balance_generation": "power-101",
        "torque_nm": [1.0, 1.2],
        "power_balance_torque_nm": [1.0, 1.2],
        "speed_rad_s": [100.0, 100.0],
        "power_balance_speed_rad_s": [100.0, 100.0],
        "mechanical_output_w": [100.0, 120.0],
        "power_balance_mechanical_output_w": [100.0, 120.0],
        "harmonic_window_samples": [20, 120],
        "loss_harmonic_window_samples": [20, 120],
        "time_average_window_s": [0.02, 0.12],
        "loss_time_average_window_s": [0.02, 0.12],
        "iron_loss_w": [5.0, 6.0],
        "power_balance_iron_loss_w": [5.0, 6.0],
        "copper_loss_w": [3.0, 4.0],
        "power_balance_copper_loss_w": [3.0, 4.0],
        "mechanical_loss_w": [2.0, 2.0],
        "power_balance_mechanical_loss_w": [2.0, 2.0],
        "electrical_input_w": [110.0, 132.0],
        "power_balance_electrical_input_w": [110.0, 132.0],
        "power_balance_sha256": "1" * 64,
        "reported_power_balance_sha256": "1" * 64,
    }
    identity[
        "skew_slice_weight_rotor_angle_phase_periodicity_generation_identity"
    ] = {
        "skew_generation": "skew-101",
        "weight_skew_generation": "skew-101",
        "angle_skew_generation": "skew-101",
        "phase_skew_generation": "skew-101",
        "periodicity_skew_generation": "skew-101",
        "solve_skew_generation": "skew-101",
        "result_skew_generation": "skew-101",
        "slice_ids": [1, 2, 3],
        "result_slice_ids": [1, 2, 3],
        "quadrature_weights": [0.25, 0.5, 0.25],
        "result_quadrature_weights": [0.25, 0.5, 0.25],
        "rotor_angles_deg": [-5.0, 0.0, 5.0],
        "result_rotor_angles_deg": [-5.0, 0.0, 5.0],
        "current_phase_ids": ["abc@-5", "abc@0", "abc@5"],
        "result_current_phase_ids": ["abc@-5", "abc@0", "abc@5"],
        "periodic_map_ids": ["p1", "p2", "p3"],
        "result_periodic_map_ids": ["p1", "p2", "p3"],
        "slice_solve_sha256": ["2" * 64, "3" * 64, "4" * 64],
        "result_slice_solve_sha256": ["2" * 64, "3" * 64, "4" * 64],
        "slice_torque_nm": [0.8, 1.0, 1.2],
        "result_slice_torque_nm": [0.8, 1.0, 1.2],
        "weighted_torque_nm": 1.0,
        "reported_weighted_torque_nm": 1.0,
        "skew_result_sha256": "5" * 64,
        "reported_skew_result_sha256": "5" * 64,
    }
    return payload


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


def _payload_v28():
    payload = _payload_v27()
    identity = payload["artifact_identity"]
    identity[
        "pwm_current_harmonic_time_electrical_angle_torque_loss_mesh_result_generation_identity"
    ] = {
        "pwm_generation": "pwm-observables-151",
        "current_pwm_generation": "pwm-observables-151",
        "time_pwm_generation": "pwm-observables-151",
        "angle_pwm_generation": "pwm-observables-151",
        "torque_pwm_generation": "pwm-observables-151",
        "loss_pwm_generation": "pwm-observables-151",
        "mesh_pwm_generation": "pwm-observables-151",
        "result_pwm_generation": "pwm-observables-151",
        "harmonic_orders": [1, 5, 7, 11],
        "result_harmonic_orders": [1, 5, 7, 11],
        "current_harmonic_a": [100.0, 8.0, 5.0, 2.0],
        "result_current_harmonic_a": [100.0, 8.0, 5.0, 2.0],
        "current_phase_deg": [0.0, -20.0, 15.0, -5.0],
        "result_current_phase_deg": [0.0, -20.0, 15.0, -5.0],
        "time_s": [0.020, 0.021, 0.022, 0.023],
        "result_time_s": [0.020, 0.021, 0.022, 0.023],
        "electrical_angle_deg": [0.0, 72.0, 144.0, 216.0],
        "result_electrical_angle_deg": [0.0, 72.0, 144.0, 216.0],
        "pole_pairs": 4,
        "result_pole_pairs": 4,
        "torque_window_s": [0.020, 0.023],
        "result_torque_window_s": [0.020, 0.023],
        "torque_average_nm": 42.0,
        "result_torque_average_nm": 42.0,
        "loss_average_w": 350.0,
        "result_loss_average_w": 350.0,
        "mesh_sha256": "1" * 64,
        "result_mesh_sha256": "1" * 64,
        "result_sha256": "2" * 64,
        "accepted_result_sha256": "2" * 64,
    }
    identity[
        "skew_slice_angle_weight_frame_interpolation_torque_ripple_mesh_generation_identity"
    ] = {
        "skew_generation": "skew-average-151",
        "angle_skew_generation": "skew-average-151",
        "weight_skew_generation": "skew-average-151",
        "frame_skew_generation": "skew-average-151",
        "interpolation_skew_generation": "skew-average-151",
        "torque_skew_generation": "skew-average-151",
        "mesh_skew_generation": "skew-average-151",
        "result_skew_generation": "skew-average-151",
        "slice_angles_deg": [-3.0, 0.0, 3.0],
        "result_slice_angles_deg": [-3.0, 0.0, 3.0],
        "slice_weights": [0.25, 0.5, 0.25],
        "result_slice_weights": [0.25, 0.5, 0.25],
        "rotor_frame": "mechanical-ccw",
        "result_rotor_frame": "mechanical-ccw",
        "interpolation_rule": "periodic-cubic",
        "result_interpolation_rule": "periodic-cubic",
        "slice_torque_nm": [40.0, 44.0, 42.0],
        "result_slice_torque_nm": [40.0, 44.0, 42.0],
        "torque_average_nm": 42.5,
        "result_torque_average_nm": 42.5,
        "torque_ripple_nm": 4.0,
        "result_torque_ripple_nm": 4.0,
        "mesh_sha256": "3" * 64,
        "result_mesh_sha256": "3" * 64,
        "result_sha256": "4" * 64,
        "accepted_result_sha256": "4" * 64,
    }
    return payload


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


def _payload_v30():
    payload = _payload_v29()
    identity = payload["artifact_identity"]
    generation = "ipm-dq-171"
    identity["ipm_dq_inductance_current_angle_park_saturation_flux_derivative_reciprocity_mesh_result_identity"] = {
        "dq_generation": generation, "current_dq_generation": generation,
        "frame_dq_generation": generation, "saturation_dq_generation": generation,
        "flux_dq_generation": generation, "derivative_dq_generation": generation,
        "reciprocity_dq_generation": generation, "mesh_dq_generation": generation,
        "result_dq_generation": generation,
        "current_magnitude_a": 100.0, "result_current_magnitude_a": 100.0,
        "current_angle_electrical_deg": 30.0, "result_current_angle_electrical_deg": 30.0,
        "park_frame": "rotor_d_aligned_ccw_power_invariant",
        "result_park_frame": "rotor_d_aligned_ccw_power_invariant",
        "saturation_operating_point_a": [86.6025403784, 50.0],
        "result_saturation_operating_point_a": [86.6025403784, 50.0],
        "flux_linkage_derivative_h": [[0.003, 0.0002], [0.0002, 0.006]],
        "result_flux_linkage_derivative_h": [[0.003, 0.0002], [0.0002, 0.006]],
        "reciprocity_tolerance_h": 1.0e-9, "result_reciprocity_tolerance_h": 1.0e-9,
        "mesh_sha256": "1" * 64, "result_mesh_sha256": "1" * 64,
        "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64,
    }
    generation = "srm-coenergy-171"
    torque = (2.1 - 1.9) / math.radians(2.0)
    identity["srm_torque_current_position_coenergy_periodicity_phase_sequence_mesh_result_identity"] = {
        "srm_generation": generation, "current_srm_generation": generation,
        "position_srm_generation": generation, "coenergy_srm_generation": generation,
        "periodicity_srm_generation": generation, "phase_srm_generation": generation,
        "mesh_srm_generation": generation, "result_srm_generation": generation,
        "current_a": [0.0, 25.0, 50.0], "result_current_a": [0.0, 25.0, 50.0],
        "rotor_position_mechanical_deg": [-1.0, 0.0, 1.0],
        "result_rotor_position_mechanical_deg": [-1.0, 0.0, 1.0],
        "coenergy_j_at_50a": [1.9, 2.0, 2.1], "result_coenergy_j_at_50a": [1.9, 2.0, 2.1],
        "torque_nm_at_50a": torque, "result_torque_nm_at_50a": torque,
        "sector_period_mechanical_deg": 30.0, "result_sector_period_mechanical_deg": 30.0,
        "phase_sequence": ["A", "B", "C"], "result_phase_sequence": ["A", "B", "C"],
        "mesh_sha256": "3" * 64, "result_mesh_sha256": "3" * 64,
        "result_sha256": "4" * 64, "accepted_result_sha256": "4" * 64,
    }
    return payload


def _payload_v31():
    payload = _payload_v30(); identity = payload["artifact_identity"]
    generation = "pwm-loss-181"
    identity["pwm_iron_loss_sampling_sideband_angle_alias_volume_energy_result_identity"] = {
        "loss_generation": generation,
        **{key: generation for key in ("sampling_loss_generation", "sideband_loss_generation", "angle_loss_generation", "alias_loss_generation", "volume_loss_generation", "energy_loss_generation", "mesh_loss_generation", "result_loss_generation")},
        "sample_period_s": 2.5e-6, "result_sample_period_s": 2.5e-6,
        "samples_per_fundamental_cycle": 4000, "result_samples_per_fundamental_cycle": 4000,
        "carrier_frequency_hz": 10000.0, "result_carrier_frequency_hz": 10000.0,
        "fundamental_frequency_hz": 100.0, "result_fundamental_frequency_hz": 100.0,
        "carrier_sidebands_hz": [9900.0, 10100.0], "result_carrier_sidebands_hz": [9900.0, 10100.0],
        "pole_pairs": 4, "result_pole_pairs": 4,
        "mechanical_angle_deg": [0.0, 5.0, 10.0], "result_mechanical_angle_deg": [0.0, 5.0, 10.0],
        "electrical_angle_deg": [0.0, 20.0, 40.0], "result_electrical_angle_deg": [0.0, 20.0, 40.0],
        "alias_filter": "nyquist_guard_and_sideband_keep", "result_alias_filter": "nyquist_guard_and_sideband_keep",
        "gross_volume_m3": 0.001, "result_gross_volume_m3": 0.001,
        "active_volume_m3": 0.00095, "result_active_volume_m3": 0.00095,
        "stacking_factor": 0.95, "result_stacking_factor": 0.95,
        "mean_iron_loss_w": 50.0, "result_mean_iron_loss_w": 50.0,
        "cycle_energy_j": 0.5, "result_cycle_energy_j": 0.5,
        "mesh_sha256": "1" * 64, "result_mesh_sha256": "1" * 64,
        "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64,
    }
    generation = "skew-torque-181"
    identity["skew_slice_torque_weight_axial_phase_periodicity_ripple_mesh_result_identity"] = {
        "skew_generation": generation,
        **{key: generation for key in ("weight_skew_generation", "axial_skew_generation", "phase_skew_generation", "periodicity_skew_generation", "torque_skew_generation", "ripple_skew_generation", "mesh_skew_generation", "result_skew_generation")},
        "axial_locations_m": [-0.01, 0.0, 0.01], "result_axial_locations_m": [-0.01, 0.0, 0.01],
        "slice_weights": [0.25, 0.5, 0.25], "result_slice_weights": [0.25, 0.5, 0.25],
        "electrical_phase_offsets_deg": [-2.0, 0.0, 2.0], "result_electrical_phase_offsets_deg": [-2.0, 0.0, 2.0],
        "periodic_wrap_electrical_deg": 360.0, "result_periodic_wrap_electrical_deg": 360.0,
        "slice_mean_torque_nm": [10.0, 10.2, 10.0], "result_slice_mean_torque_nm": [10.0, 10.2, 10.0],
        "skew_mean_torque_nm": 10.1, "result_skew_mean_torque_nm": 10.1,
        "slice_ripple_harmonics_nm": [{"6": 0.3}, {"6": 0.2}, {"6": 0.3}],
        "result_slice_ripple_harmonics_nm": [{"6": 0.3}, {"6": 0.2}, {"6": 0.3}],
        "skew_ripple_harmonics_nm": {"6": 0.25}, "result_skew_ripple_harmonics_nm": {"6": 0.25},
        "mesh_sha256": "3" * 64, "result_mesh_sha256": "3" * 64,
        "result_sha256": "4" * 64, "accepted_result_sha256": "4" * 64,
    }
    return payload


def _payload_v32():
    payload = _payload_v31()
    identity = payload["artifact_identity"]
    generation = "ipm-demag-191"
    identity[
        "ipm_demagnetization_knee_temperature_current_angle_region_fraction_mesh_result_identity"
    ] = {
        "demag_generation": generation,
        **{
            key: generation
            for key in (
                "knee_demag_generation",
                "temperature_demag_generation",
                "current_demag_generation",
                "angle_demag_generation",
                "region_demag_generation",
                "fraction_demag_generation",
                "mesh_demag_generation",
                "result_demag_generation",
            )
        },
        "knee_criterion": "b_parallel_below_temperature_knee",
        "result_knee_criterion": "b_parallel_below_temperature_knee",
        "magnet_temperature_c": 120.0,
        "result_magnet_temperature_c": 120.0,
        "phase_current_rms_a": 200.0,
        "result_phase_current_rms_a": 200.0,
        "current_angle_electrical_deg": 135.0,
        "result_current_angle_electrical_deg": 135.0,
        "irreversible_region_labels": ["magnet_1/edge", "magnet_2/edge"],
        "result_irreversible_region_labels": ["magnet_1/edge", "magnet_2/edge"],
        "demagnetized_fraction": 0.015,
        "result_demagnetized_fraction": 0.015,
        "operating_point_owner": "ipm/case-191/120C/200A/135deg",
        "result_operating_point_owner": "ipm/case-191/120C/200A/135deg",
        "mesh_sha256": "1" * 64,
        "result_mesh_sha256": "1" * 64,
        "result_sha256": "2" * 64,
        "accepted_result_sha256": "2" * 64,
    }
    generation = "synrm-dq-191"
    rows = [
        {"id_a": -100.0, "iq_a": 100.0, "psi_d_wb": 0.18, "psi_q_wb": 0.07},
        {"id_a": -50.0, "iq_a": 150.0, "psi_d_wb": 0.21, "psi_q_wb": 0.09},
        {"id_a": 0.0, "iq_a": 180.0, "psi_d_wb": 0.24, "psi_q_wb": 0.11},
    ]
    identity[
        "synrm_dq_map_angle_saturation_cross_coupling_mtpa_torque_mesh_result_identity"
    ] = {
        "map_generation": generation,
        **{
            key: generation
            for key in (
                "angle_map_generation",
                "saturation_map_generation",
                "cross_coupling_map_generation",
                "mtpa_map_generation",
                "torque_map_generation",
                "mesh_map_generation",
                "result_map_generation",
            )
        },
        "electrical_angle_deg": [0.0, 30.0, 60.0, 90.0],
        "result_electrical_angle_deg": [0.0, 30.0, 60.0, 90.0],
        "saturation_branch": "nonlinear_forward",
        "result_saturation_branch": "nonlinear_forward",
        "flux_map_rows": rows,
        "result_flux_map_rows": [dict(row) for row in rows],
        "dpsi_d_diq_h": [1.5e-4, 1.7e-4, 1.9e-4],
        "result_dpsi_d_diq_h": [1.5e-4, 1.7e-4, 1.9e-4],
        "dpsi_q_did_h": [1.5e-4, 1.7e-4, 1.9e-4],
        "result_dpsi_q_did_h": [1.5e-4, 1.7e-4, 1.9e-4],
        "mtpa_row_indices": [0, 1, 2],
        "result_mtpa_row_indices": [0, 1, 2],
        "pole_pairs": 2,
        "result_pole_pairs": 2,
        "torque_reconstruction": "1.5*p*(psi_d*iq-psi_q*id)",
        "result_torque_reconstruction": "1.5*p*(psi_d*iq-psi_q*id)",
        "torque_nm": [75.0, 81.0, 86.4],
        "result_torque_nm": [75.0, 81.0, 86.4],
        "mesh_sha256": "3" * 64,
        "result_mesh_sha256": "3" * 64,
        "map_sha256": "4" * 64,
        "accepted_map_sha256": "4" * 64,
    }
    return payload


def _payload_v33():
    payload = _payload_v32()
    identity = payload["artifact_identity"]
    generation = "srm-commutation-201"
    identity[
        "srm_commutation_phase_dwell_chop_overlap_coenergy_torque_loss_angle_mesh_result_identity"
    ] = {
        "srm_generation": generation,
        **{
            key: generation
            for key in (
                "phase_generation",
                "dwell_generation",
                "chop_generation",
                "overlap_generation",
                "coenergy_generation",
                "torque_generation",
                "loss_generation",
                "angle_generation",
                "mesh_generation",
                "result_generation",
            )
        },
        "phase_sequence": ["A", "B", "C"],
        "result_phase_sequence": ["A", "B", "C"],
        "turn_on_deg": [0.0, 30.0, 60.0],
        "result_turn_on_deg": [0.0, 30.0, 60.0],
        "turn_off_deg": [20.0, 50.0, 80.0],
        "result_turn_off_deg": [20.0, 50.0, 80.0],
        "current_chop_a": 100.0,
        "result_current_chop_a": 100.0,
        "overlap_deg": 5.0,
        "result_overlap_deg": 5.0,
        "angle_grid_rad": [0.0, 0.1, 0.2],
        "result_angle_grid_rad": [0.0, 0.1, 0.2],
        "coenergy_j": [0.0, 0.5, 1.0],
        "result_coenergy_j": [0.0, 0.5, 1.0],
        "torque_nm": [5.0, 5.0, 5.0],
        "result_torque_nm": [5.0, 5.0, 5.0],
        "copper_loss_w": 120.0,
        "result_copper_loss_w": 120.0,
        "iron_loss_w": 30.0,
        "result_iron_loss_w": 30.0,
        "total_loss_w": 150.0,
        "result_total_loss_w": 150.0,
        "mesh_sha256": "1" * 64,
        "result_mesh_sha256": "1" * 64,
        "result_sha256": "2" * 64,
        "accepted_result_sha256": "2" * 64,
    }
    generation = "axial-flux-pm-201"
    identity[
        "axial_flux_pm_sector_airgap_end_effect_torque_force_surface_direction_frame_mesh_result_identity"
    ] = {
        "axial_generation": generation,
        **{
            key: generation
            for key in (
                "sector_generation",
                "airgap_generation",
                "end_effect_generation",
                "torque_generation",
                "force_generation",
                "direction_generation",
                "frame_generation",
                "mesh_generation",
                "result_generation",
            )
        },
        "sector_multiplier": 12,
        "result_sector_multiplier": 12,
        "air_gaps_m": [0.001, 0.001],
        "result_air_gaps_m": [0.001, 0.001],
        "end_effect_factor": 0.96,
        "result_end_effect_factor": 0.96,
        "surface_coordinates": [[0.0, 0.0], [10.0, 100.0], [20.0, 200.0]],
        "result_surface_coordinates": [
            [0.0, 0.0],
            [10.0, 100.0],
            [20.0, 200.0],
        ],
        "torque_surface_nm": [40.0, 55.0, 60.0],
        "result_torque_surface_nm": [40.0, 55.0, 60.0],
        "axial_force_surface_n": [0.0, 12.0, 18.0],
        "result_axial_force_surface_n": [0.0, 12.0, 18.0],
        "force_direction": "+z",
        "result_force_direction": "+z",
        "axial_frame": "rotor_global_z",
        "result_axial_frame": "rotor_global_z",
        "mesh_sha256": "3" * 64,
        "result_mesh_sha256": "3" * 64,
        "result_lineage_sha256": "4" * 64,
        "accepted_result_lineage_sha256": "4" * 64,
    }
    return payload


def _payload_v34():
    payload = _payload_v33()
    identity = payload["artifact_identity"]

    generation = "pm-demag-operating-point-211"
    identity[
        "pm_demagnetization_temperature_recoil_loadline_operating_point_knee_margin_angle_mesh_owner_result_identity"
    ] = {
        "demag_generation": generation,
        **{
            key: generation
            for key in (
                "temperature_generation",
                "recoil_generation",
                "loadline_generation",
                "operating_point_generation",
                "knee_generation",
                "margin_generation",
                "angle_generation",
                "mesh_generation",
                "owner_generation",
                "result_generation",
            )
        },
        "reference_temperature_c": 20.0,
        "result_reference_temperature_c": 20.0,
        "operating_temperature_c": 120.0,
        "result_operating_temperature_c": 120.0,
        "remanence_reference_t": 1.2,
        "result_remanence_reference_t": 1.2,
        "remanence_temperature_coefficient_per_c": -0.001,
        "result_remanence_temperature_coefficient_per_c": -0.001,
        "remanence_operating_t": 1.08,
        "result_remanence_operating_t": 1.08,
        "coercivity_reference_a_m": 900000.0,
        "result_coercivity_reference_a_m": 900000.0,
        "coercivity_temperature_coefficient_per_c": -0.002,
        "result_coercivity_temperature_coefficient_per_c": -0.002,
        "coercivity_operating_a_m": 720000.0,
        "result_coercivity_operating_a_m": 720000.0,
        "recoil_permeability_relative": 1.05,
        "result_recoil_permeability_relative": 1.05,
        "loadline_slope_t_per_a_m": 1.0e-6,
        "result_loadline_slope_t_per_a_m": 1.0e-6,
        "operating_field_a_m": -600000.0,
        "result_operating_field_a_m": -600000.0,
        "operating_flux_density_t": 0.48,
        "result_operating_flux_density_t": 0.48,
        "knee_field_a_m": -700000.0,
        "result_knee_field_a_m": -700000.0,
        "irreversible_margin_a_m": 100000.0,
        "result_irreversible_margin_a_m": 100000.0,
        "rotor_angle_rad": 0.3,
        "result_rotor_angle_rad": 0.3,
        "demag_mesh_sha256": "1" * 64,
        "result_demag_mesh_sha256": "1" * 64,
        "demag_result_owner": "motor/pm-demag-211",
        "accepted_demag_result_owner": "motor/pm-demag-211",
        "demag_result_sha256": "2" * 64,
        "accepted_demag_result_sha256": "2" * 64,
    }

    generation = "eccentricity-ump-211"
    angle_grid = [index * math.pi / 8.0 for index in range(5)]
    identity[
        "eccentricity_static_dynamic_frame_radial_force_harmonic_ump_torque_pole_periodicity_angle_owner_result_identity"
    ] = {
        "eccentricity_generation": generation,
        **{
            key: generation
            for key in (
                "static_generation",
                "dynamic_generation",
                "frame_generation",
                "harmonic_generation",
                "force_generation",
                "torque_generation",
                "periodicity_generation",
                "angle_generation",
                "owner_generation",
                "result_generation",
            )
        },
        "static_eccentricity_m": [1.0e-4, 0.0],
        "result_static_eccentricity_m": [1.0e-4, 0.0],
        "dynamic_eccentricity_amplitude_m": 5.0e-5,
        "result_dynamic_eccentricity_amplitude_m": 5.0e-5,
        "mechanical_frame": "stator_global_xy",
        "result_mechanical_frame": "stator_global_xy",
        "radial_force_harmonics_n": [[0, 0.0, 0.0], [1, 100.0, 0.0], [2, 0.0, 0.0]],
        "result_radial_force_harmonics_n": [[0, 0.0, 0.0], [1, 100.0, 0.0], [2, 0.0, 0.0]],
        "unbalanced_magnetic_pull_n": [100.0, 0.0],
        "result_unbalanced_magnetic_pull_n": [100.0, 0.0],
        "torque_nm": 20.0,
        "result_torque_nm": 20.0,
        "pole_pairs": 4,
        "result_pole_pairs": 4,
        "periodicity_angle_rad": math.pi / 2.0,
        "result_periodicity_angle_rad": math.pi / 2.0,
        "angle_grid_rad": angle_grid,
        "result_angle_grid_rad": angle_grid,
        "eccentricity_mesh_sha256": "3" * 64,
        "result_eccentricity_mesh_sha256": "3" * 64,
        "eccentricity_result_owner": "motor/eccentricity-211",
        "accepted_eccentricity_result_owner": "motor/eccentricity-211",
        "eccentricity_result_sha256": "4" * 64,
        "accepted_eccentricity_result_sha256": "4" * 64,
    }
    return payload


def _payload_v35():
    payload = _payload_v34()
    identity = payload["artifact_identity"]
    generation = "skew-slice-torque-235"
    angles = [-0.1, 0.0, 0.1]
    phases = [[harmonic, [harmonic * 4 * angle for angle in angles]] for harmonic in [1, 6]]
    identity["skew_slice_torque_angle_axial_weight_harmonic_phase_pole_periodicity_mean_ripple_mesh_owner_result_identity"] = {
        "skew_generation": generation,
        **{key: generation for key in (
            "slice_generation", "weight_generation", "phase_generation", "periodicity_generation",
            "torque_generation", "ripple_generation", "mesh_generation", "owner_generation",
            "result_generation")},
        "slice_angles_rad": angles, "result_slice_angles_rad": angles,
        "axial_weights": [0.25, 0.5, 0.25], "result_axial_weights": [0.25, 0.5, 0.25],
        "harmonic_phase_shifts_rad": phases, "result_harmonic_phase_shifts_rad": phases,
        "pole_pairs": 4, "result_pole_pairs": 4,
        "pole_periodicity_angle_rad": math.pi / 2.0, "result_pole_periodicity_angle_rad": math.pi / 2.0,
        "slice_mean_torque_nm": [47.0, 50.0, 49.0], "result_slice_mean_torque_nm": [47.0, 50.0, 49.0],
        "weighted_mean_torque_nm": 49.0, "result_weighted_mean_torque_nm": 49.0,
        "torque_ripple_spectrum_nm": [[0, 49.0], [6, 2.0]],
        "result_torque_ripple_spectrum_nm": [[0, 49.0], [6, 2.0]],
        "skew_mesh_sha256": "1" * 64, "result_skew_mesh_sha256": "1" * 64,
        "skew_result_owner": "motor/skew-235", "accepted_skew_result_owner": "motor/skew-235",
        "skew_result_sha256": "2" * 64, "accepted_skew_result_sha256": "2" * 64,
    }
    generation = "iron-loss-separation-235"
    frequency, b_peak, waveform_factor, temperature_factor = 100.0, 1.2, 1.1, 1.1
    components = [
        2.0 * frequency * b_peak**2 * waveform_factor * temperature_factor,
        0.1 * frequency**2 * b_peak**2 * waveform_factor,
        0.05 * frequency**1.5 * b_peak**1.5 * waveform_factor,
    ]
    identity["ironloss_hysteresis_eddy_excess_waveform_frequency_coeff_volume_temperature_total_owner_result_identity"] = {
        "ironloss_generation": generation,
        **{key: generation for key in (
            "component_generation", "waveform_generation", "frequency_generation", "coefficient_generation",
            "volume_generation", "temperature_generation", "total_generation", "owner_generation",
            "result_generation")},
        "b_waveform_peak_t": b_peak, "result_b_waveform_peak_t": b_peak,
        "waveform_factor": waveform_factor, "result_waveform_factor": waveform_factor,
        "frequency_hz": frequency, "result_frequency_hz": frequency,
        "material_coefficients": [2.0, 0.1, 0.05], "result_material_coefficients": [2.0, 0.1, 0.05],
        "active_volume_m3": 0.001, "result_active_volume_m3": 0.001,
        "temperature_c": 80.0, "result_temperature_c": 80.0,
        "temperature_factor": temperature_factor, "result_temperature_factor": temperature_factor,
        "loss_components_w_m3": components, "result_loss_components_w_m3": components,
        "total_iron_loss_w": sum(components) * 0.001, "result_total_iron_loss_w": sum(components) * 0.001,
        "waveform_sha256": "3" * 64, "result_waveform_sha256": "3" * 64,
        "ironloss_owner": "motor/iron-loss-235", "accepted_ironloss_owner": "motor/iron-loss-235",
        "ironloss_result_sha256": "4" * 64, "accepted_ironloss_result_sha256": "4" * 64,
    }
    return payload


def _payload_v36():
    payload = _payload_v35()
    identity = payload["artifact_identity"]
    generation = "dq-mtpa-236"
    pole_pairs = 4
    i_d, i_q = -50.0, 100.0
    l_d, l_q, flux_pm = 0.001, 0.0015, 0.075
    flux = [flux_pm + l_d * i_d, l_q * i_q]
    torque = 1.5 * pole_pairs * (flux[0] * i_q - flux[1] * i_d)
    identity["dq_flux_inductance_torque_mtpa_current_angle_speed_convention_owner_result_identity"] = {
        "dq_generation": generation,
        **{key: generation for key in (
            "park_generation", "flux_generation", "inductance_generation", "torque_generation",
            "mtpa_generation", "current_generation", "angle_generation", "speed_generation",
            "owner_generation", "result_generation")},
        "park_convention": "power_invariant_q_leads_d", "result_park_convention": "power_invariant_q_leads_d",
        "pole_pairs": pole_pairs, "result_pole_pairs": pole_pairs,
        "current_dq_a": [i_d, i_q], "result_current_dq_a": [i_d, i_q],
        "current_magnitude_a": math.hypot(i_d, i_q), "result_current_magnitude_a": math.hypot(i_d, i_q),
        "current_angle_rad": math.atan2(i_q, i_d), "result_current_angle_rad": math.atan2(i_q, i_d),
        "pm_flux_linkage_wb_turn": flux_pm, "result_pm_flux_linkage_wb_turn": flux_pm,
        "flux_linkage_dq_wb_turn": flux, "result_flux_linkage_dq_wb_turn": flux,
        "differential_inductance_dq_h": [l_d, l_q], "result_differential_inductance_dq_h": [l_d, l_q],
        "torque_nm": torque, "result_torque_nm": torque,
        "mechanical_speed_rad_s": 100.0, "result_mechanical_speed_rad_s": 100.0,
        "electrical_speed_rad_s": 400.0, "result_electrical_speed_rad_s": 400.0,
        "dq_owner": "motor/dq-236", "accepted_dq_owner": "motor/dq-236",
        "dq_result_sha256": "1" * 64, "accepted_dq_result_sha256": "1" * 64,
    }
    generation = "iron-loss-energy-236"
    frequency, flux_peak, temperature_factor = 100.0, 1.2, 1.1
    coefficients = [2.0, 0.1, 0.05]
    components = [
        coefficients[0] * frequency * flux_peak**2 * temperature_factor,
        coefficients[1] * frequency**2 * flux_peak**2,
        coefficients[2] * frequency**1.5 * flux_peak**1.5,
    ]
    regions = [["stator", 0.0007], ["rotor", 0.0003]]
    total_power = sum(components) * sum(row[1] for row in regions)
    identity["iron_loss_component_frequency_flux_region_thermal_energy_balance_owner_result_identity"] = {
        "iron_loss_generation": generation,
        **{key: generation for key in (
            "component_generation", "frequency_generation", "flux_generation", "region_generation",
            "thermal_generation", "power_generation", "energy_generation", "owner_generation",
            "result_generation")},
        "frequency_hz": frequency, "result_frequency_hz": frequency,
        "flux_peak_t": flux_peak, "result_flux_peak_t": flux_peak,
        "loss_coefficients": coefficients, "result_loss_coefficients": coefficients,
        "loss_components_w_m3": components, "result_loss_components_w_m3": components,
        "regional_volumes_m3": regions, "result_regional_volumes_m3": regions,
        "temperature_c": 80.0, "result_temperature_c": 80.0,
        "temperature_factor": temperature_factor, "result_temperature_factor": temperature_factor,
        "total_iron_loss_w": total_power, "result_total_iron_loss_w": total_power,
        "integration_duration_s": 0.2, "result_integration_duration_s": 0.2,
        "loss_energy_j": total_power * 0.2, "result_loss_energy_j": total_power * 0.2,
        "iron_loss_owner": "motor/iron-loss-236", "accepted_iron_loss_owner": "motor/iron-loss-236",
        "iron_loss_result_sha256": "2" * 64, "accepted_iron_loss_result_sha256": "2" * 64,
    }
    return payload


def _payload_v37():
    payload = _payload_v36()
    identity = payload["artifact_identity"]
    generation = "induction-power-246"
    frequency, pole_pairs, slip = 50.0, 2, 0.05
    synchronous_speed = 2.0 * math.pi * frequency / pole_pairs
    mechanical_speed = (1.0 - slip) * synchronous_speed
    airgap_power, input_power = 1000.0, 1200.0
    torque = airgap_power / synchronous_speed
    rotor_loss = slip * airgap_power
    mechanical_output = (1.0 - slip) * airgap_power
    identity[
        "induction_motor_slip_synchronous_speed_airgap_power_torque_rotor_loss_mechanical_output_efficiency_owner_result_identity"
    ] = {
        "induction_generation": generation,
        **{key: generation for key in (
            "frequency_generation", "speed_generation", "slip_generation",
            "airgap_generation", "torque_generation", "rotor_loss_generation",
            "mechanical_generation", "efficiency_generation", "owner_generation",
            "result_generation")},
        "electrical_frequency_hz": frequency, "result_electrical_frequency_hz": frequency,
        "pole_pairs": pole_pairs, "result_pole_pairs": pole_pairs,
        "slip": slip, "result_slip": slip,
        "synchronous_speed_rad_s": synchronous_speed,
        "result_synchronous_speed_rad_s": synchronous_speed,
        "mechanical_speed_rad_s": mechanical_speed,
        "result_mechanical_speed_rad_s": mechanical_speed,
        "airgap_power_w": airgap_power, "result_airgap_power_w": airgap_power,
        "electromagnetic_torque_nm": torque,
        "result_electromagnetic_torque_nm": torque,
        "rotor_copper_loss_w": rotor_loss,
        "result_rotor_copper_loss_w": rotor_loss,
        "mechanical_output_w": mechanical_output,
        "result_mechanical_output_w": mechanical_output,
        "input_power_w": input_power, "result_input_power_w": input_power,
        "efficiency": mechanical_output / input_power,
        "result_efficiency": mechanical_output / input_power,
        "motor_owner": "motor/induction-246",
        "accepted_motor_owner": "motor/induction-246",
        "motor_result_sha256": "1" * 64,
        "accepted_motor_result_sha256": "1" * 64,
    }
    generation = "axial-flux-246"
    sector_factor = 12
    sector_torque = [2.0, 2.1, 1.9]
    full_torque = [sector_factor * item for item in sector_torque]
    average_torque = sum(full_torque) / len(full_torque)
    ripple = (max(full_torque) - min(full_torque)) / average_torque
    phases = [0.0, -2.0 * math.pi / 3.0, 2.0 * math.pi / 3.0]
    identity[
        "axial_flux_motor_sector_periodicity_dual_airgap_axial_flux_torque_ripple_backemf_frame_mesh_owner_result_identity"
    ] = {
        "axial_flux_generation": generation,
        **{key: generation for key in (
            "sector_generation", "airgap_generation", "flux_generation",
            "torque_generation", "ripple_generation", "backemf_generation",
            "frame_generation", "mesh_generation", "owner_generation",
            "result_generation")},
        "sector_factor": sector_factor, "result_sector_factor": sector_factor,
        "sector_angle_rad": 2.0 * math.pi / sector_factor,
        "result_sector_angle_rad": 2.0 * math.pi / sector_factor,
        "dual_airgap_m": [0.001, 0.001],
        "result_dual_airgap_m": [0.001, 0.001],
        "sector_axial_flux_per_gap_wb": [0.002, 0.002],
        "result_sector_axial_flux_per_gap_wb": [0.002, 0.002],
        "sector_torque_samples_nm": sector_torque,
        "result_sector_torque_samples_nm": sector_torque,
        "full_machine_torque_samples_nm": full_torque,
        "result_full_machine_torque_samples_nm": full_torque,
        "average_torque_nm": average_torque,
        "result_average_torque_nm": average_torque,
        "torque_ripple_ratio": ripple, "result_torque_ripple_ratio": ripple,
        "backemf_phase_angles_rad": phases,
        "result_backemf_phase_angles_rad": phases,
        "coordinate_frame": "cylindrical_z_axial",
        "result_coordinate_frame": "cylindrical_z_axial",
        "mesh_owner": "mesh/axial-flux-246",
        "accepted_mesh_owner": "mesh/axial-flux-246",
        "axial_flux_result_sha256": "2" * 64,
        "accepted_axial_flux_result_sha256": "2" * 64,
    }
    return payload


def _payload_v38():
    payload = _payload_v37()
    identity = payload["artifact_identity"]
    generation = "wound-field-258"
    field_current = 5.0
    excitation_inductance = 0.04
    excitation_flux = field_current * excitation_inductance
    torque_angle = math.pi / 6.0
    pole_pairs = 2
    stator_current = 10.0
    torque = (
        1.5
        * pole_pairs
        * excitation_flux
        * math.sqrt(2.0)
        * stator_current
        * math.sin(torque_angle)
    )
    speed = 100.0
    mechanical_power = torque * speed
    field_resistance = 1.0
    stator_resistance = 0.2
    field_loss = field_current**2 * field_resistance
    stator_loss = 3.0 * stator_current**2 * stator_resistance
    active_power = mechanical_power + field_loss + stator_loss
    line_voltage = 40.0
    apparent_power = math.sqrt(3.0) * line_voltage * stator_current
    identity[
        "wound_field_synchronous_excitation_flux_torque_angle_powerfactor_field_stator_loss_mechanical_energy_mesh_owner_result_identity"
    ] = {
        "wound_field_generation": generation,
        **{
            key: generation
            for key in (
                "excitation_generation", "flux_generation",
                "torque_generation", "powerfactor_generation",
                "field_loss_generation", "stator_loss_generation",
                "mechanical_generation", "energy_generation",
                "mesh_generation", "owner_generation", "result_generation",
            )
        },
        "field_current_a": field_current,
        "result_field_current_a": field_current,
        "excitation_inductance_h": excitation_inductance,
        "result_excitation_inductance_h": excitation_inductance,
        "excitation_flux_linkage_wb_turn": excitation_flux,
        "result_excitation_flux_linkage_wb_turn": excitation_flux,
        "torque_angle_rad": torque_angle,
        "result_torque_angle_rad": torque_angle,
        "pole_pairs": pole_pairs,
        "result_pole_pairs": pole_pairs,
        "stator_current_rms_a": stator_current,
        "result_stator_current_rms_a": stator_current,
        "electromagnetic_torque_nm": torque,
        "result_electromagnetic_torque_nm": torque,
        "mechanical_speed_rad_s": speed,
        "result_mechanical_speed_rad_s": speed,
        "mechanical_power_w": mechanical_power,
        "result_mechanical_power_w": mechanical_power,
        "field_resistance_ohm": field_resistance,
        "result_field_resistance_ohm": field_resistance,
        "field_copper_loss_w": field_loss,
        "result_field_copper_loss_w": field_loss,
        "stator_phase_resistance_ohm": stator_resistance,
        "result_stator_phase_resistance_ohm": stator_resistance,
        "stator_copper_loss_w": stator_loss,
        "result_stator_copper_loss_w": stator_loss,
        "line_voltage_rms_v": line_voltage,
        "result_line_voltage_rms_v": line_voltage,
        "apparent_power_va": apparent_power,
        "result_apparent_power_va": apparent_power,
        "active_input_power_w": active_power,
        "result_active_input_power_w": active_power,
        "power_factor": active_power / apparent_power,
        "result_power_factor": active_power / apparent_power,
        "energy_balance_residual_w": 0.0,
        "result_energy_balance_residual_w": 0.0,
        "energy_tolerance_w": 1.0e-9,
        "result_energy_tolerance_w": 1.0e-9,
        "mesh_owner": "mesh:wound-field-258",
        "accepted_mesh_owner": "mesh:wound-field-258",
        "motor_result_sha256": "1" * 64,
        "accepted_motor_result_sha256": "1" * 64,
    }

    generation = "flux-switching-258"
    torque_samples = [10.0, 11.0, 9.0, 10.0]
    average_torque = sum(torque_samples) / len(torque_samples)
    torque_ripple = (max(torque_samples) - min(torque_samples)) / average_torque
    phase_angles = [0.0, -2.0 * math.pi / 3.0, 2.0 * math.pi / 3.0]
    identity[
        "flux_switching_pm_slot_pole_polarity_phase_harmonic_backemf_torque_ripple_periodicity_mesh_owner_result_identity"
    ] = {
        "flux_switching_generation": generation,
        **{
            key: generation
            for key in (
                "slot_pole_generation", "polarity_generation",
                "phase_generation", "harmonic_generation",
                "backemf_generation", "torque_generation",
                "periodicity_generation", "mesh_generation",
                "owner_generation", "result_generation",
            )
        },
        "slot_count": 12,
        "result_slot_count": 12,
        "pole_count": 10,
        "result_pole_count": 10,
        "magnet_polarity_sequence": [1, -1] * 5,
        "result_magnet_polarity_sequence": [1, -1] * 5,
        "phase_sequence": "ABC",
        "result_phase_sequence": "ABC",
        "working_harmonic_order": 5,
        "result_working_harmonic_order": 5,
        "backemf_phase_angles_rad": phase_angles,
        "result_backemf_phase_angles_rad": phase_angles,
        "torque_samples_nm": torque_samples,
        "result_torque_samples_nm": torque_samples,
        "average_torque_nm": average_torque,
        "result_average_torque_nm": average_torque,
        "torque_ripple_ratio": torque_ripple,
        "result_torque_ripple_ratio": torque_ripple,
        "periodic_multiplier": 2,
        "result_periodic_multiplier": 2,
        "sector_slot_count": 6,
        "result_sector_slot_count": 6,
        "sector_pole_count": 5,
        "result_sector_pole_count": 5,
        "mesh_owner": "mesh:flux-switching-258",
        "accepted_mesh_owner": "mesh:flux-switching-258",
        "flux_switching_result_sha256": "2" * 64,
        "accepted_flux_switching_result_sha256": "2" * 64,
    }
    return payload


_SKEW = "skewed_rotor_slice_angle_phase_weight_torque_ripple_power_model_owner_result_identity"


_DEMAG = "pm_irreversible_demag_temperature_recoil_knee_operating_flux_torque_mesh_owner_result_identity"


def _payload_v39():
    payload = _payload_v38()
    identity = payload["artifact_identity"]
    generation = "skewed-rotor-271"
    angles = [-5.0, 0.0, 5.0]
    pole_pairs = 2
    phases = [pole_pairs * angle for angle in angles]
    weights = [0.25, 0.5, 0.25]
    torque = [10.0, 10.4, 10.0]
    phasors = [[math.cos(math.radians(phase)), math.sin(math.radians(phase))] for phase in phases]
    mean_torque = sum(weight * item for weight, item in zip(weights, torque))
    ripple = math.hypot(
        sum(weight * pair[0] for weight, pair in zip(weights, phasors)),
        sum(weight * pair[1] for weight, pair in zip(weights, phasors)),
    )
    identity[_SKEW] = {
        "skew_generation": generation,
        **{key: generation for key in ("slice_generation", "phase_generation", "weight_generation", "torque_generation", "ripple_generation", "power_generation", "owner_generation", "result_generation")},
        "pole_pairs": pole_pairs,
        "result_pole_pairs": pole_pairs,
        "slice_angles_mechanical_deg": angles,
        "result_slice_angles_mechanical_deg": angles,
        "slice_phase_offsets_electrical_deg": phases,
        "result_slice_phase_offsets_electrical_deg": phases,
        "axial_weights": weights,
        "result_axial_weights": weights,
        "slice_mean_torque_nm": torque,
        "result_slice_mean_torque_nm": torque,
        "weighted_mean_torque_nm": mean_torque,
        "result_weighted_mean_torque_nm": mean_torque,
        "slice_ripple_phasor": phasors,
        "result_slice_ripple_phasor": phasors,
        "weighted_ripple_residual": ripple,
        "result_weighted_ripple_residual": ripple,
        "mechanical_speed_rad_s": 100.0,
        "result_mechanical_speed_rad_s": 100.0,
        "mechanical_power_w": mean_torque * 100.0,
        "result_mechanical_power_w": mean_torque * 100.0,
        "model_owner": "motor:skewed-rotor-271",
        "accepted_model_owner": "motor:skewed-rotor-271",
        "skew_result_sha256": "1" * 64,
        "accepted_skew_result_sha256": "1" * 64,
    }

    generation = "pm-demag-271"
    temperature = 120.0
    br_reference = 1.2
    coefficient = -1.1e-3
    br_temperature = br_reference * (1.0 + coefficient * (temperature - 20.0))
    recoil_mu = 1.05
    operating_h = -8.0e5
    operating_b = br_temperature + 4.0e-7 * math.pi * recoil_mu * operating_h
    loss = 0.08
    identity[_DEMAG] = {
        "demag_generation": generation,
        **{key: generation for key in ("temperature_generation", "recoil_generation", "knee_generation", "operating_generation", "remanence_generation", "flux_generation", "torque_generation", "mesh_generation", "owner_generation", "result_generation")},
        "reference_temperature_c": 20.0,
        "result_reference_temperature_c": 20.0,
        "magnet_temperature_c": temperature,
        "result_magnet_temperature_c": temperature,
        "remanence_reference_t": br_reference,
        "result_remanence_reference_t": br_reference,
        "remanence_temperature_coefficient_per_k": coefficient,
        "result_remanence_temperature_coefficient_per_k": coefficient,
        "temperature_adjusted_remanence_t": br_temperature,
        "result_temperature_adjusted_remanence_t": br_temperature,
        "recoil_relative_permeability": recoil_mu,
        "result_recoil_relative_permeability": recoil_mu,
        "operating_h_a_per_m": operating_h,
        "result_operating_h_a_per_m": operating_h,
        "operating_b_t": operating_b,
        "result_operating_b_t": operating_b,
        "knee_h_a_per_m": -7.0e5,
        "result_knee_h_a_per_m": -7.0e5,
        "irreversible_region": True,
        "result_irreversible_region": True,
        "remanence_loss_fraction": loss,
        "result_remanence_loss_fraction": loss,
        "airgap_flux_before_wb": 1.0e-2,
        "result_airgap_flux_before_wb": 1.0e-2,
        "airgap_flux_after_wb": 1.0e-2 * (1.0 - loss),
        "result_airgap_flux_after_wb": 1.0e-2 * (1.0 - loss),
        "torque_before_nm": 10.0,
        "result_torque_before_nm": 10.0,
        "torque_after_nm": 10.0 * (1.0 - loss),
        "result_torque_after_nm": 10.0 * (1.0 - loss),
        "mesh_owner": "mesh:pm-demag-271",
        "accepted_mesh_owner": "mesh:pm-demag-271",
        "demag_result_sha256": "2" * 64,
        "accepted_demag_result_sha256": "2" * 64,
    }
    return payload


_INDUCTION_V40 = "induction_cage_slip_rotorbar_current_loss_endring_airgap_torque_power_model_mesh_result_identity"


_FIELDWEAKENING = "ipm_fieldweakening_dq_flux_voltage_limit_current_angle_speed_torque_power_model_result_identity"


def _payload_v40():
    payload = _payload_v39()
    identity = payload["artifact_identity"]
    generation = "induction-cage-280"
    synchronous_speed = 157.07963267948966
    rotor_speed = 141.3716694115407
    slip = (synchronous_speed - rotor_speed) / synchronous_speed
    bar_count = 24
    bar_current = 80.0
    bar_resistance = 1.5e-4
    ring_current = 260.0
    ring_resistance = 3.0e-5
    bar_loss = bar_count * bar_current**2 * bar_resistance
    ring_loss = 2.0 * ring_current**2 * ring_resistance
    rotor_loss = bar_loss + ring_loss
    airgap_power = rotor_loss / slip
    torque = airgap_power / synchronous_speed
    mechanical_power = torque * rotor_speed
    identity[_INDUCTION_V40] = {
        "induction_generation": generation,
        **{key: generation for key in ("slip_generation", "rotorbar_generation", "endring_generation", "loss_generation", "airgap_generation", "torque_generation", "power_generation", "model_generation", "mesh_generation", "result_generation")},
        "synchronous_mechanical_speed_rad_s": synchronous_speed,
        "result_synchronous_mechanical_speed_rad_s": synchronous_speed,
        "rotor_mechanical_speed_rad_s": rotor_speed,
        "result_rotor_mechanical_speed_rad_s": rotor_speed,
        "slip": slip,
        "result_slip": slip,
        "rotor_bar_count": bar_count,
        "result_rotor_bar_count": bar_count,
        "rotor_bar_current_rms_a": bar_current,
        "result_rotor_bar_current_rms_a": bar_current,
        "rotor_bar_resistance_ohm": bar_resistance,
        "result_rotor_bar_resistance_ohm": bar_resistance,
        "endring_current_rms_a": ring_current,
        "result_endring_current_rms_a": ring_current,
        "endring_segment_resistance_ohm": ring_resistance,
        "result_endring_segment_resistance_ohm": ring_resistance,
        "rotor_bar_loss_w": bar_loss,
        "result_rotor_bar_loss_w": bar_loss,
        "endring_loss_w": ring_loss,
        "result_endring_loss_w": ring_loss,
        "rotor_copper_loss_w": rotor_loss,
        "result_rotor_copper_loss_w": rotor_loss,
        "airgap_power_w": airgap_power,
        "result_airgap_power_w": airgap_power,
        "electromagnetic_torque_nm": torque,
        "result_electromagnetic_torque_nm": torque,
        "mechanical_power_w": mechanical_power,
        "result_mechanical_power_w": mechanical_power,
        "model_owner": "motor:induction-cage-280",
        "accepted_model_owner": "motor:induction-cage-280",
        "mesh_owner": "mesh:induction-cage-280",
        "accepted_mesh_owner": "mesh:induction-cage-280",
        "induction_result_sha256": "5" * 64,
        "accepted_induction_result_sha256": "5" * 64,
    }

    generation = "ipm-fieldweakening-280"
    pole_pairs = 4
    resistance = 0.05
    ld = 5.0e-3
    lq = 9.0e-3
    magnet_flux = 0.1
    current_d = -8.0
    current_q = 12.0
    electrical_speed = 400.0
    flux_d = magnet_flux + ld * current_d
    flux_q = lq * current_q
    voltage_d = resistance * current_d - electrical_speed * flux_q
    voltage_q = resistance * current_q + electrical_speed * flux_d
    current_magnitude = math.hypot(current_d, current_q)
    voltage_magnitude = math.hypot(voltage_d, voltage_q)
    current_angle = math.degrees(math.atan2(-current_d, current_q))
    torque = 1.5 * pole_pairs * (flux_d * current_q - flux_q * current_d)
    mechanical_speed = electrical_speed / pole_pairs
    copper_loss = 1.5 * resistance * current_magnitude**2
    electrical_power = 1.5 * (voltage_d * current_d + voltage_q * current_q)
    mechanical_power = torque * mechanical_speed
    values = {
        "phase_resistance_ohm": resistance, "ld_h": ld, "lq_h": lq,
        "magnet_flux_wb": magnet_flux, "current_d_a": current_d,
        "current_q_a": current_q, "flux_d_wb": flux_d, "flux_q_wb": flux_q,
        "electrical_speed_rad_s": electrical_speed,
        "mechanical_speed_rad_s": mechanical_speed, "voltage_d_v": voltage_d,
        "voltage_q_v": voltage_q, "current_magnitude_a": current_magnitude,
        "current_limit_a": 20.0, "voltage_magnitude_v": voltage_magnitude,
        "voltage_limit_v": 60.0, "current_angle_deg": current_angle,
        "electromagnetic_torque_nm": torque, "copper_loss_w": copper_loss,
        "electrical_power_w": electrical_power, "mechanical_power_w": mechanical_power,
    }
    identity[_FIELDWEAKENING] = {
        "fieldweakening_generation": generation,
        **{key: generation for key in ("flux_generation", "voltage_generation", "current_generation", "angle_generation", "speed_generation", "torque_generation", "power_generation", "model_generation", "result_generation")},
        "pole_pairs": pole_pairs,
        "result_pole_pairs": pole_pairs,
        **values,
        **{f"result_{key}": value for key, value in values.items()},
        "model_owner": "motor:ipm-fieldweakening-280",
        "accepted_model_owner": "motor:ipm-fieldweakening-280",
        "fieldweakening_result_sha256": "6" * 64,
        "accepted_fieldweakening_result_sha256": "6" * 64,
    }
    return payload


_SRM = "srm_inductance_position_current_coenergy_torque_ripple_power_model_result_generation_identity"


_AXIAL = "axial_flux_sector_periodicity_skew_end_effect_flux_torque_loss_power_mesh_result_generation_identity"


def _derivative(values: list[float], coordinates: list[float]) -> list[float]:
    result = []
    for index in range(len(values)):
        if index == 0:
            left, right = 0, 1
        elif index == len(values) - 1:
            left, right = index - 1, index
        else:
            left, right = index - 1, index + 1
        result.append(
            (values[right] - values[left])
            / (coordinates[right] - coordinates[left])
        )
    return result


def _payload_v41():
    payload = _payload_v40()
    identity = payload["artifact_identity"]

    generation = "srm-map-724"
    positions = [0.0, 0.1, 0.2, 0.3, 0.4]
    currents = [5.0, 10.0]
    inductance = [
        [0.0105, 0.0125, 0.0150, 0.0170, 0.0185],
        [0.0100, 0.0120, 0.0145, 0.0165, 0.0180],
    ]
    selected_current = 10.0
    coenergy = [0.5 * item * selected_current**2 for item in inductance[1]]
    torque = _derivative(coenergy, positions)
    average_torque = sum(torque) / len(torque)
    phase_count = 3
    phase_resistance = 0.05
    copper_loss = phase_count * phase_resistance * selected_current**2
    mechanical_speed = 100.0
    values = {
        "rotor_position_rad": positions,
        "current_samples_a": currents,
        "inductance_h_by_current": inductance,
        "selected_current_a": selected_current,
        "coenergy_j": coenergy,
        "torque_nm": torque,
        "average_torque_nm": average_torque,
        "torque_ripple_nm": max(torque) - min(torque),
        "phase_count": phase_count,
        "phase_resistance_ohm": phase_resistance,
        "copper_loss_w": copper_loss,
        "mechanical_speed_rad_s": mechanical_speed,
        "mechanical_power_w": average_torque * mechanical_speed,
        "electrical_power_w": average_torque * mechanical_speed + copper_loss,
    }
    identity[_SRM] = {
        "srm_generation": generation,
        **{
            key: generation
            for key in (
                "position_generation", "current_generation",
                "inductance_generation", "coenergy_generation",
                "torque_generation", "ripple_generation", "loss_generation",
                "power_generation", "model_generation", "result_generation",
            )
        },
        **values,
        **{f"result_{key}": value for key, value in values.items()},
        "model_owner": "motor:srm-map-724",
        "accepted_model_owner": "motor:srm-map-724",
        "srm_result_sha256": "5" * 64,
        "accepted_srm_result_sha256": "5" * 64,
    }

    generation = "axial-flux-724"
    sector_count = 8
    pole_pairs = 4
    skew_angle = 2.0
    skew_argument = pole_pairs * math.radians(skew_angle) / 2.0
    skew_factor = math.sin(skew_argument) / skew_argument
    end_effect = 0.95
    uncorrected_flux = 0.5
    corrected_flux = uncorrected_flux * skew_factor * end_effect
    full_torque = 26.2
    current_q = full_torque / (1.5 * pole_pairs * corrected_flux)
    phase_resistance = 0.1
    copper_loss = 3.0 * phase_resistance * current_q**2
    iron_loss = 50.0
    mechanical_speed = 100.0
    values = {
        "sector_count": sector_count,
        "sector_angle_deg": 360.0 / sector_count,
        "periodicity": "periodic",
        "pole_pairs": pole_pairs,
        "skew_angle_deg": skew_angle,
        "skew_factor": skew_factor,
        "end_effect_factor": end_effect,
        "uncorrected_airgap_flux_wb": uncorrected_flux,
        "corrected_airgap_flux_wb": corrected_flux,
        "current_q_a": current_q,
        "sector_torque_nm": full_torque / sector_count,
        "full_torque_nm": full_torque,
        "phase_resistance_ohm": phase_resistance,
        "copper_loss_w": copper_loss,
        "iron_loss_w": iron_loss,
        "mechanical_speed_rad_s": mechanical_speed,
        "mechanical_power_w": full_torque * mechanical_speed,
        "electrical_power_w": full_torque * mechanical_speed + copper_loss + iron_loss,
    }
    identity[_AXIAL] = {
        "axial_flux_generation": generation,
        **{
            key: generation
            for key in (
                "sector_generation", "periodicity_generation", "skew_generation",
                "end_effect_generation", "flux_generation", "torque_generation",
                "loss_generation", "power_generation", "mesh_generation",
                "result_generation",
            )
        },
        **values,
        **{f"result_{key}": value for key, value in values.items()},
        "mesh_owner": "mesh:axial-flux-724",
        "accepted_mesh_owner": "mesh:axial-flux-724",
        "axial_flux_result_sha256": "6" * 64,
        "accepted_axial_flux_result_sha256": "6" * 64,
    }
    return payload


_IPM = (
    "ipm_dq_current_flux_torque_voltage_powerfactor_mtpv_energy_mesh_"
    "result_generation_identity"
)


_INDUCTION = (
    "inductionmotor_slip_rotorfrequency_copperloss_torque_airgappower_"
    "mechanicalpower_efficiency_result_generation_identity"
)


def _payload_v42():
    payload = _payload_v41()
    identity = payload["artifact_identity"]
    generation = "ipm-dq-842"
    pole_pairs = 4
    current_d = -20.0
    current_q = 50.0
    pm_flux = 0.1
    inductance_d = 1.0e-3
    inductance_q = 2.0e-3
    flux_d = pm_flux + inductance_d * current_d
    flux_q = inductance_q * current_q
    torque = 1.5 * pole_pairs * (flux_d * current_q - flux_q * current_d)
    resistance = 0.05
    electrical_speed = 1000.0
    voltage_d = resistance * current_d - electrical_speed * flux_q
    voltage_q = resistance * current_q + electrical_speed * flux_d
    voltage_magnitude = math.hypot(voltage_d, voltage_q)
    voltage_limit = 150.0
    active_power = 1.5 * (voltage_d * current_d + voltage_q * current_q)
    apparent_power = 1.5 * voltage_magnitude * math.hypot(current_d, current_q)
    power_factor = active_power / apparent_power
    field_energy = 0.5 * (inductance_d * current_d**2 + inductance_q * current_q**2)
    coenergy = pm_flux * current_d + field_energy
    identity[_IPM] = {
        "dq_generation": generation,
        **{
            key: generation
            for key in (
                "current_generation", "flux_generation", "torque_generation",
                "voltage_generation", "powerfactor_generation", "mtpv_generation",
                "energy_generation", "mesh_generation", "result_generation",
            )
        },
        "pole_pairs": pole_pairs,
        "result_pole_pairs": pole_pairs,
        "current_d_a": current_d,
        "result_current_d_a": current_d,
        "current_q_a": current_q,
        "result_current_q_a": current_q,
        "pm_flux_linkage_wb": pm_flux,
        "result_pm_flux_linkage_wb": pm_flux,
        "inductance_d_h": inductance_d,
        "result_inductance_d_h": inductance_d,
        "inductance_q_h": inductance_q,
        "result_inductance_q_h": inductance_q,
        "flux_d_wb": flux_d,
        "result_flux_d_wb": flux_d,
        "flux_q_wb": flux_q,
        "result_flux_q_wb": flux_q,
        "torque_nm": torque,
        "result_torque_nm": torque,
        "phase_resistance_ohm": resistance,
        "result_phase_resistance_ohm": resistance,
        "electrical_speed_rad_s": electrical_speed,
        "result_electrical_speed_rad_s": electrical_speed,
        "voltage_d_v": voltage_d,
        "result_voltage_d_v": voltage_d,
        "voltage_q_v": voltage_q,
        "result_voltage_q_v": voltage_q,
        "voltage_magnitude_v": voltage_magnitude,
        "result_voltage_magnitude_v": voltage_magnitude,
        "voltage_limit_v": voltage_limit,
        "result_voltage_limit_v": voltage_limit,
        "active_power_w": active_power,
        "result_active_power_w": active_power,
        "apparent_power_va": apparent_power,
        "result_apparent_power_va": apparent_power,
        "power_factor": power_factor,
        "result_power_factor": power_factor,
        "mtpv_branch": "negative_id_high_speed",
        "result_mtpv_branch": "negative_id_high_speed",
        "mtpv_voltage_margin_v": voltage_limit - voltage_magnitude,
        "result_mtpv_voltage_margin_v": voltage_limit - voltage_magnitude,
        "field_energy_j": field_energy,
        "result_field_energy_j": field_energy,
        "coenergy_j": coenergy,
        "result_coenergy_j": coenergy,
        "mesh_owner": "mesh:ipm-dq-842",
        "accepted_mesh_owner": "mesh:ipm-dq-842",
        "ipm_result_sha256": "1" * 64,
        "accepted_ipm_result_sha256": "1" * 64,
    }

    generation = "induction-power-842"
    pole_pairs = 2
    supply_frequency = 50.0
    synchronous_speed = 2.0 * math.pi * supply_frequency / pole_pairs
    rotor_speed = 150.0
    slip = (synchronous_speed - rotor_speed) / synchronous_speed
    rotor_frequency = slip * supply_frequency
    torque = 20.0
    airgap_power = torque * synchronous_speed
    rotor_loss = slip * airgap_power
    converted = (1.0 - slip) * airgap_power
    mechanical_loss = 30.0
    mechanical_power = converted - mechanical_loss
    input_power = airgap_power + 100.0 + 50.0
    efficiency = mechanical_power / input_power
    identity[_INDUCTION] = {
        "induction_generation": generation,
        **{
            key: generation
            for key in (
                "slip_generation", "frequency_generation", "loss_generation",
                "torque_generation", "airgap_power_generation",
                "mechanical_power_generation", "efficiency_generation", "result_generation",
            )
        },
        "pole_pairs": pole_pairs,
        "result_pole_pairs": pole_pairs,
        "supply_frequency_hz": supply_frequency,
        "result_supply_frequency_hz": supply_frequency,
        "synchronous_speed_rad_s": synchronous_speed,
        "result_synchronous_speed_rad_s": synchronous_speed,
        "rotor_speed_rad_s": rotor_speed,
        "result_rotor_speed_rad_s": rotor_speed,
        "slip": slip,
        "result_slip": slip,
        "rotor_electrical_frequency_hz": rotor_frequency,
        "result_rotor_electrical_frequency_hz": rotor_frequency,
        "torque_nm": torque,
        "result_torque_nm": torque,
        "airgap_power_w": airgap_power,
        "result_airgap_power_w": airgap_power,
        "stator_copper_loss_w": 100.0,
        "result_stator_copper_loss_w": 100.0,
        "rotor_copper_loss_w": rotor_loss,
        "result_rotor_copper_loss_w": rotor_loss,
        "core_loss_w": 50.0,
        "result_core_loss_w": 50.0,
        "mechanical_loss_w": mechanical_loss,
        "result_mechanical_loss_w": mechanical_loss,
        "converted_power_w": converted,
        "result_converted_power_w": converted,
        "mechanical_power_w": mechanical_power,
        "result_mechanical_power_w": mechanical_power,
        "input_power_w": input_power,
        "result_input_power_w": input_power,
        "efficiency": efficiency,
        "result_efficiency": efficiency,
        "motor_result_sha256": "2" * 64,
        "accepted_motor_result_sha256": "2" * 64,
    }
    return payload


_PMSM_FORCE = "pmsm_radialforce_space_time_harmonics_torque_nvh_power_energy_generation_identity"


_WOUNDFIELD = "woundfield_synchronous_excitation_flux_torque_copperloss_efficiency_energy_generation_identity"


def _payload_v43():
    payload = _payload_v42()
    generation = "pmsm-force-nvh-843"
    payload["artifact_identity"][_PMSM_FORCE] = {
        "pmsm_generation": generation,
        **{key: generation for key in (
            "radial_force_generation", "space_harmonic_generation", "time_harmonic_generation",
            "torque_generation", "nvh_generation", "power_generation", "energy_generation",
            "mesh_generation", "result_generation")},
        "radial_force_space_orders": [6, 12], "result_radial_force_space_orders": [6, 12],
        "radial_force_time_orders": [6, 12], "result_radial_force_time_orders": [6, 12],
        "force_phase_deg": 30.0, "result_force_phase_deg": 30.0,
        "torque_ripple_rms_nm": 0.8, "result_torque_ripple_rms_nm": 0.8,
        "modal_excitation_n": 0.3, "result_modal_excitation_n": 0.3,
        "electromagnetic_power_w": 1000.0, "result_electromagnetic_power_w": 1000.0,
        "mechanical_power_w": 950.0, "result_mechanical_power_w": 950.0,
        "energy_closure_residual": 1.0e-6, "result_energy_closure_residual": 1.0e-6,
        "mesh_owner": "mesh:pmsm-force-nvh-843", "result_mesh_owner": "mesh:pmsm-force-nvh-843",
        "pmsm_result_sha256": "5" * 64, "accepted_pmsm_result_sha256": "5" * 64,
    }
    generation = "woundfield-sync-843"
    payload["artifact_identity"][_WOUNDFIELD] = {
        "machine_generation": generation,
        **{key: generation for key in (
            "field_excitation_generation", "flux_generation", "torque_generation",
            "copper_loss_generation", "efficiency_generation", "energy_generation",
            "mesh_generation", "result_generation")},
        "field_current_a": 5.0, "result_field_current_a": 5.0,
        "flux_linkage_wb": 0.2, "result_flux_linkage_wb": 0.2,
        "torque_nm": 20.0, "result_torque_nm": 20.0,
        "stator_copper_loss_w": 100.0, "result_stator_copper_loss_w": 100.0,
        "field_copper_loss_w": 20.0, "result_field_copper_loss_w": 20.0,
        "input_power_w": 2000.0, "result_input_power_w": 2000.0,
        "mechanical_power_w": 1600.0, "result_mechanical_power_w": 1600.0,
        "efficiency": 0.8, "result_efficiency": 0.8,
        "mesh_owner": "mesh:woundfield-sync-843", "result_mesh_owner": "mesh:woundfield-sync-843",
        "machine_result_sha256": "6" * 64, "accepted_machine_result_sha256": "6" * 64,
    }
    return payload


def _identity_v45():
    generation = "test-845"
    return {
        "v45_public_ipmsm_torque_ripple_radial_force_modal_power_efficiency_energy_mesh_owner_mismatch": {
            "generation": generation, **{key: generation for key in ("torque_generation", "radial_force_generation", "modal_generation", "power_generation", "efficiency_generation", "energy_generation", "mesh_generation", "result_generation")},
            "torque_ripple_rms_nm": 0.8, "result_torque_ripple_rms_nm": 0.8, "radial_force_space_orders": [6, 12], "result_radial_force_space_orders": [6, 12], "modal_excitation_n": 0.3, "result_modal_excitation_n": 0.3, "electromagnetic_power_w": 1000.0, "result_electromagnetic_power_w": 1000.0, "mechanical_power_w": 950.0, "result_mechanical_power_w": 950.0, "efficiency": 0.95, "result_efficiency": 0.95, "energy_closure_residual": 1e-10, "result_energy_closure_residual": 1e-10, "mesh_owner": "mesh:test", "result_mesh_owner": "mesh:test", "result_sha256": "a" * 64, "accepted_result_sha256": "a" * 64,
        },
        "v45_public_induction_machine_slip_rotor_loss_torque_heat_power_energy_result_mismatch": {
            "generation": generation, **{key: generation for key in ("slip_generation", "rotor_loss_generation", "torque_generation", "heat_generation", "power_generation", "energy_generation", "result_generation")},
            "slip": 0.05, "result_slip": 0.05, "rotor_copper_loss_w": 50.0, "result_rotor_copper_loss_w": 50.0, "torque_nm": 20.0, "result_torque_nm": 20.0, "heat_loss_w": 50.0, "result_heat_loss_w": 50.0, "electrical_power_w": 1000.0, "result_electrical_power_w": 1000.0, "mechanical_power_w": 900.0, "result_mechanical_power_w": 900.0, "energy_closure_residual": 1e-10, "result_energy_closure_residual": 1e-10, "result_sha256": "b" * 64, "accepted_result_sha256": "b" * 64,
        },
    }


def _identity_v47() -> dict[str, object]:
    dq_generation = "dq-v47"
    window_generation = "window-v47"
    return {
        DQ_V47: {
            "generation": dq_generation,
            **{
                key: dq_generation
                for key in (
                    "dq_generation",
                    "phase_order_generation",
                    "electrical_angle_generation",
                    "pole_pair_generation",
                    "transform_generation",
                    "result_generation",
                )
            },
            "phase_order": ["A", "B", "C"],
            "result_phase_order": ["A", "B", "C"],
            "electrical_angle_origin_deg": 0.0,
            "result_electrical_angle_origin_deg": 0.0,
            "pole_pairs": 4,
            "result_pole_pairs": 4,
            "transform_identity": "power_invariant_park",
            "result_transform_identity": "power_invariant_park",
            "angle_direction": "electrical_ccw",
            "result_angle_direction": "electrical_ccw",
            "result_sha256": "1" * 64,
            "accepted_result_sha256": "1" * 64,
        },
        WINDOW: {
            "generation": window_generation,
            **{
                key: window_generation
                for key in (
                    "integration_window_generation",
                    "parameter_row_generation",
                    "torque_generation",
                    "loss_generation",
                    "result_generation",
                )
            },
            "integration_window_s": [0.01, 0.02],
            "result_integration_window_s": [0.01, 0.02],
            "parameter_row_key": "speed=3000rpm,current=100A",
            "result_parameter_row_key": "speed=3000rpm,current=100A",
            "torque_mean_nm": 12.5,
            "result_torque_mean_nm": 12.5,
            "loss_total_w": 345.0,
            "result_loss_total_w": 345.0,
            "study_owner": "study:motor-v47",
            "result_study_owner": "study:motor-v47",
            "result_sha256": "2" * 64,
            "accepted_result_sha256": "2" * 64,
        },
    }


def _identity_v48() -> dict[str, object]:
    skew_generation = "skew-aggregation-v48-901"
    pwm_generation = "pwm-timeline-v48-901"
    times = [0.0, 5.0e-6, 10.0e-6, 15.0e-6]
    states = [[1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 1, 1]]
    currents = [[10.0, -5.0, -5.0], [9.8, -4.7, -5.1], [9.5, -4.5, -5.0], [9.2, -4.2, -5.0]]
    voltages = [[300.0, -150.0, -150.0], [150.0, 150.0, -300.0], [-300.0, 150.0, 150.0], [-150.0, 300.0, -150.0]]
    losses = [2.0, 2.1, 1.9, 2.2]
    return {
        SKEW_V48: {
            "generation": skew_generation,
            "slice_generation": skew_generation,
            "angle_generation": skew_generation,
            "harmonic_generation": skew_generation,
            "phase_generation": skew_generation,
            "result_generation": skew_generation,
            "slice_ids": ["slice:0", "slice:1", "slice:2"],
            "result_slice_ids": ["slice:0", "slice:1", "slice:2"],
            "slice_weights": [0.25, 0.50, 0.25],
            "result_slice_weights": [0.25, 0.50, 0.25],
            "rotor_angles_deg": [-5.0, 0.0, 5.0],
            "result_rotor_angles_deg": [-5.0, 0.0, 5.0],
            "torque_harmonic_phasors_nm": [[1.0, 0.0], [0.8, 0.1], [0.6, -0.1]],
            "result_torque_harmonic_phasors_nm": [[1.0, 0.0], [0.8, 0.1], [0.6, -0.1]],
            "phase_origins_deg": [0.0, 0.0, 0.0],
            "result_phase_origins_deg": [0.0, 0.0, 0.0],
            "machine_state_owner": "machine-state:skew-v48-901",
            "result_machine_state_owner": "machine-state:skew-v48-901",
            "result_sha256": "5" * 64,
            "accepted_result_sha256": "5" * 64,
        },
        PWM: {
            "generation": pwm_generation,
            "carrier_generation": pwm_generation,
            "control_generation": pwm_generation,
            "switch_generation": pwm_generation,
            "electrical_generation": pwm_generation,
            "loss_generation": pwm_generation,
            "result_generation": pwm_generation,
            "sample_times_s": times,
            "result_sample_times_s": times,
            "carrier_frequency_hz": 20000.0,
            "result_carrier_frequency_hz": 20000.0,
            "control_sample_divider": 2,
            "result_control_sample_divider": 2,
            "switch_states": states,
            "result_switch_states": states,
            "phase_current_a": currents,
            "result_phase_current_a": currents,
            "phase_voltage_v": voltages,
            "result_phase_voltage_v": voltages,
            "loss_w": losses,
            "result_loss_w": losses,
            "timeline_owner": "timeline:pwm-v48-901",
            "result_timeline_owner": "timeline:pwm-v48-901",
            "result_sha256": "6" * 64,
            "accepted_result_sha256": "6" * 64,
        },
    }


def _identity_v49() -> dict[str, object]:
    demag_generation = "demag-operating-point-v49-901"
    iron_generation = "iron-loss-spectrum-v49-901"
    state = {"magnet:north": 0.998, "magnet:south": 0.997}
    harmonics = [
        {"order": 1, "frequency_hz": 400.0, "hysteresis_w": 8.0, "eddy_w": 3.0, "excess_w": 1.0},
        {"order": 3, "frequency_hz": 1200.0, "hysteresis_w": 1.5, "eddy_w": 1.2, "excess_w": 0.3},
    ]
    coefficients = {"hysteresis": 1.0, "eddy": 0.020, "excess": 0.004}
    return {
        DEMAG_V49: {
            "generation": demag_generation,
            "temperature_generation": demag_generation,
            "current_generation": demag_generation,
            "angle_generation": demag_generation,
            "state_generation": demag_generation,
            "result_generation": demag_generation,
            "temperature_c": 120.0,
            "result_temperature_c": 120.0,
            "phase_current_a": [80.0, -40.0, -40.0],
            "result_phase_current_a": [80.0, -40.0, -40.0],
            "rotor_angle_electrical_deg": 90.0,
            "result_rotor_angle_electrical_deg": 90.0,
            "irreversible_magnet_state": state,
            "result_irreversible_magnet_state": state,
            "operating_point_owner": "operating-point:demag-v49-901",
            "result_operating_point_owner": "operating-point:demag-v49-901",
            "result_sha256": "a" * 64,
            "accepted_result_sha256": "a" * 64,
        },
        IRON_V49: {
            "generation": iron_generation,
            "time_generation": iron_generation,
            "frequency_generation": iron_generation,
            "harmonic_generation": iron_generation,
            "coefficient_generation": iron_generation,
            "result_generation": iron_generation,
            "time_window_s": [0.0, 0.01],
            "result_time_window_s": [0.0, 0.01],
            "harmonic_rows": harmonics,
            "result_harmonic_rows": harmonics,
            "loss_coefficients": coefficients,
            "result_loss_coefficients": coefficients,
            "loss_owner": "loss-table:iron-v49-901",
            "result_loss_owner": "loss-table:iron-v49-901",
            "result_sha256": "b" * 64,
            "accepted_result_sha256": "b" * 64,
        },
    }


def _identity_v50() -> dict[str, object]:
    dq_generation = "dq-operating-point-v50-901"
    thermal_generation = "electrothermal-v50-901"
    currents = {"id_a": -35.0, "iq_a": 82.0, "phase_order": "uvw"}
    saliency = {"ld_h": 0.0018, "lq_h": 0.0032, "psi_pm_wb": 0.092}
    loss = {"stator_iron_w": 42.0, "rotor_iron_w": 11.0, "copper_w": 68.0}
    convection = [{"boundary": "housing", "h_w_m2k": 18.0, "ambient_c": 25.0}]
    temperatures = {"winding_c": 96.0, "magnet_c": 78.0, "housing_c": 54.0}
    materials = {"copper": "cu:v3", "steel": "steel:v7", "magnet": "pm:v5"}
    return {
        DQ_V50: {
            "generation": dq_generation, "angle_generation": dq_generation, "current_generation": dq_generation,
            "saliency_generation": dq_generation, "operating_point_generation": dq_generation, "result_generation": dq_generation,
            "park_angle_electrical_deg": 37.5, "result_park_angle_electrical_deg": 37.5,
            "dq_currents": currents, "result_dq_currents": currents,
            "saliency_parameters": saliency, "result_saliency_parameters": saliency,
            "operating_point_id": "operating-point:dq-v50-901", "result_operating_point_id": "operating-point:dq-v50-901",
            "result_owner": "dq-result:motor-v50-901", "accepted_result_owner": "dq-result:motor-v50-901",
            "result_sha256": "1" * 64, "accepted_result_sha256": "1" * 64,
        },
        THERMAL: {
            "generation": thermal_generation, "loss_generation": thermal_generation, "boundary_generation": thermal_generation,
            "temperature_generation": thermal_generation, "material_generation": thermal_generation, "result_generation": thermal_generation,
            "loss_map": loss, "replayed_loss_map": loss,
            "convection_boundaries": convection, "replayed_convection_boundaries": convection,
            "temperature_map": temperatures, "replayed_temperature_map": temperatures,
            "material_revisions": materials, "replayed_material_revisions": materials,
            "thermal_owner": "thermal-result:motor-v50-901", "replayed_thermal_owner": "thermal-result:motor-v50-901",
            "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64,
        },
    }


def _identity_v51() -> dict[str, object]:
    generation = "motor-public-v51"
    result = {"result_sha256": "f" * 64, "accepted_result_sha256": "f" * 64}
    mechanical = [index * 11.25 for index in range(9)]
    electrical = [value * 4.0 for value in mechanical]
    torque = [10.0, 10.4, 10.0, 9.6, 10.0, 10.4, 10.0, 9.6, 10.0]
    harmonics = [{"order": 0, "amplitude_nm": 10.0}, {"order": 2, "amplitude_nm": 0.4}]
    return {
        TORQUE_V51: {
            "generation": generation, "rotor_generation": generation, "angle_generation": generation,
            "period_generation": generation, "window_generation": generation, "fft_generation": generation,
            "owner_generation": generation, "result_generation": generation, "pole_pairs": 4, "result_pole_pairs": 4,
            "rotor_angles_mechanical_deg": mechanical, "result_rotor_angles_mechanical_deg": mechanical,
            "rotor_angles_electrical_deg": electrical, "result_rotor_angles_electrical_deg": electrical,
            "mechanical_period_deg": 90.0, "result_mechanical_period_deg": 90.0, "electrical_period_deg": 360.0,
            "result_electrical_period_deg": 360.0, "sample_window_mechanical_deg": [0.0, 90.0],
            "result_sample_window_mechanical_deg": [0.0, 90.0], "torque_samples_nm": torque,
            "result_torque_samples_nm": torque, "fft_harmonics": harmonics, "result_fft_harmonics": harmonics,
            "torque_owner": "torque:motor-v51", "result_torque_owner": "torque:motor-v51", **result,
        },
        WINDING: {
            "generation": generation, "temperature_generation": generation, "resistance_generation": generation,
            "length_generation": generation, "fill_generation": generation, "loss_generation": generation,
            "owner_generation": generation, "result_generation": generation, "reference_temperature_c": 20.0,
            "winding_temperature_c": 120.0, "result_winding_temperature_c": 120.0,
            "copper_temperature_coefficient_per_k": 0.00393, "resistance_reference_ohm": 0.08,
            "resistance_at_temperature_ohm": 0.11144, "result_resistance_at_temperature_ohm": 0.11144,
            "active_length_m": 0.42, "end_turn_length_m": 0.18, "result_end_turn_length_m": 0.18,
            "slot_fill_factor": 0.48, "result_slot_fill_factor": 0.48, "current_rms_a": 20.0,
            "copper_loss_w": 44.576, "result_copper_loss_w": 44.576, "winding_owner": "winding:phase-u-v51",
            "result_winding_owner": "winding:phase-u-v51", **result,
        },
    }


def _generations(generation: str, fields: tuple[str, ...]) -> dict[str, str]:
    return {"generation": generation, **{field: generation for field in fields}}


def _identity_v52() -> dict[str, object]:
    result = {"result_sha256": "f" * 64, "accepted_result_sha256": "f" * 64}
    harmonics = [{"order": 1, "amplitude_t": 1.2}, {"order": 3, "amplitude_t": 0.08}]
    coefficients = {"hysteresis": 1.7, "classical_eddy": 0.012, "excess": 0.08}
    angles = [index * 0.5 for index in range(13)]
    torque = [0.12, 0.10, 0.04, -0.04, -0.10, -0.12, -0.10, -0.04, 0.04, 0.10, 0.12, 0.10, 0.12]
    return {
        IRON_LOSS: {
            **_generations("iron-loss-v52", ("window_generation", "harmonic_generation", "frequency_generation", "coefficient_generation", "owner_generation", "result_generation")),
            "fft_window": "hann_periodic", "result_fft_window": "hann_periodic",
            "sample_count": 1024, "result_sample_count": 1024,
            "harmonics": harmonics, "result_harmonics": harmonics,
            "rotation_frequency_hz": 50.0, "result_rotation_frequency_hz": 50.0,
            "pole_pairs": 4, "result_pole_pairs": 4,
            "electrical_frequency_hz": 200.0, "result_electrical_frequency_hz": 200.0,
            "loss_coefficients": coefficients, "result_loss_coefficients": coefficients,
            "waveform_owner": "waveform:iron-loss-v52", "result_waveform_owner": "waveform:iron-loss-v52",
            **result,
        },
        COGGING: {
            **_generations("cogging-v52", ("period_generation", "sampling_generation", "phase_generation", "owner_generation", "result_generation")),
            "slot_count": 12, "result_slot_count": 12, "pole_count": 10, "result_pole_count": 10,
            "cogging_period_mechanical_deg": 6.0, "result_cogging_period_mechanical_deg": 6.0,
            "sample_angles_mechanical_deg": angles, "result_sample_angles_mechanical_deg": angles,
            "torque_samples_nm": torque, "result_torque_samples_nm": torque,
            "phase_alignment": "slot_center_to_pole_center", "result_phase_alignment": "slot_center_to_pole_center",
            "torque_owner": "torque:cogging-v52", "result_torque_owner": "torque:cogging-v52", **result,
        },
    }


def _identity_v53():
    weights = [0.25, 0.5, 0.25]; angles = [-5.0, 0.0, 5.0]
    harmonics = [{"order": 1, "amplitude_nm": 12.0, "phase_deg": 0.0}, {"order": 6, "amplitude_nm": 0.35, "phase_deg": 25.0}]
    skew = {**_generations("skew-v53", ("slice_generation", "angle_generation", "harmonic_generation", "owner_generation", "result_generation")), "slice_weights": weights, "result_slice_weights": weights, "skew_angles_mechanical_deg": angles, "result_skew_angles_mechanical_deg": angles, "harmonic_torque": harmonics, "result_harmonic_torque": harmonics, "rotor_owner": "rotor:v53", "result_rotor_owner": "rotor:v53", "result_sha256": "1" * 64, "accepted_result_sha256": "1" * 64}
    operating = {"b_t": 0.62, "h_a_per_m": -420000.0}; recoil = {"relative_permeability": 1.05, "coercivity_a_per_m": 900000.0}
    demag = {**_generations("demag-v53", ("operating_generation", "temperature_generation", "recoil_generation", "irreversible_generation", "owner_generation", "result_generation")), "operating_point": operating, "result_operating_point": operating, "temperature_c": 140.0, "result_temperature_c": 140.0, "recoil_line": recoil, "result_recoil_line": recoil, "irreversible_demag_fraction": 0.03, "result_irreversible_demag_fraction": 0.03, "irreversible_state": "partially_demagnetized", "result_irreversible_state": "partially_demagnetized", "magnet_owner": "magnet:v53", "result_magnet_owner": "magnet:v53", "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64}
    return {SKEW_V53: skew, DEMAG_V53: demag}


def _payload_v54():
    mechanical = [0.0, 5.0, 10.0]; electrical = [0.0, 20.0, 40.0]; harmonics = [{"order": 6, "amplitude_nm": 0.35, "phase_electrical_deg": 25.0}]
    torque = {**_generations("torque-v54", ("harmonic_generation", "mechanical_generation", "electrical_generation", "polepair_generation", "owner_generation", "result_generation")), "pole_pairs": 4, "result_pole_pairs": 4, "mechanical_angles_deg": mechanical, "result_mechanical_angles_deg": mechanical, "electrical_angles_deg": electrical, "result_electrical_angles_deg": electrical, "torque_harmonics": harmonics, "result_torque_harmonics": harmonics, "result_owner": "result:v54", "accepted_result_owner": "result:v54", "result_sha256": "1" * 64, "accepted_result_sha256": "1" * 64}
    knee = {"b_t": 0.24, "h_a_per_m": -640000.0, "criterion": "operating_point_below_knee"}; current = [120.0, -60.0, -60.0]
    demag = {**_generations("demag-v54", ("knee_generation", "temperature_generation", "current_generation", "recovery_generation", "owner_generation", "result_generation")), "knee_criterion": knee, "result_knee_criterion": knee, "temperature_c": 160.0, "result_temperature_c": 160.0, "current_vector_abc_a": current, "result_current_vector_abc_a": current, "irreversible_demag_fraction": 0.05, "result_irreversible_demag_fraction": 0.05, "post_recovery_remanence_fraction": 0.95, "result_post_recovery_remanence_fraction": 0.95, "recovery_state": "partially_demagnetized", "result_recovery_state": "partially_demagnetized", "magnet_owner": "magnet:v54", "result_magnet_owner": "magnet:v54", "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64}
    return {TORQUE_V54: torque, DEMAG_V54: demag}


def _payload_v55():
    from math import atan2, degrees

    generation = lambda name, fields: {"generation": name, **{field: name for field in fields}}
    current = {"d": -20.0, "q": 80.0}; inductance = {"d": 1.0e-3, "q": 1.5e-3}; pm_flux = 0.08
    flux = {"d": pm_flux + inductance["d"] * current["d"], "q": inductance["q"] * current["q"]}; pole_pairs = 4
    angle = degrees(atan2(current["q"], current["d"])) % 360.0
    saliency = 1.5 * pole_pairs * (inductance["d"] - inductance["q"]) * current["d"] * current["q"]
    torque = 1.5 * pole_pairs * (flux["d"] * current["q"] - flux["q"] * current["d"])
    dq = {**generation("dq-v55", ("flux_generation", "current_generation", "angle_generation", "saliency_generation", "torque_generation", "owner_generation", "result_generation")), "pole_pairs": pole_pairs, "result_pole_pairs": pole_pairs, "current_dq_a": current, "result_current_dq_a": current, "inductance_dq_h": inductance, "result_inductance_dq_h": inductance, "pm_flux_linkage_wb": pm_flux, "result_pm_flux_linkage_wb": pm_flux, "flux_linkage_dq_wb": flux, "result_flux_linkage_dq_wb": flux, "current_electrical_angle_deg": angle, "result_current_electrical_angle_deg": angle, "saliency_torque_nm": saliency, "result_saliency_torque_nm": saliency, "torque_nm": torque, "result_torque_nm": torque, "result_owner": "result:dq-v55", "accepted_result_owner": "result:dq-v55", "result_sha256": "5" * 64, "accepted_result_sha256": "5" * 64}
    components = {"hysteresis_w": 12.0, "eddy_w": 8.0, "excess_w": 2.0}; waveform = [0.0, 1.4, 0.0, -1.4, 0.0]
    iron = {**generation("iron-v55", ("component_generation", "frequency_generation", "waveform_generation", "material_generation", "owner_generation", "result_generation")), "loss_components": components, "result_loss_components": components, "total_iron_loss_w": 22.0, "result_total_iron_loss_w": 22.0, "frequency_hz": 400.0, "result_frequency_hz": 400.0, "flux_density_waveform_t": waveform, "result_flux_density_waveform_t": waveform, "peak_flux_density_t": 1.4, "result_peak_flux_density_t": 1.4, "material_revision": "steel-v55-r7", "result_material_revision": "steel-v55-r7", "material_owner": "material:steel-v55", "result_material_owner": "material:steel-v55", "result_sha256": "6" * 64, "accepted_result_sha256": "6" * 64}
    return {DQ_V55: dq, IRON_V55: iron}


def _identity_v56() -> dict[str, object]:
    generation = "motor-public-v56-test"; generations = lambda fields: {field: generation for field in fields}
    speed = 3000.0; torque = 20.0; output = torque * speed * 2.0 * math.pi / 60.0
    losses = {"copper_w": 240.0, "iron_w": 110.0, "mechanical_w": 50.0}; input_power = output + sum(losses.values())
    frequency = 50.0; poles = 4; synchronous = 120.0 * frequency / poles; rotor = 1440.0; slip = (synchronous - rotor) / synchronous
    return {
        MAP: {"generation": generation, **generations(("speed_generation", "torque_generation", "input_generation", "output_generation", "loss_generation", "efficiency_generation", "owner_generation", "result_generation")), "speed_rpm": speed, "result_speed_rpm": speed, "torque_nm": torque, "result_torque_nm": torque, "input_power_w": input_power, "result_input_power_w": input_power, "output_power_w": output, "result_output_power_w": output, "loss_components_w": losses, "result_loss_components_w": losses, "efficiency": output / input_power, "result_efficiency": output / input_power, "result_owner": "result:map-v56", "accepted_result_owner": "result:map-v56", "result_sha256": "1" * 64, "accepted_result_sha256": "1" * 64},
        INDUCTION: {"generation": generation, **generations(("frequency_generation", "pole_generation", "speed_generation", "slip_generation", "rotorfrequency_generation", "torque_generation", "owner_generation", "result_generation")), "supply_frequency_hz": frequency, "result_supply_frequency_hz": frequency, "pole_count": poles, "result_pole_count": poles, "synchronous_speed_rpm": synchronous, "result_synchronous_speed_rpm": synchronous, "rotor_speed_rpm": rotor, "result_rotor_speed_rpm": rotor, "slip": slip, "result_slip": slip, "rotor_frequency_hz": slip * frequency, "result_rotor_frequency_hz": slip * frequency, "torque_nm": 48.0, "result_torque_nm": 48.0, "torque_state": "motoring", "result_torque_state": "motoring", "result_owner": "result:induction-v56", "accepted_result_owner": "result:induction-v56", "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64},
    }
