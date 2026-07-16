from __future__ import annotations

import copy
import json
import math

from radia_mcp.radia_ngsolve.pwm_controlled_motor_loss_gate import (
    pwm_controlled_motor_loss_gate,
)
from radia_mcp.radia_ngsolve.server import pwm_controlled_motor_loss_gate as mcp_gate


def _payload() -> dict:
    count = 101
    time_s = [index * 1.0e-4 for index in range(count)]
    phase_currents = []
    for index in range(count):
        angle = 2.0 * math.pi * index / (count - 1)
        phase_currents.append(
            [
                10.0 * math.sin(angle),
                10.0 * math.sin(angle - 2.0 * math.pi / 3.0),
                10.0 * math.sin(angle + 2.0 * math.pi / 3.0),
            ]
        )
    power_components = [[10.0 + index, 2.0 - 0.5 * index] for index in range(count)]
    reported_power = [sum(row) for row in power_components]
    eddy_bins = [[3.0, 1.0], [2.0, 0.5], [1.0, 0.5], [0.0, 0.0], [0.0, 0.0]]
    hysteresis = [[2.0, 0.5]] + [[0.0, 0.0] for _ in range(4)]
    iron = [[5.0, 1.5]] + [[0.0, 0.0] for _ in range(4)]
    return {
        "time_series": {
            "time_s": time_s,
            "angle_deg": [3600.0 * value for value in time_s],
            "speed_rpm": [600.0] * count,
            "torque_nm": [20.0 + math.sin(index / 5.0) for index in range(count)],
            "phase_currents_a": phase_currents,
            "control_command_a": [[-10.0, 20.0] for _ in range(count)],
            "control_feedback_a": [[-9.8, 19.6] for _ in range(count)],
            "power_components_w": power_components,
            "reported_total_power_w": reported_power,
        },
        "loss_spectrum": {
            "frequency_hz": [0.0, 100.0, 200.0, 300.0, 400.0],
            "eddy_components_w": eddy_bins,
            "hysteresis_components_w": hysteresis,
            "iron_components_w": iron,
            "reported_eddy_total_w": [sum(row) for row in eddy_bins],
            "reported_hysteresis_total_w": [sum(row) for row in hysteresis],
            "reported_iron_total_w": [sum(row) for row in iron],
        },
        "ratio_diagnostic_policy": "exclude_ratio_when_denominator_current_is_below_floor",
    }


def _with_artifact_identity(payload: dict) -> dict:
    payload["artifact_identity"] = {
        "torque_cycle_generation": "periodic-cycle-8",
        "loss_cycle_generation": "periodic-cycle-8",
        "waveform_segments": [
            {
                "segment_id": "initial",
                "phase_origin_deg": 0.0,
                "start_time_s": 0.0,
                "end_time_s": 0.0049,
            },
            {
                "segment_id": "restart",
                "phase_origin_deg": 0.0,
                "start_time_s": 0.005,
                "end_time_s": 0.01,
            },
        ],
        "angle_convention": {
            "torque_angle_basis": "mechanical",
            "dq_current_angle_basis": "electrical",
            "joined_angle_basis": "mechanical",
            "pole_pairs": 4,
            "dq_to_joined_basis_transform_applied": True,
        },
        "loss_normalization": {
            "copper_loss_scope": "total_machine",
            "iron_loss_scope": "total_machine",
            "magnet_loss_scope": "total_machine",
            "phase_count": 3,
            "per_phase_to_total_applied": True,
        },
        "dq_phase_order": {
            "winding_connection_phase_order": ["U", "V", "W"],
            "current_table_phase_order": ["U", "V", "W"],
            "abc_to_dq_input_phase_order": ["U", "V", "W"],
            "phase_order_generation": "phase-order-12",
        },
        "torque_ripple_aggregation": {
            "cycle_generation": "periodic-cycle-8",
            "cycle_start_sample": 0,
            "cycle_end_sample_exclusive": 100,
            "exported_sample_count": 101,
            "aggregation_sample_count": 100,
            "repeated_cycle_endpoint_present": True,
            "repeated_endpoint_removed_before_aggregation": True,
        },
        "efficiency_average_window": {
            "input_window_start_sample": 75,
            "input_window_end_sample_exclusive": 100,
            "output_window_start_sample": 75,
            "output_window_end_sample_exclusive": 100,
            "periodic_cycle_generation": "periodic-cycle-8",
            "input_power_cycle_generation": "periodic-cycle-8",
            "output_power_cycle_generation": "periodic-cycle-8",
            "startup_samples_excluded": True,
        },
        "demag_temperature_material_state": {
            "magnet_temperature_c": 120.0,
            "recoil_curve_temperature_c": 120.0,
            "knee_curve_temperature_c": 120.0,
            "magnet_state_generation": "magnet-state-13",
            "recoil_curve_state_generation": "magnet-state-13",
            "knee_curve_state_generation": "magnet-state-13",
        },
        "torque_angle_basis_identity": {
            "pole_pairs": 4,
            "reference_angle_basis": "mechanical",
            "candidate_angle_basis": "mechanical",
            "waveform_alignment_basis": "mechanical",
            "reference_angle_grid_generation": "angle-grid-14",
            "candidate_angle_grid_generation": "angle-grid-14",
            "reference_to_electrical_scale": 4.0,
            "candidate_to_electrical_scale": 4.0,
        },
        "loss_harmonic_rotor_window_identity": {
            "window_angle_basis": "mechanical",
            "window_start_deg": 0.0,
            "window_end_deg": 90.0,
            "pole_pairs": 4,
            "expected_electrical_span_deg": 360.0,
            "rotor_position_generation": "rotor-position-14",
            "sample_rotor_position_generations": ["rotor-position-14"] * 8,
            "loss_solve_generation": "loss-solve-14",
            "harmonic_transform_solve_generation": "loss-solve-14",
        },
        "iron_loss_coefficient_frequency_basis_identity": {
            "waveform_frequency_basis_hz": [50.0, 150.0, 250.0],
            "hysteresis_coefficient_frequency_basis_hz": [50.0, 150.0, 250.0],
            "eddy_coefficient_frequency_basis_hz": [50.0, 150.0, 250.0],
            "coefficient_set_generation": "iron-loss-coefficients-15",
            "hysteresis_coefficient_generation": "iron-loss-coefficients-15",
            "eddy_coefficient_generation": "iron-loss-coefficients-15",
            "waveform_solve_generation": "loss-solve-15",
            "loss_result_solve_generation": "loss-solve-15",
        },
        "dq_current_phase_convention_identity": {
            "source_phase_order": "U-V-W",
            "dq_transform_phase_order": "U-V-W",
            "source_q_axis_lead": "q_leads_d_positive_electrical",
            "result_q_axis_lead": "q_leads_d_positive_electrical",
            "source_current_angle_deg_electrical": 30.0,
            "result_current_angle_deg_electrical": 30.0,
            "electrical_angle_zero_axis": "phase_U_positive_peak",
            "result_electrical_angle_zero_axis": "phase_U_positive_peak",
            "current_command_generation": "current-command-15",
            "result_generation": "current-command-15",
        },
        "torque_average_period_angle_basis_identity": {
            "pole_pairs": 4,
            "torque_sample_angle_basis": "electrical",
            "integration_window_angle_basis": "electrical",
            "reported_window_angle_basis": "electrical",
            "window_start_deg": 0.0,
            "window_end_deg": 360.0,
            "reported_window_span_deg": 360.0,
            "equivalent_mechanical_span_deg": 90.0,
            "sample_period_generation": "torque-period-16",
            "integration_period_generation": "torque-period-16",
            "result_period_generation": "torque-period-16",
        },
        "lamination_stacking_factor_loss_conductivity_identity": {
            "material_generation": "lamination-material-16",
            "stacking_factor_material_generation": "lamination-material-16",
            "loss_material_generation": "lamination-material-16",
            "stacking_factor": 0.95,
            "geometric_lamination_volume_m3": 1.0e-3,
            "effective_magnetic_volume_m3": 9.5e-4,
            "volume_stacking_factor_application_count": 1,
            "conductivity_basis": "lamination_effective_cross_section",
            "eddy_loss_conductivity_basis": "lamination_effective_cross_section",
            "conductivity_stacking_factor_application_count": 1,
            "solve_generation": "loss-solve-16",
            "loss_result_solve_generation": "loss-solve-16",
        },
        "dq_park_transform_power_invariant_scaling_identity": {
            "solve_generation": "dq-solve-17",
            "dq_voltage_result_generation": "dq-solve-17",
            "dq_current_result_generation": "dq-solve-17",
            "voltage_park_transform_basis": "power_invariant",
            "current_park_transform_basis": "power_invariant",
            "reported_power_basis": "power_invariant",
            "abc_instantaneous_power_w": 500.0,
            "dq0_instantaneous_power_w": 500.0,
            "dq0_power_scale_to_abc": 1.0,
            "power_scale_application_count": 1,
            "park_transform_sha256": "1" * 64,
            "power_closure_transform_sha256": "1" * 64,
        },
        "demag_recoil_temperature_operating_point_identity": {
            "material_generation": "magnet-material-17",
            "recoil_curve_material_generation": "magnet-material-17",
            "field_solution_material_generation": "magnet-material-17",
            "temperature_state_generation": "temperature-17",
            "recoil_curve_temperature_state_generation": "temperature-17",
            "operating_point_temperature_state_generation": "temperature-17",
            "magnet_temperature_k": 393.15,
            "recoil_curve_temperature_k": 393.15,
            "operating_point_temperature_k": 393.15,
            "solve_generation": "demag-solve-17",
            "operating_point_solve_generation": "demag-solve-17",
            "recoil_curve_sha256": "2" * 64,
            "operating_point_recoil_curve_sha256": "2" * 64,
        },
    }
    return payload


def test_pwm_motor_loss_gate_accepts_balanced_control_and_summary_row_semantics() -> None:
    result = pwm_controlled_motor_loss_gate(_payload())
    assert result["status"] == "ok"
    assert result["metrics"]["three_phase_kcl_global_relative_error"] < 1.0e-12
    assert result["checks"]["eddy_aggregate_matches_harmonic_bin_sum"] is True


def test_pwm_motor_loss_gate_rejects_control_tracking_failure() -> None:
    payload = copy.deepcopy(_payload())
    for row in payload["time_series"]["control_feedback_a"]:
        row[1] = 10.0
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert result["checks"]["tail_current_control_tracks_commands"] is False


def test_pwm_motor_loss_gate_rejects_wrong_iron_summary() -> None:
    payload = copy.deepcopy(_payload())
    payload["loss_spectrum"]["iron_components_w"][0][0] += 1.0
    payload["loss_spectrum"]["reported_iron_total_w"][0] += 1.0
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert result["checks"]["aggregate_iron_equals_eddy_plus_hysteresis"] is False


def test_pwm_motor_loss_gate_rejects_unguarded_ratio_outputs() -> None:
    payload = _payload()
    payload["ratio_diagnostic_policy"] = "accept_all_resistance_and_inductance_ratios"
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert result["checks"]["near_zero_current_ratio_outputs_are_diagnostic_only"] is False


def test_pwm_motor_loss_mcp_tool_dispatches() -> None:
    result = json.loads(mcp_gate(_payload()))
    assert result["status"] == "ok"


def test_generalization_v7_public_loss_component_resampling_alias() -> None:
    payload = _payload()
    time_s = payload["time_series"]["time_s"]
    payload["time_series"]["component_time_s"] = [
        time_s.copy(),
        [value + 5.0e-5 for value in time_s],
    ]
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert result["checks"]["power_component_time_axes_match_common_axis_knotwise"] is False


def test_accepts_bound_cycle_generation_and_restart_phase_origin() -> None:
    result = pwm_controlled_motor_loss_gate(_with_artifact_identity(_payload()))
    assert result["status"] == "ok"
    assert result["warnings"] == []


def test_v8_public_torque_loss_cycle_generation_mix() -> None:
    payload = _with_artifact_identity(_payload())
    payload["artifact_identity"]["loss_cycle_generation"] = "periodic-cycle-7"
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert result["checks"]["torque_and_loss_share_periodic_cycle_generation"] is False


def test_v8_public_restart_phase_origin_shift() -> None:
    payload = _with_artifact_identity(_payload())
    payload["artifact_identity"]["waveform_segments"][1]["phase_origin_deg"] = 30.0
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert result["checks"]["restart_segments_preserve_phase_origin"] is False


def test_v9_public_electrical_mechanical_angle_convention_mix() -> None:
    payload = _with_artifact_identity(_payload())
    payload["artifact_identity"]["angle_convention"].update(
        {
            "joined_angle_basis": "electrical",
            "dq_to_joined_basis_transform_applied": False,
        }
    )
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"]["torque_and_dq_tables_share_transformed_angle_basis"]
        is False
    )


def test_v9_public_per_phase_total_copper_loss_mix() -> None:
    payload = _with_artifact_identity(_payload())
    payload["artifact_identity"]["loss_normalization"].update(
        {
            "copper_loss_scope": "per_phase",
            "per_phase_to_total_applied": False,
        }
    )
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert result["checks"]["loss_components_share_total_machine_scope"] is False


def test_v10_public_dq_transform_phase_order_mismatch() -> None:
    payload = _with_artifact_identity(_payload())
    payload["artifact_identity"]["dq_phase_order"][
        "abc_to_dq_input_phase_order"
    ] = ["U", "W", "V"]
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"]["abc_to_dq_phase_order_matches_winding_connection"]
        is False
    )


def test_v10_public_torque_ripple_duplicate_cycle_endpoint() -> None:
    payload = _with_artifact_identity(_payload())
    payload["artifact_identity"]["torque_ripple_aggregation"].update(
        {
            "aggregation_sample_count": 101,
            "repeated_endpoint_removed_before_aggregation": False,
        }
    )
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "torque_ripple_aggregation_excludes_repeated_cycle_endpoint"
        ]
        is False
    )


def test_v11_public_motor_efficiency_average_window_mismatch() -> None:
    payload = _with_artifact_identity(_payload())
    payload["artifact_identity"]["efficiency_average_window"].update(
        {
            "input_window_start_sample": 0,
            "input_power_cycle_generation": "startup-cycle-8",
            "startup_samples_excluded": False,
        }
    )
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "efficiency_input_and_output_share_periodic_average_window"
        ]
        is False
    )


def test_v11_public_demag_margin_temperature_material_state_mismatch() -> None:
    payload = _with_artifact_identity(_payload())
    payload["artifact_identity"]["demag_temperature_material_state"].update(
        {
            "recoil_curve_temperature_c": 20.0,
            "knee_curve_temperature_c": 20.0,
            "recoil_curve_state_generation": "magnet-state-12",
            "knee_curve_state_generation": "magnet-state-12",
        }
    )
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "demag_margin_uses_current_temperature_material_state"
        ]
        is False
    )


def test_v12_public_motor_torque_electrical_mechanical_angle_mismatch() -> None:
    payload = _with_artifact_identity(_payload())
    payload["artifact_identity"]["torque_angle_basis_identity"].update(
        {
            "candidate_angle_basis": "electrical",
            "candidate_to_electrical_scale": 1.0,
        }
    )
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"]["motor_torque_waveforms_share_mechanical_angle_basis"]
        is False
    )


def test_v12_public_loss_harmonic_window_rotor_position_generation_mismatch() -> None:
    payload = _with_artifact_identity(_payload())
    payload["artifact_identity"]["loss_harmonic_rotor_window_identity"][
        "sample_rotor_position_generations"
    ][3] = "rotor-position-13"
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "loss_harmonic_window_uses_one_rotor_position_generation"
        ]
        is False
    )


def test_v13_public_iron_loss_coefficient_frequency_basis_mismatch() -> None:
    payload = _with_artifact_identity(_payload())
    payload["artifact_identity"][
        "iron_loss_coefficient_frequency_basis_identity"
    ].update(
        {
            "hysteresis_coefficient_frequency_basis_hz": [60.0, 180.0, 300.0],
            "hysteresis_coefficient_generation": "iron-loss-coefficients-14",
        }
    )
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"]["iron_loss_coefficients_share_waveform_frequency_basis"]
        is False
    )


def test_v13_public_dq_current_angle_phase_convention_mismatch() -> None:
    payload = _with_artifact_identity(_payload())
    payload["artifact_identity"]["dq_current_phase_convention_identity"].update(
        {
            "dq_transform_phase_order": "U-W-V",
            "result_q_axis_lead": "q_lags_d_positive_electrical",
        }
    )
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"]["dq_currents_share_phase_order_and_q_axis_convention"]
        is False
    )


def test_v14_public_torque_average_electrical_mechanical_period_mismatch() -> None:
    payload = _with_artifact_identity(_payload())
    payload["artifact_identity"]["torque_average_period_angle_basis_identity"].update(
        {
            "reported_window_angle_basis": "mechanical",
            "equivalent_mechanical_span_deg": 360.0,
            "result_period_generation": "torque-period-15",
        }
    )
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "torque_average_uses_one_electrical_or_mechanical_period_basis"
        ]
        is False
    )


def test_v14_public_lamination_stacking_factor_loss_conductivity_basis_mismatch() -> None:
    payload = _with_artifact_identity(_payload())
    payload["artifact_identity"][
        "lamination_stacking_factor_loss_conductivity_identity"
    ].update(
        {
            "eddy_loss_conductivity_basis": "bulk_material_cross_section",
            "conductivity_stacking_factor_application_count": 0,
            "loss_result_solve_generation": "loss-solve-15",
        }
    )
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "lamination_loss_uses_consistent_stacking_factor_and_conductivity_basis"
        ]
        is False
    )


def test_v15_public_dq_park_transform_power_invariant_scaling_mismatch() -> None:
    payload = _with_artifact_identity(_payload())
    payload["artifact_identity"][
        "dq_park_transform_power_invariant_scaling_identity"
    ].update(
        {
            "current_park_transform_basis": "amplitude_invariant",
            "dq0_power_scale_to_abc": 1.5,
            "dq0_instantaneous_power_w": 750.0,
            "power_closure_transform_sha256": "5" * 64,
        }
    )
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "dq_power_uses_one_power_invariant_park_transform_scaling"
        ]
        is False
    )


def test_v15_public_demag_recoil_temperature_operating_point_generation_mismatch() -> None:
    payload = _with_artifact_identity(_payload())
    payload["artifact_identity"][
        "demag_recoil_temperature_operating_point_identity"
    ].update(
        {
            "recoil_curve_temperature_state_generation": "temperature-16",
            "recoil_curve_temperature_k": 353.15,
            "operating_point_solve_generation": "demag-solve-16",
            "operating_point_recoil_curve_sha256": "5" * 64,
        }
    )
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "demag_operating_point_uses_current_recoil_temperature_state"
        ]
        is False
    )
