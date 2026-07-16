"""Solver-neutral gates for PWM current control and harmonic-loss result tables."""

from __future__ import annotations

import math
from typing import Any


def _vector(value: Any, name: str, *, minimum: int = 1) -> list[float]:
    try:
        parsed = [float(item) for item in value]
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a numeric sequence") from exc
    if len(parsed) < minimum or not all(math.isfinite(item) for item in parsed):
        raise ValueError(f"{name} must contain at least {minimum} finite values")
    return parsed


def _matrix(
    value: Any,
    name: str,
    *,
    rows: int,
    columns: int | None = None,
) -> list[list[float]]:
    try:
        parsed = [[float(item) for item in row] for row in value]
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a numeric matrix") from exc
    if len(parsed) != rows or not parsed:
        raise ValueError(f"{name} must have {rows} rows")
    width = len(parsed[0])
    if width == 0 or (columns is not None and width != columns):
        raise ValueError(f"{name} has an invalid column count")
    if any(len(row) != width for row in parsed):
        raise ValueError(f"{name} rows must have equal width")
    if not all(math.isfinite(item) for row in parsed for item in row):
        raise ValueError(f"{name} must contain finite values")
    return parsed


def _relative(residual: float, scale: float) -> float:
    return abs(residual) / max(abs(scale), 1.0e-30)


def pwm_controlled_motor_loss_gate(
    payload: dict[str, Any],
    *,
    max_three_phase_kcl_relative_error: float = 1.0e-3,
    max_tail_control_tracking_rms_relative_error: float = 0.05,
    max_angle_speed_integral_relative_error: float = 2.0e-4,
    max_power_sum_relative_error: float = 1.0e-10,
    max_loss_identity_relative_error: float = 1.0e-10,
    max_frequency_step_relative_span: float = 1.0e-10,
) -> dict[str, Any]:
    """Validate a PWM-controlled motor time series and harmonic-loss package.

    The loss spectrum uses a common table convention: row zero stores aggregate
    values, while positive-frequency rows store harmonic bins. Eddy loss is
    reconstructable from those bins; hysteresis and combined iron loss may be
    aggregate-only. Ratio-derived R/L values are excluded when current is near
    zero because they are ill-conditioned rather than physical outliers.
    """

    if not isinstance(payload, dict):
        raise ValueError("payload must be an object")
    tolerances = [
        max_three_phase_kcl_relative_error,
        max_tail_control_tracking_rms_relative_error,
        max_angle_speed_integral_relative_error,
        max_power_sum_relative_error,
        max_loss_identity_relative_error,
        max_frequency_step_relative_span,
    ]
    if any(not math.isfinite(float(value)) or float(value) < 0.0 for value in tolerances):
        raise ValueError("tolerances must be finite and nonnegative")

    time_series = payload.get("time_series")
    loss_spectrum = payload.get("loss_spectrum")
    if not isinstance(time_series, dict) or not isinstance(loss_spectrum, dict):
        raise ValueError("time_series and loss_spectrum objects are required")

    identity_value = payload.get("artifact_identity")
    identity_present = isinstance(identity_value, dict)
    cycle_generation_ok = True
    restart_phase_origin_ok = True
    angle_convention_ok = True
    loss_normalization_ok = True
    dq_phase_order_ok = True
    torque_ripple_aggregation_ok = True
    efficiency_average_window_ok = True
    demag_temperature_material_state_ok = True
    torque_angle_basis_identity_ok = True
    loss_harmonic_rotor_window_identity_ok = True
    iron_loss_coefficient_frequency_basis_identity_ok = True
    dq_current_phase_convention_identity_ok = True
    torque_average_period_angle_basis_identity_ok = True
    lamination_stacking_factor_loss_conductivity_identity_ok = True
    dq_park_transform_power_invariant_scaling_identity_ok = True
    demag_recoil_temperature_operating_point_identity_ok = True
    torque_ripple_fft_angle_endpoint_generation_identity_ok = True
    iron_loss_spatial_harmonic_mesh_volume_identity_ok = True
    dq_torque_map_park_transform_angle_sign_identity_ok = True
    efficiency_map_power_averaging_window_identity_ok = True
    if identity_value is not None and not identity_present:
        cycle_generation_ok = False
        restart_phase_origin_ok = False
        angle_convention_ok = False
        loss_normalization_ok = False
        dq_phase_order_ok = False
        torque_ripple_aggregation_ok = False
        efficiency_average_window_ok = False
        demag_temperature_material_state_ok = False
        torque_angle_basis_identity_ok = False
        loss_harmonic_rotor_window_identity_ok = False
        iron_loss_coefficient_frequency_basis_identity_ok = False
        dq_current_phase_convention_identity_ok = False
        torque_average_period_angle_basis_identity_ok = False
        lamination_stacking_factor_loss_conductivity_identity_ok = False
        dq_park_transform_power_invariant_scaling_identity_ok = False
        demag_recoil_temperature_operating_point_identity_ok = False
        torque_ripple_fft_angle_endpoint_generation_identity_ok = False
        iron_loss_spatial_harmonic_mesh_volume_identity_ok = False
        dq_torque_map_park_transform_angle_sign_identity_ok = False
        efficiency_map_power_averaging_window_identity_ok = False
    elif identity_present:
        torque_generation = str(identity_value.get("torque_cycle_generation", ""))
        loss_generation = str(identity_value.get("loss_cycle_generation", ""))
        cycle_generation_ok = (
            bool(torque_generation) and torque_generation == loss_generation
        )
        segments = identity_value.get("waveform_segments")
        if not isinstance(segments, list) or not segments:
            restart_phase_origin_ok = False
        else:
            phase_origins = []
            previous_end = -math.inf
            for segment in segments:
                if not isinstance(segment, dict):
                    restart_phase_origin_ok = False
                    break
                try:
                    phase_origin = float(segment["phase_origin_deg"])
                    start = float(segment["start_time_s"])
                    end = float(segment["end_time_s"])
                except (KeyError, TypeError, ValueError):
                    restart_phase_origin_ok = False
                    break
                if not all(math.isfinite(value) for value in (phase_origin, start, end)):
                    restart_phase_origin_ok = False
                    break
                if not str(segment.get("segment_id", "")) or end < start or start < previous_end:
                    restart_phase_origin_ok = False
                    break
                phase_origins.append(phase_origin)
                previous_end = end
            restart_phase_origin_ok = (
                restart_phase_origin_ok
                and bool(phase_origins)
                and len(set(phase_origins)) == 1
            )
        angle_convention = identity_value.get("angle_convention")
        if angle_convention is not None:
            angle_convention_ok = (
                isinstance(angle_convention, dict)
                and angle_convention.get("torque_angle_basis") == "mechanical"
                and angle_convention.get("dq_current_angle_basis") == "electrical"
                and angle_convention.get("joined_angle_basis") == "mechanical"
                and isinstance(angle_convention.get("pole_pairs"), int)
                and angle_convention["pole_pairs"] > 0
                and angle_convention.get("dq_to_joined_basis_transform_applied")
                is True
            )
        loss_normalization = identity_value.get("loss_normalization")
        if loss_normalization is not None:
            loss_normalization_ok = (
                isinstance(loss_normalization, dict)
                and loss_normalization.get("copper_loss_scope") == "total_machine"
                and loss_normalization.get("iron_loss_scope") == "total_machine"
                and loss_normalization.get("magnet_loss_scope") == "total_machine"
                and isinstance(loss_normalization.get("phase_count"), int)
                and loss_normalization["phase_count"] >= 2
                and loss_normalization.get("per_phase_to_total_applied") is True
            )
        phase_order = identity_value.get("dq_phase_order")
        if phase_order is not None:
            winding_order = (
                phase_order.get("winding_connection_phase_order")
                if isinstance(phase_order, dict)
                else None
            )
            current_order = (
                phase_order.get("current_table_phase_order")
                if isinstance(phase_order, dict)
                else None
            )
            transform_order = (
                phase_order.get("abc_to_dq_input_phase_order")
                if isinstance(phase_order, dict)
                else None
            )
            dq_phase_order_ok = (
                isinstance(phase_order, dict)
                and isinstance(winding_order, list)
                and len(winding_order) == 3
                and all(isinstance(name, str) and name for name in winding_order)
                and len(set(winding_order)) == 3
                and current_order == winding_order
                and transform_order == winding_order
                and bool(phase_order.get("phase_order_generation"))
            )
        aggregation = identity_value.get("torque_ripple_aggregation")
        if aggregation is not None:
            try:
                start = int(aggregation["cycle_start_sample"])
                end = int(aggregation["cycle_end_sample_exclusive"])
                exported_count = int(aggregation["exported_sample_count"])
                aggregation_count = int(aggregation["aggregation_sample_count"])
            except (KeyError, TypeError, ValueError):
                start = end = exported_count = aggregation_count = -1
            unique_cycle_count = end - start
            torque_ripple_aggregation_ok = (
                isinstance(aggregation, dict)
                and start >= 0
                and unique_cycle_count >= 3
                and exported_count == unique_cycle_count + 1
                and aggregation_count == unique_cycle_count
                and aggregation.get("repeated_cycle_endpoint_present") is True
                and aggregation.get("repeated_endpoint_removed_before_aggregation")
                is True
                and aggregation.get("cycle_generation") == torque_generation
            )
        efficiency_window = identity_value.get("efficiency_average_window")
        if efficiency_window is not None:
            try:
                input_start = int(efficiency_window["input_window_start_sample"])
                input_end = int(
                    efficiency_window["input_window_end_sample_exclusive"]
                )
                output_start = int(
                    efficiency_window["output_window_start_sample"]
                )
                output_end = int(
                    efficiency_window["output_window_end_sample_exclusive"]
                )
            except (KeyError, TypeError, ValueError):
                input_start = input_end = output_start = output_end = -1
            efficiency_average_window_ok = (
                isinstance(efficiency_window, dict)
                and input_start >= 0
                and input_end > input_start
                and input_start == output_start
                and input_end == output_end
                and efficiency_window.get("periodic_cycle_generation")
                == torque_generation
                and efficiency_window.get("input_power_cycle_generation")
                == torque_generation
                and efficiency_window.get("output_power_cycle_generation")
                == torque_generation
                and efficiency_window.get("startup_samples_excluded") is True
            )
        demag_state = identity_value.get("demag_temperature_material_state")
        if demag_state is not None:
            try:
                magnet_temperature = float(demag_state["magnet_temperature_c"])
                recoil_temperature = float(
                    demag_state["recoil_curve_temperature_c"]
                )
                knee_temperature = float(demag_state["knee_curve_temperature_c"])
            except (KeyError, TypeError, ValueError):
                magnet_temperature = recoil_temperature = knee_temperature = math.nan
            state_generation = str(
                demag_state.get("magnet_state_generation", "")
            ) if isinstance(demag_state, dict) else ""
            demag_temperature_material_state_ok = (
                isinstance(demag_state, dict)
                and all(
                    math.isfinite(value)
                    for value in (
                        magnet_temperature,
                        recoil_temperature,
                        knee_temperature,
                    )
                )
                and recoil_temperature == magnet_temperature
                and knee_temperature == magnet_temperature
                and bool(state_generation)
                and demag_state.get("recoil_curve_state_generation")
                == state_generation
                and demag_state.get("knee_curve_state_generation")
                == state_generation
            )

        torque_basis_value = identity_value.get("torque_angle_basis_identity")
        if torque_basis_value is not None:
            torque_basis = (
                torque_basis_value if isinstance(torque_basis_value, dict) else {}
            )
            pole_pairs = torque_basis.get("pole_pairs")
            reference_scale = torque_basis.get("reference_to_electrical_scale")
            candidate_scale = torque_basis.get("candidate_to_electrical_scale")
            reference_grid = str(
                torque_basis.get("reference_angle_grid_generation", "")
            )
            candidate_grid = str(
                torque_basis.get("candidate_angle_grid_generation", "")
            )
            torque_angle_basis_identity_ok = (
                isinstance(pole_pairs, int)
                and not isinstance(pole_pairs, bool)
                and pole_pairs > 0
                and torque_basis.get("reference_angle_basis") == "mechanical"
                and torque_basis.get("candidate_angle_basis") == "mechanical"
                and torque_basis.get("waveform_alignment_basis") == "mechanical"
                and bool(reference_grid)
                and reference_grid == candidate_grid
                and isinstance(reference_scale, (int, float))
                and not isinstance(reference_scale, bool)
                and isinstance(candidate_scale, (int, float))
                and not isinstance(candidate_scale, bool)
                and math.isfinite(float(reference_scale))
                and math.isfinite(float(candidate_scale))
                and math.isclose(float(reference_scale), float(pole_pairs))
                and math.isclose(float(candidate_scale), float(pole_pairs))
            )

        loss_window_value = identity_value.get(
            "loss_harmonic_rotor_window_identity"
        )
        if loss_window_value is not None:
            loss_window = (
                loss_window_value if isinstance(loss_window_value, dict) else {}
            )
            sample_generations = loss_window.get(
                "sample_rotor_position_generations"
            )
            pole_pairs = loss_window.get("pole_pairs")
            try:
                start_deg = float(loss_window["window_start_deg"])
                end_deg = float(loss_window["window_end_deg"])
                expected_span_deg = float(
                    loss_window["expected_electrical_span_deg"]
                )
            except (KeyError, TypeError, ValueError):
                start_deg = end_deg = expected_span_deg = math.nan
            rotor_generation = str(
                loss_window.get("rotor_position_generation", "")
            )
            loss_solve_generation = str(
                loss_window.get("loss_solve_generation", "")
            )
            loss_harmonic_rotor_window_identity_ok = (
                loss_window.get("window_angle_basis") == "mechanical"
                and isinstance(pole_pairs, int)
                and not isinstance(pole_pairs, bool)
                and pole_pairs > 0
                and all(
                    math.isfinite(value)
                    for value in (start_deg, end_deg, expected_span_deg)
                )
                and end_deg > start_deg
                and expected_span_deg > 0.0
                and math.isclose(
                    float(pole_pairs) * (end_deg - start_deg),
                    expected_span_deg,
                    rel_tol=1.0e-12,
                    abs_tol=1.0e-12,
                )
                and bool(rotor_generation)
                and isinstance(sample_generations, list)
                and bool(sample_generations)
                and all(
                    generation == rotor_generation
                    for generation in sample_generations
                )
                and bool(loss_solve_generation)
                and loss_window.get("harmonic_transform_solve_generation")
                == loss_solve_generation
            )

        coefficient_basis_value = identity_value.get(
            "iron_loss_coefficient_frequency_basis_identity"
        )
        if coefficient_basis_value is not None:
            coefficient_basis = (
                coefficient_basis_value
                if isinstance(coefficient_basis_value, dict)
                else {}
            )
            try:
                waveform_frequencies = [
                    float(value)
                    for value in coefficient_basis["waveform_frequency_basis_hz"]
                ]
                hysteresis_frequencies = [
                    float(value)
                    for value in coefficient_basis[
                        "hysteresis_coefficient_frequency_basis_hz"
                    ]
                ]
                eddy_frequencies = [
                    float(value)
                    for value in coefficient_basis[
                        "eddy_coefficient_frequency_basis_hz"
                    ]
                ]
            except (KeyError, TypeError, ValueError):
                waveform_frequencies = []
                hysteresis_frequencies = []
                eddy_frequencies = []
            coefficient_generation = str(
                coefficient_basis.get("coefficient_set_generation", "")
            )
            waveform_solve_generation = str(
                coefficient_basis.get("waveform_solve_generation", "")
            )
            iron_loss_coefficient_frequency_basis_identity_ok = (
                bool(waveform_frequencies)
                and all(
                    math.isfinite(value) and value > 0.0
                    for value in waveform_frequencies
                    + hysteresis_frequencies
                    + eddy_frequencies
                )
                and all(
                    right > left
                    for left, right in zip(
                        waveform_frequencies, waveform_frequencies[1:]
                    )
                )
                and hysteresis_frequencies == waveform_frequencies
                and eddy_frequencies == waveform_frequencies
                and bool(coefficient_generation)
                and coefficient_basis.get("hysteresis_coefficient_generation")
                == coefficient_generation
                and coefficient_basis.get("eddy_coefficient_generation")
                == coefficient_generation
                and bool(waveform_solve_generation)
                and coefficient_basis.get("loss_result_solve_generation")
                == waveform_solve_generation
            )

        dq_convention_value = identity_value.get(
            "dq_current_phase_convention_identity"
        )
        if dq_convention_value is not None:
            dq_convention = (
                dq_convention_value if isinstance(dq_convention_value, dict) else {}
            )
            try:
                source_angle = float(
                    dq_convention["source_current_angle_deg_electrical"]
                )
                result_angle = float(
                    dq_convention["result_current_angle_deg_electrical"]
                )
            except (KeyError, TypeError, ValueError):
                source_angle = result_angle = math.nan
            source_phase_order = str(dq_convention.get("source_phase_order", ""))
            source_q_axis_lead = str(dq_convention.get("source_q_axis_lead", ""))
            zero_axis = str(dq_convention.get("electrical_angle_zero_axis", ""))
            current_generation = str(
                dq_convention.get("current_command_generation", "")
            )
            dq_current_phase_convention_identity_ok = (
                source_phase_order
                in {"U-V-W", "U-W-V", "V-U-W", "V-W-U", "W-U-V", "W-V-U"}
                and dq_convention.get("dq_transform_phase_order")
                == source_phase_order
                and source_q_axis_lead
                in {
                    "q_leads_d_positive_electrical",
                    "q_lags_d_positive_electrical",
                }
                and dq_convention.get("result_q_axis_lead") == source_q_axis_lead
                and math.isfinite(source_angle)
                and math.isfinite(result_angle)
                and math.isclose(
                    source_angle, result_angle, rel_tol=0.0, abs_tol=1.0e-12
                )
                and bool(zero_axis)
                and dq_convention.get("result_electrical_angle_zero_axis")
                == zero_axis
                and bool(current_generation)
                and dq_convention.get("result_generation") == current_generation
            )

        torque_period_value = identity_value.get(
            "torque_average_period_angle_basis_identity"
        )
        if torque_period_value is not None:
            torque_period = (
                torque_period_value
                if isinstance(torque_period_value, dict)
                else {}
            )
            try:
                pole_pairs = int(torque_period.get("pole_pairs"))
                start_deg = float(torque_period.get("window_start_deg"))
                end_deg = float(torque_period.get("window_end_deg"))
                reported_span_deg = float(
                    torque_period.get("reported_window_span_deg")
                )
                equivalent_mechanical_span_deg = float(
                    torque_period.get("equivalent_mechanical_span_deg")
                )
            except (TypeError, ValueError):
                pole_pairs = 0
                start_deg = end_deg = reported_span_deg = math.nan
                equivalent_mechanical_span_deg = math.nan
            sample_basis = str(
                torque_period.get("torque_sample_angle_basis", "")
            )
            period_generation = str(
                torque_period.get("sample_period_generation", "")
            )
            span_deg = end_deg - start_deg
            expected_mechanical_span = (
                span_deg / pole_pairs
                if sample_basis == "electrical" and pole_pairs > 0
                else span_deg
            )
            torque_average_period_angle_basis_identity_ok = (
                pole_pairs > 0
                and sample_basis in {"electrical", "mechanical"}
                and torque_period.get("integration_window_angle_basis")
                == sample_basis
                and torque_period.get("reported_window_angle_basis")
                == sample_basis
                and all(
                    math.isfinite(value)
                    for value in (
                        start_deg,
                        end_deg,
                        reported_span_deg,
                        equivalent_mechanical_span_deg,
                    )
                )
                and span_deg > 0.0
                and math.isclose(
                    reported_span_deg, span_deg, rel_tol=0.0, abs_tol=1.0e-12
                )
                and math.isclose(
                    equivalent_mechanical_span_deg,
                    expected_mechanical_span,
                    rel_tol=0.0,
                    abs_tol=1.0e-12,
                )
                and bool(period_generation)
                and torque_period.get("integration_period_generation")
                == period_generation
                and torque_period.get("result_period_generation")
                == period_generation
            )

        lamination_value = identity_value.get(
            "lamination_stacking_factor_loss_conductivity_identity"
        )
        if lamination_value is not None:
            lamination = lamination_value if isinstance(lamination_value, dict) else {}
            try:
                stacking_factor = float(lamination.get("stacking_factor"))
                geometric_volume = float(
                    lamination.get("geometric_lamination_volume_m3")
                )
                effective_volume = float(
                    lamination.get("effective_magnetic_volume_m3")
                )
                volume_application_count = int(
                    lamination.get("volume_stacking_factor_application_count")
                )
                conductivity_application_count = int(
                    lamination.get(
                        "conductivity_stacking_factor_application_count"
                    )
                )
            except (TypeError, ValueError):
                stacking_factor = geometric_volume = effective_volume = math.nan
                volume_application_count = conductivity_application_count = -1
            material_generation = str(
                lamination.get("material_generation", "")
            )
            solve_generation = str(lamination.get("solve_generation", ""))
            conductivity_basis = str(lamination.get("conductivity_basis", ""))
            lamination_stacking_factor_loss_conductivity_identity_ok = (
                math.isfinite(stacking_factor)
                and 0.0 < stacking_factor <= 1.0
                and math.isfinite(geometric_volume)
                and geometric_volume > 0.0
                and math.isfinite(effective_volume)
                and math.isclose(
                    effective_volume,
                    geometric_volume * stacking_factor,
                    rel_tol=1.0e-12,
                    abs_tol=1.0e-15,
                )
                and volume_application_count == 1
                and conductivity_application_count == 1
                and conductivity_basis == "lamination_effective_cross_section"
                and lamination.get("eddy_loss_conductivity_basis")
                == conductivity_basis
                and bool(material_generation)
                and lamination.get("stacking_factor_material_generation")
                == material_generation
                and lamination.get("loss_material_generation")
                == material_generation
                and bool(solve_generation)
                and lamination.get("loss_result_solve_generation")
                == solve_generation
            )

        dq_power_value = identity_value.get(
            "dq_park_transform_power_invariant_scaling_identity"
        )
        if dq_power_value is not None:
            dq_power = dq_power_value if isinstance(dq_power_value, dict) else {}
            try:
                abc_power = float(dq_power.get("abc_instantaneous_power_w"))
                dq0_power = float(dq_power.get("dq0_instantaneous_power_w"))
                power_scale = float(dq_power.get("dq0_power_scale_to_abc"))
                scale_count = int(dq_power.get("power_scale_application_count"))
            except (TypeError, ValueError):
                abc_power = dq0_power = power_scale = math.nan
                scale_count = -1
            solve_generation = str(dq_power.get("solve_generation", ""))
            transform_digest = str(
                dq_power.get("park_transform_sha256", "")
            ).lower()
            dq_park_transform_power_invariant_scaling_identity_ok = (
                bool(solve_generation)
                and dq_power.get("dq_voltage_result_generation")
                == solve_generation
                and dq_power.get("dq_current_result_generation")
                == solve_generation
                and dq_power.get("voltage_park_transform_basis")
                == "power_invariant"
                and dq_power.get("current_park_transform_basis")
                == "power_invariant"
                and dq_power.get("reported_power_basis") == "power_invariant"
                and all(
                    math.isfinite(value)
                    for value in (abc_power, dq0_power, power_scale)
                )
                and math.isclose(power_scale, 1.0, rel_tol=0.0, abs_tol=0.0)
                and scale_count == 1
                and math.isclose(
                    dq0_power * power_scale,
                    abc_power,
                    rel_tol=1.0e-12,
                    abs_tol=1.0e-12,
                )
                and len(transform_digest) == 64
                and all(
                    character in "0123456789abcdef"
                    for character in transform_digest
                )
                and str(
                    dq_power.get("power_closure_transform_sha256", "")
                ).lower()
                == transform_digest
            )

        demag_value = identity_value.get(
            "demag_recoil_temperature_operating_point_identity"
        )
        if demag_value is not None:
            demag = demag_value if isinstance(demag_value, dict) else {}
            try:
                magnet_temperature = float(demag.get("magnet_temperature_k"))
                recoil_temperature = float(
                    demag.get("recoil_curve_temperature_k")
                )
                operating_temperature = float(
                    demag.get("operating_point_temperature_k")
                )
            except (TypeError, ValueError):
                magnet_temperature = recoil_temperature = math.nan
                operating_temperature = math.nan
            material_generation = str(demag.get("material_generation", ""))
            temperature_generation = str(
                demag.get("temperature_state_generation", "")
            )
            solve_generation = str(demag.get("solve_generation", ""))
            recoil_digest = str(demag.get("recoil_curve_sha256", "")).lower()
            demag_recoil_temperature_operating_point_identity_ok = (
                bool(material_generation)
                and demag.get("recoil_curve_material_generation")
                == material_generation
                and demag.get("field_solution_material_generation")
                == material_generation
                and bool(temperature_generation)
                and demag.get("recoil_curve_temperature_state_generation")
                == temperature_generation
                and demag.get("operating_point_temperature_state_generation")
                == temperature_generation
                and all(
                    math.isfinite(value)
                    for value in (
                        magnet_temperature,
                        recoil_temperature,
                        operating_temperature,
                    )
                )
                and magnet_temperature > 0.0
                and math.isclose(
                    recoil_temperature,
                    magnet_temperature,
                    rel_tol=0.0,
                    abs_tol=1.0e-12,
                )
                and math.isclose(
                    operating_temperature,
                    magnet_temperature,
                    rel_tol=0.0,
                    abs_tol=1.0e-12,
                )
                and bool(solve_generation)
                and demag.get("operating_point_solve_generation")
                == solve_generation
                and len(recoil_digest) == 64
                and all(
                    character in "0123456789abcdef" for character in recoil_digest
                )
                and str(
                    demag.get("operating_point_recoil_curve_sha256", "")
                ).lower()
                == recoil_digest
            )

        torque_fft_value = identity_value.get(
            "torque_ripple_fft_angle_endpoint_generation_identity"
        )
        if torque_fft_value is not None:
            torque_fft = torque_fft_value if isinstance(torque_fft_value, dict) else {}
            solve_generation = str(torque_fft.get("torque_solve_generation", ""))
            resampling_generation = str(
                torque_fft.get("angle_resampling_generation", "")
            )
            angle_digest = str(torque_fft.get("resampled_angle_sha256", "")).lower()
            try:
                pole_pairs = int(torque_fft.get("pole_pairs"))
                start_deg = float(torque_fft.get("periodic_start_deg"))
                end_deg = float(torque_fft.get("periodic_end_deg"))
            except (TypeError, ValueError):
                pole_pairs = 0
                start_deg = end_deg = math.nan
            torque_ripple_fft_angle_endpoint_generation_identity_ok = (
                bool(solve_generation)
                and torque_fft.get("torque_sample_generation") == solve_generation
                and bool(resampling_generation)
                and torque_fft.get("fft_endpoint_resampling_generation")
                == resampling_generation
                and torque_fft.get("sample_angle_basis") == "electrical"
                and torque_fft.get("fft_endpoint_angle_basis") == "electrical"
                and pole_pairs > 0
                and math.isfinite(start_deg)
                and math.isfinite(end_deg)
                and math.isclose(
                    end_deg - start_deg, 360.0, rel_tol=0.0, abs_tol=1.0e-12
                )
                and torque_fft.get("duplicate_endpoint_removed") is True
                and len(angle_digest) == 64
                and all(character in "0123456789abcdef" for character in angle_digest)
                and str(torque_fft.get("fft_angle_sha256", "")).lower()
                == angle_digest
            )

        loss_volume_value = identity_value.get(
            "iron_loss_spatial_harmonic_mesh_volume_identity"
        )
        if loss_volume_value is not None:
            loss_volume = (
                loss_volume_value if isinstance(loss_volume_value, dict) else {}
            )
            solve_generation = str(loss_volume.get("loss_solve_generation", ""))
            mesh_generation = str(loss_volume.get("field_mesh_generation", ""))
            topology_generation = str(
                loss_volume.get("adaptive_mesh_topology_generation", "")
            )
            loss_ids = loss_volume.get("loss_element_ids")
            volume_ids = loss_volume.get("volume_weight_element_ids")
            volume_digest = str(
                loss_volume.get("element_volume_sha256", "")
            ).lower()
            iron_loss_spatial_harmonic_mesh_volume_identity_ok = (
                bool(solve_generation)
                and loss_volume.get("spatial_harmonic_result_generation")
                == solve_generation
                and bool(mesh_generation)
                and loss_volume.get("loss_density_mesh_generation")
                == mesh_generation
                and loss_volume.get("element_volume_mesh_generation")
                == mesh_generation
                and bool(topology_generation)
                and loss_volume.get("volume_weight_topology_generation")
                == topology_generation
                and isinstance(loss_ids, list)
                and bool(loss_ids)
                and all(
                    isinstance(value, int) and not isinstance(value, bool)
                    for value in loss_ids
                )
                and len(set(loss_ids)) == len(loss_ids)
                and volume_ids == loss_ids
                and len(volume_digest) == 64
                and all(character in "0123456789abcdef" for character in volume_digest)
                and str(loss_volume.get("loss_volume_weight_sha256", "")).lower()
                == volume_digest
            )

        dq_torque_map_value = identity_value.get(
            "dq_torque_map_park_transform_angle_sign_identity"
        )
        if dq_torque_map_value is not None:
            dq_torque_map = (
                dq_torque_map_value
                if isinstance(dq_torque_map_value, dict)
                else {}
            )
            map_generation = str(dq_torque_map.get("torque_map_generation", ""))
            park_generation = str(
                dq_torque_map.get("park_transform_generation", "")
            )
            park_digest = str(dq_torque_map.get("park_matrix_sha256", "")).lower()
            q_axis_sign = dq_torque_map.get("q_axis_sign")
            dq_torque_map_park_transform_angle_sign_identity_ok = (
                bool(map_generation)
                and dq_torque_map.get("park_transform_torque_map_generation")
                == map_generation
                and dq_torque_map.get("dq_current_map_torque_map_generation")
                == map_generation
                and bool(park_generation)
                and dq_torque_map.get("torque_map_park_transform_generation")
                == park_generation
                and dq_torque_map.get("electrical_angle_origin") == "rotor_d_axis"
                and dq_torque_map.get("torque_map_electrical_angle_origin")
                == dq_torque_map.get("electrical_angle_origin")
                and q_axis_sign in (-1, 1)
                and dq_torque_map.get("torque_map_q_axis_sign") == q_axis_sign
                and len(park_digest) == 64
                and all(character in "0123456789abcdef" for character in park_digest)
                and str(
                    dq_torque_map.get("torque_map_park_matrix_sha256", "")
                ).lower()
                == park_digest
            )

        efficiency_window_value = identity_value.get(
            "efficiency_map_power_averaging_window_identity"
        )
        if efficiency_window_value is not None:
            efficiency_window = (
                efficiency_window_value
                if isinstance(efficiency_window_value, dict)
                else {}
            )
            map_generation = str(
                efficiency_window.get("efficiency_map_generation", "")
            )
            cycle_generation = str(
                efficiency_window.get("steady_cycle_generation", "")
            )
            window_digest = str(
                efficiency_window.get("power_window_sha256", "")
            ).lower()
            try:
                start = int(efficiency_window.get("window_start_sample"))
                end = int(efficiency_window.get("window_end_sample"))
                electrical_window = [
                    int(value)
                    for value in efficiency_window.get("electrical_power_window", [])
                ]
                mechanical_window = [
                    int(value)
                    for value in efficiency_window.get("mechanical_power_window", [])
                ]
            except (TypeError, ValueError):
                start = end = -1
                electrical_window = []
                mechanical_window = []
            efficiency_map_power_averaging_window_identity_ok = (
                bool(map_generation)
                and efficiency_window.get("electrical_input_power_map_generation")
                == map_generation
                and efficiency_window.get("mechanical_output_power_map_generation")
                == map_generation
                and bool(cycle_generation)
                and efficiency_window.get("electrical_power_window_cycle_generation")
                == cycle_generation
                and efficiency_window.get("mechanical_power_window_cycle_generation")
                == cycle_generation
                and efficiency_window.get("loss_power_window_cycle_generation")
                == cycle_generation
                and 0 <= start < end
                and electrical_window == [start, end]
                and mechanical_window == [start, end]
                and len(window_digest) == 64
                and all(
                    character in "0123456789abcdef" for character in window_digest
                )
                and str(
                    efficiency_window.get("efficiency_power_window_sha256", "")
                ).lower()
                == window_digest
            )

    time_s = _vector(time_series.get("time_s"), "time_series.time_s", minimum=5)
    count = len(time_s)
    angle_deg = _vector(time_series.get("angle_deg"), "time_series.angle_deg", minimum=count)
    speed_rpm = _vector(time_series.get("speed_rpm"), "time_series.speed_rpm", minimum=count)
    torque_nm = _vector(time_series.get("torque_nm"), "time_series.torque_nm", minimum=count)
    total_power = _vector(
        time_series.get("reported_total_power_w"),
        "time_series.reported_total_power_w",
        minimum=count,
    )
    if any(len(values) != count for values in (angle_deg, speed_rpm, torque_nm, total_power)):
        raise ValueError("time-series vectors must have equal lengths")
    phase_currents = _matrix(
        time_series.get("phase_currents_a"),
        "time_series.phase_currents_a",
        rows=count,
        columns=3,
    )
    commands = _matrix(
        time_series.get("control_command_a"),
        "time_series.control_command_a",
        rows=count,
        columns=2,
    )
    feedback = _matrix(
        time_series.get("control_feedback_a"),
        "time_series.control_feedback_a",
        rows=count,
        columns=2,
    )
    power_components = _matrix(
        time_series.get("power_components_w"),
        "time_series.power_components_w",
        rows=count,
    )
    component_time_value = time_series.get("component_time_s")
    component_time_axes: list[list[float]] | None = None
    component_time_alignment_error = 0.0
    if component_time_value is not None:
        component_time_axes = _matrix(
            component_time_value,
            "time_series.component_time_s",
            rows=len(power_components[0]),
            columns=count,
        )
        component_time_alignment_error = max(
            abs(component_time_axes[component][index] - time_s[index])
            for component in range(len(component_time_axes))
            for index in range(count)
        )

    frequencies = _vector(
        loss_spectrum.get("frequency_hz"),
        "loss_spectrum.frequency_hz",
        minimum=3,
    )
    loss_count = len(frequencies)
    eddy = _matrix(
        loss_spectrum.get("eddy_components_w"),
        "loss_spectrum.eddy_components_w",
        rows=loss_count,
    )
    hysteresis = _matrix(
        loss_spectrum.get("hysteresis_components_w"),
        "loss_spectrum.hysteresis_components_w",
        rows=loss_count,
        columns=len(eddy[0]),
    )
    iron = _matrix(
        loss_spectrum.get("iron_components_w"),
        "loss_spectrum.iron_components_w",
        rows=loss_count,
        columns=len(eddy[0]),
    )
    reported_eddy = _vector(
        loss_spectrum.get("reported_eddy_total_w"),
        "loss_spectrum.reported_eddy_total_w",
        minimum=loss_count,
    )
    reported_hysteresis = _vector(
        loss_spectrum.get("reported_hysteresis_total_w"),
        "loss_spectrum.reported_hysteresis_total_w",
        minimum=loss_count,
    )
    reported_iron = _vector(
        loss_spectrum.get("reported_iron_total_w"),
        "loss_spectrum.reported_iron_total_w",
        minimum=loss_count,
    )
    if any(len(values) != loss_count for values in (reported_eddy, reported_hysteresis, reported_iron)):
        raise ValueError("loss total vectors must match the frequency axis")

    time_deltas = [right - left for left, right in zip(time_s, time_s[1:])]
    frequency_deltas = [right - left for left, right in zip(frequencies, frequencies[1:])]
    phase_scale = max(abs(value) for row in phase_currents for value in row)
    kcl_relative = max(abs(sum(row)) for row in phase_currents) / max(phase_scale, 1.0e-30)

    tail_start = max(1, int(0.75 * count))
    tracking_errors = []
    for component in range(2):
        command_rms = math.sqrt(
            sum(commands[index][component] ** 2 for index in range(tail_start, count))
            / (count - tail_start)
        )
        error_rms = math.sqrt(
            sum(
                (feedback[index][component] - commands[index][component]) ** 2
                for index in range(tail_start, count)
            )
            / (count - tail_start)
        )
        tracking_errors.append(error_rms / max(command_rms, 1.0e-30))

    integrated_angle = [angle_deg[0]]
    for index, delta_t in enumerate(time_deltas):
        integrated_angle.append(
            integrated_angle[-1]
            + 3.0 * (speed_rpm[index] + speed_rpm[index + 1]) * delta_t
        )
    angle_span = max(angle_deg) - min(angle_deg)
    angle_integral_relative = max(
        abs(actual - expected) for actual, expected in zip(angle_deg, integrated_angle)
    ) / max(angle_span, 1.0e-30)

    power_scale = max(abs(value) for value in total_power)
    power_sum_relative = max(
        abs(sum(components) - reported)
        for components, reported in zip(power_components, total_power)
    ) / max(power_scale, 1.0e-30)

    loss_scale = max(
        abs(value)
        for matrix in (eddy, hysteresis, iron)
        for row in matrix
        for value in row
    )
    loss_total_residual = max(
        abs(sum(matrix[index]) - reported[index])
        for matrix, reported in (
            (eddy, reported_eddy),
            (hysteresis, reported_hysteresis),
            (iron, reported_iron),
        )
        for index in range(loss_count)
    ) / max(loss_scale, 1.0e-30)
    eddy_reconstruction = max(
        [
            abs(eddy[0][component] - sum(row[component] for row in eddy[1:]))
            for component in range(len(eddy[0]))
        ]
        + [abs(reported_eddy[0] - sum(reported_eddy[1:]))]
    ) / max(loss_scale, 1.0e-30)
    aggregate_iron_decomposition = max(
        [
            abs(iron[0][component] - eddy[0][component] - hysteresis[0][component])
            for component in range(len(eddy[0]))
        ]
        + [abs(reported_iron[0] - reported_eddy[0] - reported_hysteresis[0])]
    ) / max(loss_scale, 1.0e-30)
    aggregate_only_bin_residual = max(
        [abs(value) for row in hysteresis[1:] for value in row]
        + [abs(value) for row in iron[1:] for value in row]
        + [abs(value) for value in reported_hysteresis[1:]]
        + [abs(value) for value in reported_iron[1:]]
    ) / max(loss_scale, 1.0e-30)
    frequency_step_span = (
        (max(frequency_deltas) - min(frequency_deltas))
        / max(abs(sum(frequency_deltas) / len(frequency_deltas)), 1.0e-30)
    )

    checks = {
        "time_axis_strictly_increases": all(delta > 0.0 for delta in time_deltas),
        "frequency_axis_has_zero_summary_then_uniform_positive_bins": frequencies[0] == 0.0
        and all(delta > 0.0 for delta in frequency_deltas)
        and frequency_step_span <= max_frequency_step_relative_span,
        "three_phase_current_satisfies_kcl": kcl_relative
        <= max_three_phase_kcl_relative_error,
        "tail_current_control_tracks_commands": max(tracking_errors)
        <= max_tail_control_tracking_rms_relative_error,
        "angle_matches_integrated_speed": angle_integral_relative
        <= max_angle_speed_integral_relative_error,
        "circuit_power_total_matches_component_sum": power_sum_relative
        <= max_power_sum_relative_error,
        "power_component_time_axes_match_common_axis_knotwise": (
            component_time_axes is None or component_time_alignment_error <= 1.0e-12
        ),
        "loss_total_columns_match_component_sums": loss_total_residual
        <= max_loss_identity_relative_error,
        "eddy_aggregate_matches_harmonic_bin_sum": eddy_reconstruction
        <= max_loss_identity_relative_error,
        "aggregate_iron_equals_eddy_plus_hysteresis": aggregate_iron_decomposition
        <= max_loss_identity_relative_error,
        "hysteresis_and_iron_are_aggregate_only": aggregate_only_bin_residual
        <= max_loss_identity_relative_error,
        "near_zero_current_ratio_outputs_are_diagnostic_only": payload.get(
            "ratio_diagnostic_policy"
        )
        == "exclude_ratio_when_denominator_current_is_below_floor",
        "torque_and_loss_share_periodic_cycle_generation": cycle_generation_ok,
        "restart_segments_preserve_phase_origin": restart_phase_origin_ok,
        "torque_and_dq_tables_share_transformed_angle_basis": angle_convention_ok,
        "loss_components_share_total_machine_scope": loss_normalization_ok,
        "abc_to_dq_phase_order_matches_winding_connection": dq_phase_order_ok,
        "torque_ripple_aggregation_excludes_repeated_cycle_endpoint": (
            torque_ripple_aggregation_ok
        ),
        "efficiency_input_and_output_share_periodic_average_window": (
            efficiency_average_window_ok
        ),
        "demag_margin_uses_current_temperature_material_state": (
            demag_temperature_material_state_ok
        ),
        "motor_torque_waveforms_share_mechanical_angle_basis": (
            torque_angle_basis_identity_ok
        ),
        "loss_harmonic_window_uses_one_rotor_position_generation": (
            loss_harmonic_rotor_window_identity_ok
        ),
        "iron_loss_coefficients_share_waveform_frequency_basis": (
            iron_loss_coefficient_frequency_basis_identity_ok
        ),
        "dq_currents_share_phase_order_and_q_axis_convention": (
            dq_current_phase_convention_identity_ok
        ),
        "torque_average_uses_one_electrical_or_mechanical_period_basis": (
            torque_average_period_angle_basis_identity_ok
        ),
        "lamination_loss_uses_consistent_stacking_factor_and_conductivity_basis": (
            lamination_stacking_factor_loss_conductivity_identity_ok
        ),
        "dq_power_uses_one_power_invariant_park_transform_scaling": (
            dq_park_transform_power_invariant_scaling_identity_ok
        ),
        "demag_operating_point_uses_current_recoil_temperature_state": (
            demag_recoil_temperature_operating_point_identity_ok
        ),
        "torque_ripple_fft_uses_current_electrical_angle_resampling": (
            torque_ripple_fft_angle_endpoint_generation_identity_ok
        ),
        "iron_loss_spatial_harmonics_use_current_mesh_volume_weights": (
            iron_loss_spatial_harmonic_mesh_volume_identity_ok
        ),
        "dq_torque_map_uses_current_park_angle_and_q_axis_sign": (
            dq_torque_map_park_transform_angle_sign_identity_ok
        ),
        "efficiency_map_powers_share_current_steady_cycle_window": (
            efficiency_map_power_averaging_window_identity_ok
        ),
    }
    tail_torque = torque_nm[tail_start:]
    tail_torque_mean = sum(tail_torque) / len(tail_torque)
    tail_torque_ripple = (max(tail_torque) - min(tail_torque)) / max(
        abs(tail_torque_mean), 1.0e-30
    )
    return {
        "policy": "pwm_controlled_motor_loss_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "warnings": [] if identity_present else ["artifact_identity_not_recorded"],
        "metrics": {
            "time_row_count": count,
            "loss_row_count": loss_count,
            "loss_component_count": len(eddy[0]),
            "three_phase_kcl_global_relative_error": kcl_relative,
            "tail_control_tracking_rms_relative_errors": tracking_errors,
            "angle_speed_integral_relative_error": angle_integral_relative,
            "circuit_power_sum_global_relative_error": power_sum_relative,
            "power_component_time_axis_maximum_absolute_error_s": component_time_alignment_error,
            "loss_total_column_relative_error": loss_total_residual,
            "eddy_harmonic_reconstruction_relative_error": eddy_reconstruction,
            "aggregate_iron_decomposition_relative_error": aggregate_iron_decomposition,
            "aggregate_only_bin_relative_residual": aggregate_only_bin_residual,
            "frequency_step_relative_span": frequency_step_span,
            "tail_torque_mean_nm": tail_torque_mean,
            "tail_torque_peak_to_peak_relative": tail_torque_ripple,
        },
        "lesson": (
            "Gate PWM motor results with three-phase KCL, tail current-command tracking, "
            "integrated speed/angle, and power-column closure. In harmonic-loss tables, "
            "bind every time-domain power/loss component to the common time axis knot by knot; "
            "matching only an integrated loss cannot prove sample identity. "
            "treat row zero as the aggregate summary: eddy aggregate reconstructs from "
            "positive-frequency bins, while hysteresis and combined iron loss can remain "
            "aggregate-only. Bind abc-to-dq input ordering to the winding connection and "
            "remove a repeated periodic endpoint before mean/ripple aggregation. Exclude R/L "
            "ratios when their denominator current is near zero."
        ),
    }
