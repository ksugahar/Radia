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


def _sample_index(value: Any) -> int:
    if isinstance(value, bool):
        raise ValueError("sample indices must be integers")
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError("sample indices must be integers") from exc
    if not math.isfinite(parsed) or not parsed.is_integer():
        raise ValueError("sample indices must be integers")
    return int(parsed)


def _loss_power_balance_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, dict):
        return False
    generation = str(value.get("power_balance_generation", "")).strip()
    try:
        torque = [float(item) for item in value.get("torque_nm", [])]
        result_torque = [
            float(item) for item in value.get("power_balance_torque_nm", [])
        ]
        speed = [float(item) for item in value.get("speed_rad_s", [])]
        result_speed = [
            float(item) for item in value.get("power_balance_speed_rad_s", [])
        ]
        output = [float(item) for item in value.get("mechanical_output_w", [])]
        result_output = [
            float(item)
            for item in value.get("power_balance_mechanical_output_w", [])
        ]
        time_window = [
            float(item) for item in value.get("time_average_window_s", [])
        ]
        result_time_window = [
            float(item) for item in value.get("loss_time_average_window_s", [])
        ]
        iron = [float(item) for item in value.get("iron_loss_w", [])]
        result_iron = [
            float(item) for item in value.get("power_balance_iron_loss_w", [])
        ]
        copper = [float(item) for item in value.get("copper_loss_w", [])]
        result_copper = [
            float(item) for item in value.get("power_balance_copper_loss_w", [])
        ]
        mechanical_loss = [
            float(item) for item in value.get("mechanical_loss_w", [])
        ]
        result_mechanical_loss = [
            float(item)
            for item in value.get("power_balance_mechanical_loss_w", [])
        ]
        electrical_input = [
            float(item) for item in value.get("electrical_input_w", [])
        ]
        result_electrical_input = [
            float(item)
            for item in value.get("power_balance_electrical_input_w", [])
        ]
        harmonic_window = [
            int(item) for item in value.get("harmonic_window_samples", [])
        ]
        result_harmonic_window = [
            int(item) for item in value.get("loss_harmonic_window_samples", [])
        ]
    except (TypeError, ValueError):
        return False
    digest = str(value.get("power_balance_sha256", "")).lower()
    count = len(torque)
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "torque_speed_power_balance_generation",
                "harmonic_window_power_balance_generation",
                "time_average_power_balance_generation",
                "iron_loss_power_balance_generation",
                "copper_loss_power_balance_generation",
                "mechanical_loss_power_balance_generation",
                "result_power_balance_generation",
            )
        )
        and count >= 2
        and all(
            len(items) == count
            for items in (
                speed,
                output,
                iron,
                copper,
                mechanical_loss,
                electrical_input,
            )
        )
        and all(
            math.isfinite(item)
            for items in (
                torque,
                speed,
                output,
                iron,
                copper,
                mechanical_loss,
                electrical_input,
            )
            for item in items
        )
        and result_torque == torque
        and result_speed == speed
        and result_output == output
        and result_iron == iron
        and result_copper == copper
        and result_mechanical_loss == mechanical_loss
        and result_electrical_input == electrical_input
        and all(
            math.isclose(
                output[index],
                torque[index] * speed[index],
                rel_tol=1.0e-12,
                abs_tol=1.0e-12,
            )
            for index in range(count)
        )
        and all(
            math.isclose(
                electrical_input[index],
                output[index]
                + iron[index]
                + copper[index]
                + mechanical_loss[index],
                rel_tol=1.0e-12,
                abs_tol=1.0e-12,
            )
            for index in range(count)
        )
        and len(harmonic_window) == 2
        and 0 <= harmonic_window[0] < harmonic_window[1]
        and result_harmonic_window == harmonic_window
        and len(time_window) == 2
        and all(math.isfinite(item) for item in time_window)
        and 0.0 <= time_window[0] < time_window[1]
        and result_time_window == time_window
        and _valid_sha256(digest)
        and value.get("reported_power_balance_sha256") == digest
    )


def _skew_slice_quadrature_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, dict):
        return False
    generation = str(value.get("skew_generation", "")).strip()
    try:
        slice_ids = [int(item) for item in value.get("slice_ids", [])]
        result_slice_ids = [int(item) for item in value.get("result_slice_ids", [])]
        weights = [float(item) for item in value.get("quadrature_weights", [])]
        result_weights = [
            float(item) for item in value.get("result_quadrature_weights", [])
        ]
        angles = [float(item) for item in value.get("rotor_angles_deg", [])]
        result_angles = [
            float(item) for item in value.get("result_rotor_angles_deg", [])
        ]
        torque = [float(item) for item in value.get("slice_torque_nm", [])]
        result_torque = [
            float(item) for item in value.get("result_slice_torque_nm", [])
        ]
        weighted_torque = float(value.get("weighted_torque_nm"))
        reported_weighted_torque = float(value.get("reported_weighted_torque_nm"))
    except (TypeError, ValueError):
        return False
    phases = [str(item).strip() for item in value.get("current_phase_ids", [])]
    result_phases = [
        str(item).strip() for item in value.get("result_current_phase_ids", [])
    ]
    periodic_maps = [
        str(item).strip() for item in value.get("periodic_map_ids", [])
    ]
    result_periodic_maps = [
        str(item).strip() for item in value.get("result_periodic_map_ids", [])
    ]
    solve_digests = [
        str(item).lower() for item in value.get("slice_solve_sha256", [])
    ]
    result_solve_digests = [
        str(item).lower() for item in value.get("result_slice_solve_sha256", [])
    ]
    digest = str(value.get("skew_result_sha256", "")).lower()
    count = len(slice_ids)
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "weight_skew_generation",
                "angle_skew_generation",
                "phase_skew_generation",
                "periodicity_skew_generation",
                "solve_skew_generation",
                "result_skew_generation",
            )
        )
        and count >= 2
        and all(item > 0 for item in slice_ids)
        and len(set(slice_ids)) == count
        and result_slice_ids == slice_ids
        and all(len(items) == count for items in (weights, angles, phases, periodic_maps, solve_digests, torque))
        and all(math.isfinite(item) and item >= 0.0 for item in weights)
        and math.isclose(sum(weights), 1.0, rel_tol=0.0, abs_tol=1.0e-12)
        and result_weights == weights
        and all(math.isfinite(item) for item in angles)
        and all(left < right for left, right in zip(angles, angles[1:]))
        and result_angles == angles
        and all(phases)
        and result_phases == phases
        and all(periodic_maps)
        and result_periodic_maps == periodic_maps
        and all(_valid_sha256(item) for item in solve_digests)
        and result_solve_digests == solve_digests
        and all(math.isfinite(item) for item in torque)
        and result_torque == torque
        and math.isfinite(weighted_torque)
        and math.isclose(
            weighted_torque,
            sum(weight * item for weight, item in zip(weights, torque)),
            rel_tol=1.0e-12,
            abs_tol=1.0e-12,
        )
        and math.isclose(
            reported_weighted_torque,
            weighted_torque,
            rel_tol=0.0,
            abs_tol=1.0e-12,
        )
        and _valid_sha256(digest)
        and value.get("reported_skew_result_sha256") == digest
    )


def _torque_map_interpolation_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, dict):
        return False
    generation = str(value.get("map_generation", "")).strip()
    method = str(value.get("interpolation_method", "")).strip()
    try:
        current = [float(item) for item in value.get("current_axis_a", [])]
        result_current = [float(item) for item in value.get("result_current_axis_a", [])]
        angle = [float(item) for item in value.get("electrical_angle_axis_deg", [])]
        result_angle = [float(item) for item in value.get("result_electrical_angle_axis_deg", [])]
        temperature = [float(item) for item in value.get("temperature_axis_c", [])]
        result_temperature = [float(item) for item in value.get("result_temperature_axis_c", [])]
        speed = [float(item) for item in value.get("speed_axis_rpm", [])]
        result_speed = [float(item) for item in value.get("result_speed_axis_rpm", [])]
        period = float(value.get("angle_period_deg"))
        result_period = float(value.get("result_angle_period_deg"))
        query = [float(item) for item in value.get("query_point", [])]
        result_query = [float(item) for item in value.get("result_query_point", [])]
        torque = float(value.get("interpolated_torque_nm"))
        result_torque = float(value.get("result_interpolated_torque_nm"))
    except (TypeError, ValueError):
        return False
    tensor_digest = str(value.get("torque_tensor_sha256", "")).lower()
    result_digest = str(value.get("result_sha256", "")).lower()
    axes = (current, angle, temperature, speed)
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "current_map_generation",
                "angle_map_generation",
                "temperature_map_generation",
                "speed_map_generation",
                "interpolation_map_generation",
                "query_map_generation",
                "result_map_generation",
            )
        )
        and all(len(axis) >= 2 for axis in axes)
        and all(all(math.isfinite(item) for item in axis) for axis in axes)
        and all(all(left < right for left, right in zip(axis, axis[1:])) for axis in axes)
        and result_current == current
        and result_angle == angle
        and result_temperature == temperature
        and result_speed == speed
        and math.isfinite(period)
        and period > 0.0
        and math.isclose(result_period, period, rel_tol=0.0, abs_tol=1.0e-12)
        and method == "multilinear_periodic_angle"
        and value.get("result_interpolation_method") == method
        and _valid_sha256(tensor_digest)
        and value.get("result_torque_tensor_sha256") == tensor_digest
        and len(query) == 4
        and all(math.isfinite(item) for item in query)
        and result_query == query
        and math.isfinite(torque)
        and math.isclose(result_torque, torque, rel_tol=0.0, abs_tol=1.0e-12)
        and _valid_sha256(result_digest)
        and value.get("accepted_result_sha256") == result_digest
    )


def _demagnetization_margin_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, dict):
        return False
    generation = str(value.get("demag_generation", "")).strip()
    magnet_ids = [str(item).strip() for item in value.get("magnet_ids", [])]
    result_magnet_ids = [str(item).strip() for item in value.get("result_magnet_ids", [])]
    try:
        temperature = float(value.get("temperature_c"))
        result_temperature = float(value.get("result_temperature_c"))
        coercivity = float(value.get("coercivity_a_m"))
        result_coercivity = float(value.get("result_coercivity_a_m"))
        recoil = float(value.get("recoil_relative_permeability"))
        result_recoil = float(value.get("result_recoil_relative_permeability"))
        operating_h = float(value.get("minimum_operating_h_a_m"))
        result_operating_h = float(value.get("result_minimum_operating_h_a_m"))
        margin = float(value.get("demagnetization_margin_a_m"))
        result_margin = float(value.get("result_demagnetization_margin_a_m"))
    except (TypeError, ValueError):
        return False
    curve_digest = str(value.get("material_curve_sha256", "")).lower()
    field_digest = str(value.get("operating_point_field_sha256", "")).lower()
    result_digest = str(value.get("result_sha256", "")).lower()
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "material_demag_generation",
                "temperature_demag_generation",
                "recoil_demag_generation",
                "operating_point_demag_generation",
                "margin_demag_generation",
                "result_demag_generation",
            )
        )
        and bool(magnet_ids)
        and all(magnet_ids)
        and len(set(magnet_ids)) == len(magnet_ids)
        and result_magnet_ids == magnet_ids
        and all(math.isfinite(item) for item in (temperature, coercivity, recoil, operating_h, margin))
        and coercivity > 0.0
        and recoil > 0.0
        and margin >= 0.0
        and math.isclose(result_temperature, temperature, rel_tol=0.0, abs_tol=1.0e-12)
        and math.isclose(result_coercivity, coercivity, rel_tol=0.0, abs_tol=1.0e-9)
        and math.isclose(result_recoil, recoil, rel_tol=0.0, abs_tol=1.0e-12)
        and math.isclose(result_operating_h, operating_h, rel_tol=0.0, abs_tol=1.0e-9)
        and math.isclose(result_margin, margin, rel_tol=0.0, abs_tol=1.0e-9)
        and _valid_sha256(curve_digest)
        and value.get("result_material_curve_sha256") == curve_digest
        and _valid_sha256(field_digest)
        and value.get("result_operating_point_field_sha256") == field_digest
        and _valid_sha256(result_digest)
        and value.get("accepted_result_sha256") == result_digest
    )


def _iron_loss_component_harmonic_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, dict):
        return False
    generation = str(value.get("loss_generation", "")).strip()
    try:
        hysteresis = float(value.get("hysteresis_loss_w"))
        result_hysteresis = float(value.get("result_hysteresis_loss_w"))
        eddy = float(value.get("eddy_loss_w"))
        result_eddy = float(value.get("result_eddy_loss_w"))
        excess = float(value.get("excess_loss_w"))
        result_excess = float(value.get("result_excess_loss_w"))
        total = float(value.get("total_iron_loss_w"))
        result_total = float(value.get("result_total_iron_loss_w"))
        orders = [int(item) for item in value.get("harmonic_orders", [])]
        result_orders = [int(item) for item in value.get("result_harmonic_orders", [])]
        frequencies = [float(item) for item in value.get("harmonic_frequencies_hz", [])]
        result_frequencies = [float(item) for item in value.get("result_harmonic_frequencies_hz", [])]
        volume = float(value.get("integration_volume_m3"))
        result_volume = float(value.get("result_integration_volume_m3"))
    except (TypeError, ValueError):
        return False
    return (
        bool(generation)
        and all(value.get(key) == generation for key in (
            "component_loss_generation", "harmonic_loss_generation", "frequency_loss_generation",
            "material_loss_generation", "volume_loss_generation", "result_loss_generation"))
        and all(math.isfinite(item) and item >= 0.0 for item in (hysteresis, eddy, excess, total))
        and result_hysteresis == hysteresis and result_eddy == eddy and result_excess == excess
        and math.isclose(total, hysteresis + eddy + excess, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(result_total, total, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and bool(orders) and all(item > 0 for item in orders) and len(set(orders)) == len(orders)
        and result_orders == orders and len(frequencies) == len(orders)
        and all(math.isfinite(item) and item > 0.0 for item in frequencies)
        and all(left < right for left, right in zip(frequencies, frequencies[1:]))
        and result_frequencies == frequencies
        and volume > 0.0 and result_volume == volume
        and _valid_sha256(value.get("material_law_sha256"))
        and value.get("result_material_law_sha256") == value.get("material_law_sha256")
        and _valid_sha256(value.get("mesh_sha256"))
        and value.get("result_mesh_sha256") == value.get("mesh_sha256")
        and _valid_sha256(value.get("result_sha256"))
        and value.get("accepted_result_sha256") == value.get("result_sha256")
    )


def _skew_slice_phase_angle_mesh_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, dict):
        return False
    generation = str(value.get("skew_generation", "")).strip()
    try:
        phases = [float(item) for item in value.get("slice_phase_deg", [])]
        result_phases = [float(item) for item in value.get("result_slice_phase_deg", [])]
        angles = [float(item) for item in value.get("mechanical_angle_deg", [])]
        result_angles = [float(item) for item in value.get("result_mechanical_angle_deg", [])]
        weights = [float(item) for item in value.get("slice_weights", [])]
        result_weights = [float(item) for item in value.get("result_slice_weights", [])]
        periodicity = int(value.get("periodicity"))
        result_periodicity = int(value.get("result_periodicity"))
        torque = [float(item) for item in value.get("slice_torque_nm", [])]
        result_torque = [float(item) for item in value.get("result_slice_torque_nm", [])]
        average = float(value.get("skew_averaged_torque_nm"))
        result_average = float(value.get("result_skew_averaged_torque_nm"))
    except (TypeError, ValueError):
        return False
    meshes = [str(item).lower() for item in value.get("slice_mesh_sha256", [])]
    result_meshes = [str(item).lower() for item in value.get("result_slice_mesh_sha256", [])]
    count = len(phases)
    return (
        bool(generation)
        and all(value.get(key) == generation for key in (
            "phase_skew_generation", "angle_skew_generation", "weight_skew_generation",
            "periodicity_skew_generation", "mesh_skew_generation", "result_skew_generation"))
        and count >= 2 and all(len(items) == count for items in (angles, weights, meshes, torque))
        and all(math.isfinite(item) for item in phases + angles + weights + torque)
        and result_phases == phases and result_angles == angles and result_weights == weights
        and all(item >= 0.0 for item in weights)
        and math.isclose(sum(weights), 1.0, rel_tol=0.0, abs_tol=1.0e-12)
        and periodicity > 0 and result_periodicity == periodicity
        and all(_valid_sha256(item) for item in meshes) and result_meshes == meshes
        and result_torque == torque and math.isfinite(average)
        and math.isclose(average, sum(weight * item for weight, item in zip(weights, torque)),
                         rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(result_average, average, rel_tol=0.0, abs_tol=1.0e-12)
        and _valid_sha256(value.get("result_sha256"))
        and value.get("accepted_result_sha256") == value.get("result_sha256")
    )


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
    skew_slice_torque_mechanical_phase_offset_generation_identity_ok = True
    winding_phase_belt_slot_numbering_sequence_generation_identity_ok = True
    iron_loss_harmonic_coefficient_unit_basis_identity_ok = True
    demag_temperature_phase_operating_point_identity_ok = True
    skew_slice_angle_weight_periodicity_identity_ok = True
    incremental_inductance_perturbation_phase_state_identity_ok = True
    dq_transform_rotor_angle_phase_order_generation_identity_ok = True
    iron_loss_frequency_harmonic_material_curve_generation_identity_ok = True
    motion_skew_force_harmonic_generation_identity_ok = True
    irreversible_demag_state_generation_identity_ok = True
    winding_current_torque_identity_ok = True
    demag_knee_operating_identity_ok = True
    loss_power_balance_identity_ok = True
    skew_slice_quadrature_identity_ok = True
    torque_map_interpolation_identity_ok = True
    demagnetization_margin_identity_ok = True
    iron_loss_component_harmonic_identity_ok = True
    skew_slice_phase_angle_mesh_identity_ok = True
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
        skew_slice_torque_mechanical_phase_offset_generation_identity_ok = False
        winding_phase_belt_slot_numbering_sequence_generation_identity_ok = False
        iron_loss_harmonic_coefficient_unit_basis_identity_ok = False
        demag_temperature_phase_operating_point_identity_ok = False
        skew_slice_angle_weight_periodicity_identity_ok = False
        incremental_inductance_perturbation_phase_state_identity_ok = False
        dq_transform_rotor_angle_phase_order_generation_identity_ok = False
        iron_loss_frequency_harmonic_material_curve_generation_identity_ok = False
        motion_skew_force_harmonic_generation_identity_ok = False
        irreversible_demag_state_generation_identity_ok = False
        winding_current_torque_identity_ok = False
        demag_knee_operating_identity_ok = False
        loss_power_balance_identity_ok = False
        skew_slice_quadrature_identity_ok = False
        torque_map_interpolation_identity_ok = False
        demagnetization_margin_identity_ok = False
        iron_loss_component_harmonic_identity_ok = False
        skew_slice_phase_angle_mesh_identity_ok = False
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
            canonical_park = identity_value.get(
                "dq_park_transform_power_invariant_scaling_identity"
            )
            canonical_dq = identity_value.get("dq_current_phase_convention_identity")
            if not isinstance(canonical_park, dict) or not isinstance(
                canonical_dq, dict
            ):
                canonical_park = {}
                canonical_dq = {}
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
                and type(q_axis_sign) is int
                and q_axis_sign in (-1, 1)
                and dq_torque_map.get("torque_map_q_axis_sign") == q_axis_sign
                and canonical_dq.get("source_q_axis_lead")
                in {
                    "q_leads_d_positive_electrical",
                    "q_lags_d_positive_electrical",
                }
                and dq_torque_map.get("q_axis_convention")
                == canonical_dq.get("source_q_axis_lead")
                and dq_torque_map.get("torque_map_q_axis_convention")
                == dq_torque_map.get("q_axis_convention")
                and q_axis_sign
                == (
                    1
                    if canonical_dq.get("source_q_axis_lead")
                    == "q_leads_d_positive_electrical"
                    else -1
                )
                and len(park_digest) == 64
                and all(character in "0123456789abcdef" for character in park_digest)
                and str(
                    dq_torque_map.get("torque_map_park_matrix_sha256", "")
                ).lower()
                == park_digest
                and park_digest
                == str(canonical_park.get("park_transform_sha256", "")).lower()
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
            reference_window = identity_value.get("efficiency_average_window")
            try:
                start = _sample_index(
                    efficiency_window.get("window_start_sample")
                )
                end = _sample_index(efficiency_window.get("window_end_sample"))
                electrical_window = [
                    _sample_index(value)
                    for value in efficiency_window.get("electrical_power_window", [])
                ]
                mechanical_window = [
                    _sample_index(value)
                    for value in efficiency_window.get("mechanical_power_window", [])
                ]
                loss_window = [
                    _sample_index(value)
                    for value in efficiency_window.get("loss_power_window", [])
                ]
                reference_start = _sample_index(
                    reference_window.get("input_window_start_sample")
                )
                reference_end = _sample_index(
                    reference_window.get("input_window_end_sample_exclusive")
                )
                sample_count = len(time_series.get("time_s", []))
            except (TypeError, ValueError):
                start = end = -1
                electrical_window = []
                mechanical_window = []
                loss_window = []
                reference_start = reference_end = -1
                sample_count = -1
            efficiency_map_power_averaging_window_identity_ok = (
                isinstance(reference_window, dict)
                and bool(map_generation)
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
                and end <= sample_count
                and [start, end] == [reference_start, reference_end]
                and electrical_window == [start, end]
                and mechanical_window == [start, end]
                and loss_window == [start, end]
                and len(window_digest) == 64
                and all(
                    character in "0123456789abcdef" for character in window_digest
                )
                and str(
                    efficiency_window.get("efficiency_power_window_sha256", "")
                ).lower()
                == window_digest
            )

        skew_slice_value = identity_value.get(
            "skew_slice_torque_mechanical_phase_offset_generation_identity"
        )
        if skew_slice_value is not None:
            skew_slice = skew_slice_value if isinstance(skew_slice_value, dict) else {}
            average_generation = str(
                skew_slice.get("torque_average_generation", "")
            )
            geometry_generation = str(
                skew_slice.get("slice_geometry_generation", "")
            )
            slice_ids = skew_slice.get("slice_ids")
            phase_slice_ids = skew_slice.get("phase_offset_slice_ids")
            offsets = skew_slice.get("mechanical_phase_offsets_deg")
            applied_offsets = skew_slice.get("applied_mechanical_phase_offsets_deg")
            waveform_digests = skew_slice.get("slice_torque_waveform_sha256")
            averaged_digests = skew_slice.get(
                "averaged_slice_torque_waveform_sha256"
            )
            phase_digest = str(
                skew_slice.get("phase_offset_map_sha256", "")
            ).lower()
            skew_slice_torque_mechanical_phase_offset_generation_identity_ok = (
                bool(average_generation)
                and skew_slice.get("slice_torque_average_generation")
                == average_generation
                and bool(geometry_generation)
                and skew_slice.get("phase_offset_slice_geometry_generation")
                == geometry_generation
                and skew_slice.get("torque_waveform_slice_geometry_generation")
                == geometry_generation
                and isinstance(slice_ids, list)
                and bool(slice_ids)
                and all(
                    isinstance(value, int) and not isinstance(value, bool)
                    for value in slice_ids
                )
                and len(set(slice_ids)) == len(slice_ids)
                and phase_slice_ids == slice_ids
                and isinstance(offsets, list)
                and len(offsets) == len(slice_ids)
                and all(
                    isinstance(value, (int, float))
                    and not isinstance(value, bool)
                    and math.isfinite(float(value))
                    for value in offsets
                )
                and applied_offsets == offsets
                and isinstance(waveform_digests, list)
                and len(waveform_digests) == len(slice_ids)
                and all(
                    len(str(value)) == 64
                    and all(
                        character in "0123456789abcdef"
                        for character in str(value).lower()
                    )
                    for value in waveform_digests
                )
                and averaged_digests == waveform_digests
                and len(phase_digest) == 64
                and all(
                    character in "0123456789abcdef" for character in phase_digest
                )
                and str(
                    skew_slice.get("applied_phase_offset_map_sha256", "")
                ).lower()
                == phase_digest
            )

        winding_value = identity_value.get(
            "winding_phase_belt_slot_numbering_sequence_generation_identity"
        )
        if winding_value is not None:
            winding = winding_value if isinstance(winding_value, dict) else {}
            winding_generation = str(winding.get("winding_generation", ""))
            numbering_generation = str(
                winding.get("slot_numbering_generation", "")
            )
            slot_numbers = winding.get("slot_numbers")
            phase_sequence = winding.get("phase_sequence")
            slot_map_digest = str(
                winding.get("slot_phase_map_sha256", "")
            ).lower()
            winding_phase_belt_slot_numbering_sequence_generation_identity_ok = (
                bool(winding_generation)
                and winding.get("phase_belt_winding_generation")
                == winding_generation
                and winding.get("mmf_harmonic_winding_generation")
                == winding_generation
                and bool(numbering_generation)
                and winding.get("phase_belt_slot_numbering_generation")
                == numbering_generation
                and winding.get("mmf_slot_numbering_generation")
                == numbering_generation
                and isinstance(slot_numbers, list)
                and bool(slot_numbers)
                and all(
                    isinstance(value, int)
                    and not isinstance(value, bool)
                    and value > 0
                    for value in slot_numbers
                )
                and len(set(slot_numbers)) == len(slot_numbers)
                and winding.get("phase_belt_slot_numbers") == slot_numbers
                and isinstance(phase_sequence, list)
                and len(phase_sequence) == len(slot_numbers)
                and all(
                    isinstance(value, str)
                    and len(value) == 2
                    and value[0] in "ABC"
                    and value[1] in "+-"
                    for value in phase_sequence
                )
                and winding.get("mmf_phase_sequence") == phase_sequence
                and len(slot_map_digest) == 64
                and all(
                    character in "0123456789abcdef"
                    for character in slot_map_digest
                )
                and str(winding.get("mmf_slot_phase_map_sha256", "")).lower()
                == slot_map_digest
            )

        harmonic_loss_value = identity_value.get(
            "iron_loss_harmonic_frequency_coefficient_unit_basis_identity"
        )
        if harmonic_loss_value is not None:
            harmonic_loss = (
                harmonic_loss_value if isinstance(harmonic_loss_value, dict) else {}
            )
            loss_generation = str(harmonic_loss.get("loss_generation", ""))
            digest = str(harmonic_loss.get("loss_basis_sha256", "")).lower()
            try:
                frequencies = [
                    float(value)
                    for value in harmonic_loss.get("harmonic_frequencies_hz", [])
                ]
                evaluated_frequencies = [
                    float(value)
                    for value in harmonic_loss.get(
                        "evaluated_harmonic_frequencies_hz", []
                    )
                ]
                coefficients = [
                    float(value)
                    for value in harmonic_loss.get("loss_coefficients", [])
                ]
                evaluated_coefficients = [
                    float(value)
                    for value in harmonic_loss.get(
                        "evaluated_loss_coefficients", []
                    )
                ]
            except (TypeError, ValueError):
                frequencies = evaluated_frequencies = []
                coefficients = evaluated_coefficients = []
            iron_loss_harmonic_coefficient_unit_basis_identity_ok = (
                bool(loss_generation)
                and harmonic_loss.get("harmonic_loss_generation") == loss_generation
                and harmonic_loss.get("coefficient_loss_generation") == loss_generation
                and harmonic_loss.get("frequency_unit") == "Hz"
                and harmonic_loss.get("coefficient_frequency_unit") == "Hz"
                and harmonic_loss.get("flux_density_unit") == "T"
                and harmonic_loss.get("coefficient_flux_density_unit") == "T"
                and bool(frequencies)
                and all(
                    math.isfinite(value) and value > 0.0 for value in frequencies
                )
                and all(
                    right > left for left, right in zip(frequencies, frequencies[1:])
                )
                and evaluated_frequencies == frequencies
                and len(coefficients) == len(evaluated_coefficients)
                and bool(coefficients)
                and all(math.isfinite(value) and value >= 0.0 for value in coefficients)
                and evaluated_coefficients == coefficients
                and len(digest) == 64
                and all(character in "0123456789abcdef" for character in digest)
                and str(
                    harmonic_loss.get("evaluated_loss_basis_sha256", "")
                ).lower()
                == digest
            )

        demag_operating_value = identity_value.get(
            "demagnetization_temperature_current_phase_operating_point_identity"
        )
        if demag_operating_value is not None:
            demag_operating = (
                demag_operating_value
                if isinstance(demag_operating_value, dict)
                else {}
            )
            generation = str(
                demag_operating.get("operating_point_generation", "")
            )
            digest = str(
                demag_operating.get("operating_point_sha256", "")
            ).lower()
            try:
                temperature = float(demag_operating.get("magnet_temperature_c"))
                margin_temperature = float(
                    demag_operating.get("demag_margin_temperature_c")
                )
                phase = float(demag_operating.get("current_phase_deg"))
                margin_phase = float(
                    demag_operating.get("demag_margin_current_phase_deg")
                )
            except (TypeError, ValueError):
                temperature = margin_temperature = phase = margin_phase = math.nan
            demag_temperature_phase_operating_point_identity_ok = (
                bool(generation)
                and demag_operating.get("temperature_operating_point_generation")
                == generation
                and demag_operating.get("current_phase_operating_point_generation")
                == generation
                and demag_operating.get("demag_margin_operating_point_generation")
                == generation
                and math.isfinite(temperature)
                and margin_temperature == temperature
                and math.isfinite(phase)
                and margin_phase == phase
                and len(digest) == 64
                and all(character in "0123456789abcdef" for character in digest)
                and str(
                    demag_operating.get("demag_margin_operating_point_sha256", "")
                ).lower()
                == digest
            )

        skew_average = identity_value.get(
            "skew_slice_torque_angle_weight_periodicity_generation_identity"
        )
        if skew_average is not None:
            skew_average = skew_average if isinstance(skew_average, dict) else {}
            generation = str(skew_average.get("skew_generation", "")).strip()
            try:
                slice_ids = [int(value) for value in skew_average.get("slice_ids", [])]
                torque_slice_ids = [
                    int(value) for value in skew_average.get("torque_slice_ids", [])
                ]
                angles = [
                    float(value) for value in skew_average.get("slice_angles_deg", [])
                ]
                torque_angles = [
                    float(value)
                    for value in skew_average.get("torque_slice_angles_deg", [])
                ]
                weights = [
                    float(value) for value in skew_average.get("quadrature_weights", [])
                ]
                torque_weights = [
                    float(value)
                    for value in skew_average.get("torque_quadrature_weights", [])
                ]
                wrap = float(skew_average.get("periodic_wrap_deg"))
                torque_wrap = float(skew_average.get("torque_periodic_wrap_deg"))
            except (TypeError, ValueError):
                slice_ids = torque_slice_ids = []
                angles = torque_angles = weights = torque_weights = []
                wrap = torque_wrap = math.nan
            digest = str(
                skew_average.get("skew_average_table_sha256", "")
            ).lower()
            skew_slice_angle_weight_periodicity_identity_ok = (
                bool(generation)
                and all(
                    skew_average.get(key) == generation
                    for key in (
                        "torque_skew_generation",
                        "angle_skew_generation",
                        "weight_skew_generation",
                        "periodicity_skew_generation",
                    )
                )
                and bool(slice_ids)
                and len(set(slice_ids)) == len(slice_ids)
                and torque_slice_ids == slice_ids
                and len(angles) == len(slice_ids)
                and all(math.isfinite(value) for value in angles)
                and torque_angles == angles
                and len(weights) == len(slice_ids)
                and all(math.isfinite(value) and value >= 0.0 for value in weights)
                and math.isclose(sum(weights), 1.0, rel_tol=1.0e-12, abs_tol=1.0e-12)
                and torque_weights == weights
                and math.isfinite(wrap)
                and wrap > 0.0
                and torque_wrap == wrap
                and len(digest) == 64
                and all(character in "0123456789abcdef" for character in digest)
                and str(
                    skew_average.get("torque_skew_average_table_sha256", "")
                ).lower()
                == digest
            )

        incremental = identity_value.get(
            "incremental_inductance_current_perturbation_phase_state_generation_identity"
        )
        if incremental is not None:
            incremental = incremental if isinstance(incremental, dict) else {}
            operating_generation = str(
                incremental.get("operating_point_generation", "")
            ).strip()
            base_solve = str(incremental.get("base_solve_generation", "")).strip()
            perturbation_solves = [
                str(value)
                for value in incremental.get("perturbation_solve_generations", [])
            ]
            matrix_solves = [
                str(value)
                for value in incremental.get(
                    "matrix_perturbation_solve_generations", []
                )
            ]
            phases = [str(value) for value in incremental.get("phase_names", [])]
            matrix_phases = [
                str(value) for value in incremental.get("matrix_phase_names", [])
            ]
            try:
                currents = [
                    [float(value) for value in row]
                    for row in incremental.get("perturbation_currents_a", [])
                ]
                matrix_currents = [
                    [float(value) for value in row]
                    for row in incremental.get(
                        "matrix_perturbation_currents_a", []
                    )
                ]
            except (TypeError, ValueError):
                currents = matrix_currents = []
            digest = str(
                incremental.get("incremental_inductance_table_sha256", "")
            ).lower()
            incremental_inductance_perturbation_phase_state_identity_ok = (
                bool(operating_generation)
                and all(
                    incremental.get(key) == operating_generation
                    for key in (
                        "matrix_operating_point_generation",
                        "perturbation_operating_point_generation",
                        "phase_state_operating_point_generation",
                    )
                )
                and bool(base_solve)
                and incremental.get("matrix_base_solve_generation") == base_solve
                and bool(perturbation_solves)
                and all(perturbation_solves)
                and len(set(perturbation_solves)) == len(perturbation_solves)
                and matrix_solves == perturbation_solves
                and bool(phases)
                and len(set(phases)) == len(phases)
                and all(phase.strip() for phase in phases)
                and len(phases) == len(perturbation_solves)
                and matrix_phases == phases
                and len(currents) == len(phases)
                and all(
                    len(row) == len(phases)
                    and all(math.isfinite(value) for value in row)
                    for row in currents
                )
                and matrix_currents == currents
                and len(digest) == 64
                and all(character in "0123456789abcdef" for character in digest)
                and str(
                    incremental.get(
                        "resolved_incremental_inductance_table_sha256", ""
                    )
                ).lower()
                == digest
            )

        dq_transform = identity_value.get(
            "dq_transform_rotor_angle_phase_order_generation_identity"
        )
        if dq_transform is not None:
            dq_transform = dq_transform if isinstance(dq_transform, dict) else {}
            generation = str(
                dq_transform.get("operating_point_generation", "")
            ).strip()
            phase_order_value = dq_transform.get("phase_order")
            dq_phase_order_value = dq_transform.get("dq_phase_order")
            phase_order = (
                [str(value) for value in phase_order_value]
                if isinstance(phase_order_value, list)
                else []
            )
            dq_phase_order = (
                [str(value) for value in dq_phase_order_value]
                if isinstance(dq_phase_order_value, list)
                else []
            )
            try:
                rotor_angle = float(
                    dq_transform.get("rotor_mechanical_angle_deg")
                )
                dq_rotor_angle = float(
                    dq_transform.get("dq_rotor_mechanical_angle_deg")
                )
                electrical_offset = float(
                    dq_transform.get("electrical_offset_deg")
                )
                dq_electrical_offset = float(
                    dq_transform.get("dq_electrical_offset_deg")
                )
                phase_values = [
                    float(value) for value in dq_transform.get("phase_values", [])
                ]
                dq_phase_values = [
                    float(value)
                    for value in dq_transform.get("dq_source_phase_values", [])
                ]
            except (TypeError, ValueError):
                rotor_angle = dq_rotor_angle = math.nan
                electrical_offset = dq_electrical_offset = math.nan
                phase_values = dq_phase_values = []
            pole_pairs = dq_transform.get("pole_pairs")
            digest = str(dq_transform.get("dq_transform_table_sha256", "")).lower()
            dq_transform_rotor_angle_phase_order_generation_identity_ok = (
                bool(generation)
                and all(
                    dq_transform.get(key) == generation
                    for key in (
                        "rotor_angle_operating_point_generation",
                        "electrical_offset_operating_point_generation",
                        "phase_order_operating_point_generation",
                        "dq_result_operating_point_generation",
                    )
                )
                and math.isfinite(rotor_angle)
                and dq_rotor_angle == rotor_angle
                and isinstance(pole_pairs, int)
                and not isinstance(pole_pairs, bool)
                and pole_pairs > 0
                and dq_transform.get("dq_pole_pairs") == pole_pairs
                and math.isfinite(electrical_offset)
                and dq_electrical_offset == electrical_offset
                and bool(phase_order)
                and all(value.strip() for value in phase_order)
                and len(set(phase_order)) == len(phase_order)
                and dq_phase_order == phase_order
                and len(phase_values) == len(phase_order)
                and all(math.isfinite(value) for value in phase_values)
                and dq_phase_values == phase_values
                and len(digest) == 64
                and all(character in "0123456789abcdef" for character in digest)
                and str(
                    dq_transform.get("resolved_dq_transform_table_sha256", "")
                ).lower()
                == digest
            )

        iron_loss_inputs = identity_value.get(
            "iron_loss_frequency_harmonic_material_curve_generation_identity"
        )
        if iron_loss_inputs is not None:
            iron_loss_inputs = (
                iron_loss_inputs if isinstance(iron_loss_inputs, dict) else {}
            )
            generation = str(
                iron_loss_inputs.get("loss_study_generation", "")
            ).strip()
            curve_ids_value = iron_loss_inputs.get("material_curve_ids")
            loss_curve_ids_value = iron_loss_inputs.get("loss_material_curve_ids")
            curve_ids = (
                [str(value) for value in curve_ids_value]
                if isinstance(curve_ids_value, list)
                else []
            )
            loss_curve_ids = (
                [str(value) for value in loss_curve_ids_value]
                if isinstance(loss_curve_ids_value, list)
                else []
            )
            try:
                frequency_hz = float(
                    iron_loss_inputs.get("fundamental_frequency_hz")
                )
                loss_frequency_hz = float(
                    iron_loss_inputs.get("loss_frequency_hz")
                )
                harmonic_orders = [
                    int(value)
                    for value in iron_loss_inputs.get("harmonic_orders", [])
                ]
                loss_harmonic_orders = [
                    int(value)
                    for value in iron_loss_inputs.get("loss_harmonic_orders", [])
                ]
                amplitudes = [
                    float(value)
                    for value in iron_loss_inputs.get("harmonic_amplitudes_t", [])
                ]
                loss_amplitudes = [
                    float(value)
                    for value in iron_loss_inputs.get(
                        "loss_harmonic_amplitudes_t", []
                    )
                ]
            except (TypeError, ValueError):
                frequency_hz = loss_frequency_hz = math.nan
                harmonic_orders = loss_harmonic_orders = []
                amplitudes = loss_amplitudes = []
            digest = str(
                iron_loss_inputs.get("loss_input_table_sha256", "")
            ).lower()
            iron_loss_frequency_harmonic_material_curve_generation_identity_ok = (
                bool(generation)
                and all(
                    iron_loss_inputs.get(key) == generation
                    for key in (
                        "frequency_loss_study_generation",
                        "harmonic_spectrum_loss_study_generation",
                        "material_curve_loss_study_generation",
                        "loss_result_study_generation",
                    )
                )
                and math.isfinite(frequency_hz)
                and frequency_hz > 0.0
                and loss_frequency_hz == frequency_hz
                and bool(harmonic_orders)
                and all(value > 0 for value in harmonic_orders)
                and harmonic_orders == sorted(harmonic_orders)
                and len(set(harmonic_orders)) == len(harmonic_orders)
                and loss_harmonic_orders == harmonic_orders
                and len(amplitudes) == len(harmonic_orders)
                and all(math.isfinite(value) and value >= 0.0 for value in amplitudes)
                and loss_amplitudes == amplitudes
                and bool(curve_ids)
                and all(value.strip() for value in curve_ids)
                and len(set(curve_ids)) == len(curve_ids)
                and loss_curve_ids == curve_ids
                and len(digest) == 64
                and all(character in "0123456789abcdef" for character in digest)
                and str(
                    iron_loss_inputs.get("resolved_loss_input_table_sha256", "")
                ).lower()
                == digest
            )

        motion_force = identity_value.get(
            "motion_skew_force_harmonic_time_angle_phase_generation_identity"
        )
        if motion_force is not None:
            motion_force = motion_force if isinstance(motion_force, dict) else {}
            generation = str(
                motion_force.get("motion_study_generation", "")
            ).strip()
            try:
                time_values = [float(value) for value in motion_force.get("time_s", [])]
                force_time_values = [
                    float(value) for value in motion_force.get("force_time_s", [])
                ]
                angle_values = [
                    float(value)
                    for value in motion_force.get("mechanical_angle_deg", [])
                ]
                force_angle_values = [
                    float(value)
                    for value in motion_force.get("force_mechanical_angle_deg", [])
                ]
                skew_angles = [
                    float(value)
                    for value in motion_force.get("skew_slice_angles_deg", [])
                ]
                force_skew_angles = [
                    float(value)
                    for value in motion_force.get("force_skew_slice_angles_deg", [])
                ]
                weights = [
                    float(value) for value in motion_force.get("slice_weights", [])
                ]
                force_weights = [
                    float(value)
                    for value in motion_force.get("force_slice_weights", [])
                ]
                phase_reference = float(motion_force.get("phase_reference_deg"))
                force_phase_reference = float(
                    motion_force.get("force_phase_reference_deg")
                )
                harmonic_orders = [
                    int(value) for value in motion_force.get("harmonic_orders", [])
                ]
                force_harmonic_orders = [
                    int(value)
                    for value in motion_force.get("force_harmonic_orders", [])
                ]
                force_harmonics = [
                    float(value)
                    for value in motion_force.get("force_harmonics_n", [])
                ]
                reported_force_harmonics = [
                    float(value)
                    for value in motion_force.get("reported_force_harmonics_n", [])
                ]
            except (TypeError, ValueError):
                time_values = force_time_values = []
                angle_values = force_angle_values = []
                skew_angles = force_skew_angles = []
                weights = force_weights = []
                phase_reference = force_phase_reference = math.nan
                harmonic_orders = force_harmonic_orders = []
                force_harmonics = reported_force_harmonics = []
            digest = str(
                motion_force.get("force_harmonic_table_sha256", "")
            ).lower()
            motion_skew_force_harmonic_generation_identity_ok = (
                bool(generation)
                and all(
                    motion_force.get(key) == generation
                    for key in (
                        "time_motion_study_generation",
                        "angle_motion_study_generation",
                        "skew_motion_study_generation",
                        "phase_motion_study_generation",
                        "force_result_motion_study_generation",
                    )
                )
                and len(time_values) >= 2
                and len(angle_values) == len(time_values)
                and all(math.isfinite(value) for value in time_values + angle_values)
                and all(
                    right > left for left, right in zip(time_values, time_values[1:])
                )
                and force_time_values == time_values
                and force_angle_values == angle_values
                and bool(skew_angles)
                and len(weights) == len(skew_angles)
                and all(math.isfinite(value) for value in skew_angles + weights)
                and all(value >= 0.0 for value in weights)
                and math.isclose(sum(weights), 1.0, rel_tol=0.0, abs_tol=1.0e-12)
                and force_skew_angles == skew_angles
                and force_weights == weights
                and math.isfinite(phase_reference)
                and force_phase_reference == phase_reference
                and bool(harmonic_orders)
                and all(value > 0 for value in harmonic_orders)
                and harmonic_orders == sorted(set(harmonic_orders))
                and force_harmonic_orders == harmonic_orders
                and len(force_harmonics) == len(harmonic_orders)
                and all(math.isfinite(value) for value in force_harmonics)
                and reported_force_harmonics == force_harmonics
                and len(digest) == 64
                and all(character in "0123456789abcdef" for character in digest)
                and str(
                    motion_force.get("resolved_force_harmonic_table_sha256", "")
                ).lower()
                == digest
            )

        demag_state = identity_value.get(
            "ipm_irreversible_demag_recoil_temperature_operating_generation_identity"
        )
        if demag_state is not None:
            demag_state = demag_state if isinstance(demag_state, dict) else {}
            generation = str(
                demag_state.get("demag_study_generation", "")
            ).strip()
            operating_point = str(
                demag_state.get("operating_point_id", "")
            ).strip()
            try:
                temperature = float(demag_state.get("temperature_c"))
                result_temperature = float(demag_state.get("result_temperature_c"))
                orientations = [
                    [float(component) for component in vector]
                    for vector in demag_state.get("magnet_orientation_vectors", [])
                ]
                result_orientations = [
                    [float(component) for component in vector]
                    for vector in demag_state.get(
                        "result_magnet_orientation_vectors", []
                    )
                ]
                margins = [
                    float(value)
                    for value in demag_state.get("demag_margin_a_per_m", [])
                ]
                result_margins = [
                    float(value)
                    for value in demag_state.get(
                        "reported_demag_margin_a_per_m", []
                    )
                ]
            except (TypeError, ValueError):
                temperature = result_temperature = math.nan
                orientations = result_orientations = []
                margins = result_margins = []
            recoil_digest = str(
                demag_state.get("recoil_curve_sha256", "")
            ).lower()
            state_digest = str(demag_state.get("magnet_state_sha256", "")).lower()
            irreversible_demag_state_generation_identity_ok = (
                bool(generation)
                and all(
                    demag_state.get(key) == generation
                    for key in (
                        "recoil_curve_demag_study_generation",
                        "temperature_demag_study_generation",
                        "operating_point_demag_study_generation",
                        "magnet_orientation_demag_study_generation",
                        "result_demag_study_generation",
                    )
                )
                and math.isfinite(temperature)
                and result_temperature == temperature
                and bool(operating_point)
                and demag_state.get("result_operating_point_id") == operating_point
                and bool(orientations)
                and all(
                    len(vector) in {2, 3}
                    and all(math.isfinite(component) for component in vector)
                    and any(component != 0.0 for component in vector)
                    for vector in orientations
                )
                and result_orientations == orientations
                and len(recoil_digest) == 64
                and all(
                    character in "0123456789abcdef" for character in recoil_digest
                )
                and str(demag_state.get("result_recoil_curve_sha256", "")).lower()
                == recoil_digest
                and len(state_digest) == 64
                and all(
                    character in "0123456789abcdef" for character in state_digest
                )
                and str(demag_state.get("result_magnet_state_sha256", "")).lower()
                == state_digest
                and bool(margins)
                and all(math.isfinite(value) for value in margins)
                and result_margins == margins
            )

        winding_current_torque_identity_ok = _winding_current_torque_identity_ok(
            identity_value.get(
                "winding_current_phase_circuit_sequence_torque_generation_identity"
            )
        )
        demag_knee_operating_identity_ok = _demag_knee_operating_identity_ok(
            identity_value.get(
                "demagnetization_knee_temperature_recoil_operating_generation_identity"
            )
        )
        loss_power_balance_identity_ok = _loss_power_balance_identity_ok(
            identity_value.get(
                "loss_torque_speed_power_balance_harmonic_window_generation_identity"
            )
        )
        skew_slice_quadrature_identity_ok = _skew_slice_quadrature_identity_ok(
            identity_value.get(
                "skew_slice_weight_rotor_angle_phase_periodicity_generation_identity"
            )
        )
        torque_map_interpolation_identity_ok = _torque_map_interpolation_identity_ok(
            identity_value.get(
                "torque_map_current_angle_temperature_speed_interpolation_generation_identity"
            )
        )
        demagnetization_margin_identity_ok = _demagnetization_margin_identity_ok(
            identity_value.get(
                "demagnetization_margin_operating_point_temperature_recoil_generation_identity"
            )
        )
        iron_loss_component_harmonic_identity_ok = _iron_loss_component_harmonic_identity_ok(
            identity_value.get(
                "iron_loss_hysteresis_eddy_excess_harmonic_frequency_material_volume_generation_identity"
            )
        )
        skew_slice_phase_angle_mesh_identity_ok = _skew_slice_phase_angle_mesh_identity_ok(
            identity_value.get(
                "skew_slice_torque_phase_angle_weight_periodicity_mesh_generation_identity"
            )
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
        "skew_torque_average_uses_current_slice_phase_offset_map": (
            skew_slice_torque_mechanical_phase_offset_generation_identity_ok
        ),
        "winding_phase_belt_uses_current_slot_numbering_and_sequence": (
            winding_phase_belt_slot_numbering_sequence_generation_identity_ok
        ),
        "iron_loss_harmonics_and_coefficients_share_frequency_flux_units": (
            iron_loss_harmonic_coefficient_unit_basis_identity_ok
        ),
        "demag_margin_uses_current_temperature_and_current_phase_state": (
            demag_temperature_phase_operating_point_identity_ok
        ),
        "skew_torque_uses_current_slice_angles_weights_and_periodicity": (
            skew_slice_angle_weight_periodicity_identity_ok
        ),
        "incremental_inductance_uses_current_perturbation_phase_and_state": (
            incremental_inductance_perturbation_phase_state_identity_ok
        ),
        "dq_transform_uses_current_rotor_angle_offset_and_phase_order": (
            dq_transform_rotor_angle_phase_order_generation_identity_ok
        ),
        "iron_loss_uses_current_frequency_harmonics_and_material_curves": (
            iron_loss_frequency_harmonic_material_curve_generation_identity_ok
        ),
        "force_harmonics_use_current_motion_skew_time_angle_and_phase": (
            motion_skew_force_harmonic_generation_identity_ok
        ),
        "irreversible_demag_uses_current_recoil_temperature_operating_state": (
            irreversible_demag_state_generation_identity_ok
        ),
        "motor_torque_uses_current_winding_phase_circuit_sequence_and_angles": (
            winding_current_torque_identity_ok
        ),
        "demag_margin_uses_current_knee_temperature_recoil_and_operating_state": (
            demag_knee_operating_identity_ok
        ),
        "loss_power_balance_uses_current_torque_speed_windows_and_loss_components": (
            loss_power_balance_identity_ok
        ),
        "skew_slice_result_uses_current_weights_angles_phases_and_periodicity": (
            skew_slice_quadrature_identity_ok
        ),
        "torque_map_uses_current_axes_interpolation_query_and_result_generation": (
            torque_map_interpolation_identity_ok
        ),
        "demagnetization_margin_uses_current_temperature_recoil_material_and_operating_point": (
            demagnetization_margin_identity_ok
        ),
        "iron_loss_uses_current_components_harmonics_frequency_material_volume_and_result": (
            iron_loss_component_harmonic_identity_ok
        ),
        "skew_torque_uses_current_slice_phases_angles_weights_periodicity_meshes_and_result": (
            skew_slice_phase_angle_mesh_identity_ok
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
