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


def _rotating_sector_torque_frame_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, dict):
        return False
    generation = str(value.get("sector_generation", "")).strip()
    try:
        pole_pairs = int(value.get("pole_pairs"))
        result_pole_pairs = int(value.get("result_pole_pairs"))
        sector_angle = float(value.get("sector_angle_deg"))
        result_sector_angle = float(value.get("result_sector_angle_deg"))
        periodic_phase = float(value.get("periodic_phase_deg"))
        result_periodic_phase = float(value.get("result_periodic_phase_deg"))
        pairs = [[int(item) for item in pair] for pair in value.get("periodic_pair_ids", [])]
        result_pairs = [
            [int(item) for item in pair]
            for pair in value.get("result_periodic_pair_ids", [])
        ]
        orientations = [int(item) for item in value.get("periodic_pair_orientation", [])]
        result_orientations = [
            int(item) for item in value.get("result_periodic_pair_orientation", [])
        ]
        skew = [float(item) for item in value.get("skew_slice_deg", [])]
        result_skew = [float(item) for item in value.get("result_skew_slice_deg", [])]
        weights = [float(item) for item in value.get("skew_slice_weights", [])]
        result_weights = [
            float(item) for item in value.get("result_skew_slice_weights", [])
        ]
        angles = [float(item) for item in value.get("rotor_mechanical_angle_deg", [])]
        result_angles = [
            float(item) for item in value.get("result_rotor_mechanical_angle_deg", [])
        ]
        torque = [float(item) for item in value.get("slice_torque_nm", [])]
        result_torque = [float(item) for item in value.get("result_slice_torque_nm", [])]
        average = float(value.get("torque_average_nm"))
        result_average = float(value.get("result_torque_average_nm"))
    except (TypeError, ValueError):
        return False
    count = len(skew)
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "pole_pair_sector_generation",
                "periodic_sector_generation",
                "skew_sector_generation",
                "rotor_frame_sector_generation",
                "torque_sector_generation",
                "mesh_sector_generation",
                "result_sector_generation",
            )
        )
        and pole_pairs > 0
        and result_pole_pairs == pole_pairs
        and math.isfinite(sector_angle)
        and math.isclose(
            sector_angle, 360.0 / (2.0 * pole_pairs), rel_tol=0.0, abs_tol=1.0e-12
        )
        and result_sector_angle == sector_angle
        and periodic_phase in {0.0, 180.0}
        and result_periodic_phase == periodic_phase
        and bool(pairs)
        and all(len(pair) == 2 and pair[0] > 0 and pair[1] > 0 for pair in pairs)
        and len({tuple(pair) for pair in pairs}) == len(pairs)
        and result_pairs == pairs
        and len(orientations) == len(pairs)
        and all(item in {-1, 1} for item in orientations)
        and result_orientations == orientations
        and count >= 2
        and len(weights) == len(torque) == count
        and all(math.isfinite(item) for item in skew + weights + torque)
        and all(item >= 0.0 for item in weights)
        and math.isclose(sum(weights), 1.0, rel_tol=0.0, abs_tol=1.0e-12)
        and result_skew == skew
        and result_weights == weights
        and bool(angles)
        and all(math.isfinite(item) for item in angles)
        and result_angles == angles
        and value.get("torque_frame") == "rotor-mechanical-ccw"
        and value.get("result_torque_frame") == "rotor-mechanical-ccw"
        and result_torque == torque
        and math.isclose(
            average,
            sum(weight * item for weight, item in zip(weights, torque)),
            rel_tol=1.0e-12,
            abs_tol=1.0e-12,
        )
        and math.isclose(result_average, average, rel_tol=0.0, abs_tol=1.0e-12)
        and _valid_sha256(value.get("sector_mesh_sha256"))
        and value.get("result_sector_mesh_sha256") == value.get("sector_mesh_sha256")
        and _valid_sha256(value.get("torque_result_sha256"))
        and value.get("accepted_torque_result_sha256") == value.get("torque_result_sha256")
    )


def _iron_loss_decomposition_element_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, dict):
        return False
    generation = str(value.get("decomposition_generation", "")).strip()
    try:
        temperature = float(value.get("material_temperature_c"))
        result_temperature = float(value.get("result_material_temperature_c"))
        orders = [int(item) for item in value.get("harmonic_orders", [])]
        result_orders = [int(item) for item in value.get("result_harmonic_orders", [])]
        frequencies = [float(item) for item in value.get("frequency_hz", [])]
        result_frequencies = [float(item) for item in value.get("result_frequency_hz", [])]
        element_ids = [int(item) for item in value.get("element_ids", [])]
        result_element_ids = [int(item) for item in value.get("result_element_ids", [])]
        volumes = [float(item) for item in value.get("element_volume_m3", [])]
        result_volumes = [float(item) for item in value.get("result_element_volume_m3", [])]
        volume = float(value.get("integration_volume_m3"))
        result_volume = float(value.get("result_integration_volume_m3"))
    except (TypeError, ValueError):
        return False
    components = value.get("harmonic_loss_w")
    result_components = value.get("result_harmonic_loss_w")
    component_rows_ok = isinstance(components, dict) and set(components) == {
        "hysteresis",
        "eddy",
        "excess",
    }
    if component_rows_ok:
        try:
            component_rows = {
                key: [float(item) for item in components[key]] for key in components
            }
            result_component_rows = {
                key: [float(item) for item in result_components[key]]
                for key in result_components
            }
        except (TypeError, ValueError):
            component_rows_ok = False
            component_rows = result_component_rows = {}
    else:
        component_rows = result_component_rows = {}
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "model_decomposition_generation",
                "temperature_decomposition_generation",
                "frequency_decomposition_generation",
                "volume_decomposition_generation",
                "material_decomposition_generation",
                "result_decomposition_generation",
            )
        )
        and value.get("loss_model") == "bertotti-three-term"
        and value.get("result_loss_model") == "bertotti-three-term"
        and math.isfinite(temperature)
        and result_temperature == temperature
        and bool(orders)
        and all(item > 0 for item in orders)
        and len(set(orders)) == len(orders)
        and result_orders == orders
        and len(frequencies) == len(orders)
        and all(math.isfinite(item) and item > 0.0 for item in frequencies)
        and all(
            math.isclose(frequency, frequencies[0] * order / orders[0], rel_tol=1.0e-12)
            for order, frequency in zip(orders, frequencies)
        )
        and result_frequencies == frequencies
        and component_rows_ok
        and all(
            len(row) == len(orders)
            and all(math.isfinite(item) and item >= 0.0 for item in row)
            for row in component_rows.values()
        )
        and result_component_rows == component_rows
        and bool(element_ids)
        and all(item > 0 for item in element_ids)
        and len(set(element_ids)) == len(element_ids)
        and result_element_ids == element_ids
        and len(volumes) == len(element_ids)
        and all(math.isfinite(item) and item > 0.0 for item in volumes)
        and result_volumes == volumes
        and math.isclose(volume, sum(volumes), rel_tol=1.0e-12, abs_tol=1.0e-15)
        and result_volume == volume
        and _valid_sha256(value.get("material_state_sha256"))
        and value.get("result_material_state_sha256") == value.get("material_state_sha256")
        and _valid_sha256(value.get("mesh_sha256"))
        and value.get("result_mesh_sha256") == value.get("mesh_sha256")
        and _valid_sha256(value.get("loss_result_sha256"))
        and value.get("accepted_loss_result_sha256") == value.get("loss_result_sha256")
    )


def _pwm_observable_alignment_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, dict):
        return False
    generation = str(value.get("pwm_generation", "")).strip()
    try:
        orders = [int(item) for item in value.get("harmonic_orders", [])]
        result_orders = [int(item) for item in value.get("result_harmonic_orders", [])]
        currents = [float(item) for item in value.get("current_harmonic_a", [])]
        result_currents = [
            float(item) for item in value.get("result_current_harmonic_a", [])
        ]
        phases = [float(item) for item in value.get("current_phase_deg", [])]
        result_phases = [
            float(item) for item in value.get("result_current_phase_deg", [])
        ]
        times = [float(item) for item in value.get("time_s", [])]
        result_times = [float(item) for item in value.get("result_time_s", [])]
        angles = [float(item) for item in value.get("electrical_angle_deg", [])]
        result_angles = [
            float(item) for item in value.get("result_electrical_angle_deg", [])
        ]
        pole_pairs = int(value.get("pole_pairs"))
        result_pole_pairs = int(value.get("result_pole_pairs"))
        window = [float(item) for item in value.get("torque_window_s", [])]
        result_window = [
            float(item) for item in value.get("result_torque_window_s", [])
        ]
        torque_average = float(value.get("torque_average_nm"))
        result_torque_average = float(value.get("result_torque_average_nm"))
        loss_average = float(value.get("loss_average_w"))
        result_loss_average = float(value.get("result_loss_average_w"))
    except (TypeError, ValueError):
        return False
    harmonic_count = len(orders)
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "current_pwm_generation",
                "time_pwm_generation",
                "angle_pwm_generation",
                "torque_pwm_generation",
                "loss_pwm_generation",
                "mesh_pwm_generation",
                "result_pwm_generation",
            )
        )
        and harmonic_count > 0
        and len(set(orders)) == harmonic_count
        and all(item > 0 for item in orders)
        and len(currents) == len(phases) == harmonic_count
        and all(math.isfinite(item) and item >= 0.0 for item in currents)
        and all(math.isfinite(item) for item in phases)
        and result_orders == orders
        and result_currents == currents
        and result_phases == phases
        and len(times) >= 2
        and len(angles) == len(times)
        and all(math.isfinite(item) for item in times + angles)
        and all(right > left for left, right in zip(times, times[1:]))
        and result_times == times
        and result_angles == angles
        and pole_pairs > 0
        and result_pole_pairs == pole_pairs
        and len(window) == 2
        and math.isclose(window[0], times[0], rel_tol=0.0, abs_tol=1.0e-15)
        and math.isclose(window[1], times[-1], rel_tol=0.0, abs_tol=1.0e-15)
        and result_window == window
        and math.isfinite(torque_average)
        and math.isclose(
            result_torque_average, torque_average, rel_tol=0.0, abs_tol=1.0e-12
        )
        and math.isfinite(loss_average)
        and loss_average >= 0.0
        and math.isclose(
            result_loss_average, loss_average, rel_tol=0.0, abs_tol=1.0e-12
        )
        and _valid_sha256(value.get("mesh_sha256"))
        and value.get("result_mesh_sha256") == value.get("mesh_sha256")
        and _valid_sha256(value.get("result_sha256"))
        and value.get("accepted_result_sha256") == value.get("result_sha256")
    )


def _skew_slice_average_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, dict):
        return False
    generation = str(value.get("skew_generation", "")).strip()
    try:
        angles = [float(item) for item in value.get("slice_angles_deg", [])]
        result_angles = [
            float(item) for item in value.get("result_slice_angles_deg", [])
        ]
        weights = [float(item) for item in value.get("slice_weights", [])]
        result_weights = [float(item) for item in value.get("result_slice_weights", [])]
        torque = [float(item) for item in value.get("slice_torque_nm", [])]
        result_torque = [float(item) for item in value.get("result_slice_torque_nm", [])]
        average = float(value.get("torque_average_nm"))
        result_average = float(value.get("result_torque_average_nm"))
        ripple = float(value.get("torque_ripple_nm"))
        result_ripple = float(value.get("result_torque_ripple_nm"))
    except (TypeError, ValueError):
        return False
    count = len(angles)
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "angle_skew_generation",
                "weight_skew_generation",
                "frame_skew_generation",
                "interpolation_skew_generation",
                "torque_skew_generation",
                "mesh_skew_generation",
                "result_skew_generation",
            )
        )
        and count >= 3
        and len(weights) == len(torque) == count
        and all(math.isfinite(item) for item in angles + weights + torque)
        and all(item > 0.0 for item in weights)
        and math.isclose(sum(weights), 1.0, rel_tol=0.0, abs_tol=1.0e-12)
        and result_angles == angles
        and result_weights == weights
        and value.get("rotor_frame") == "mechanical-ccw"
        and value.get("result_rotor_frame") == "mechanical-ccw"
        and value.get("interpolation_rule") == "periodic-cubic"
        and value.get("result_interpolation_rule") == "periodic-cubic"
        and result_torque == torque
        and math.isclose(
            average,
            sum(weight * item for weight, item in zip(weights, torque)),
            rel_tol=1.0e-12,
            abs_tol=1.0e-12,
        )
        and math.isclose(result_average, average, rel_tol=0.0, abs_tol=1.0e-12)
        and math.isclose(ripple, max(torque) - min(torque), rel_tol=0.0, abs_tol=1.0e-12)
        and math.isclose(result_ripple, ripple, rel_tol=0.0, abs_tol=1.0e-12)
        and _valid_sha256(value.get("mesh_sha256"))
        and value.get("result_mesh_sha256") == value.get("mesh_sha256")
        and _valid_sha256(value.get("result_sha256"))
        and value.get("accepted_result_sha256") == value.get("result_sha256")
    )


def _iron_loss_component_volume_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, dict):
        return False
    generation = str(value.get("iron_loss_generation", "")).strip()
    try:
        frequency = float(value.get("frequency_hz"))
        result_frequency = float(value.get("result_frequency_hz"))
        orders = [int(item) for item in value.get("harmonic_orders", [])]
        result_orders = [int(item) for item in value.get("result_harmonic_orders", [])]
        flux = [float(item) for item in value.get("flux_density_harmonic_t", [])]
        result_flux = [
            float(item) for item in value.get("result_flux_density_harmonic_t", [])
        ]
        hysteresis = [
            float(item) for item in value.get("hysteresis_component_w", [])
        ]
        result_hysteresis = [
            float(item) for item in value.get("result_hysteresis_component_w", [])
        ]
        eddy = [float(item) for item in value.get("eddy_component_w", [])]
        result_eddy = [
            float(item) for item in value.get("result_eddy_component_w", [])
        ]
        anomalous = [
            float(item) for item in value.get("anomalous_component_w", [])
        ]
        result_anomalous = [
            float(item) for item in value.get("result_anomalous_component_w", [])
        ]
        total = float(value.get("total_iron_loss_w"))
        result_total = float(value.get("result_total_iron_loss_w"))
        element_ids = [int(item) for item in value.get("element_ids", [])]
        result_element_ids = [
            int(item) for item in value.get("result_element_ids", [])
        ]
        volumes = [float(item) for item in value.get("element_volumes_m3", [])]
        result_volumes = [
            float(item) for item in value.get("result_element_volumes_m3", [])
        ]
    except (TypeError, ValueError):
        return False
    count = len(orders)
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "hysteresis_iron_loss_generation",
                "eddy_iron_loss_generation",
                "anomalous_iron_loss_generation",
                "frequency_iron_loss_generation",
                "harmonic_iron_loss_generation",
                "material_iron_loss_generation",
                "volume_iron_loss_generation",
                "mesh_iron_loss_generation",
                "result_iron_loss_generation",
            )
        )
        and math.isfinite(frequency)
        and frequency > 0.0
        and math.isclose(result_frequency, frequency, rel_tol=0.0, abs_tol=1.0e-12)
        and count > 0
        and len(set(orders)) == count
        and all(item > 0 for item in orders)
        and all(left < right for left, right in zip(orders, orders[1:]))
        and result_orders == orders
        and all(len(items) == count for items in (flux, hysteresis, eddy, anomalous))
        and all(
            math.isfinite(item) and item >= 0.0
            for items in (flux, hysteresis, eddy, anomalous)
            for item in items
        )
        and result_flux == flux
        and result_hysteresis == hysteresis
        and result_eddy == eddy
        and result_anomalous == anomalous
        and math.isfinite(total)
        and total >= 0.0
        and math.isclose(
            total,
            sum(hysteresis) + sum(eddy) + sum(anomalous),
            rel_tol=1.0e-12,
            abs_tol=1.0e-12,
        )
        and math.isclose(result_total, total, rel_tol=0.0, abs_tol=1.0e-12)
        and bool(element_ids)
        and len(set(element_ids)) == len(element_ids)
        and all(item > 0 for item in element_ids)
        and len(volumes) == len(element_ids)
        and all(math.isfinite(item) and item > 0.0 for item in volumes)
        and result_element_ids == element_ids
        and result_volumes == volumes
        and _valid_sha256(value.get("material_coefficients_sha256"))
        and value.get("result_material_coefficients_sha256")
        == value.get("material_coefficients_sha256")
        and _valid_sha256(value.get("mesh_sha256"))
        and value.get("result_mesh_sha256") == value.get("mesh_sha256")
        and _valid_sha256(value.get("result_sha256"))
        and value.get("accepted_result_sha256") == value.get("result_sha256")
    )


def _induction_power_frame_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, dict):
        return False
    generation = str(value.get("induction_generation", "")).strip()
    try:
        stator_frequency = float(value.get("stator_frequency_hz"))
        result_stator_frequency = float(value.get("result_stator_frequency_hz"))
        slip = float(value.get("slip"))
        result_slip = float(value.get("result_slip"))
        rotor_frequency = float(value.get("rotor_frequency_hz"))
        result_rotor_frequency = float(value.get("result_rotor_frequency_hz"))
        pole_pairs = int(value.get("pole_pairs"))
        result_pole_pairs = int(value.get("result_pole_pairs"))
        rotor_current = [
            float(item) for item in value.get("rotor_current_rms_a", [])
        ]
        result_rotor_current = [
            float(item) for item in value.get("result_rotor_current_rms_a", [])
        ]
        torque = float(value.get("torque_nm"))
        result_torque = float(value.get("result_torque_nm"))
        speed = float(value.get("mechanical_speed_rad_s"))
        result_speed = float(value.get("result_mechanical_speed_rad_s"))
        output = float(value.get("mechanical_output_w"))
        result_output = float(value.get("result_mechanical_output_w"))
        stator_copper = float(value.get("stator_copper_loss_w"))
        result_stator_copper = float(value.get("result_stator_copper_loss_w"))
        rotor_copper = float(value.get("rotor_copper_loss_w"))
        result_rotor_copper = float(value.get("result_rotor_copper_loss_w"))
        iron = float(value.get("iron_loss_w"))
        result_iron = float(value.get("result_iron_loss_w"))
        mechanical_loss = float(value.get("mechanical_loss_w"))
        result_mechanical_loss = float(value.get("result_mechanical_loss_w"))
        electrical_input = float(value.get("electrical_input_w"))
        result_electrical_input = float(value.get("result_electrical_input_w"))
    except (TypeError, ValueError):
        return False
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "stator_frequency_induction_generation",
                "slip_induction_generation",
                "rotor_frequency_induction_generation",
                "rotor_current_induction_generation",
                "frame_induction_generation",
                "torque_induction_generation",
                "power_induction_generation",
                "loss_induction_generation",
                "result_induction_generation",
            )
        )
        and math.isfinite(stator_frequency)
        and stator_frequency > 0.0
        and math.isclose(
            result_stator_frequency, stator_frequency, rel_tol=0.0, abs_tol=1.0e-12
        )
        and math.isfinite(slip)
        and 0.0 <= slip < 1.0
        and math.isclose(result_slip, slip, rel_tol=0.0, abs_tol=1.0e-12)
        and math.isclose(
            rotor_frequency,
            slip * stator_frequency,
            rel_tol=1.0e-12,
            abs_tol=1.0e-12,
        )
        and math.isclose(
            result_rotor_frequency, rotor_frequency, rel_tol=0.0, abs_tol=1.0e-12
        )
        and pole_pairs > 0
        and result_pole_pairs == pole_pairs
        and len(rotor_current) == 3
        and all(math.isfinite(item) and item >= 0.0 for item in rotor_current)
        and result_rotor_current == rotor_current
        and value.get("reference_frame") == "stator-mechanical-ccw"
        and value.get("result_reference_frame") == value.get("reference_frame")
        and all(
            math.isfinite(item)
            for item in (
                torque,
                speed,
                output,
                stator_copper,
                rotor_copper,
                iron,
                mechanical_loss,
                electrical_input,
            )
        )
        and all(item >= 0.0 for item in (stator_copper, rotor_copper, iron, mechanical_loss))
        and math.isclose(
            speed,
            (1.0 - slip) * 2.0 * math.pi * stator_frequency / pole_pairs,
            rel_tol=1.0e-12,
            abs_tol=1.0e-12,
        )
        and math.isclose(output, torque * speed, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(
            electrical_input,
            output + stator_copper + rotor_copper + iron + mechanical_loss,
            rel_tol=1.0e-12,
            abs_tol=1.0e-12,
        )
        and math.isclose(result_torque, torque, rel_tol=0.0, abs_tol=1.0e-12)
        and math.isclose(result_speed, speed, rel_tol=0.0, abs_tol=1.0e-12)
        and math.isclose(result_output, output, rel_tol=0.0, abs_tol=1.0e-12)
        and math.isclose(result_stator_copper, stator_copper, rel_tol=0.0, abs_tol=1.0e-12)
        and math.isclose(result_rotor_copper, rotor_copper, rel_tol=0.0, abs_tol=1.0e-12)
        and math.isclose(result_iron, iron, rel_tol=0.0, abs_tol=1.0e-12)
        and math.isclose(
            result_mechanical_loss, mechanical_loss, rel_tol=0.0, abs_tol=1.0e-12
        )
        and math.isclose(
            result_electrical_input, electrical_input, rel_tol=0.0, abs_tol=1.0e-12
        )
        and _valid_sha256(value.get("result_sha256"))
        and value.get("accepted_result_sha256") == value.get("result_sha256")
    )


def _ipm_dq_inductance_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, dict):
        return False
    generation = str(value.get("dq_generation", "")).strip()
    try:
        current = float(value.get("current_magnitude_a"))
        angle = float(value.get("current_angle_electrical_deg"))
        operating = [float(item) for item in value.get("saturation_operating_point_a", [])]
        result_operating = [float(item) for item in value.get("result_saturation_operating_point_a", [])]
        matrix = [[float(item) for item in row] for row in value.get("flux_linkage_derivative_h", [])]
        result_matrix = [[float(item) for item in row] for row in value.get("result_flux_linkage_derivative_h", [])]
        tolerance = float(value.get("reciprocity_tolerance_h"))
    except (TypeError, ValueError):
        return False
    return (
        bool(generation)
        and all(value.get(key) == generation for key in (
            "current_dq_generation", "frame_dq_generation", "saturation_dq_generation",
            "flux_dq_generation", "derivative_dq_generation", "reciprocity_dq_generation",
            "mesh_dq_generation", "result_dq_generation"))
        and math.isfinite(current) and current > 0.0
        and value.get("result_current_magnitude_a") == current
        and math.isfinite(angle) and value.get("result_current_angle_electrical_deg") == angle
        and value.get("park_frame") == "rotor_d_aligned_ccw_power_invariant"
        and value.get("result_park_frame") == value.get("park_frame")
        and len(operating) == 2 and all(math.isfinite(item) for item in operating)
        and result_operating == operating
        and len(matrix) == 2 and all(len(row) == 2 for row in matrix)
        and all(math.isfinite(item) for row in matrix for item in row)
        and matrix[0][0] > 0.0 and matrix[1][1] > 0.0
        and result_matrix == matrix
        and math.isfinite(tolerance) and tolerance >= 0.0
        and value.get("result_reciprocity_tolerance_h") == tolerance
        and abs(matrix[0][1] - matrix[1][0]) <= tolerance
        and _valid_sha256(value.get("mesh_sha256"))
        and value.get("result_mesh_sha256") == value.get("mesh_sha256")
        and _valid_sha256(value.get("result_sha256"))
        and value.get("accepted_result_sha256") == value.get("result_sha256")
    )


def _srm_coenergy_torque_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, dict):
        return False
    generation = str(value.get("srm_generation", "")).strip()
    try:
        currents = [float(item) for item in value.get("current_a", [])]
        positions = [float(item) for item in value.get("rotor_position_mechanical_deg", [])]
        coenergy = [float(item) for item in value.get("coenergy_j_at_50a", [])]
        torque = float(value.get("torque_nm_at_50a"))
        period = float(value.get("sector_period_mechanical_deg"))
    except (TypeError, ValueError):
        return False
    derivative = (
        (coenergy[2] - coenergy[0]) / math.radians(positions[2] - positions[0])
        if len(coenergy) == 3 and len(positions) == 3 and positions[2] != positions[0]
        else math.nan
    )
    return (
        bool(generation)
        and all(value.get(key) == generation for key in (
            "current_srm_generation", "position_srm_generation", "coenergy_srm_generation",
            "periodicity_srm_generation", "phase_srm_generation", "mesh_srm_generation",
            "result_srm_generation"))
        and len(currents) == 3 and all(math.isfinite(item) and item >= 0.0 for item in currents)
        and currents[0] < currents[1] < currents[2]
        and value.get("result_current_a") == currents
        and len(positions) == 3 and positions[0] < positions[1] < positions[2]
        and math.isclose(positions[1], 0.0, abs_tol=1.0e-15)
        and math.isclose(positions[2], -positions[0], rel_tol=1.0e-12)
        and value.get("result_rotor_position_mechanical_deg") == positions
        and len(coenergy) == 3 and all(math.isfinite(item) for item in coenergy)
        and value.get("result_coenergy_j_at_50a") == coenergy
        and math.isfinite(torque) and math.isclose(torque, derivative, rel_tol=1.0e-10)
        and value.get("result_torque_nm_at_50a") == torque
        and math.isfinite(period) and period > 0.0
        and value.get("result_sector_period_mechanical_deg") == period
        and value.get("phase_sequence") == ["A", "B", "C"]
        and value.get("result_phase_sequence") == value.get("phase_sequence")
        and _valid_sha256(value.get("mesh_sha256"))
        and value.get("result_mesh_sha256") == value.get("mesh_sha256")
        and _valid_sha256(value.get("result_sha256"))
        and value.get("accepted_result_sha256") == value.get("result_sha256")
    )


def _pwm_sampling_loss_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, dict):
        return False
    generation = str(value.get("loss_generation", "")).strip()
    try:
        dt = float(value.get("sample_period_s"))
        result_dt = float(value.get("result_sample_period_s"))
        samples = int(value.get("samples_per_fundamental_cycle"))
        result_samples = int(value.get("result_samples_per_fundamental_cycle"))
        carrier = float(value.get("carrier_frequency_hz"))
        result_carrier = float(value.get("result_carrier_frequency_hz"))
        fundamental = float(value.get("fundamental_frequency_hz"))
        result_fundamental = float(value.get("result_fundamental_frequency_hz"))
        sidebands = [float(item) for item in value.get("carrier_sidebands_hz", [])]
        result_sidebands = [float(item) for item in value.get("result_carrier_sidebands_hz", [])]
        pole_pairs = int(value.get("pole_pairs"))
        result_pole_pairs = int(value.get("result_pole_pairs"))
        mechanical = [float(item) for item in value.get("mechanical_angle_deg", [])]
        result_mechanical = [float(item) for item in value.get("result_mechanical_angle_deg", [])]
        electrical = [float(item) for item in value.get("electrical_angle_deg", [])]
        result_electrical = [float(item) for item in value.get("result_electrical_angle_deg", [])]
        gross = float(value.get("gross_volume_m3")); result_gross = float(value.get("result_gross_volume_m3"))
        active = float(value.get("active_volume_m3")); result_active = float(value.get("result_active_volume_m3"))
        stacking = float(value.get("stacking_factor")); result_stacking = float(value.get("result_stacking_factor"))
        loss = float(value.get("mean_iron_loss_w")); result_loss = float(value.get("result_mean_iron_loss_w"))
        energy = float(value.get("cycle_energy_j")); result_energy = float(value.get("result_cycle_energy_j"))
    except (TypeError, ValueError):
        return False
    return (
        bool(generation)
        and all(value.get(key) == generation for key in (
            "sampling_loss_generation", "sideband_loss_generation", "angle_loss_generation",
            "alias_loss_generation", "volume_loss_generation", "energy_loss_generation",
            "mesh_loss_generation", "result_loss_generation"))
        and dt > 0.0 and result_dt == dt and carrier > fundamental > 0.0
        and result_carrier == carrier and result_fundamental == fundamental
        and samples >= 4 and result_samples == samples
        and math.isclose(samples * dt, 1.0 / fundamental, rel_tol=1.0e-12)
        and sidebands == [carrier - fundamental, carrier + fundamental]
        and result_sidebands == sidebands and 0.5 / dt > max(sidebands)
        and pole_pairs > 0 and result_pole_pairs == pole_pairs
        and len(mechanical) == len(electrical) >= 3 and result_mechanical == mechanical
        and result_electrical == electrical
        and all(math.isclose(e, pole_pairs * m, rel_tol=0.0, abs_tol=1.0e-12) for m, e in zip(mechanical, electrical))
        and value.get("alias_filter") == "nyquist_guard_and_sideband_keep"
        and value.get("result_alias_filter") == value.get("alias_filter")
        and gross > 0.0 and result_gross == gross and 0.0 < stacking <= 1.0
        and result_stacking == stacking and math.isclose(active, gross * stacking, rel_tol=1.0e-12)
        and result_active == active and loss >= 0.0 and result_loss == loss
        and math.isclose(energy, loss / fundamental, rel_tol=1.0e-12)
        and result_energy == energy and _valid_sha256(value.get("mesh_sha256"))
        and value.get("result_mesh_sha256") == value.get("mesh_sha256")
        and _valid_sha256(value.get("result_sha256"))
        and value.get("accepted_result_sha256") == value.get("result_sha256")
    )


def _skew_slice_torque_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, dict):
        return False
    generation = str(value.get("skew_generation", "")).strip()
    try:
        axial = [float(item) for item in value.get("axial_locations_m", [])]
        result_axial = [float(item) for item in value.get("result_axial_locations_m", [])]
        weights = [float(item) for item in value.get("slice_weights", [])]
        result_weights = [float(item) for item in value.get("result_slice_weights", [])]
        offsets = [float(item) for item in value.get("electrical_phase_offsets_deg", [])]
        result_offsets = [float(item) for item in value.get("result_electrical_phase_offsets_deg", [])]
        wrap = float(value.get("periodic_wrap_electrical_deg")); result_wrap = float(value.get("result_periodic_wrap_electrical_deg"))
        torques = [float(item) for item in value.get("slice_mean_torque_nm", [])]
        result_torques = [float(item) for item in value.get("result_slice_mean_torque_nm", [])]
        mean = float(value.get("skew_mean_torque_nm")); result_mean = float(value.get("result_skew_mean_torque_nm"))
    except (TypeError, ValueError):
        return False
    ripple = value.get("slice_ripple_harmonics_nm")
    result_ripple = value.get("result_slice_ripple_harmonics_nm")
    aggregate = value.get("skew_ripple_harmonics_nm")
    result_aggregate = value.get("result_skew_ripple_harmonics_nm")
    return (
        bool(generation)
        and all(value.get(key) == generation for key in (
            "weight_skew_generation", "axial_skew_generation", "phase_skew_generation",
            "periodicity_skew_generation", "torque_skew_generation", "ripple_skew_generation",
            "mesh_skew_generation", "result_skew_generation"))
        and len(axial) == len(weights) == len(offsets) == len(torques) >= 3
        and all(right > left for left, right in zip(axial, axial[1:])) and result_axial == axial
        and all(item >= 0.0 for item in weights) and math.isclose(sum(weights), 1.0, abs_tol=1.0e-12)
        and result_weights == weights and result_offsets == offsets
        and math.isclose(offsets[0], -offsets[-1], abs_tol=1.0e-12)
        and wrap == 360.0 and result_wrap == wrap and result_torques == torques
        and math.isclose(mean, sum(w * t for w, t in zip(weights, torques)), rel_tol=1.0e-12)
        and result_mean == mean and isinstance(ripple, list) and len(ripple) == len(weights)
        and result_ripple == ripple and isinstance(aggregate, dict) and result_aggregate == aggregate
        and set(aggregate) == {"6"}
        and math.isclose(float(aggregate["6"]), sum(w * float(row["6"]) for w, row in zip(weights, ripple)), rel_tol=1.0e-12)
        and _valid_sha256(value.get("mesh_sha256")) and value.get("result_mesh_sha256") == value.get("mesh_sha256")
        and _valid_sha256(value.get("result_sha256")) and value.get("accepted_result_sha256") == value.get("result_sha256")
    )


def _ipm_demagnetization_closure_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, dict):
        return False
    generation = str(value.get("demag_generation", "")).strip()
    try:
        temperature = float(value.get("magnet_temperature_c"))
        result_temperature = float(value.get("result_magnet_temperature_c"))
        current = float(value.get("phase_current_rms_a"))
        result_current = float(value.get("result_phase_current_rms_a"))
        angle = float(value.get("current_angle_electrical_deg"))
        result_angle = float(value.get("result_current_angle_electrical_deg"))
        fraction = float(value.get("demagnetized_fraction"))
        result_fraction = float(value.get("result_demagnetized_fraction"))
    except (TypeError, ValueError):
        return False
    regions = [
        str(item).strip() for item in value.get("irreversible_region_labels", [])
    ]
    result_regions = [
        str(item).strip()
        for item in value.get("result_irreversible_region_labels", [])
    ]
    return (
        bool(generation)
        and all(
            value.get(key) == generation
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
        )
        and value.get("knee_criterion") == "b_parallel_below_temperature_knee"
        and value.get("result_knee_criterion") == value.get("knee_criterion")
        and math.isfinite(temperature)
        and math.isclose(result_temperature, temperature, rel_tol=0.0, abs_tol=1.0e-12)
        and math.isfinite(current)
        and current >= 0.0
        and math.isclose(result_current, current, rel_tol=0.0, abs_tol=1.0e-12)
        and math.isfinite(angle)
        and math.isclose(result_angle, angle, rel_tol=0.0, abs_tol=1.0e-12)
        and bool(regions)
        and all(regions)
        and len(set(regions)) == len(regions)
        and result_regions == regions
        and math.isfinite(fraction)
        and 0.0 <= fraction <= 1.0
        and math.isclose(result_fraction, fraction, rel_tol=0.0, abs_tol=1.0e-15)
        and bool(str(value.get("operating_point_owner", "")).strip())
        and value.get("result_operating_point_owner")
        == value.get("operating_point_owner")
        and _valid_sha256(value.get("mesh_sha256"))
        and value.get("result_mesh_sha256") == value.get("mesh_sha256")
        and _valid_sha256(value.get("result_sha256"))
        and value.get("accepted_result_sha256") == value.get("result_sha256")
    )


def _synrm_dq_map_closure_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, dict):
        return False
    generation = str(value.get("map_generation", "")).strip()
    rows = value.get("flux_map_rows")
    result_rows = value.get("result_flux_map_rows")
    if not isinstance(rows, list) or not isinstance(result_rows, list):
        return False
    try:
        angles = [float(item) for item in value.get("electrical_angle_deg", [])]
        result_angles = [
            float(item) for item in value.get("result_electrical_angle_deg", [])
        ]
        normalized_rows = [
            {
                key: float(row[key])
                for key in ("id_a", "iq_a", "psi_d_wb", "psi_q_wb")
            }
            for row in rows
        ]
        normalized_result_rows = [
            {
                key: float(row[key])
                for key in ("id_a", "iq_a", "psi_d_wb", "psi_q_wb")
            }
            for row in result_rows
        ]
        dpsi_d_diq = [float(item) for item in value.get("dpsi_d_diq_h", [])]
        result_dpsi_d_diq = [
            float(item) for item in value.get("result_dpsi_d_diq_h", [])
        ]
        dpsi_q_did = [float(item) for item in value.get("dpsi_q_did_h", [])]
        result_dpsi_q_did = [
            float(item) for item in value.get("result_dpsi_q_did_h", [])
        ]
        mtpa = [int(item) for item in value.get("mtpa_row_indices", [])]
        result_mtpa = [int(item) for item in value.get("result_mtpa_row_indices", [])]
        pole_pairs = int(value.get("pole_pairs"))
        result_pole_pairs = int(value.get("result_pole_pairs"))
        torque = [float(item) for item in value.get("torque_nm", [])]
        result_torque = [float(item) for item in value.get("result_torque_nm", [])]
    except (KeyError, TypeError, ValueError):
        return False
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "angle_map_generation",
                "saturation_map_generation",
                "cross_coupling_map_generation",
                "mtpa_map_generation",
                "torque_map_generation",
                "mesh_map_generation",
                "result_map_generation",
            )
        )
        and len(angles) >= 3
        and all(math.isfinite(item) for item in angles)
        and all(left < right for left, right in zip(angles, angles[1:]))
        and result_angles == angles
        and value.get("saturation_branch") == "nonlinear_forward"
        and value.get("result_saturation_branch") == value.get("saturation_branch")
        and len(normalized_rows) >= 3
        and all(math.isfinite(item) for row in normalized_rows for item in row.values())
        and normalized_result_rows == normalized_rows
        and len(dpsi_d_diq) == len(dpsi_q_did) == len(normalized_rows)
        and all(math.isfinite(item) for item in dpsi_d_diq + dpsi_q_did)
        and all(
            math.isclose(left, right, rel_tol=1.0e-12, abs_tol=1.0e-15)
            for left, right in zip(dpsi_d_diq, dpsi_q_did)
        )
        and result_dpsi_d_diq == dpsi_d_diq
        and result_dpsi_q_did == dpsi_q_did
        and bool(mtpa)
        and len(set(mtpa)) == len(mtpa)
        and all(0 <= item < len(normalized_rows) for item in mtpa)
        and result_mtpa == mtpa
        and pole_pairs > 0
        and result_pole_pairs == pole_pairs
        and value.get("torque_reconstruction")
        == "1.5*p*(psi_d*iq-psi_q*id)"
        and value.get("result_torque_reconstruction")
        == value.get("torque_reconstruction")
        and len(torque) == len(normalized_rows)
        and all(math.isfinite(item) for item in torque)
        and result_torque == torque
        and _valid_sha256(value.get("mesh_sha256"))
        and value.get("result_mesh_sha256") == value.get("mesh_sha256")
        and _valid_sha256(value.get("map_sha256"))
        and value.get("accepted_map_sha256") == value.get("map_sha256")
    )


def _srm_commutation_closure_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, dict):
        return False
    generation = str(value.get("srm_generation", "")).strip()
    try:
        turn_on = [float(item) for item in value.get("turn_on_deg", [])]
        result_turn_on = [
            float(item) for item in value.get("result_turn_on_deg", [])
        ]
        turn_off = [float(item) for item in value.get("turn_off_deg", [])]
        result_turn_off = [
            float(item) for item in value.get("result_turn_off_deg", [])
        ]
        chop = float(value.get("current_chop_a"))
        result_chop = float(value.get("result_current_chop_a"))
        overlap = float(value.get("overlap_deg"))
        result_overlap = float(value.get("result_overlap_deg"))
        angles = [float(item) for item in value.get("angle_grid_rad", [])]
        result_angles = [
            float(item) for item in value.get("result_angle_grid_rad", [])
        ]
        coenergy = [float(item) for item in value.get("coenergy_j", [])]
        result_coenergy = [
            float(item) for item in value.get("result_coenergy_j", [])
        ]
        torque = [float(item) for item in value.get("torque_nm", [])]
        result_torque = [
            float(item) for item in value.get("result_torque_nm", [])
        ]
        copper_loss = float(value.get("copper_loss_w"))
        result_copper_loss = float(value.get("result_copper_loss_w"))
        iron_loss = float(value.get("iron_loss_w"))
        result_iron_loss = float(value.get("result_iron_loss_w"))
        total_loss = float(value.get("total_loss_w"))
        result_total_loss = float(value.get("result_total_loss_w"))
    except (TypeError, ValueError):
        return False
    phases = [str(item).strip() for item in value.get("phase_sequence", [])]
    result_phases = [
        str(item).strip() for item in value.get("result_phase_sequence", [])
    ]
    derivatives: list[float] = []
    if len(angles) == len(coenergy) and len(angles) >= 3:
        for index in range(len(angles)):
            left = max(index - 1, 0)
            right = min(index + 1, len(angles) - 1)
            width = angles[right] - angles[left]
            derivatives.append(
                (coenergy[right] - coenergy[left]) / width
                if width > 0.0
                else math.nan
            )
    dwell = [off - on for on, off in zip(turn_on, turn_off)]
    return (
        bool(generation)
        and all(
            value.get(key) == generation
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
        )
        and len(phases) >= 2
        and all(phases)
        and len(set(phases)) == len(phases)
        and result_phases == phases
        and len(turn_on) == len(turn_off) == len(phases)
        and all(math.isfinite(item) for item in turn_on + turn_off)
        and all(item > 0.0 for item in dwell)
        and result_turn_on == turn_on
        and result_turn_off == turn_off
        and math.isfinite(chop)
        and chop > 0.0
        and math.isclose(result_chop, chop, rel_tol=0.0, abs_tol=1.0e-12)
        and math.isfinite(overlap)
        and 0.0 <= overlap < min(dwell)
        and math.isclose(result_overlap, overlap, rel_tol=0.0, abs_tol=1.0e-12)
        and len(angles) >= 3
        and all(math.isfinite(item) for item in angles)
        and all(left < right for left, right in zip(angles, angles[1:]))
        and result_angles == angles
        and len(coenergy) == len(angles)
        and all(math.isfinite(item) for item in coenergy)
        and result_coenergy == coenergy
        and len(torque) == len(angles) == len(derivatives)
        and all(math.isfinite(item) for item in torque + derivatives)
        and all(
            math.isclose(observed, expected, rel_tol=1.0e-10, abs_tol=1.0e-12)
            for observed, expected in zip(torque, derivatives)
        )
        and result_torque == torque
        and all(
            math.isfinite(item) and item >= 0.0
            for item in (copper_loss, iron_loss, total_loss)
        )
        and math.isclose(
            total_loss,
            copper_loss + iron_loss,
            rel_tol=1.0e-12,
            abs_tol=1.0e-12,
        )
        and math.isclose(
            result_copper_loss, copper_loss, rel_tol=0.0, abs_tol=1.0e-12
        )
        and math.isclose(
            result_iron_loss, iron_loss, rel_tol=0.0, abs_tol=1.0e-12
        )
        and math.isclose(
            result_total_loss, total_loss, rel_tol=0.0, abs_tol=1.0e-12
        )
        and _valid_sha256(value.get("mesh_sha256"))
        and value.get("result_mesh_sha256") == value.get("mesh_sha256")
        and _valid_sha256(value.get("result_sha256"))
        and value.get("accepted_result_sha256") == value.get("result_sha256")
    )


def _axial_flux_pm_closure_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, dict):
        return False
    generation = str(value.get("axial_generation", "")).strip()
    try:
        sector_multiplier = int(value.get("sector_multiplier"))
        result_sector_multiplier = int(value.get("result_sector_multiplier"))
        air_gaps = [float(item) for item in value.get("air_gaps_m", [])]
        result_air_gaps = [
            float(item) for item in value.get("result_air_gaps_m", [])
        ]
        end_effect = float(value.get("end_effect_factor"))
        result_end_effect = float(value.get("result_end_effect_factor"))
        coordinates = [
            [float(item) for item in row]
            for row in value.get("surface_coordinates", [])
        ]
        result_coordinates = [
            [float(item) for item in row]
            for row in value.get("result_surface_coordinates", [])
        ]
        torque = [float(item) for item in value.get("torque_surface_nm", [])]
        result_torque = [
            float(item) for item in value.get("result_torque_surface_nm", [])
        ]
        force = [
            float(item) for item in value.get("axial_force_surface_n", [])
        ]
        result_force = [
            float(item) for item in value.get("result_axial_force_surface_n", [])
        ]
    except (TypeError, ValueError):
        return False
    coordinate_keys = [tuple(row) for row in coordinates]
    return (
        bool(generation)
        and all(
            value.get(key) == generation
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
        )
        and sector_multiplier > 0
        and result_sector_multiplier == sector_multiplier
        and len(air_gaps) == 2
        and all(math.isfinite(item) and item > 0.0 for item in air_gaps)
        and result_air_gaps == air_gaps
        and math.isfinite(end_effect)
        and 0.0 < end_effect <= 1.0
        and math.isclose(result_end_effect, end_effect, rel_tol=0.0, abs_tol=1.0e-12)
        and len(coordinates) >= 3
        and all(len(row) == 2 for row in coordinates)
        and all(math.isfinite(item) for row in coordinates for item in row)
        and len(set(coordinate_keys)) == len(coordinate_keys)
        and result_coordinates == coordinates
        and len(torque) == len(force) == len(coordinates)
        and all(math.isfinite(item) for item in torque + force)
        and result_torque == torque
        and result_force == force
        and value.get("force_direction") in {"+z", "-z"}
        and value.get("result_force_direction") == value.get("force_direction")
        and value.get("axial_frame") == "rotor_global_z"
        and value.get("result_axial_frame") == value.get("axial_frame")
        and _valid_sha256(value.get("mesh_sha256"))
        and value.get("result_mesh_sha256") == value.get("mesh_sha256")
        and _valid_sha256(value.get("result_lineage_sha256"))
        and value.get("accepted_result_lineage_sha256")
        == value.get("result_lineage_sha256")
    )


def _pm_demagnetization_operating_point_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, dict):
        return False
    generation = str(value.get("demag_generation", "")).strip()
    fields = (
        "reference_temperature_c", "operating_temperature_c", "remanence_reference_t",
        "remanence_temperature_coefficient_per_c", "remanence_operating_t",
        "coercivity_reference_a_m", "coercivity_temperature_coefficient_per_c",
        "coercivity_operating_a_m", "recoil_permeability_relative",
        "loadline_slope_t_per_a_m", "operating_field_a_m", "operating_flux_density_t",
        "knee_field_a_m", "irreversible_margin_a_m", "rotor_angle_rad",
        "demag_mesh_sha256",
    )
    try:
        numbers = {field: float(value.get(field)) for field in fields[:-1]}
    except (TypeError, ValueError):
        return False
    delta_temperature = numbers["operating_temperature_c"] - numbers["reference_temperature_c"]
    expected_remanence = numbers["remanence_reference_t"] * (
        1.0 + numbers["remanence_temperature_coefficient_per_c"] * delta_temperature
    )
    expected_coercivity = numbers["coercivity_reference_a_m"] * (
        1.0 + numbers["coercivity_temperature_coefficient_per_c"] * delta_temperature
    )
    expected_flux_density = numbers["remanence_operating_t"] + (
        numbers["loadline_slope_t_per_a_m"] * numbers["operating_field_a_m"]
    )
    expected_margin = numbers["operating_field_a_m"] - numbers["knee_field_a_m"]
    return (
        bool(generation)
        and all(value.get(key) == generation for key in (
            "temperature_generation", "recoil_generation", "loadline_generation",
            "operating_point_generation", "knee_generation", "margin_generation",
            "angle_generation", "mesh_generation", "owner_generation", "result_generation"))
        and all(math.isfinite(number) for number in numbers.values())
        and numbers["operating_temperature_c"] >= numbers["reference_temperature_c"]
        and numbers["remanence_reference_t"] > 0.0
        and numbers["coercivity_reference_a_m"] > 0.0
        and numbers["recoil_permeability_relative"] > 0.0
        and numbers["loadline_slope_t_per_a_m"] > 0.0
        and math.isclose(numbers["remanence_operating_t"], expected_remanence, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(numbers["coercivity_operating_a_m"], expected_coercivity, rel_tol=1.0e-12, abs_tol=1.0e-6)
        and numbers["coercivity_operating_a_m"] > 0.0
        and -numbers["coercivity_operating_a_m"] <= numbers["knee_field_a_m"] < 0.0
        and math.isclose(numbers["operating_flux_density_t"], expected_flux_density, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(numbers["irreversible_margin_a_m"], expected_margin, rel_tol=1.0e-12, abs_tol=1.0e-6)
        and numbers["irreversible_margin_a_m"] > 0.0
        and all(value.get(f"result_{field}") == value.get(field) for field in fields)
        and _valid_sha256(value.get("demag_mesh_sha256"))
        and bool(str(value.get("demag_result_owner", "")).strip())
        and value.get("accepted_demag_result_owner") == value.get("demag_result_owner")
        and _valid_sha256(value.get("demag_result_sha256"))
        and value.get("accepted_demag_result_sha256") == value.get("demag_result_sha256")
    )


def _eccentricity_ump_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, dict):
        return False
    generation = str(value.get("eccentricity_generation", "")).strip()
    try:
        static = [float(item) for item in value.get("static_eccentricity_m", [])]
        dynamic = float(value.get("dynamic_eccentricity_amplitude_m"))
        harmonics = [[int(row[0]), float(row[1]), float(row[2])] for row in value.get("radial_force_harmonics_n", [])]
        ump = [float(item) for item in value.get("unbalanced_magnetic_pull_n", [])]
        torque = float(value.get("torque_nm"))
        pole_pairs = int(value.get("pole_pairs"))
        periodicity = float(value.get("periodicity_angle_rad"))
        angles = [float(item) for item in value.get("angle_grid_rad", [])]
    except (IndexError, TypeError, ValueError):
        return False
    first_harmonic = next((row[1:] for row in harmonics if row[0] == 1), None)
    static_norm = math.hypot(*static) if len(static) == 2 else 0.0
    ump_norm = math.hypot(*ump) if len(ump) == 2 else 0.0
    aligned = (
        static_norm > 0.0
        and ump_norm > 0.0
        and abs(static[0] * ump[1] - static[1] * ump[0]) <= 1.0e-12 * static_norm * ump_norm
        and static[0] * ump[0] + static[1] * ump[1] > 0.0
    )
    mirrored = (
        "static_eccentricity_m", "dynamic_eccentricity_amplitude_m", "mechanical_frame",
        "radial_force_harmonics_n", "unbalanced_magnetic_pull_n", "torque_nm",
        "pole_pairs", "periodicity_angle_rad", "angle_grid_rad", "eccentricity_mesh_sha256",
    )
    return (
        bool(generation)
        and all(value.get(key) == generation for key in (
            "static_generation", "dynamic_generation", "frame_generation",
            "harmonic_generation", "force_generation", "torque_generation",
            "periodicity_generation", "angle_generation", "owner_generation", "result_generation"))
        and len(static) == 2 and all(math.isfinite(item) for item in static)
        and math.isfinite(dynamic) and dynamic >= 0.0
        and value.get("mechanical_frame") == "stator_global_xy"
        and len(harmonics) >= 2
        and len({row[0] for row in harmonics}) == len(harmonics)
        and all(row[0] >= 0 and all(math.isfinite(item) for item in row[1:]) for row in harmonics)
        and len(ump) == 2 and all(math.isfinite(item) for item in ump)
        and first_harmonic is not None
        and all(math.isclose(item, expected, rel_tol=1.0e-12, abs_tol=1.0e-12) for item, expected in zip(ump, first_harmonic))
        and aligned
        and math.isfinite(torque)
        and pole_pairs > 0
        and math.isclose(periodicity, 2.0 * math.pi / pole_pairs, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and len(angles) >= 3 and all(math.isfinite(item) for item in angles)
        and all(left < right for left, right in zip(angles, angles[1:]))
        and math.isclose(angles[0], 0.0, rel_tol=0.0, abs_tol=1.0e-12)
        and math.isclose(angles[-1], periodicity, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and all(value.get(f"result_{field}") == value.get(field) for field in mirrored)
        and _valid_sha256(value.get("eccentricity_mesh_sha256"))
        and bool(str(value.get("eccentricity_result_owner", "")).strip())
        and value.get("accepted_eccentricity_result_owner") == value.get("eccentricity_result_owner")
        and _valid_sha256(value.get("eccentricity_result_sha256"))
        and value.get("accepted_eccentricity_result_sha256") == value.get("eccentricity_result_sha256")
    )


def _skew_harmonic_torque_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, dict):
        return False
    generation = str(value.get("skew_generation", "")).strip()
    try:
        angles = [float(item) for item in value.get("slice_angles_rad", [])]
        weights = [float(item) for item in value.get("axial_weights", [])]
        phases = [
            [int(row[0]), [float(item) for item in row[1]]]
            for row in value.get("harmonic_phase_shifts_rad", [])
        ]
        pole_pairs = int(value.get("pole_pairs"))
        periodicity = float(value.get("pole_periodicity_angle_rad"))
        torques = [float(item) for item in value.get("slice_mean_torque_nm", [])]
        mean_torque = float(value.get("weighted_mean_torque_nm"))
        ripple = [[int(row[0]), float(row[1])] for row in value.get("torque_ripple_spectrum_nm", [])]
    except (IndexError, TypeError, ValueError):
        return False
    mirrored = (
        "slice_angles_rad", "axial_weights", "harmonic_phase_shifts_rad", "pole_pairs",
        "pole_periodicity_angle_rad", "slice_mean_torque_nm", "weighted_mean_torque_nm",
        "torque_ripple_spectrum_nm", "skew_mesh_sha256",
    )
    return (
        bool(generation)
        and all(value.get(key) == generation for key in (
            "slice_generation", "weight_generation", "phase_generation",
            "periodicity_generation", "torque_generation", "ripple_generation",
            "mesh_generation", "owner_generation", "result_generation"))
        and len(angles) == len(weights) == len(torques) >= 3
        and all(math.isfinite(item) for item in angles + weights + torques + [periodicity, mean_torque])
        and all(left < right for left, right in zip(angles, angles[1:]))
        and all(item >= 0.0 for item in weights)
        and math.isclose(sum(weights), 1.0, rel_tol=0.0, abs_tol=1.0e-12)
        and pole_pairs > 0
        and math.isclose(periodicity, 2.0 * math.pi / pole_pairs, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and bool(phases)
        and len({row[0] for row in phases}) == len(phases)
        and all(
            harmonic > 0
            and len(row_phases) == len(angles)
            and all(
                math.isclose(phase, harmonic * pole_pairs * angle, rel_tol=1.0e-12, abs_tol=1.0e-12)
                for phase, angle in zip(row_phases, angles)
            )
            for harmonic, row_phases in phases
        )
        and math.isclose(mean_torque, sum(weight * torque for weight, torque in zip(weights, torques)), rel_tol=1.0e-12, abs_tol=1.0e-12)
        and bool(ripple) and ripple[0][0] == 0
        and len({row[0] for row in ripple}) == len(ripple)
        and math.isclose(ripple[0][1], mean_torque, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and all(row[0] >= 0 and math.isfinite(row[1]) and row[1] >= 0.0 for row in ripple)
        and all(value.get(f"result_{field}") == value.get(field) for field in mirrored)
        and _valid_sha256(value.get("skew_mesh_sha256"))
        and bool(str(value.get("skew_result_owner", "")).strip())
        and value.get("accepted_skew_result_owner") == value.get("skew_result_owner")
        and _valid_sha256(value.get("skew_result_sha256"))
        and value.get("accepted_skew_result_sha256") == value.get("skew_result_sha256")
    )


def _ironloss_separation_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, dict):
        return False
    generation = str(value.get("ironloss_generation", "")).strip()
    try:
        b_peak = float(value.get("b_waveform_peak_t"))
        waveform_factor = float(value.get("waveform_factor"))
        frequency = float(value.get("frequency_hz"))
        coefficients = [float(item) for item in value.get("material_coefficients", [])]
        volume = float(value.get("active_volume_m3"))
        temperature = float(value.get("temperature_c"))
        temperature_factor = float(value.get("temperature_factor"))
        components = [float(item) for item in value.get("loss_components_w_m3", [])]
        total = float(value.get("total_iron_loss_w"))
    except (TypeError, ValueError):
        return False
    if len(coefficients) != 3 or len(components) != 3:
        return False
    expected = [
        coefficients[0] * frequency * b_peak**2 * waveform_factor * temperature_factor,
        coefficients[1] * frequency**2 * b_peak**2 * waveform_factor,
        coefficients[2] * frequency**1.5 * b_peak**1.5 * waveform_factor,
    ]
    mirrored = (
        "b_waveform_peak_t", "waveform_factor", "frequency_hz", "material_coefficients",
        "active_volume_m3", "temperature_c", "temperature_factor", "loss_components_w_m3",
        "total_iron_loss_w", "waveform_sha256",
    )
    return (
        bool(generation)
        and all(value.get(key) == generation for key in (
            "component_generation", "waveform_generation", "frequency_generation",
            "coefficient_generation", "volume_generation", "temperature_generation",
            "total_generation", "owner_generation", "result_generation"))
        and all(math.isfinite(item) for item in [b_peak, waveform_factor, frequency, volume, temperature, temperature_factor, total] + coefficients + components)
        and b_peak > 0.0 and waveform_factor > 0.0 and frequency > 0.0 and volume > 0.0
        and temperature_factor > 0.0 and coefficients[0] > 0.0
        and coefficients[1] >= 0.0 and coefficients[2] >= 0.0
        and all(item >= 0.0 for item in components)
        and all(math.isclose(item, target, rel_tol=1.0e-12, abs_tol=1.0e-12) for item, target in zip(components, expected))
        and math.isclose(total, sum(components) * volume, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and all(value.get(f"result_{field}") == value.get(field) for field in mirrored)
        and _valid_sha256(value.get("waveform_sha256"))
        and bool(str(value.get("ironloss_owner", "")).strip())
        and value.get("accepted_ironloss_owner") == value.get("ironloss_owner")
        and _valid_sha256(value.get("ironloss_result_sha256"))
        and value.get("accepted_ironloss_result_sha256") == value.get("ironloss_result_sha256")
    )


def _dq_mtpa_torque_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, dict):
        return False
    generation = str(value.get("dq_generation", "")).strip()
    try:
        pole_pairs = int(value.get("pole_pairs"))
        current = [float(item) for item in value.get("current_dq_a", [])]
        current_magnitude = float(value.get("current_magnitude_a"))
        current_angle = float(value.get("current_angle_rad"))
        pm_flux = float(value.get("pm_flux_linkage_wb_turn"))
        flux = [float(item) for item in value.get("flux_linkage_dq_wb_turn", [])]
        inductance = [float(item) for item in value.get("differential_inductance_dq_h", [])]
        torque = float(value.get("torque_nm"))
        mechanical_speed = float(value.get("mechanical_speed_rad_s"))
        electrical_speed = float(value.get("electrical_speed_rad_s"))
    except (TypeError, ValueError):
        return False
    if len(current) != 2 or len(flux) != 2 or len(inductance) != 2:
        return False
    values = current + flux + inductance + [current_magnitude, current_angle, pm_flux, torque, mechanical_speed, electrical_speed]
    if not all(math.isfinite(item) for item in values):
        return False
    i_d, i_q = current
    l_d, l_q = inductance
    expected_flux = [pm_flux + l_d * i_d, l_q * i_q]
    expected_torque = 1.5 * pole_pairs * (flux[0] * i_q - flux[1] * i_d)
    mtpa_residual = pm_flux * i_d + (l_d - l_q) * (i_d**2 - i_q**2)
    mirrored = (
        "park_convention", "pole_pairs", "current_dq_a", "current_magnitude_a",
        "current_angle_rad", "pm_flux_linkage_wb_turn", "flux_linkage_dq_wb_turn",
        "differential_inductance_dq_h", "torque_nm", "mechanical_speed_rad_s",
        "electrical_speed_rad_s",
    )
    return (
        bool(generation)
        and all(value.get(key) == generation for key in (
            "park_generation", "flux_generation", "inductance_generation", "torque_generation",
            "mtpa_generation", "current_generation", "angle_generation", "speed_generation",
            "owner_generation", "result_generation"))
        and value.get("park_convention") == "power_invariant_q_leads_d"
        and pole_pairs > 0 and l_d > 0.0 and l_q > 0.0 and pm_flux > 0.0
        and math.isclose(current_magnitude, math.hypot(i_d, i_q), rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isclose(current_angle, math.atan2(i_q, i_d), rel_tol=1.0e-12, abs_tol=1.0e-15)
        and all(math.isclose(item, expected, rel_tol=1.0e-12, abs_tol=1.0e-15) for item, expected in zip(flux, expected_flux))
        and math.isclose(torque, expected_torque, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and abs(mtpa_residual) <= 1.0e-12 * max(abs(pm_flux * i_d), 1.0)
        and mechanical_speed >= 0.0
        and math.isclose(electrical_speed, pole_pairs * mechanical_speed, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and all(value.get(f"result_{field}") == value.get(field) for field in mirrored)
        and bool(str(value.get("dq_owner", "")).strip())
        and value.get("accepted_dq_owner") == value.get("dq_owner")
        and _valid_sha256(value.get("dq_result_sha256"))
        and value.get("accepted_dq_result_sha256") == value.get("dq_result_sha256")
    )


def _iron_loss_energy_balance_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, dict):
        return False
    generation = str(value.get("iron_loss_generation", "")).strip()
    try:
        frequency = float(value.get("frequency_hz"))
        flux_peak = float(value.get("flux_peak_t"))
        coefficients = [float(item) for item in value.get("loss_coefficients", [])]
        components = [float(item) for item in value.get("loss_components_w_m3", [])]
        regions = [[str(row[0]), float(row[1])] for row in value.get("regional_volumes_m3", [])]
        temperature = float(value.get("temperature_c"))
        temperature_factor = float(value.get("temperature_factor"))
        total_power = float(value.get("total_iron_loss_w"))
        duration = float(value.get("integration_duration_s"))
        energy = float(value.get("loss_energy_j"))
    except (IndexError, TypeError, ValueError):
        return False
    if len(coefficients) != 3 or len(components) != 3 or not regions:
        return False
    values = coefficients + components + [row[1] for row in regions] + [frequency, flux_peak, temperature, temperature_factor, total_power, duration, energy]
    if not all(math.isfinite(item) for item in values):
        return False
    expected_components = [
        coefficients[0] * frequency * flux_peak**2 * temperature_factor,
        coefficients[1] * frequency**2 * flux_peak**2,
        coefficients[2] * frequency**1.5 * flux_peak**1.5,
    ]
    expected_power = sum(expected_components) * sum(row[1] for row in regions)
    mirrored = (
        "frequency_hz", "flux_peak_t", "loss_coefficients", "loss_components_w_m3",
        "regional_volumes_m3", "temperature_c", "temperature_factor",
        "total_iron_loss_w", "integration_duration_s", "loss_energy_j",
    )
    return (
        bool(generation)
        and all(value.get(key) == generation for key in (
            "component_generation", "frequency_generation", "flux_generation", "region_generation",
            "thermal_generation", "power_generation", "energy_generation", "owner_generation",
            "result_generation"))
        and frequency > 0.0 and flux_peak > 0.0 and temperature_factor > 0.0
        and coefficients[0] > 0.0 and coefficients[1] >= 0.0 and coefficients[2] >= 0.0
        and len({row[0] for row in regions}) == len(regions) and all(row[0] and row[1] > 0.0 for row in regions)
        and all(math.isclose(item, expected, rel_tol=1.0e-12, abs_tol=1.0e-12) for item, expected in zip(components, expected_components))
        and math.isclose(total_power, expected_power, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and duration > 0.0 and math.isclose(energy, total_power * duration, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and all(value.get(f"result_{field}") == value.get(field) for field in mirrored)
        and bool(str(value.get("iron_loss_owner", "")).strip())
        and value.get("accepted_iron_loss_owner") == value.get("iron_loss_owner")
        and _valid_sha256(value.get("iron_loss_result_sha256"))
        and value.get("accepted_iron_loss_result_sha256") == value.get("iron_loss_result_sha256")
    )


def _induction_motor_power_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, dict):
        return False
    generation = str(value.get("induction_generation", "")).strip()
    try:
        frequency = float(value.get("electrical_frequency_hz"))
        pole_pairs = int(value.get("pole_pairs"))
        slip = float(value.get("slip"))
        synchronous_speed = float(value.get("synchronous_speed_rad_s"))
        mechanical_speed = float(value.get("mechanical_speed_rad_s"))
        airgap_power = float(value.get("airgap_power_w"))
        torque = float(value.get("electromagnetic_torque_nm"))
        rotor_loss = float(value.get("rotor_copper_loss_w"))
        mechanical_output = float(value.get("mechanical_output_w"))
        input_power = float(value.get("input_power_w"))
        efficiency = float(value.get("efficiency"))
    except (TypeError, ValueError):
        return False
    fields = (
        "electrical_frequency_hz", "pole_pairs", "slip",
        "synchronous_speed_rad_s", "mechanical_speed_rad_s",
        "airgap_power_w", "electromagnetic_torque_nm",
        "rotor_copper_loss_w", "mechanical_output_w", "input_power_w",
        "efficiency",
    )
    numbers = (
        frequency, slip, synchronous_speed, mechanical_speed, airgap_power,
        torque, rotor_loss, mechanical_output, input_power, efficiency,
    )
    return (
        bool(generation)
        and all(value.get(key) == generation for key in (
            "frequency_generation", "speed_generation", "slip_generation",
            "airgap_generation", "torque_generation", "rotor_loss_generation",
            "mechanical_generation", "efficiency_generation",
            "owner_generation", "result_generation"))
        and all(math.isfinite(item) for item in numbers)
        and frequency > 0.0
        and pole_pairs > 0
        and 0.0 < slip < 1.0
        and airgap_power > 0.0
        and input_power >= airgap_power
        and math.isclose(
            synchronous_speed,
            2.0 * math.pi * frequency / pole_pairs,
            rel_tol=1.0e-12,
            abs_tol=1.0e-12,
        )
        and math.isclose(
            mechanical_speed,
            (1.0 - slip) * synchronous_speed,
            rel_tol=1.0e-12,
            abs_tol=1.0e-12,
        )
        and math.isclose(
            torque,
            airgap_power / synchronous_speed,
            rel_tol=1.0e-12,
            abs_tol=1.0e-12,
        )
        and math.isclose(
            rotor_loss, slip * airgap_power,
            rel_tol=1.0e-12, abs_tol=1.0e-12
        )
        and math.isclose(
            mechanical_output, (1.0 - slip) * airgap_power,
            rel_tol=1.0e-12, abs_tol=1.0e-12
        )
        and math.isclose(
            mechanical_output, torque * mechanical_speed,
            rel_tol=1.0e-12, abs_tol=1.0e-12
        )
        and 0.0 < efficiency <= 1.0
        and math.isclose(
            efficiency, mechanical_output / input_power,
            rel_tol=1.0e-12, abs_tol=1.0e-15
        )
        and all(value.get(f"result_{field}") == value.get(field) for field in fields)
        and bool(str(value.get("motor_owner", "")).strip())
        and value.get("accepted_motor_owner") == value.get("motor_owner")
        and _valid_sha256(value.get("motor_result_sha256"))
        and value.get("accepted_motor_result_sha256")
        == value.get("motor_result_sha256")
    )


def _axial_flux_periodicity_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, dict):
        return False
    generation = str(value.get("axial_flux_generation", "")).strip()
    try:
        sector_factor = int(value.get("sector_factor"))
        sector_angle = float(value.get("sector_angle_rad"))
        airgaps = [float(item) for item in value.get("dual_airgap_m", [])]
        flux = [
            float(item) for item in value.get("sector_axial_flux_per_gap_wb", [])
        ]
        sector_torque = [
            float(item) for item in value.get("sector_torque_samples_nm", [])
        ]
        full_torque = [
            float(item) for item in value.get("full_machine_torque_samples_nm", [])
        ]
        average_torque = float(value.get("average_torque_nm"))
        torque_ripple = float(value.get("torque_ripple_ratio"))
        phase_angles = [
            float(item) for item in value.get("backemf_phase_angles_rad", [])
        ]
    except (TypeError, ValueError):
        return False
    numbers = (
        [sector_angle, average_torque, torque_ripple]
        + airgaps + flux + sector_torque + full_torque + phase_angles
    )
    mirrored = (
        "sector_factor", "sector_angle_rad", "dual_airgap_m",
        "sector_axial_flux_per_gap_wb", "sector_torque_samples_nm",
        "full_machine_torque_samples_nm", "average_torque_nm",
        "torque_ripple_ratio", "backemf_phase_angles_rad", "coordinate_frame",
    )
    if (
        not all(math.isfinite(item) for item in numbers)
        or len(airgaps) != 2
        or len(flux) != 2
        or len(sector_torque) < 3
        or len(full_torque) != len(sector_torque)
        or len(phase_angles) != 3
    ):
        return False
    expected_full = [sector_factor * item for item in sector_torque]
    expected_average = sum(expected_full) / len(expected_full)
    phase_vector = [
        sum(math.cos(angle) for angle in phase_angles),
        sum(math.sin(angle) for angle in phase_angles),
    ]
    return (
        bool(generation)
        and all(value.get(key) == generation for key in (
            "sector_generation", "airgap_generation", "flux_generation",
            "torque_generation", "ripple_generation", "backemf_generation",
            "frame_generation", "mesh_generation", "owner_generation",
            "result_generation"))
        and sector_factor >= 2
        and math.isclose(
            sector_angle, 2.0 * math.pi / sector_factor,
            rel_tol=1.0e-12, abs_tol=1.0e-12
        )
        and all(item > 0.0 for item in airgaps)
        and all(item > 0.0 for item in flux)
        and all(math.isclose(item, expected, rel_tol=1.0e-12, abs_tol=1.0e-12)
                for item, expected in zip(full_torque, expected_full))
        and average_torque > 0.0
        and math.isclose(
            average_torque, expected_average,
            rel_tol=1.0e-12, abs_tol=1.0e-12
        )
        and math.isclose(
            torque_ripple,
            (max(full_torque) - min(full_torque)) / average_torque,
            rel_tol=1.0e-12, abs_tol=1.0e-12
        )
        and all(abs(item) <= 1.0e-12 for item in phase_vector)
        and value.get("coordinate_frame") == "cylindrical_z_axial"
        and all(value.get(f"result_{field}") == value.get(field) for field in mirrored)
        and bool(str(value.get("mesh_owner", "")).strip())
        and value.get("accepted_mesh_owner") == value.get("mesh_owner")
        and _valid_sha256(value.get("axial_flux_result_sha256"))
        and value.get("accepted_axial_flux_result_sha256")
        == value.get("axial_flux_result_sha256")
    )


def _wound_field_synchronous_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, dict):
        return False
    generation = str(value.get("wound_field_generation", "")).strip()
    owner = str(value.get("mesh_owner", "")).strip()
    try:
        field_current = float(value.get("field_current_a"))
        excitation_inductance = float(value.get("excitation_inductance_h"))
        excitation_flux = float(value.get("excitation_flux_linkage_wb_turn"))
        torque_angle = float(value.get("torque_angle_rad"))
        pole_pairs = int(value.get("pole_pairs"))
        stator_current = float(value.get("stator_current_rms_a"))
        torque = float(value.get("electromagnetic_torque_nm"))
        speed = float(value.get("mechanical_speed_rad_s"))
        mechanical_power = float(value.get("mechanical_power_w"))
        field_resistance = float(value.get("field_resistance_ohm"))
        field_loss = float(value.get("field_copper_loss_w"))
        stator_resistance = float(value.get("stator_phase_resistance_ohm"))
        stator_loss = float(value.get("stator_copper_loss_w"))
        line_voltage = float(value.get("line_voltage_rms_v"))
        apparent_power = float(value.get("apparent_power_va"))
        active_power = float(value.get("active_input_power_w"))
        power_factor = float(value.get("power_factor"))
        residual = float(value.get("energy_balance_residual_w"))
        tolerance = float(value.get("energy_tolerance_w"))
    except (TypeError, ValueError):
        return False
    numbers = (
        field_current, excitation_inductance, excitation_flux, torque_angle,
        stator_current, torque, speed, mechanical_power, field_resistance,
        field_loss, stator_resistance, stator_loss, line_voltage,
        apparent_power, active_power, power_factor, residual, tolerance,
    )
    if not all(math.isfinite(item) for item in numbers):
        return False
    expected_flux = field_current * excitation_inductance
    expected_torque = (
        1.5
        * pole_pairs
        * excitation_flux
        * math.sqrt(2.0)
        * stator_current
        * math.sin(torque_angle)
    )
    expected_mechanical_power = torque * speed
    expected_field_loss = field_current**2 * field_resistance
    expected_stator_loss = 3.0 * stator_current**2 * stator_resistance
    expected_apparent_power = math.sqrt(3.0) * line_voltage * stator_current
    expected_active_power = mechanical_power + field_loss + stator_loss
    expected_residual = active_power - expected_active_power
    mirrored = (
        "field_current_a", "excitation_inductance_h",
        "excitation_flux_linkage_wb_turn", "torque_angle_rad",
        "pole_pairs", "stator_current_rms_a", "electromagnetic_torque_nm",
        "mechanical_speed_rad_s", "mechanical_power_w",
        "field_resistance_ohm", "field_copper_loss_w",
        "stator_phase_resistance_ohm", "stator_copper_loss_w",
        "line_voltage_rms_v", "apparent_power_va",
        "active_input_power_w", "power_factor",
        "energy_balance_residual_w", "energy_tolerance_w",
    )
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "excitation_generation", "flux_generation",
                "torque_generation", "powerfactor_generation",
                "field_loss_generation", "stator_loss_generation",
                "mechanical_generation", "energy_generation",
                "mesh_generation", "owner_generation", "result_generation",
            )
        )
        and all(value.get(f"result_{field}") == value.get(field) for field in mirrored)
        and field_current > 0.0
        and excitation_inductance > 0.0
        and excitation_flux > 0.0
        and 0.0 < torque_angle < math.pi / 2.0
        and pole_pairs > 0
        and stator_current > 0.0
        and speed > 0.0
        and field_resistance > 0.0
        and stator_resistance > 0.0
        and line_voltage > 0.0
        and math.isclose(excitation_flux, expected_flux, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isclose(torque, expected_torque, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and torque > 0.0
        and math.isclose(mechanical_power, expected_mechanical_power, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(field_loss, expected_field_loss, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(stator_loss, expected_stator_loss, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(apparent_power, expected_apparent_power, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(active_power, expected_active_power, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and 0.0 < power_factor <= 1.0
        and math.isclose(power_factor, active_power / apparent_power, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isclose(residual, expected_residual, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and tolerance >= 0.0
        and abs(residual) <= tolerance
        and bool(owner)
        and value.get("accepted_mesh_owner") == owner
        and _valid_sha256(value.get("motor_result_sha256"))
        and value.get("accepted_motor_result_sha256")
        == value.get("motor_result_sha256")
    )


def _flux_switching_pm_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, dict):
        return False
    generation = str(value.get("flux_switching_generation", "")).strip()
    owner = str(value.get("mesh_owner", "")).strip()
    try:
        slots = int(value.get("slot_count"))
        poles = int(value.get("pole_count"))
        polarity = [int(item) for item in value.get("magnet_polarity_sequence", [])]
        harmonic = int(value.get("working_harmonic_order"))
        phase_angles = [float(item) for item in value.get("backemf_phase_angles_rad", [])]
        torque_samples = [float(item) for item in value.get("torque_samples_nm", [])]
        average_torque = float(value.get("average_torque_nm"))
        torque_ripple = float(value.get("torque_ripple_ratio"))
        multiplier = int(value.get("periodic_multiplier"))
        sector_slots = int(value.get("sector_slot_count"))
        sector_poles = int(value.get("sector_pole_count"))
    except (TypeError, ValueError):
        return False
    numbers = phase_angles + torque_samples + [average_torque, torque_ripple]
    if (
        not all(math.isfinite(item) for item in numbers)
        or len(phase_angles) != 3
        or len(torque_samples) < 3
    ):
        return False
    expected_average = sum(torque_samples) / len(torque_samples)
    expected_ripple = (
        (max(torque_samples) - min(torque_samples)) / expected_average
        if expected_average > 0.0
        else math.nan
    )
    phase_vector = [
        sum(math.cos(angle) for angle in phase_angles),
        sum(math.sin(angle) for angle in phase_angles),
    ]
    mirrored = (
        "slot_count", "pole_count", "magnet_polarity_sequence",
        "phase_sequence", "working_harmonic_order",
        "backemf_phase_angles_rad", "torque_samples_nm",
        "average_torque_nm", "torque_ripple_ratio", "periodic_multiplier",
        "sector_slot_count", "sector_pole_count",
    )
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "slot_pole_generation", "polarity_generation",
                "phase_generation", "harmonic_generation",
                "backemf_generation", "torque_generation",
                "periodicity_generation", "mesh_generation",
                "owner_generation", "result_generation",
            )
        )
        and all(value.get(f"result_{field}") == value.get(field) for field in mirrored)
        and slots > 0
        and poles >= 2
        and poles % 2 == 0
        and len(polarity) == poles
        and all(item in {-1, 1} for item in polarity)
        and all(left * right == -1 for left, right in zip(polarity, polarity[1:]))
        and value.get("phase_sequence") == "ABC"
        and harmonic == poles // 2
        and all(abs(item) <= 1.0e-12 for item in phase_vector)
        and expected_average > 0.0
        and math.isclose(average_torque, expected_average, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(torque_ripple, expected_ripple, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and multiplier == math.gcd(slots, poles)
        and multiplier >= 2
        and sector_slots * multiplier == slots
        and sector_poles * multiplier == poles
        and bool(owner)
        and value.get("accepted_mesh_owner") == owner
        and _valid_sha256(value.get("flux_switching_result_sha256"))
        and value.get("accepted_flux_switching_result_sha256")
        == value.get("flux_switching_result_sha256")
    )


def _skewed_rotor_slice_closure_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, dict):
        return False
    generation = str(value.get("skew_generation", "")).strip()
    try:
        pole_pairs = int(value.get("pole_pairs"))
        angles = [float(item) for item in value.get("slice_angles_mechanical_deg", [])]
        phases = [float(item) for item in value.get("slice_phase_offsets_electrical_deg", [])]
        weights = [float(item) for item in value.get("axial_weights", [])]
        torque = [float(item) for item in value.get("slice_mean_torque_nm", [])]
        phasors = [[float(item) for item in row] for row in value.get("slice_ripple_phasor", [])]
        mean_torque = float(value.get("weighted_mean_torque_nm"))
        ripple = float(value.get("weighted_ripple_residual"))
        speed = float(value.get("mechanical_speed_rad_s"))
        power = float(value.get("mechanical_power_w"))
    except (TypeError, ValueError):
        return False
    count = len(angles)
    expected_mean = sum(weight * item for weight, item in zip(weights, torque))
    expected_ripple = (
        math.hypot(
            sum(weight * pair[0] for weight, pair in zip(weights, phasors)),
            sum(weight * pair[1] for weight, pair in zip(weights, phasors)),
        )
        if count and len(phasors) == count and all(len(pair) == 2 for pair in phasors)
        else math.nan
    )
    mirrored = (
        "pole_pairs", "slice_angles_mechanical_deg",
        "slice_phase_offsets_electrical_deg", "axial_weights",
        "slice_mean_torque_nm", "slice_ripple_phasor",
        "weighted_mean_torque_nm", "weighted_ripple_residual",
        "mechanical_speed_rad_s", "mechanical_power_w",
    )
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "slice_generation", "phase_generation", "weight_generation",
                "torque_generation", "ripple_generation", "power_generation",
                "owner_generation", "result_generation",
            )
        )
        and pole_pairs > 0
        and count >= 3
        and len(phases) == len(weights) == len(torque) == len(phasors) == count
        and all(math.isfinite(item) for item in angles + phases + weights + torque)
        and all(left < right for left, right in zip(angles, angles[1:]))
        and math.isclose(angles[0], -angles[-1], rel_tol=0.0, abs_tol=1.0e-12)
        and all(
            math.isclose(phase, pole_pairs * angle, rel_tol=1.0e-12, abs_tol=1.0e-12)
            for phase, angle in zip(phases, angles)
        )
        and all(weight > 0.0 for weight in weights)
        and math.isclose(sum(weights), 1.0, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and all(
            len(pair) == 2
            and all(math.isfinite(item) for item in pair)
            and math.isclose(math.hypot(*pair), 1.0, rel_tol=1.0e-12, abs_tol=1.0e-12)
            for pair in phasors
        )
        and math.isfinite(mean_torque) and mean_torque > 0.0
        and math.isclose(mean_torque, expected_mean, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isfinite(ripple) and 0.0 <= ripple < 1.0
        and math.isclose(ripple, expected_ripple, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isfinite(speed) and speed > 0.0
        and math.isfinite(power) and power > 0.0
        and math.isclose(power, mean_torque * speed, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and all(value.get(f"result_{field}") == value.get(field) for field in mirrored)
        and str(value.get("model_owner", "")).startswith("motor:")
        and value.get("accepted_model_owner") == value.get("model_owner")
        and _valid_sha256(value.get("skew_result_sha256"))
        and value.get("accepted_skew_result_sha256") == value.get("skew_result_sha256")
    )


def _pm_irreversible_demag_closure_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, dict):
        return False
    generation = str(value.get("demag_generation", "")).strip()
    fields = (
        "reference_temperature_c", "magnet_temperature_c",
        "remanence_reference_t", "remanence_temperature_coefficient_per_k",
        "temperature_adjusted_remanence_t", "recoil_relative_permeability",
        "operating_h_a_per_m", "operating_b_t", "knee_h_a_per_m",
        "remanence_loss_fraction", "airgap_flux_before_wb",
        "airgap_flux_after_wb", "torque_before_nm", "torque_after_nm",
    )
    try:
        numbers = {field: float(value.get(field)) for field in fields}
    except (TypeError, ValueError):
        return False
    expected_br = numbers["remanence_reference_t"] * (
        1.0
        + numbers["remanence_temperature_coefficient_per_k"]
        * (numbers["magnet_temperature_c"] - numbers["reference_temperature_c"])
    )
    expected_b = expected_br + 4.0e-7 * math.pi * numbers["recoil_relative_permeability"] * numbers["operating_h_a_per_m"]
    expected_irreversible = numbers["operating_h_a_per_m"] < numbers["knee_h_a_per_m"]
    retained = 1.0 - numbers["remanence_loss_fraction"]
    mirrored = fields + ("irreversible_region",)
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "temperature_generation", "recoil_generation", "knee_generation",
                "operating_generation", "remanence_generation", "flux_generation",
                "torque_generation", "mesh_generation", "owner_generation",
                "result_generation",
            )
        )
        and all(math.isfinite(number) for number in numbers.values())
        and numbers["magnet_temperature_c"] >= numbers["reference_temperature_c"]
        and numbers["remanence_reference_t"] > 0.0
        and numbers["remanence_temperature_coefficient_per_k"] < 0.0
        and math.isclose(numbers["temperature_adjusted_remanence_t"], expected_br, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and numbers["recoil_relative_permeability"] > 0.0
        and numbers["operating_h_a_per_m"] < 0.0
        and numbers["operating_b_t"] > 0.0
        and math.isclose(numbers["operating_b_t"], expected_b, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and bool(value.get("irreversible_region")) is expected_irreversible
        and expected_irreversible
        and 0.0 < numbers["remanence_loss_fraction"] < 1.0
        and numbers["airgap_flux_before_wb"] > 0.0
        and numbers["torque_before_nm"] > 0.0
        and math.isclose(numbers["airgap_flux_after_wb"], retained * numbers["airgap_flux_before_wb"], rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isclose(numbers["torque_after_nm"], retained * numbers["torque_before_nm"], rel_tol=1.0e-12, abs_tol=1.0e-12)
        and all(value.get(f"result_{field}") == value.get(field) for field in mirrored)
        and str(value.get("mesh_owner", "")).startswith("mesh:")
        and value.get("accepted_mesh_owner") == value.get("mesh_owner")
        and _valid_sha256(value.get("demag_result_sha256"))
        and value.get("accepted_demag_result_sha256") == value.get("demag_result_sha256")
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
    rotating_sector_torque_frame_identity_ok = True
    iron_loss_decomposition_element_identity_ok = True
    pwm_observable_alignment_identity_ok = True
    skew_slice_average_identity_ok = True
    iron_loss_component_volume_identity_ok = True
    induction_power_frame_identity_ok = True
    ipm_dq_inductance_identity_ok = True
    srm_coenergy_torque_identity_ok = True
    pwm_sampling_loss_identity_ok = True
    skew_slice_torque_v31_identity_ok = True
    ipm_demagnetization_closure_identity_ok = True
    synrm_dq_map_closure_identity_ok = True
    srm_commutation_closure_identity_ok = True
    axial_flux_pm_closure_identity_ok = True
    pm_demagnetization_identity_ok = True
    eccentricity_ump_identity_ok = True
    skew_harmonic_torque_identity_ok = True
    ironloss_separation_identity_ok = True
    dq_mtpa_torque_identity_ok = True
    iron_loss_energy_balance_identity_ok = True
    induction_motor_power_identity_ok = True
    axial_flux_periodicity_identity_ok = True
    wound_field_synchronous_identity_ok = True
    flux_switching_pm_identity_ok = True
    skewed_rotor_slice_closure_identity_ok = True
    pm_irreversible_demag_closure_identity_ok = True
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
        rotating_sector_torque_frame_identity_ok = False
        iron_loss_decomposition_element_identity_ok = False
        pwm_observable_alignment_identity_ok = False
        skew_slice_average_identity_ok = False
        iron_loss_component_volume_identity_ok = False
        induction_power_frame_identity_ok = False
        ipm_dq_inductance_identity_ok = False
        srm_coenergy_torque_identity_ok = False
        pwm_sampling_loss_identity_ok = False
        skew_slice_torque_v31_identity_ok = False
        ipm_demagnetization_closure_identity_ok = False
        synrm_dq_map_closure_identity_ok = False
        srm_commutation_closure_identity_ok = False
        axial_flux_pm_closure_identity_ok = False
        pm_demagnetization_identity_ok = False
        eccentricity_ump_identity_ok = False
        skew_harmonic_torque_identity_ok = False
        ironloss_separation_identity_ok = False
        dq_mtpa_torque_identity_ok = False
        iron_loss_energy_balance_identity_ok = False
        induction_motor_power_identity_ok = False
        axial_flux_periodicity_identity_ok = False
        wound_field_synchronous_identity_ok = False
        flux_switching_pm_identity_ok = False
        skewed_rotor_slice_closure_identity_ok = False
        pm_irreversible_demag_closure_identity_ok = False
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
        rotating_sector_torque_frame_identity_ok = _rotating_sector_torque_frame_identity_ok(
            identity_value.get(
                "rotating_sector_pole_pair_periodic_phase_skew_slice_torque_frame_generation_identity"
            )
        )
        iron_loss_decomposition_element_identity_ok = _iron_loss_decomposition_element_identity_ok(
            identity_value.get(
                "iron_loss_harmonic_decomposition_model_temperature_frequency_element_volume_result_generation_identity"
            )
        )
        pwm_observable_alignment_identity_ok = _pwm_observable_alignment_identity_ok(
            identity_value.get(
                "pwm_current_harmonic_time_electrical_angle_torque_loss_mesh_result_generation_identity"
            )
        )
        skew_slice_average_identity_ok = _skew_slice_average_identity_ok(
            identity_value.get(
                "skew_slice_angle_weight_frame_interpolation_torque_ripple_mesh_generation_identity"
            )
        )
        iron_loss_component_volume_identity_ok = _iron_loss_component_volume_identity_ok(
            identity_value.get(
                "iron_loss_component_harmonic_frequency_volume_generation_identity"
            )
        )
        induction_power_frame_identity_ok = _induction_power_frame_identity_ok(
            identity_value.get(
                "induction_slip_rotor_current_torque_power_frame_generation_identity"
            )
        )
        ipm_dq_inductance_identity_ok = _ipm_dq_inductance_identity_ok(
            identity_value.get(
                "ipm_dq_inductance_current_angle_park_saturation_flux_derivative_reciprocity_mesh_result_identity"
            )
        )
        srm_coenergy_torque_identity_ok = _srm_coenergy_torque_identity_ok(
            identity_value.get(
                "srm_torque_current_position_coenergy_periodicity_phase_sequence_mesh_result_identity"
            )
        )
        pwm_sampling_loss_identity_ok = _pwm_sampling_loss_identity_ok(
            identity_value.get(
                "pwm_iron_loss_sampling_sideband_angle_alias_volume_energy_result_identity"
            )
        )
        skew_slice_torque_v31_identity_ok = _skew_slice_torque_identity_ok(
            identity_value.get(
                "skew_slice_torque_weight_axial_phase_periodicity_ripple_mesh_result_identity"
            )
        )
        ipm_demagnetization_closure_identity_ok = (
            _ipm_demagnetization_closure_identity_ok(
                identity_value.get(
                    "ipm_demagnetization_knee_temperature_current_angle_region_fraction_mesh_result_identity"
                )
            )
        )
        synrm_dq_map_closure_identity_ok = _synrm_dq_map_closure_identity_ok(
            identity_value.get(
                "synrm_dq_map_angle_saturation_cross_coupling_mtpa_torque_mesh_result_identity"
            )
        )
        srm_commutation_closure_identity_ok = _srm_commutation_closure_identity_ok(
            identity_value.get(
                "srm_commutation_phase_dwell_chop_overlap_coenergy_torque_loss_angle_mesh_result_identity"
            )
        )
        axial_flux_pm_closure_identity_ok = _axial_flux_pm_closure_identity_ok(
            identity_value.get(
                "axial_flux_pm_sector_airgap_end_effect_torque_force_surface_direction_frame_mesh_result_identity"
            )
        )
        pm_demagnetization_identity_ok = _pm_demagnetization_operating_point_identity_ok(
            identity_value.get(
                "pm_demagnetization_temperature_recoil_loadline_operating_point_knee_margin_angle_mesh_owner_result_identity"
            )
        )
        eccentricity_ump_identity_ok = _eccentricity_ump_identity_ok(
            identity_value.get(
                "eccentricity_static_dynamic_frame_radial_force_harmonic_ump_torque_pole_periodicity_angle_owner_result_identity"
            )
        )
        skew_harmonic_torque_identity_ok = _skew_harmonic_torque_identity_ok(
            identity_value.get(
                "skew_slice_torque_angle_axial_weight_harmonic_phase_pole_periodicity_mean_ripple_mesh_owner_result_identity"
            )
        )
        ironloss_separation_identity_ok = _ironloss_separation_identity_ok(
            identity_value.get(
                "ironloss_hysteresis_eddy_excess_waveform_frequency_coeff_volume_temperature_total_owner_result_identity"
            )
        )
        dq_mtpa_torque_identity_ok = _dq_mtpa_torque_identity_ok(
            identity_value.get(
                "dq_flux_inductance_torque_mtpa_current_angle_speed_convention_owner_result_identity"
            )
        )
        iron_loss_energy_balance_identity_ok = _iron_loss_energy_balance_identity_ok(
            identity_value.get(
                "iron_loss_component_frequency_flux_region_thermal_energy_balance_owner_result_identity"
            )
        )
        induction_motor_power_identity_ok = _induction_motor_power_identity_ok(
            identity_value.get(
                "induction_motor_slip_synchronous_speed_airgap_power_torque_rotor_loss_mechanical_output_efficiency_owner_result_identity"
            )
        )
        axial_flux_periodicity_identity_ok = _axial_flux_periodicity_identity_ok(
            identity_value.get(
                "axial_flux_motor_sector_periodicity_dual_airgap_axial_flux_torque_ripple_backemf_frame_mesh_owner_result_identity"
            )
        )
        wound_field_synchronous_identity_ok = _wound_field_synchronous_identity_ok(
            identity_value.get(
                "wound_field_synchronous_excitation_flux_torque_angle_powerfactor_field_stator_loss_mechanical_energy_mesh_owner_result_identity"
            )
        )
        flux_switching_pm_identity_ok = _flux_switching_pm_identity_ok(
            identity_value.get(
                "flux_switching_pm_slot_pole_polarity_phase_harmonic_backemf_torque_ripple_periodicity_mesh_owner_result_identity"
            )
        )
        skewed_rotor_slice_closure_identity_ok = _skewed_rotor_slice_closure_identity_ok(
            identity_value.get(
                "skewed_rotor_slice_angle_phase_weight_torque_ripple_power_model_owner_result_identity"
            )
        )
        pm_irreversible_demag_closure_identity_ok = _pm_irreversible_demag_closure_identity_ok(
            identity_value.get(
                "pm_irreversible_demag_temperature_recoil_knee_operating_flux_torque_mesh_owner_result_identity"
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
        "rotating_sector_torque_uses_current_pole_pairs_periodic_phase_skew_frame_mesh_and_result": (
            rotating_sector_torque_frame_identity_ok
        ),
        "iron_loss_decomposition_uses_current_model_temperature_frequency_elements_volume_and_result": (
            iron_loss_decomposition_element_identity_ok
        ),
        "pwm_observables_use_current_harmonics_time_angle_torque_loss_mesh_and_result": (
            pwm_observable_alignment_identity_ok
        ),
        "skew_average_uses_current_angles_weights_frame_interpolation_torque_mesh_and_result": (
            skew_slice_average_identity_ok
        ),
        "iron_loss_components_use_current_frequency_harmonics_material_volumes_mesh_and_result": (
            iron_loss_component_volume_identity_ok
        ),
        "induction_motor_uses_current_slip_rotor_frequency_current_frame_torque_and_power_balance": (
            induction_power_frame_identity_ok
        ),
        "ipm_dq_inductance_uses_current_angle_park_frame_saturation_derivatives_reciprocity_mesh_and_result": (
            ipm_dq_inductance_identity_ok
        ),
        "srm_torque_uses_current_positions_coenergy_periodicity_phase_sequence_mesh_and_result": (
            srm_coenergy_torque_identity_ok
        ),
        "pwm_iron_loss_uses_current_sampling_sidebands_angles_alias_volume_energy_mesh_and_result": (
            pwm_sampling_loss_identity_ok
        ),
        "skew_torque_uses_current_slice_weights_axial_phase_periodicity_ripple_mesh_and_result": (
            skew_slice_torque_v31_identity_ok
        ),
        "ipm_demagnetization_uses_current_knee_temperature_current_angle_regions_fraction_mesh_owner_and_result": (
            ipm_demagnetization_closure_identity_ok
        ),
        "synrm_dq_map_uses_current_angles_saturation_cross_coupling_mtpa_torque_mesh_and_result": (
            synrm_dq_map_closure_identity_ok
        ),
        "srm_commutation_uses_current_phases_dwell_chop_overlap_coenergy_torque_loss_mesh_and_result": (
            srm_commutation_closure_identity_ok
        ),
        "axial_flux_pm_uses_current_sector_airgaps_end_effect_torque_force_surface_frame_mesh_and_result": (
            axial_flux_pm_closure_identity_ok
        ),
        "pm_demagnetization_uses_current_temperature_recoil_loadline_knee_margin_angle_mesh_owner_and_result": (
            pm_demagnetization_identity_ok
        ),
        "eccentricity_ump_uses_current_static_dynamic_frame_harmonics_force_torque_periodicity_angles_owner_and_result": (
            eccentricity_ump_identity_ok
        ),
        "skew_slice_torque_closes_angles_axial_weights_harmonic_phases_pole_periodicity_mean_ripple_mesh_owner_and_result": (
            skew_harmonic_torque_identity_ok
        ),
        "iron_loss_closes_hysteresis_eddy_excess_waveform_frequency_coefficients_volume_temperature_total_owner_and_result": (
            ironloss_separation_identity_ok
        ),
        "dq_map_closes_park_flux_inductance_torque_mtpa_current_angle_speed_owner_and_result": (
            dq_mtpa_torque_identity_ok
        ),
        "iron_loss_closes_components_frequency_flux_regions_thermal_power_energy_owner_and_result": (
            iron_loss_energy_balance_identity_ok
        ),
        "induction_motor_closes_slip_speed_airgap_power_torque_rotor_loss_output_efficiency_owner_and_result": (
            induction_motor_power_identity_ok
        ),
        "axial_flux_motor_closes_sector_dual_airgap_flux_torque_ripple_backemf_frame_mesh_owner_and_result": (
            axial_flux_periodicity_identity_ok
        ),
        "wound_field_motor_closes_excitation_torque_powerfactor_copper_losses_mechanical_energy_mesh_owner_and_result": (
            wound_field_synchronous_identity_ok
        ),
        "flux_switching_pm_closes_slot_pole_polarity_harmonic_backemf_torque_periodicity_mesh_owner_and_result": (
            flux_switching_pm_identity_ok
        ),
        "skewed_rotor_closes_slice_angles_phase_offsets_weights_torque_ripple_power_owner_and_result": (
            skewed_rotor_slice_closure_identity_ok
        ),
        "pm_demagnetization_closes_temperature_recoil_knee_irreversible_loss_flux_torque_mesh_and_result": (
            pm_irreversible_demag_closure_identity_ok
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
