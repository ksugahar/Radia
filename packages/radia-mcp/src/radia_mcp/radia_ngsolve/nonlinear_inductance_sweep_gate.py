"""Solver-neutral gate for nonlinear two-winding inductance sweeps."""

from __future__ import annotations

import math
from collections import defaultdict
from collections.abc import Mapping, Sequence
from typing import Any


def _valid_sha256(value: object) -> bool:
    digest = str(value or "").lower()
    return len(digest) == 64 and all(
        character in "0123456789abcdef" for character in digest
    )


def _relative_error(actual: float, expected: float) -> float:
    return abs(actual - expected) / max(abs(actual), abs(expected), 1.0e-300)


def _vector(value: Any, name: str) -> list[float]:
    if (
        not isinstance(value, Sequence)
        or isinstance(value, (str, bytes))
        or len(value) != 2
    ):
        raise ValueError(f"{name} must contain two values")
    result = [float(item) for item in value]
    if not all(math.isfinite(item) for item in result):
        raise ValueError(f"{name} must contain finite values")
    return result


def _matrix(value: Any, name: str) -> list[list[float]]:
    if (
        not isinstance(value, Sequence)
        or isinstance(value, (str, bytes))
        or len(value) != 2
        or any(
            not isinstance(row, Sequence)
            or isinstance(row, (str, bytes))
            or len(row) != 2
            for row in value
        )
    ):
        raise ValueError(f"{name} must be a 2x2 matrix")
    result = [[float(item) for item in row] for row in value]
    if not all(math.isfinite(item) for row in result for item in row):
        raise ValueError(f"{name} must contain finite values")
    return result


def _matrix_metrics(matrix: list[list[float]]) -> dict[str, float]:
    l11, m12 = matrix[0]
    m21, l22 = matrix[1]
    mutual = 0.5 * (m12 + m21)
    diagonal_product = l11 * l22
    return {
        "symmetry_relative_error": _relative_error(m12, m21),
        "determinant_H2": diagonal_product - mutual * mutual,
        "diagonal_product_H2": diagonal_product,
        "l11_H": l11,
        "l22_H": l22,
    }


def _flatten_replay_values(row: Mapping[str, Any]) -> list[float]:
    result = []
    for name in ("apparent_inductance_H", "incremental_inductance_H"):
        result.extend(item for matrix_row in row[name] for item in matrix_row)
    result.extend(row["current_A"])
    result.extend(row["flux_linkage_Vs"])
    result.extend((row["energy_J"], row["coenergy_J"]))
    return result


def _result_metadata_run_ids_are_consistent(raw: Mapping[str, Any]) -> bool:
    metadata = raw.get("result_metadata")
    if metadata is None:
        return True
    if not isinstance(metadata, Mapping) or not metadata:
        return False
    run_ids = []
    for row in metadata.values():
        if not isinstance(row, Mapping):
            return False
        run_id = row.get("run_id")
        if isinstance(run_id, bool) or not isinstance(run_id, int) or run_id < 0:
            return False
        run_ids.append(run_id)
    return len(set(run_ids)) == 1


def _matrix_operating_point_identity_matches(
    raw: Mapping[str, Any], currents: list[float]
) -> tuple[bool, bool]:
    identity_names = (
        "operating_point_id",
        "apparent_matrix_operating_point_id",
        "incremental_matrix_operating_point_id",
    )
    current_names = (
        "apparent_matrix_current_A",
        "incremental_matrix_current_A",
    )
    identity_present = any(name in raw for name in identity_names)
    current_present = any(name in raw for name in current_names)
    identity_ok = True
    current_ok = True
    if identity_present:
        values = [str(raw.get(name, "")).strip() for name in identity_names]
        identity_ok = all(values) and len(set(values)) == 1
    if current_present:
        if not all(name in raw for name in current_names):
            current_ok = False
        else:
            matrix_currents = [
                _vector(raw[name], name) for name in current_names
            ]
            current_ok = all(
                all(
                    _relative_error(actual, expected) <= 1.0e-12
                    for actual, expected in zip(value, currents)
                )
                for value in matrix_currents
            )
    return identity_ok, current_ok


def _artifact_units_are_consistent(raw: Mapping[str, Any]) -> bool:
    expected = {
        "current": "A",
        "flux_linkage": "Vs",
        "inductance": "H",
        "energy": "J",
        "coenergy": "J",
    }
    reported = raw.get("reported_units")
    artifact = raw.get("artifact_units")
    if reported is None and artifact is None:
        return True
    if not isinstance(reported, Mapping) or not isinstance(artifact, Mapping):
        return False
    return all(
        reported.get(name) == unit and artifact.get(name) == unit
        for name, unit in expected.items()
    )


def _matrix_sweep_generations_match(raw: Mapping[str, Any]) -> bool:
    names = (
        "solve_sweep_generation",
        "apparent_matrix_sweep_generation",
        "incremental_matrix_sweep_generation",
    )
    if not any(name in raw for name in names):
        return True
    values = [str(raw.get(name, "")).strip() for name in names]
    return all(values) and len(set(values)) == 1


def _matrix_port_orders_match(raw: Mapping[str, Any]) -> bool:
    order = raw.get("matrix_port_order")
    if order is None:
        return True
    if not isinstance(order, Mapping):
        return False
    names = (
        "run_current",
        "flux_linkage",
        "apparent_rows",
        "apparent_columns",
        "incremental_rows",
        "incremental_columns",
    )
    values = [order.get(name) for name in names]
    return (
        all(isinstance(value, list) for value in values)
        and all(value == ["primary", "secondary"] for value in values)
    )


def _energy_loss_basis_is_si(raw: Mapping[str, Any]) -> bool:
    basis = raw.get("energy_loss_basis")
    if basis is None:
        return True
    return (
        isinstance(basis, Mapping)
        and basis.get("stored_energy_unit") == "J"
        and basis.get("coenergy_unit") == "J"
        and basis.get("loss_series_unit") == "J"
        and basis.get("loss_series_scale_to_J") == 1.0
        and basis.get("shared_accumulation_basis") == "J"
    )


def _sparameter_reference_impedance_is_bound(raw: Mapping[str, Any]) -> bool:
    identity = raw.get("sparameter_reference_impedance")
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    ports = identity.get("port_order")
    solver_values = identity.get("solver_reference_impedance_ohm_complex")
    comparison_values = identity.get(
        "comparison_reference_impedance_ohm_complex"
    )
    if (
        not isinstance(ports, list)
        or not ports
        or len(set(ports)) != len(ports)
        or not isinstance(solver_values, list)
        or not isinstance(comparison_values, list)
        or not (len(solver_values) == len(comparison_values) == len(ports))
    ):
        return False
    try:
        solver = [complex(float(value[0]), float(value[1])) for value in solver_values]
        comparison = [
            complex(float(value[0]), float(value[1]))
            for value in comparison_values
        ]
    except (IndexError, TypeError, ValueError):
        return False
    if not all(
        math.isfinite(value.real) and math.isfinite(value.imag)
        for value in solver + comparison
    ):
        return False
    same_reference = solver == comparison
    return (
        bool(identity.get("reference_impedance_generation"))
        and (
            same_reference
            and identity.get("renormalization_applied") is False
            or not same_reference
            and identity.get("renormalization_applied") is True
            and identity.get("renormalized_port_order") == ports
            and bool(identity.get("renormalization_generation"))
        )
    )


def _frequency_axis_unit_is_bound(raw: Mapping[str, Any]) -> bool:
    identity = raw.get("frequency_axis_identity")
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    unit = str(identity.get("numeric_axis_unit", ""))
    scales = {"Hz": 1.0, "kHz": 1.0e3, "MHz": 1.0e6, "GHz": 1.0e9}
    try:
        scale = float(identity.get("scale_to_hz"))
    except (TypeError, ValueError):
        return False
    return (
        unit in scales
        and identity.get("metadata_axis_unit") == unit
        and math.isclose(scale, scales[unit], rel_tol=0.0, abs_tol=0.0)
        and identity.get("normalized_axis_unit") == "Hz"
        and identity.get("normalization_applied_once") is True
        and bool(identity.get("frequency_axis_generation"))
    )


def _sparameter_reference_planes_are_bound(raw: Mapping[str, Any]) -> bool:
    identity = raw.get("sparameter_reference_plane_identity")
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    ports = identity.get("port_order")
    original = identity.get("original_reference_plane_ids")
    target = identity.get("target_reference_plane_ids")
    compared = identity.get("compared_port_mode_reference_plane_ids")
    if not (
        isinstance(ports, list)
        and bool(ports)
        and len(set(ports)) == len(ports)
        and isinstance(original, list)
        and isinstance(target, list)
        and isinstance(compared, list)
        and len(original) == len(target) == len(compared) == len(ports)
        and all(str(value).strip() for value in original + target + compared)
    ):
        return False
    matrix_order = raw.get("matrix_port_order")
    run_order = (
        matrix_order.get("run_current")
        if isinstance(matrix_order, Mapping)
        else ports
    )
    generation = str(identity.get("deembedding_generation", "")).strip()
    return (
        run_order == ports
        and all(left != right for left, right in zip(original, target))
        and compared == target
        and identity.get("deembedding_applied") is True
        and bool(generation)
        and identity.get("sparameter_generation") == generation
    )


def _energy_q_frequency_sample_is_bound(raw: Mapping[str, Any]) -> bool:
    identity = raw.get("energy_q_frequency_identity")
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        frequencies = [
            float(identity[name])
            for name in (
                "q_frequency_hz",
                "stored_energy_frequency_hz",
                "loss_frequency_hz",
            )
        ]
    except (KeyError, TypeError, ValueError):
        return False
    sample = str(identity.get("adaptive_sample_id", "")).strip()
    return (
        all(math.isfinite(value) and value > 0.0 for value in frequencies)
        and len(set(frequencies)) == 1
        and bool(sample)
        and identity.get("stored_energy_sample_id") == sample
        and identity.get("loss_sample_id") == sample
    )


def _mixed_mode_sparameter_basis_matches_port_order(
    raw: Mapping[str, Any],
) -> bool:
    identity = raw.get("mixed_mode_sparameter_basis_identity")
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    port_order = identity.get("single_ended_port_order")
    generation = str(identity.get("port_order_generation", "")).strip()
    digest = str(identity.get("basis_matrix_sha256", "")).strip()
    return (
        isinstance(port_order, list)
        and len(port_order) == 4
        and all(isinstance(port, str) and port for port in port_order)
        and len(set(port_order)) == len(port_order)
        and identity.get("sparameter_port_order") == port_order
        and identity.get("basis_matrix_port_order") == port_order
        and bool(generation)
        and identity.get("basis_matrix_port_order_generation") == generation
        and len(digest) == 64
        and all(character in "0123456789abcdef" for character in digest)
    )


def _farfield_gain_power_frequency_sample_is_bound(raw: Mapping[str, Any]) -> bool:
    identity = raw.get("farfield_realized_gain_power_frequency_identity")
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        farfield_frequency = float(identity["farfield_frequency_hz"])
        accepted_frequency = float(identity["accepted_power_frequency_hz"])
    except (KeyError, TypeError, ValueError):
        return False
    sample_id = str(identity.get("farfield_adaptive_sample_id", "")).strip()
    result_generation = str(
        identity.get("farfield_result_generation", "")
    ).strip()
    return (
        math.isfinite(farfield_frequency)
        and farfield_frequency > 0.0
        and accepted_frequency == farfield_frequency
        and bool(sample_id)
        and identity.get("accepted_power_adaptive_sample_id") == sample_id
        and bool(result_generation)
        and identity.get("accepted_power_result_generation") == result_generation
    )


def _field_monitor_interpolation_matches_mesh_pass(raw: Mapping[str, Any]) -> bool:
    identity = raw.get("field_monitor_interpolation_mesh_pass_identity")
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    adaptive_pass = str(identity.get("active_adaptive_pass_id", "")).strip()
    mesh_generation = str(identity.get("active_mesh_generation", "")).strip()
    weight_digest = str(identity.get("interpolation_weight_sha256", "")).lower()
    integral_digest = str(identity.get("integral_weight_sha256", "")).lower()
    return (
        bool(adaptive_pass)
        and identity.get("field_monitor_adaptive_pass_id") == adaptive_pass
        and identity.get("interpolation_weight_adaptive_pass_id") == adaptive_pass
        and identity.get("integral_adaptive_pass_id") == adaptive_pass
        and bool(mesh_generation)
        and identity.get("field_monitor_mesh_generation") == mesh_generation
        and identity.get("interpolation_weight_mesh_generation") == mesh_generation
        and len(weight_digest) == len(integral_digest) == 64
        and all(
            character in "0123456789abcdef"
            for character in weight_digest + integral_digest
        )
        and integral_digest == weight_digest
    )


def _port_deembed_reference_plane_unit_is_bound(raw: Mapping[str, Any]) -> bool:
    identity = raw.get("port_deembed_reference_plane_unit_identity")
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    scales = {"m": 1.0, "cm": 1.0e-2, "mm": 1.0e-3, "um": 1.0e-6}
    model_unit = str(identity.get("model_length_unit", ""))
    offset_unit = str(identity.get("reference_plane_offset_unit", ""))
    result_unit = str(identity.get("result_reference_plane_offset_unit", ""))
    try:
        model_scale = float(identity["model_length_scale_to_m"])
        offset = float(identity["reference_plane_offset_numeric"])
        offset_scale = float(identity["reference_plane_offset_scale_to_m"])
        result_offset = float(identity["result_reference_plane_offset_numeric"])
        result_scale = float(identity["result_reference_plane_offset_scale_to_m"])
    except (KeyError, TypeError, ValueError):
        return False
    generation = str(identity.get("port_setup_generation", "")).strip()
    return (
        model_unit in scales
        and offset_unit in scales
        and result_unit in scales
        and math.isclose(model_scale, scales[model_unit], rel_tol=0.0, abs_tol=0.0)
        and math.isclose(offset_scale, scales[offset_unit], rel_tol=0.0, abs_tol=0.0)
        and math.isclose(result_scale, scales[result_unit], rel_tol=0.0, abs_tol=0.0)
        and all(
            math.isfinite(value)
            for value in (offset, offset_scale, result_offset, result_scale)
        )
        and math.isclose(
            offset * offset_scale,
            result_offset * result_scale,
            rel_tol=1.0e-12,
            abs_tol=1.0e-15,
        )
        and bool(generation)
        and identity.get("sparameter_result_generation") == generation
    )


def _sparameter_renormalization_matches_reference_impedance(
    raw: Mapping[str, Any],
) -> bool:
    identity = raw.get("sparameter_reference_impedance_renormalization_identity")
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        reference_impedance = complex(
            float(identity["reference_impedance_real_ohm"]),
            float(identity["reference_impedance_imag_ohm"]),
        )
        result_impedance = complex(
            float(identity["result_reference_impedance_real_ohm"]),
            float(identity["result_reference_impedance_imag_ohm"]),
        )
    except (KeyError, TypeError, ValueError):
        return False
    port_generation = str(identity.get("port_setup_generation", "")).strip()
    renormalization_generation = str(
        identity.get("renormalization_generation", "")
    ).strip()
    source_digest = str(identity.get("sparameter_array_sha256", "")).lower()
    result_digest = str(identity.get("renormalized_array_sha256", "")).lower()
    return (
        identity.get("reference_impedance_basis") == "complex_ohm"
        and identity.get("result_reference_impedance_basis") == "complex_ohm"
        and math.isfinite(reference_impedance.real)
        and math.isfinite(reference_impedance.imag)
        and reference_impedance.real > 0.0
        and result_impedance == reference_impedance
        and bool(port_generation)
        and identity.get("sparameter_result_generation") == port_generation
        and identity.get("renormalization_applied") is True
        and bool(renormalization_generation)
        and identity.get("result_renormalization_generation")
        == renormalization_generation
        and len(source_digest) == len(result_digest) == 64
        and all(
            character in "0123456789abcdef"
            for character in source_digest + result_digest
        )
        and result_digest == source_digest
    )


def _farfield_ludwig_polarization_basis_is_current(
    raw: Mapping[str, Any],
) -> bool:
    identity = raw.get("farfield_ludwig_polarization_basis_identity")
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    result_generation = str(
        identity.get("farfield_result_generation", "")
    ).strip()
    frame = str(identity.get("coordinate_frame_id", "")).strip()
    ludwig_basis = str(identity.get("ludwig_basis_definition", "")).strip()
    reference_axis = str(identity.get("polarization_reference_axis", "")).strip()
    basis_digest = str(identity.get("polarization_basis_sha256", "")).lower()
    result_digest = str(
        identity.get("co_cross_polarization_basis_sha256", "")
    ).lower()
    return (
        bool(result_generation)
        and identity.get("co_cross_result_generation") == result_generation
        and bool(frame)
        and identity.get("co_cross_coordinate_frame_id") == frame
        and ludwig_basis in {"ludwig_1", "ludwig_2", "ludwig_3"}
        and identity.get("co_cross_ludwig_basis_definition") == ludwig_basis
        and bool(reference_axis)
        and identity.get("co_cross_polarization_reference_axis") == reference_axis
        and len(basis_digest) == len(result_digest) == 64
        and all(
            character in "0123456789abcdef"
            for character in basis_digest + result_digest
        )
        and result_digest == basis_digest
    )


def _sparameter_power_wave_normalization_is_current(
    raw: Mapping[str, Any],
) -> bool:
    identity = raw.get("sparameter_power_wave_normalization_identity")
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        impedance = complex(
            float(identity["reference_impedance_real_ohm"]),
            float(identity["reference_impedance_imag_ohm"]),
        )
    except (KeyError, TypeError, ValueError):
        return False
    port_generation = str(identity.get("port_setup_generation", "")).strip()
    normalization_generation = str(
        identity.get("modal_normalization_generation", "")
    ).strip()
    modal_digest = str(identity.get("modal_normalization_sha256", "")).lower()
    result_digest = str(
        identity.get("sparameter_normalization_sha256", "")
    ).lower()
    return (
        bool(port_generation)
        and identity.get("modal_result_port_setup_generation") == port_generation
        and identity.get("sparameter_result_port_setup_generation")
        == port_generation
        and identity.get("incident_modal_amplitude_normalization") == "power_wave"
        and identity.get("reflected_modal_amplitude_normalization")
        == "power_wave"
        and identity.get("sparameter_normalization") == "power_wave"
        and identity.get("reference_impedance_basis") == "complex_ohm"
        and math.isfinite(impedance.real)
        and math.isfinite(impedance.imag)
        and impedance.real > 0.0
        and bool(normalization_generation)
        and identity.get("sparameter_normalization_generation")
        == normalization_generation
        and len(modal_digest) == len(result_digest) == 64
        and all(
            character in "0123456789abcdef"
            for character in modal_digest + result_digest
        )
        and result_digest == modal_digest
    )


def _fft_window_coherent_gain_is_current(raw: Mapping[str, Any]) -> bool:
    identity = raw.get("time_domain_fft_window_coherent_gain_identity")
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        sample_count = int(identity["sample_count"])
        coherent_gain = float(identity["coherent_gain"])
        gain_correction = float(identity["fft_coherent_gain_correction"])
        application_count = int(identity["coherent_gain_application_count"])
    except (KeyError, TypeError, ValueError):
        return False
    trace_generation = str(identity.get("time_trace_generation", "")).strip()
    window_generation = str(identity.get("window_generation", "")).strip()
    window_digest = str(identity.get("window_coefficients_sha256", "")).lower()
    result_digest = str(
        identity.get("fft_window_coefficients_sha256", "")
    ).lower()
    return (
        bool(trace_generation)
        and identity.get("fft_input_trace_generation") == trace_generation
        and bool(window_generation)
        and identity.get("coherent_gain_window_generation") == window_generation
        and identity.get("fft_result_window_generation") == window_generation
        and identity.get("window_definition") == "periodic_hann"
        and sample_count >= 2
        and math.isfinite(coherent_gain)
        and math.isclose(coherent_gain, 0.5, rel_tol=0.0, abs_tol=1.0e-15)
        and math.isfinite(gain_correction)
        and math.isclose(
            gain_correction,
            1.0 / coherent_gain,
            rel_tol=1.0e-12,
            abs_tol=1.0e-15,
        )
        and application_count == 1
        and len(window_digest) == len(result_digest) == 64
        and all(
            character in "0123456789abcdef"
            for character in window_digest + result_digest
        )
        and result_digest == window_digest
    )


def _sparameter_complex_impedance_renormalization_is_current(
    raw: Mapping[str, Any],
) -> bool:
    identity = raw.get("sparameter_complex_impedance_renormalization_identity")
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        source_impedance = complex(
            float(identity["source_reference_impedance_real_ohm"]),
            float(identity["source_reference_impedance_imag_ohm"]),
        )
        target_impedance = complex(
            float(identity["target_reference_impedance_real_ohm"]),
            float(identity["target_reference_impedance_imag_ohm"]),
        )
    except (KeyError, TypeError, ValueError):
        return False
    calibration_generation = str(
        identity.get("port_calibration_generation", "")
    ).strip()
    transform_generation = str(
        identity.get("renormalization_transform_generation", "")
    ).strip()
    transform_digest = str(
        identity.get("renormalization_transform_sha256", "")
    ).lower()
    result_digest = str(
        identity.get("result_renormalization_transform_sha256", "")
    ).lower()
    return (
        bool(calibration_generation)
        and identity.get("sparameter_result_port_calibration_generation")
        == calibration_generation
        and identity.get("renormalization_port_calibration_generation")
        == calibration_generation
        and identity.get("source_reference_impedance_basis") == "complex_ohm"
        and identity.get("renormalization_reference_impedance_basis")
        == "complex_ohm"
        and all(
            math.isfinite(value)
            for value in (
                source_impedance.real,
                source_impedance.imag,
                target_impedance.real,
                target_impedance.imag,
            )
        )
        and source_impedance.real > 0.0
        and target_impedance.real > 0.0
        and bool(transform_generation)
        and identity.get("result_renormalization_transform_generation")
        == transform_generation
        and len(transform_digest) == len(result_digest) == 64
        and all(
            character in "0123456789abcdef"
            for character in transform_digest + result_digest
        )
        and result_digest == transform_digest
    )


def _farfield_polarization_basis_transform_is_current(
    raw: Mapping[str, Any],
) -> bool:
    identity = raw.get("farfield_polarization_basis_transform_identity")
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    solve_generation = str(identity.get("farfield_solve_generation", "")).strip()
    transform_generation = str(
        identity.get("basis_transform_generation", "")
    ).strip()
    transform_digest = str(identity.get("basis_transform_sha256", "")).lower()
    result_digest = str(
        identity.get("comparison_basis_transform_sha256", "")
    ).lower()
    return (
        bool(solve_generation)
        and identity.get("farfield_result_generation") == solve_generation
        and identity.get("source_polarization_basis") == "spherical_theta_phi"
        and identity.get("comparison_polarization_basis")
        in {"ludwig_1", "ludwig_2", "ludwig_3"}
        and identity.get("basis_transform_applied") is True
        and bool(transform_generation)
        and identity.get("comparison_transform_generation") == transform_generation
        and identity.get("angular_coordinate_frame") == "global_spherical"
        and identity.get("comparison_angular_coordinate_frame")
        == identity.get("angular_coordinate_frame")
        and len(transform_digest) == len(result_digest) == 64
        and all(
            character in "0123456789abcdef"
            for character in transform_digest + result_digest
        )
        and result_digest == transform_digest
    )


def _mixed_mode_sparameter_pair_order_is_current(raw: Mapping[str, Any]) -> bool:
    identity = raw.get("mixed_mode_sparameter_port_pair_order_identity")
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    basis_identity = raw.get("mixed_mode_sparameter_basis_identity")
    if not isinstance(basis_identity, Mapping):
        return False
    calibration_generation = str(
        identity.get("port_calibration_generation", "")
    ).strip()
    port_order = identity.get("single_ended_port_order")
    pair_order = identity.get("mixed_mode_pair_order")
    pair_polarity = identity.get("mixed_mode_pair_polarity")
    pairing_digest = str(identity.get("pairing_sha256", "")).lower()
    basis_digest = str(identity.get("basis_matrix_sha256", "")).lower()
    paired_ports = (
        [port for pair in pair_order for port in pair]
        if isinstance(pair_order, list)
        and all(isinstance(pair, list) for pair in pair_order)
        else []
    )
    return (
        bool(calibration_generation)
        and identity.get("single_ended_result_port_calibration_generation")
        == calibration_generation
        and identity.get("mixed_mode_pairing_port_calibration_generation")
        == calibration_generation
        and isinstance(port_order, list)
        and len(port_order) >= 2
        and all(
            isinstance(port, str) and bool(port.strip()) for port in port_order
        )
        and len(set(port_order)) == len(port_order)
        and port_order == basis_identity.get("single_ended_port_order")
        and isinstance(pair_order, list)
        and pair_order
        and all(
            isinstance(pair, list)
            and len(pair) == 2
            and all(isinstance(port, str) and bool(port.strip()) for port in pair)
            for pair in pair_order
        )
        and len(paired_ports) == len(port_order)
        and len(set(paired_ports)) == len(paired_ports)
        and set(paired_ports) == set(port_order)
        and isinstance(pair_polarity, list)
        and len(pair_polarity) == len(pair_order)
        and all(type(value) is int and value in (-1, 1) for value in pair_polarity)
        and identity.get("transform_pair_order") == pair_order
        and identity.get("transform_pair_polarity") == pair_polarity
        and len(pairing_digest) == 64
        and all(character in "0123456789abcdef" for character in pairing_digest)
        and str(identity.get("transform_pairing_sha256", "")).lower()
        == pairing_digest
        and len(basis_digest) == 64
        and all(character in "0123456789abcdef" for character in basis_digest)
        and basis_digest
        == str(basis_identity.get("basis_matrix_sha256", "")).lower()
        and str(identity.get("transform_basis_matrix_sha256", "")).lower()
        == basis_digest
    )


def _nearfield_farfield_phase_center_frame_is_current(
    raw: Mapping[str, Any],
) -> bool:
    identity = raw.get("nearfield_farfield_phase_center_coordinate_frame_identity")
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    farfield_identity = raw.get("farfield_polarization_basis_transform_identity")
    if not isinstance(farfield_identity, Mapping):
        return False
    nearfield_generation = str(
        identity.get("nearfield_result_generation", "")
    ).strip()
    frame_generation = str(
        identity.get("phase_center_coordinate_frame_generation", "")
    ).strip()
    farfield_generation = str(
        identity.get("farfield_result_generation", "")
    ).strip()
    coordinates = identity.get("phase_center_coordinates_m")
    transformed_coordinates = identity.get("farfield_phase_center_coordinates_m")
    phase_digest = str(identity.get("phase_center_sha256", "")).lower()
    try:
        coordinate_values = [float(value) for value in coordinates]
        transformed_values = [float(value) for value in transformed_coordinates]
    except (TypeError, ValueError):
        coordinate_values = []
        transformed_values = []
    return (
        bool(nearfield_generation)
        and identity.get("farfield_transform_nearfield_generation")
        == nearfield_generation
        and bool(frame_generation)
        and identity.get("farfield_phase_center_frame_generation")
        == frame_generation
        and bool(farfield_generation)
        and identity.get("phase_center_farfield_result_generation")
        == farfield_generation
        and farfield_generation
        == farfield_identity.get("farfield_result_generation")
        and identity.get("phase_center_coordinate_frame") == "global_cartesian"
        and identity.get("farfield_phase_center_coordinate_frame")
        == identity.get("phase_center_coordinate_frame")
        and identity.get("phase_center_coordinate_unit") == "m"
        and identity.get("farfield_phase_center_coordinate_unit") == "m"
        and len(coordinate_values) == 3
        and all(math.isfinite(value) for value in coordinate_values)
        and len(transformed_values) == 3
        and all(math.isfinite(value) for value in transformed_values)
        and all(
            math.isclose(source, result, rel_tol=1.0e-12, abs_tol=1.0e-12)
            for source, result in zip(coordinate_values, transformed_values)
        )
        and len(phase_digest) == 64
        and all(character in "0123456789abcdef" for character in phase_digest)
        and str(identity.get("farfield_phase_center_sha256", "")).lower()
        == phase_digest
    )


def _sparameter_deembed_reference_plane_map_is_current(
    raw: Mapping[str, Any],
) -> bool:
    identity = raw.get(
        "sparameter_deembed_reference_plane_per_port_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    sparameter_generation = str(
        identity.get("sparameter_generation", "")
    ).strip()
    port_generation = str(identity.get("port_generation", "")).strip()
    port_ids = identity.get("port_ids")
    offsets = identity.get("reference_plane_offsets_m")
    applied_offsets = identity.get("applied_reference_plane_offsets_m")
    map_digest = str(identity.get("reference_plane_map_sha256", "")).lower()
    try:
        offset_values = [float(value) for value in offsets]
        applied_values = [float(value) for value in applied_offsets]
    except (TypeError, ValueError):
        offset_values = []
        applied_values = []
    return (
        bool(sparameter_generation)
        and identity.get("deembedded_result_sparameter_generation")
        == sparameter_generation
        and bool(port_generation)
        and identity.get("reference_plane_port_generation") == port_generation
        and identity.get("deembedded_result_port_generation") == port_generation
        and isinstance(port_ids, list)
        and bool(port_ids)
        and all(isinstance(value, str) and bool(value) for value in port_ids)
        and len(set(port_ids)) == len(port_ids)
        and identity.get("reference_plane_port_ids") == port_ids
        and len(offset_values) == len(port_ids)
        and all(math.isfinite(value) for value in offset_values)
        and applied_values == offset_values
        and len(map_digest) == 64
        and all(character in "0123456789abcdef" for character in map_digest)
        and str(
            identity.get("deembedded_reference_plane_map_sha256", "")
        ).lower()
        == map_digest
    )


def _time_domain_port_signal_gate_window_is_current(
    raw: Mapping[str, Any],
) -> bool:
    identity = raw.get("time_domain_port_signal_gate_window_generation_identity")
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    signal_generation = str(identity.get("signal_generation", "")).strip()
    gate_generation = str(identity.get("gate_generation", "")).strip()
    try:
        signal_start = int(identity.get("signal_sample_start"))
        signal_end = int(identity.get("signal_sample_end"))
        gate_start = int(identity.get("gate_window_start_sample"))
        gate_end = int(identity.get("gate_window_end_sample"))
        transform_window = [
            int(value) for value in identity.get("transform_gate_window", [])
        ]
    except (TypeError, ValueError):
        signal_start = signal_end = gate_start = gate_end = -1
        transform_window = []
    normalization_basis = str(identity.get("normalization_basis", ""))
    gate_digest = str(identity.get("gate_window_sha256", "")).lower()
    return (
        bool(signal_generation)
        and identity.get("gate_window_signal_generation") == signal_generation
        and identity.get("transform_signal_generation") == signal_generation
        and bool(gate_generation)
        and identity.get("transform_gate_generation") == gate_generation
        and 0 <= signal_start <= gate_start < gate_end <= signal_end
        and transform_window == [gate_start, gate_end]
        and normalization_basis
        in {"incident_wave_peak", "incident_wave_energy", "unit_impulse"}
        and identity.get("transform_normalization_basis") == normalization_basis
        and len(gate_digest) == 64
        and all(character in "0123456789abcdef" for character in gate_digest)
        and str(identity.get("transform_gate_window_sha256", "")).lower()
        == gate_digest
    )


def _sparameter_renormalization_reference_impedance_is_current(
    raw: Mapping[str, Any],
) -> bool:
    identity = raw.get(
        "sparameter_port_renormalization_reference_impedance_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    sparameter_generation = str(
        identity.get("sparameter_generation", "")
    ).strip()
    calibration_generation = str(
        identity.get("port_calibration_generation", "")
    ).strip()
    port_ids = identity.get("port_ids")
    impedances = identity.get("reference_impedances_ohm")
    applied_impedances = identity.get("applied_reference_impedances_ohm")
    grid_digest = str(identity.get("frequency_grid_sha256", "")).lower()
    map_digest = str(
        identity.get("reference_impedance_map_sha256", "")
    ).lower()
    try:
        impedance_values = [float(value) for value in impedances]
        applied_values = [float(value) for value in applied_impedances]
    except (TypeError, ValueError):
        impedance_values = []
        applied_values = []
    return (
        bool(sparameter_generation)
        and identity.get("renormalized_result_sparameter_generation")
        == sparameter_generation
        and bool(calibration_generation)
        and identity.get(
            "reference_impedance_port_calibration_generation"
        )
        == calibration_generation
        and identity.get("renormalization_port_calibration_generation")
        == calibration_generation
        and identity.get(
            "renormalized_result_port_calibration_generation"
        )
        == calibration_generation
        and isinstance(port_ids, list)
        and bool(port_ids)
        and all(isinstance(value, str) and bool(value) for value in port_ids)
        and len(set(port_ids)) == len(port_ids)
        and identity.get("reference_impedance_port_ids") == port_ids
        and len(impedance_values) == len(port_ids)
        and all(math.isfinite(value) and value > 0.0 for value in impedance_values)
        and applied_values == impedance_values
        and len(grid_digest) == 64
        and all(character in "0123456789abcdef" for character in grid_digest)
        and str(identity.get("renormalized_frequency_grid_sha256", "")).lower()
        == grid_digest
        and len(map_digest) == 64
        and all(character in "0123456789abcdef" for character in map_digest)
        and str(
            identity.get("renormalized_reference_impedance_map_sha256", "")
        ).lower()
        == map_digest
    )


def _realized_gain_excitation_and_accepted_power_are_current(
    raw: Mapping[str, Any],
) -> bool:
    identity = raw.get(
        "realized_gain_accepted_power_port_excitation_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    gain_generation = str(
        identity.get("realized_gain_generation", "")
    ).strip()
    excitation_generation = str(
        identity.get("excitation_generation", "")
    ).strip()
    port_ids = identity.get("port_ids")
    accepted_power = identity.get("accepted_power_w")
    coefficients = identity.get("excitation_coefficients_re_im")
    table_digest = str(identity.get("excitation_table_sha256", "")).lower()
    try:
        power_values = [float(value) for value in accepted_power]
        coefficient_values = [
            [float(component) for component in value]
            for value in coefficients
        ]
    except (TypeError, ValueError):
        power_values = []
        coefficient_values = []
    return (
        bool(gain_generation)
        and identity.get("result_realized_gain_generation") == gain_generation
        and bool(excitation_generation)
        and identity.get("accepted_power_excitation_generation")
        == excitation_generation
        and identity.get("port_coefficient_excitation_generation")
        == excitation_generation
        and identity.get("realized_gain_excitation_generation")
        == excitation_generation
        and isinstance(port_ids, list)
        and bool(port_ids)
        and all(isinstance(value, str) and bool(value) for value in port_ids)
        and len(set(port_ids)) == len(port_ids)
        and identity.get("accepted_power_port_ids") == port_ids
        and identity.get("excitation_coefficient_port_ids") == port_ids
        and len(power_values) == len(port_ids)
        and all(math.isfinite(value) and value >= 0.0 for value in power_values)
        and any(value > 0.0 for value in power_values)
        and identity.get("realized_gain_accepted_power_w") == accepted_power
        and len(coefficient_values) == len(port_ids)
        and all(
            len(value) == 2 and all(math.isfinite(component) for component in value)
            for value in coefficient_values
        )
        and identity.get("realized_gain_excitation_coefficients_re_im")
        == coefficients
        and identity.get("accepted_power_unit") == "W"
        and identity.get("gain_unit") == "dBi"
        and len(table_digest) == 64
        and all(character in "0123456789abcdef" for character in table_digest)
        and str(identity.get("realized_gain_excitation_table_sha256", "")).lower()
        == table_digest
    )


def _farfield_polarization_phase_center_is_current(
    raw: Mapping[str, Any],
) -> bool:
    identity = raw.get(
        "farfield_polarization_basis_phase_center_coordinate_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    farfield_generation = str(identity.get("farfield_generation", "")).strip()
    coordinate_generation = str(
        identity.get("monitor_coordinate_generation", "")
    ).strip()
    sample_ids = identity.get("sample_ids")
    theta_digests = identity.get("theta_basis_sha256")
    phi_digests = identity.get("phi_basis_sha256")
    phase_center = identity.get("phase_center_m")
    try:
        phase_center_values = [float(value) for value in phase_center]
    except (TypeError, ValueError):
        phase_center_values = []

    def valid_digest_list(values: Any) -> bool:
        return (
            isinstance(values, list)
            and bool(values)
            and all(
                isinstance(value, str)
                and len(value) == 64
                and all(character in "0123456789abcdef" for character in value.lower())
                for value in values
            )
        )

    table_digest = str(identity.get("polarization_table_sha256", "")).lower()
    return (
        bool(farfield_generation)
        and identity.get("result_farfield_generation") == farfield_generation
        and bool(coordinate_generation)
        and identity.get("theta_basis_monitor_coordinate_generation")
        == coordinate_generation
        and identity.get("phi_basis_monitor_coordinate_generation")
        == coordinate_generation
        and identity.get("phase_center_monitor_coordinate_generation")
        == coordinate_generation
        and identity.get("result_monitor_coordinate_generation")
        == coordinate_generation
        and isinstance(sample_ids, list)
        and bool(sample_ids)
        and all(isinstance(value, int) and not isinstance(value, bool) for value in sample_ids)
        and len(set(sample_ids)) == len(sample_ids)
        and identity.get("result_sample_ids") == sample_ids
        and valid_digest_list(theta_digests)
        and len(theta_digests) == len(sample_ids)
        and identity.get("result_theta_basis_sha256") == theta_digests
        and valid_digest_list(phi_digests)
        and len(phi_digests) == len(sample_ids)
        and identity.get("result_phi_basis_sha256") == phi_digests
        and len(phase_center_values) == 3
        and all(math.isfinite(value) for value in phase_center_values)
        and identity.get("result_phase_center_m") == phase_center
        and len(table_digest) == 64
        and all(character in "0123456789abcdef" for character in table_digest)
        and str(identity.get("result_polarization_table_sha256", "")).lower()
        == table_digest
    )


def _broadband_energy_q_inputs_are_current(raw: Mapping[str, Any]) -> bool:
    identity = raw.get(
        "broadband_energy_q_port_loss_normalization_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    analysis_generation = str(identity.get("analysis_generation", "")).strip()
    grid_generation = str(
        identity.get("frequency_grid_generation", "")
    ).strip()
    excitation_generation = str(
        identity.get("excitation_generation", "")
    ).strip()
    series_keys = (
        "frequencies_hz",
        "stored_energy_j",
        "port_power_w",
        "loss_power_w",
    )
    try:
        series = {
            key: [float(value) for value in identity.get(key)] for key in series_keys
        }
    except (TypeError, ValueError):
        series = {key: [] for key in series_keys}
    frequencies = series["frequencies_hz"]
    digest = str(identity.get("energy_q_input_sha256", "")).lower()
    return (
        bool(analysis_generation)
        and all(
            identity.get(key) == analysis_generation
            for key in (
                "energy_analysis_generation",
                "port_power_analysis_generation",
                "loss_analysis_generation",
                "q_result_analysis_generation",
            )
        )
        and bool(grid_generation)
        and all(
            identity.get(key) == grid_generation
            for key in (
                "energy_frequency_grid_generation",
                "port_power_frequency_grid_generation",
                "loss_frequency_grid_generation",
                "q_frequency_grid_generation",
            )
        )
        and bool(excitation_generation)
        and all(
            identity.get(key) == excitation_generation
            for key in (
                "energy_excitation_generation",
                "port_power_excitation_generation",
                "loss_excitation_generation",
            )
        )
        and bool(frequencies)
        and all(math.isfinite(value) and value > 0.0 for value in frequencies)
        and all(left < right for left, right in zip(frequencies, frequencies[1:]))
        and all(
            len(values) == len(frequencies)
            and all(math.isfinite(value) and value >= 0.0 for value in values)
            for key, values in series.items()
            if key != "frequencies_hz"
        )
        and identity.get("q_frequencies_hz") == identity.get("frequencies_hz")
        and identity.get("q_stored_energy_j") == identity.get("stored_energy_j")
        and identity.get("q_port_power_w") == identity.get("port_power_w")
        and identity.get("q_loss_power_w") == identity.get("loss_power_w")
        and len(digest) == 64
        and all(character in "0123456789abcdef" for character in digest)
        and str(identity.get("result_energy_q_input_sha256", "")).lower()
        == digest
    )


def _mixed_mode_port_metadata_is_current(raw: Mapping[str, Any]) -> bool:
    identity = raw.get(
        "mixed_mode_pair_impedance_reference_plane_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    result_generation = str(identity.get("result_generation", "")).strip()
    port_generation = str(identity.get("port_generation", "")).strip()
    pair_ids = identity.get("pair_ids")
    polarities = identity.get("pair_polarities")
    try:
        differential = [
            float(value) for value in identity.get("differential_impedances_ohm", [])
        ]
        common = [
            float(value) for value in identity.get("common_impedances_ohm", [])
        ]
        planes = [float(value) for value in identity.get("reference_planes_m", [])]
    except (TypeError, ValueError):
        differential = common = planes = []
    digest = str(identity.get("mixed_mode_port_table_sha256", "")).lower()
    flattened_ports = (
        [port for pair in pair_ids for port in pair]
        if isinstance(pair_ids, list)
        and all(isinstance(pair, list) and len(pair) == 2 for pair in pair_ids)
        else []
    )
    return (
        bool(result_generation)
        and identity.get("decoded_result_generation") == result_generation
        and bool(port_generation)
        and all(
            identity.get(key) == port_generation
            for key in (
                "pair_map_port_generation",
                "modal_impedance_port_generation",
                "polarity_port_generation",
                "reference_plane_port_generation",
            )
        )
        and bool(flattened_ports)
        and all(isinstance(value, str) and bool(value.strip()) for value in flattened_ports)
        and len(set(flattened_ports)) == len(flattened_ports)
        and identity.get("result_pair_ids") == pair_ids
        and isinstance(polarities, list)
        and len(polarities) == len(pair_ids)
        and all(
            isinstance(pair, list)
            and len(pair) == 2
            and set(pair) == {-1, 1}
            for pair in polarities
        )
        and identity.get("result_pair_polarities") == polarities
        and len(differential) == len(pair_ids)
        and all(math.isfinite(value) and value > 0.0 for value in differential)
        and identity.get("result_differential_impedances_ohm")
        == identity.get("differential_impedances_ohm")
        and len(common) == len(pair_ids)
        and all(math.isfinite(value) and value > 0.0 for value in common)
        and identity.get("result_common_impedances_ohm")
        == identity.get("common_impedances_ohm")
        and len(planes) == len(pair_ids)
        and all(math.isfinite(value) for value in planes)
        and identity.get("result_reference_planes_m")
        == identity.get("reference_planes_m")
        and len(digest) == 64
        and all(character in "0123456789abcdef" for character in digest)
        and str(identity.get("result_mixed_mode_port_table_sha256", "")).lower()
        == digest
    )


def _time_farfield_fft_inputs_are_current(raw: Mapping[str, Any]) -> bool:
    identity = raw.get(
        "time_farfield_fft_window_phase_center_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    farfield_generation = str(identity.get("farfield_generation", "")).strip()
    monitor_generation = str(identity.get("monitor_generation", "")).strip()
    try:
        times = [float(value) for value in identity.get("time_samples_s", [])]
        fft_times = [
            float(value) for value in identity.get("fft_time_samples_s", [])
        ]
        window = [float(value) for value in identity.get("window_samples", [])]
        fft_window = [
            float(value) for value in identity.get("fft_window_samples", [])
        ]
        phase_center = [float(value) for value in identity.get("phase_center_m", [])]
        result_phase_center = [
            float(value) for value in identity.get("result_phase_center_m", [])
        ]
    except (TypeError, ValueError):
        times = fft_times = window = fft_window = []
        phase_center = result_phase_center = []
    scaling = identity.get("fft_scaling")
    digest = str(identity.get("time_farfield_input_sha256", "")).lower()
    return (
        bool(farfield_generation)
        and identity.get("result_farfield_generation") == farfield_generation
        and bool(monitor_generation)
        and all(
            identity.get(key) == monitor_generation
            for key in (
                "time_grid_monitor_generation",
                "window_monitor_generation",
                "fft_scaling_monitor_generation",
                "phase_center_monitor_generation",
            )
        )
        and len(times) >= 4
        and all(math.isfinite(value) for value in times)
        and all(left < right for left, right in zip(times, times[1:]))
        and fft_times == times
        and len(window) == len(times)
        and all(math.isfinite(value) and value >= 0.0 for value in window)
        and any(value > 0.0 for value in window)
        and fft_window == window
        and scaling in {"one_sided_amplitude", "two_sided_amplitude", "power"}
        and identity.get("result_fft_scaling") == scaling
        and len(phase_center) == 3
        and all(math.isfinite(value) for value in phase_center)
        and result_phase_center == phase_center
        and len(digest) == 64
        and all(character in "0123456789abcdef" for character in digest)
        and str(identity.get("result_time_farfield_input_sha256", "")).lower()
        == digest
    )


def _waveguide_degenerate_mode_tracking_is_current(raw: Mapping[str, Any]) -> bool:
    identity = raw.get(
        "waveguide_degenerate_mode_phase_order_overlap_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    generation = str(identity.get("sweep_generation", "")).strip()
    mode_ids = identity.get("mode_ids")
    modal_order = identity.get("modal_order")
    try:
        phases = [
            float(value) for value in identity.get("phase_reference_deg", [])
        ]
        result_phases = [
            float(value)
            for value in identity.get("result_phase_reference_deg", [])
        ]
        overlaps = [
            [float(value) for value in row]
            for row in identity.get("overlap_vectors", [])
        ]
        result_overlaps = [
            [float(value) for value in row]
            for row in identity.get("result_overlap_vectors", [])
        ]
    except (TypeError, ValueError):
        phases = result_phases = []
        overlaps = result_overlaps = []
    digest = str(identity.get("mode_tracking_table_sha256", "")).lower()
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "mesh_sweep_generation",
                "phase_sweep_generation",
                "modal_order_sweep_generation",
                "overlap_sweep_generation",
                "result_sweep_generation",
            )
        )
        and isinstance(mode_ids, list)
        and len(mode_ids) >= 2
        and all(isinstance(value, str) and bool(value.strip()) for value in mode_ids)
        and len(set(mode_ids)) == len(mode_ids)
        and identity.get("result_mode_ids") == mode_ids
        and isinstance(modal_order, list)
        and len(modal_order) == len(mode_ids)
        and all(
            isinstance(value, int) and not isinstance(value, bool) and value > 0
            for value in modal_order
        )
        and len(set(modal_order)) == len(modal_order)
        and identity.get("result_modal_order") == modal_order
        and len(phases) == len(mode_ids)
        and all(math.isfinite(value) for value in phases)
        and result_phases == phases
        and len(overlaps) == len(mode_ids)
        and all(
            len(row) == len(mode_ids)
            and all(math.isfinite(value) for value in row)
            and any(value != 0.0 for value in row)
            for row in overlaps
        )
        and result_overlaps == overlaps
        and len(digest) == 64
        and all(character in "0123456789abcdef" for character in digest)
        and str(identity.get("result_mode_tracking_table_sha256", "")).lower()
        == digest
    )


def _dispersive_causal_pole_fit_is_current(raw: Mapping[str, Any]) -> bool:
    identity = raw.get(
        "dispersive_causal_pole_fit_temperature_unit_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    generation = str(identity.get("fit_generation", "")).strip()
    convention = str(identity.get("causal_convention", "")).strip()
    unit = str(identity.get("frequency_unit", "")).strip()
    try:
        temperature = float(identity.get("temperature_c"))
        result_temperature = float(identity.get("result_temperature_c"))
        poles = [
            [float(value) for value in row]
            for row in identity.get("pole_pairs_rad_per_s", [])
        ]
        result_poles = [
            [float(value) for value in row]
            for row in identity.get("result_pole_pairs_rad_per_s", [])
        ]
        residues = [
            [float(value) for value in row]
            for row in identity.get("residues", [])
        ]
        result_residues = [
            [float(value) for value in row]
            for row in identity.get("result_residues", [])
        ]
    except (TypeError, ValueError):
        temperature = result_temperature = math.nan
        poles = result_poles = []
        residues = result_residues = []
    digest = str(identity.get("pole_fit_sha256", "")).lower()
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "pole_fit_generation",
                "causal_convention_fit_generation",
                "temperature_fit_generation",
                "frequency_unit_fit_generation",
                "field_result_fit_generation",
            )
        )
        and convention in {"exp(-iwt)", "exp(+iwt)"}
        and identity.get("result_causal_convention") == convention
        and math.isfinite(temperature)
        and result_temperature == temperature
        and unit in {"Hz", "kHz", "MHz", "GHz"}
        and identity.get("result_frequency_unit") == unit
        and len(poles) >= 2
        and all(
            len(row) == 2
            and all(math.isfinite(value) for value in row)
            and row[0] < 0.0
            for row in poles
        )
        and result_poles == poles
        and len(residues) == len(poles)
        and all(
            len(row) == 2 and all(math.isfinite(value) for value in row)
            for row in residues
        )
        and result_residues == residues
        and len(digest) == 64
        and all(character in "0123456789abcdef" for character in digest)
        and str(identity.get("result_pole_fit_sha256", "")).lower() == digest
    )


def _broadband_sparameter_inputs_are_current(raw: Mapping[str, Any]) -> bool:
    identity = raw.get(
        "broadband_adaptive_mesh_sparam_renormalization_port_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    generation = str(identity.get("sweep_generation", "")).strip()
    modes = identity.get("port_mode_ids")
    method = str(identity.get("frequency_interpolation", "")).strip()
    try:
        frequencies = [
            float(value) for value in identity.get("frequency_samples_hz", [])
        ]
        result_frequencies = [
            float(value)
            for value in identity.get("result_frequency_samples_hz", [])
        ]
        impedances = [
            [float(value) for value in row]
            for row in identity.get("renormalization_impedance_ohm", [])
        ]
        result_impedances = [
            [float(value) for value in row]
            for row in identity.get("result_renormalization_impedance_ohm", [])
        ]
    except (TypeError, ValueError):
        return False
    mesh_digest = str(identity.get("adaptive_mesh_sha256", "")).lower()
    table_digest = str(identity.get("sparameter_table_sha256", "")).lower()
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "adaptive_mesh_sweep_generation",
                "frequency_interpolation_sweep_generation",
                "port_mode_sweep_generation",
                "renormalization_sweep_generation",
                "sparameter_result_sweep_generation",
            )
        )
        and _valid_sha256(mesh_digest)
        and str(identity.get("result_adaptive_mesh_sha256", "")).lower()
        == mesh_digest
        and len(frequencies) >= 3
        and all(math.isfinite(value) and value > 0.0 for value in frequencies)
        and all(left < right for left, right in zip(frequencies, frequencies[1:]))
        and result_frequencies == frequencies
        and method in {"linear", "vector_fitting", "rational"}
        and identity.get("result_frequency_interpolation") == method
        and isinstance(modes, list)
        and bool(modes)
        and all(isinstance(value, str) and bool(value.strip()) for value in modes)
        and len(set(modes)) == len(modes)
        and identity.get("result_port_mode_ids") == modes
        and len(impedances) == len(modes)
        and all(
            len(row) == 2
            and all(math.isfinite(value) for value in row)
            and row[0] > 0.0
            for row in impedances
        )
        and result_impedances == impedances
        and _valid_sha256(table_digest)
        and str(identity.get("result_sparameter_table_sha256", "")).lower()
        == table_digest
    )


def _transient_monitor_inputs_are_current(raw: Mapping[str, Any]) -> bool:
    identity = raw.get(
        "transient_monitor_time_origin_excitation_waveform_mesh_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    generation = str(identity.get("transient_generation", "")).strip()
    frame = str(identity.get("monitor_coordinate_frame", "")).strip()
    monitor_ids = identity.get("monitor_ids")
    try:
        time_origin = float(identity.get("time_origin_s"))
        result_time_origin = float(identity.get("result_time_origin_s"))
        times = [float(value) for value in identity.get("time_samples_s", [])]
        result_times = [
            float(value) for value in identity.get("result_time_samples_s", [])
        ]
    except (TypeError, ValueError):
        return False
    waveform_digest = str(identity.get("excitation_waveform_sha256", "")).lower()
    mesh_digest = str(identity.get("mesh_sha256", "")).lower()
    table_digest = str(identity.get("monitor_field_table_sha256", "")).lower()
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "time_origin_transient_generation",
                "excitation_waveform_transient_generation",
                "monitor_frame_transient_generation",
                "mesh_transient_generation",
                "field_result_transient_generation",
            )
        )
        and math.isfinite(time_origin)
        and result_time_origin == time_origin
        and _valid_sha256(waveform_digest)
        and str(identity.get("result_excitation_waveform_sha256", "")).lower()
        == waveform_digest
        and frame in {"global_xyz", "port_local"}
        and identity.get("result_monitor_coordinate_frame") == frame
        and isinstance(monitor_ids, list)
        and bool(monitor_ids)
        and all(
            isinstance(value, int) and not isinstance(value, bool) and value > 0
            for value in monitor_ids
        )
        and len(set(monitor_ids)) == len(monitor_ids)
        and identity.get("result_monitor_ids") == monitor_ids
        and _valid_sha256(mesh_digest)
        and str(identity.get("result_mesh_sha256", "")).lower() == mesh_digest
        and len(times) >= 3
        and times[0] == time_origin
        and all(math.isfinite(value) for value in times)
        and all(left < right for left, right in zip(times, times[1:]))
        and result_times == times
        and _valid_sha256(table_digest)
        and str(identity.get("result_monitor_field_table_sha256", "")).lower()
        == table_digest
    )


def _deembedded_network_inputs_are_current(raw: Mapping[str, Any]) -> bool:
    identity = raw.get(
        "deembedding_reference_plane_phase_causality_passivity_grid_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    generation = str(identity.get("deembedding_generation", "")).strip()
    modes = identity.get("port_mode_ids")
    try:
        offsets = [
            float(value) for value in identity.get("reference_plane_offsets_m", [])
        ]
        result_offsets = [
            float(value)
            for value in identity.get("result_reference_plane_offsets_m", [])
        ]
        frequencies = [
            float(value) for value in identity.get("frequency_grid_hz", [])
        ]
        result_frequencies = [
            float(value) for value in identity.get("result_frequency_grid_hz", [])
        ]
        phases = [float(value) for value in identity.get("unwrapped_phase_rad", [])]
        result_phases = [
            float(value) for value in identity.get("result_unwrapped_phase_rad", [])
        ]
        singular_values = [
            float(value)
            for value in identity.get("passivity_max_singular_values", [])
        ]
        result_singular_values = [
            float(value)
            for value in identity.get("result_passivity_max_singular_values", [])
        ]
    except (TypeError, ValueError):
        return False
    digest = str(identity.get("deembedded_network_sha256", "")).lower()
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "reference_plane_deembedding_generation",
                "phase_deembedding_generation",
                "causality_deembedding_generation",
                "passivity_deembedding_generation",
                "frequency_grid_deembedding_generation",
                "result_deembedding_generation",
            )
        )
        and isinstance(modes, list)
        and len(modes) >= 2
        and all(isinstance(value, str) and bool(value.strip()) for value in modes)
        and len(set(modes)) == len(modes)
        and identity.get("result_port_mode_ids") == modes
        and len(offsets) == len(modes)
        and all(math.isfinite(value) for value in offsets)
        and result_offsets == offsets
        and len(frequencies) >= 3
        and all(math.isfinite(value) and value > 0.0 for value in frequencies)
        and all(left < right for left, right in zip(frequencies, frequencies[1:]))
        and result_frequencies == frequencies
        and len(phases) == len(frequencies)
        and all(math.isfinite(value) for value in phases)
        and all(abs(right - left) <= math.pi for left, right in zip(phases, phases[1:]))
        and result_phases == phases
        and identity.get("causality_check_passed") is True
        and identity.get("result_causality_check_passed") is True
        and len(singular_values) == len(frequencies)
        and all(
            math.isfinite(value) and 0.0 <= value <= 1.0
            for value in singular_values
        )
        and result_singular_values == singular_values
        and _valid_sha256(digest)
        and str(identity.get("result_deembedded_network_sha256", "")).lower()
        == digest
    )


def _field_circuit_cosim_inputs_are_current(raw: Mapping[str, Any]) -> bool:
    identity = raw.get(
        "field_circuit_cosim_port_sign_impedance_power_balance_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    generation = str(identity.get("cosim_generation", "")).strip()
    port_id = str(identity.get("port_id", "")).strip()
    current_sign = str(identity.get("current_sign_convention", "")).strip()
    voltage_reference = str(identity.get("voltage_reference", "")).strip()
    convention = str(identity.get("phasor_amplitude_convention", "")).strip()
    try:
        voltage = [float(value) for value in identity.get("port_voltage_ri_v", [])]
        result_voltage = [
            float(value) for value in identity.get("result_port_voltage_ri_v", [])
        ]
        current = [float(value) for value in identity.get("port_current_ri_a", [])]
        result_current = [
            float(value) for value in identity.get("result_port_current_ri_a", [])
        ]
        impedance = [
            float(value) for value in identity.get("port_impedance_ri_ohm", [])
        ]
        result_impedance = [
            float(value)
            for value in identity.get("result_port_impedance_ri_ohm", [])
        ]
        field_power = float(identity.get("field_absorbed_power_w"))
        circuit_power = float(identity.get("circuit_delivered_power_w"))
        residual = float(identity.get("power_balance_residual_w"))
        result_residual = float(identity.get("result_power_balance_residual_w"))
    except (TypeError, ValueError):
        return False
    if len(voltage) != 2 or len(current) != 2 or len(impedance) != 2:
        return False
    complex_voltage = complex(*voltage)
    complex_current = complex(*current)
    if abs(complex_current) <= 1.0e-300:
        return False
    expected_impedance = complex_voltage / complex_current
    expected_power = (complex_voltage * complex_current.conjugate()).real
    digest = str(identity.get("cosim_result_sha256", "")).lower()
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "field_port_cosim_generation",
                "circuit_port_cosim_generation",
                "sign_cosim_generation",
                "impedance_cosim_generation",
                "power_balance_cosim_generation",
                "result_cosim_generation",
            )
        )
        and bool(port_id)
        and identity.get("result_port_id") == port_id
        and current_sign == "positive_into_field_port"
        and identity.get("result_current_sign_convention") == current_sign
        and voltage_reference == "positive_to_negative_terminal"
        and identity.get("result_voltage_reference") == voltage_reference
        and convention == "rms"
        and identity.get("result_phasor_amplitude_convention") == convention
        and all(math.isfinite(value) for value in voltage + current + impedance)
        and result_voltage == voltage
        and result_current == current
        and math.isclose(impedance[0], expected_impedance.real, rel_tol=1.0e-12)
        and math.isclose(impedance[1], expected_impedance.imag, rel_tol=1.0e-12)
        and result_impedance == impedance
        and math.isfinite(field_power)
        and math.isclose(field_power, expected_power, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isclose(circuit_power, field_power, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isclose(residual, circuit_power - field_power, rel_tol=0.0, abs_tol=1.0e-15)
        and math.isclose(result_residual, residual, rel_tol=0.0, abs_tol=1.0e-15)
        and _valid_sha256(digest)
        and str(identity.get("reported_cosim_result_sha256", "")).lower()
        == digest
    )


def _adaptive_mesh_convergence_inputs_are_current(raw: Mapping[str, Any]) -> bool:
    identity = raw.get(
        "adaptive_mesh_pass_sparameter_energy_convergence_grid_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    generation = str(identity.get("adaptive_generation", "")).strip()
    pass_ids = identity.get("mesh_pass_ids")
    try:
        cells = [float(value) for value in identity.get("mesh_cell_counts", [])]
        result_cells = [
            float(value) for value in identity.get("result_mesh_cell_counts", [])
        ]
        frequencies = [
            float(value) for value in identity.get("frequency_grid_hz", [])
        ]
        result_frequencies = [
            float(value) for value in identity.get("result_frequency_grid_hz", [])
        ]
        s_deltas = [
            float(value) for value in identity.get("maximum_sparameter_delta", [])
        ]
        result_s_deltas = [
            float(value)
            for value in identity.get("result_maximum_sparameter_delta", [])
        ]
        energy_residuals = [
            float(value)
            for value in identity.get("stored_energy_closure_residual", [])
        ]
        result_energy_residuals = [
            float(value)
            for value in identity.get("result_stored_energy_closure_residual", [])
        ]
        s_tolerance = float(identity.get("sparameter_delta_tolerance"))
        energy_tolerance = float(identity.get("energy_closure_tolerance"))
    except (TypeError, ValueError):
        return False
    digest = str(identity.get("adaptive_result_sha256", "")).lower()
    count = len(pass_ids) if isinstance(pass_ids, list) else 0
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "mesh_pass_adaptive_generation",
                "sparameter_adaptive_generation",
                "energy_adaptive_generation",
                "frequency_grid_adaptive_generation",
                "stopping_rule_adaptive_generation",
                "result_adaptive_generation",
            )
        )
        and count >= 2
        and pass_ids == list(range(count))
        and identity.get("result_mesh_pass_ids") == pass_ids
        and len(cells) == count
        and all(math.isfinite(value) and value > 0.0 for value in cells)
        and all(left < right for left, right in zip(cells, cells[1:]))
        and result_cells == cells
        and len(frequencies) >= 3
        and all(math.isfinite(value) and value > 0.0 for value in frequencies)
        and all(left < right for left, right in zip(frequencies, frequencies[1:]))
        and result_frequencies == frequencies
        and len(s_deltas) == count
        and all(math.isfinite(value) and value >= 0.0 for value in s_deltas)
        and all(left > right for left, right in zip(s_deltas, s_deltas[1:]))
        and result_s_deltas == s_deltas
        and len(energy_residuals) == count
        and all(
            math.isfinite(value) and value >= 0.0 for value in energy_residuals
        )
        and all(
            left > right
            for left, right in zip(energy_residuals, energy_residuals[1:])
        )
        and result_energy_residuals == energy_residuals
        and math.isfinite(s_tolerance)
        and s_tolerance > 0.0
        and math.isfinite(energy_tolerance)
        and energy_tolerance > 0.0
        and s_deltas[-1] <= s_tolerance
        and energy_residuals[-1] <= energy_tolerance
        and identity.get("converged_pass_id") == pass_ids[-1]
        and identity.get("result_converged_pass_id") == pass_ids[-1]
        and _valid_sha256(digest)
        and str(identity.get("reported_adaptive_result_sha256", "")).lower()
        == digest
    )


def _eigenmode_tracking_inputs_are_current(raw: Mapping[str, Any]) -> bool:
    identity = raw.get(
        "eigenmode_tracking_phase_normalization_port_coupling_mesh_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    generation = str(identity.get("tracking_generation", "")).strip()
    mode_ids = identity.get("tracked_mode_ids")
    subspaces = identity.get("modal_subspace_sha256")
    anchors = identity.get("phase_anchor_ids")
    meshes = identity.get("mesh_sha256")
    normalization = str(identity.get("normalization", "")).strip()
    try:
        sweep = [float(value) for value in identity.get("sweep_parameters", [])]
        result_sweep = [
            float(value) for value in identity.get("result_sweep_parameters", [])
        ]
        coupling = [
            [float(value) for value in row]
            for row in identity.get("port_coupling_magnitudes", [])
        ]
        result_coupling = [
            [float(value) for value in row]
            for row in identity.get("result_port_coupling_magnitudes", [])
        ]
    except (TypeError, ValueError):
        return False
    digest = str(identity.get("eigenmode_track_sha256", "")).lower()
    modes_ok = (
        isinstance(mode_ids, list)
        and len(mode_ids) >= 2
        and all(isinstance(value, str) and bool(value) for value in mode_ids)
        and len(set(mode_ids)) == len(mode_ids)
    )
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "modal_subspace_tracking_generation",
                "phase_tracking_generation",
                "normalization_tracking_generation",
                "port_coupling_tracking_generation",
                "mesh_tracking_generation",
                "result_tracking_generation",
            )
        )
        and len(sweep) >= 3
        and all(math.isfinite(value) for value in sweep)
        and all(left < right for left, right in zip(sweep, sweep[1:]))
        and result_sweep == sweep
        and modes_ok
        and identity.get("result_tracked_mode_ids") == mode_ids
        and isinstance(subspaces, list)
        and len(subspaces) == len(sweep)
        and all(_valid_sha256(value) for value in subspaces)
        and identity.get("result_modal_subspace_sha256") == subspaces
        and isinstance(anchors, list)
        and modes_ok
        and len(anchors) == len(mode_ids)
        and all(isinstance(value, str) and bool(value) for value in anchors)
        and identity.get("result_phase_anchor_ids") == anchors
        and normalization == "stored_energy_1j"
        and identity.get("result_normalization") == normalization
        and len(coupling) == len(sweep)
        and all(
            len(row) == len(mode_ids)
            and all(math.isfinite(value) and value >= 0.0 for value in row)
            for row in coupling
        )
        and result_coupling == coupling
        and isinstance(meshes, list)
        and len(meshes) == len(sweep)
        and all(_valid_sha256(value) for value in meshes)
        and identity.get("result_mesh_sha256") == meshes
        and _valid_sha256(digest)
        and str(identity.get("reported_eigenmode_track_sha256", "")).lower()
        == digest
    )


def _port_network_inputs_are_current(raw: Mapping[str, Any]) -> bool:
    identity = raw.get(
        "port_deembedding_reference_plane_impedance_mode_normalization_smatrix_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    generation = str(identity.get("network_generation", "")).strip()
    mode_ids = identity.get("port_mode_ids")
    smatrix = identity.get("smatrix_ri")
    try:
        offsets = [float(value) for value in identity.get("reference_plane_offsets_m", [])]
        result_offsets = [
            float(value)
            for value in identity.get("result_reference_plane_offsets_m", [])
        ]
        impedances = [
            [float(value) for value in pair]
            for pair in identity.get("reference_impedance_ri_ohm", [])
        ]
        result_impedances = [
            [float(value) for value in pair]
            for pair in identity.get("result_reference_impedance_ri_ohm", [])
        ]
        frequencies = [float(value) for value in identity.get("frequency_grid_hz", [])]
        result_frequencies = [
            float(value) for value in identity.get("result_frequency_grid_hz", [])
        ]
        numeric_smatrix = [
            [[float(value) for value in pair] for pair in row]
            for row in smatrix
        ]
        numeric_result_smatrix = [
            [[float(value) for value in pair] for pair in row]
            for row in identity.get("result_smatrix_ri", [])
        ]
    except (TypeError, ValueError):
        return False
    digest = str(identity.get("smatrix_sha256", "")).lower()
    mode_count = len(mode_ids) if isinstance(mode_ids, list) else 0
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "port_mode_network_generation",
                "deembedding_network_generation",
                "reference_impedance_network_generation",
                "normalization_network_generation",
                "frequency_grid_network_generation",
                "result_network_generation",
            )
        )
        and mode_count >= 2
        and all(isinstance(value, str) and bool(value) for value in mode_ids)
        and len(set(mode_ids)) == mode_count
        and identity.get("result_port_mode_ids") == mode_ids
        and len(offsets) == mode_count
        and all(math.isfinite(value) for value in offsets)
        and result_offsets == offsets
        and len(impedances) == mode_count
        and all(
            len(pair) == 2
            and all(math.isfinite(value) for value in pair)
            and pair[0] > 0.0
            for pair in impedances
        )
        and result_impedances == impedances
        and identity.get("wave_normalization") == "power_wave"
        and identity.get("result_wave_normalization") == "power_wave"
        and len(frequencies) >= 2
        and all(math.isfinite(value) and value > 0.0 for value in frequencies)
        and all(left < right for left, right in zip(frequencies, frequencies[1:]))
        and result_frequencies == frequencies
        and len(numeric_smatrix) == mode_count
        and all(
            len(row) == mode_count
            and all(
                len(pair) == 2 and all(math.isfinite(value) for value in pair)
                for pair in row
            )
            for row in numeric_smatrix
        )
        and numeric_result_smatrix == numeric_smatrix
        and _valid_sha256(digest)
        and str(identity.get("reported_smatrix_sha256", "")).lower() == digest
    )


def _farfield_result_inputs_are_current(raw: Mapping[str, Any]) -> bool:
    identity = raw.get(
        "farfield_angular_grid_polarization_coordinate_power_normalization_mesh_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    generation = str(identity.get("farfield_generation", "")).strip()
    try:
        theta = [float(value) for value in identity.get("theta_deg", [])]
        result_theta = [float(value) for value in identity.get("result_theta_deg", [])]
        phi = [float(value) for value in identity.get("phi_deg", [])]
        result_phi = [float(value) for value in identity.get("result_phi_deg", [])]
        radiated_power = float(identity.get("radiated_power_w"))
        result_radiated_power = float(identity.get("result_radiated_power_w"))
    except (TypeError, ValueError):
        return False
    mesh_digest = str(identity.get("mesh_sha256", "")).lower()
    result_digest = str(identity.get("farfield_sha256", "")).lower()
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "angular_grid_farfield_generation",
                "polarization_farfield_generation",
                "coordinate_farfield_generation",
                "power_farfield_generation",
                "mesh_farfield_generation",
                "result_farfield_generation",
            )
        )
        and len(theta) >= 2
        and all(math.isfinite(value) and 0.0 <= value <= 180.0 for value in theta)
        and all(left < right for left, right in zip(theta, theta[1:]))
        and result_theta == theta
        and len(phi) >= 2
        and all(math.isfinite(value) and 0.0 <= value < 360.0 for value in phi)
        and all(left < right for left, right in zip(phi, phi[1:]))
        and result_phi == phi
        and identity.get("polarization_basis") == "ludwig3_co_cross"
        and identity.get("result_polarization_basis")
        == identity.get("polarization_basis")
        and identity.get("coordinate_frame") == "global_xyz_z_up"
        and identity.get("result_coordinate_frame")
        == identity.get("coordinate_frame")
        and math.isfinite(radiated_power)
        and radiated_power > 0.0
        and math.isclose(result_radiated_power, radiated_power, rel_tol=1.0e-12)
        and identity.get("field_normalization") == "sqrt_radiated_power"
        and identity.get("result_field_normalization")
        == identity.get("field_normalization")
        and _valid_sha256(mesh_digest)
        and str(identity.get("result_mesh_sha256", "")).lower() == mesh_digest
        and _valid_sha256(result_digest)
        and str(identity.get("reported_farfield_sha256", "")).lower()
        == result_digest
    )


def _time_domain_port_smatrix_inputs_are_current(raw: Mapping[str, Any]) -> bool:
    identity = raw.get(
        "time_domain_port_waveform_normalization_fft_window_grid_deembedding_smatrix_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    generation = str(identity.get("time_domain_generation", "")).strip()
    mode_ids = identity.get("port_mode_ids")
    smatrix = identity.get("smatrix_ri")
    try:
        times = [float(value) for value in identity.get("time_grid_s", [])]
        result_times = [float(value) for value in identity.get("result_time_grid_s", [])]
        waveform = [float(value) for value in identity.get("incident_waveform", [])]
        result_waveform = [
            float(value) for value in identity.get("result_incident_waveform", [])
        ]
        impedances = [
            float(value) for value in identity.get("reference_impedance_ohm", [])
        ]
        result_impedances = [
            float(value)
            for value in identity.get("result_reference_impedance_ohm", [])
        ]
        frequencies = [float(value) for value in identity.get("frequency_grid_hz", [])]
        result_frequencies = [
            float(value) for value in identity.get("result_frequency_grid_hz", [])
        ]
        offsets = [float(value) for value in identity.get("deembedding_offsets_m", [])]
        result_offsets = [
            float(value) for value in identity.get("result_deembedding_offsets_m", [])
        ]
        numeric_smatrix = [
            [[float(value) for value in pair] for pair in row] for row in smatrix
        ]
        numeric_result_smatrix = [
            [[float(value) for value in pair] for pair in row]
            for row in identity.get("result_smatrix_ri", [])
        ]
    except (TypeError, ValueError):
        return False
    mode_count = len(mode_ids) if isinstance(mode_ids, list) else 0
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "waveform_time_domain_generation",
                "normalization_time_domain_generation",
                "fft_time_domain_generation",
                "grid_time_domain_generation",
                "deembedding_time_domain_generation",
                "smatrix_time_domain_generation",
                "result_time_domain_generation",
            )
        )
        and mode_count >= 2
        and all(isinstance(value, str) and value.strip() for value in mode_ids)
        and len(set(mode_ids)) == mode_count
        and identity.get("result_port_mode_ids") == mode_ids
        and len(times) >= 4
        and all(math.isfinite(value) and value >= 0.0 for value in times)
        and all(left < right for left, right in zip(times, times[1:]))
        and result_times == times
        and len(waveform) == len(times)
        and all(math.isfinite(value) for value in waveform)
        and any(value != 0.0 for value in waveform)
        and result_waveform == waveform
        and identity.get("wave_normalization") == "power-wave"
        and identity.get("result_wave_normalization") == "power-wave"
        and len(impedances) == mode_count
        and all(math.isfinite(value) and value > 0.0 for value in impedances)
        and result_impedances == impedances
        and identity.get("fft_window") == "tukey-0.2"
        and identity.get("result_fft_window") == "tukey-0.2"
        and len(frequencies) >= 2
        and all(math.isfinite(value) and value > 0.0 for value in frequencies)
        and all(left < right for left, right in zip(frequencies, frequencies[1:]))
        and result_frequencies == frequencies
        and len(offsets) == mode_count
        and all(math.isfinite(value) for value in offsets)
        and result_offsets == offsets
        and len(numeric_smatrix) == mode_count
        and all(
            len(row) == mode_count
            and all(
                len(pair) == 2 and all(math.isfinite(value) for value in pair)
                for pair in row
            )
            for row in numeric_smatrix
        )
        and numeric_result_smatrix == numeric_smatrix
        and _valid_sha256(identity.get("time_result_sha256"))
        and identity.get("accepted_time_result_sha256")
        == identity.get("time_result_sha256")
        and _valid_sha256(identity.get("smatrix_sha256"))
        and identity.get("accepted_smatrix_sha256") == identity.get("smatrix_sha256")
    )


def _huygens_near_far_inputs_are_current(raw: Mapping[str, Any]) -> bool:
    identity = raw.get(
        "huygens_box_orientation_phase_center_frequency_mesh_near_far_transform_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    generation = str(identity.get("huygens_generation", "")).strip()
    face_ids = identity.get("box_face_ids")
    try:
        signs = [int(value) for value in identity.get("outward_orientation_sign", [])]
        result_signs = [
            int(value) for value in identity.get("result_outward_orientation_sign", [])
        ]
        phase_center = [float(value) for value in identity.get("phase_center_m", [])]
        result_phase_center = [
            float(value) for value in identity.get("result_phase_center_m", [])
        ]
        frequency = float(identity.get("frequency_hz"))
        result_frequency = float(identity.get("result_frequency_hz"))
    except (TypeError, ValueError):
        return False
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "orientation_huygens_generation",
                "phase_center_huygens_generation",
                "frequency_huygens_generation",
                "mesh_huygens_generation",
                "transform_huygens_generation",
                "result_huygens_generation",
            )
        )
        and face_ids == ["-x", "+x", "-y", "+y", "-z", "+z"]
        and identity.get("result_box_face_ids") == face_ids
        and signs == [-1, 1, -1, 1, -1, 1]
        and result_signs == signs
        and len(phase_center) == 3
        and all(math.isfinite(value) for value in phase_center)
        and result_phase_center == phase_center
        and math.isfinite(frequency)
        and frequency > 0.0
        and result_frequency == frequency
        and identity.get("near_far_transform") == "equivalent-current-near-to-far"
        and identity.get("result_near_far_transform")
        == "equivalent-current-near-to-far"
        and identity.get("encloses_all_sources") is True
        and identity.get("result_encloses_all_sources") is True
        and _valid_sha256(identity.get("enclosing_mesh_sha256"))
        and identity.get("result_enclosing_mesh_sha256")
        == identity.get("enclosing_mesh_sha256")
        and _valid_sha256(identity.get("near_field_sha256"))
        and identity.get("accepted_near_field_sha256")
        == identity.get("near_field_sha256")
        and _valid_sha256(identity.get("far_field_sha256"))
        and identity.get("accepted_far_field_sha256")
        == identity.get("far_field_sha256")
    )


def _waveguide_port_mode_inputs_are_current(raw: Mapping[str, Any]) -> bool:
    identity = raw.get(
        "waveguide_port_mode_cutoff_normalization_reference_plane_mesh_field_result_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    generation = str(identity.get("port_generation", "")).strip()
    port_id = str(identity.get("port_id", "")).strip()
    mode_id = str(identity.get("mode_id", "")).strip()
    try:
        cutoff = float(identity.get("cutoff_frequency_hz"))
        result_cutoff = float(identity.get("result_cutoff_frequency_hz"))
        frequency = float(identity.get("evaluation_frequency_hz"))
        result_frequency = float(identity.get("result_evaluation_frequency_hz"))
        reference_plane = float(identity.get("reference_plane_m"))
        result_reference_plane = float(identity.get("result_reference_plane_m"))
    except (TypeError, ValueError):
        return False
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "mode_port_generation",
                "cutoff_port_generation",
                "normalization_port_generation",
                "reference_plane_port_generation",
                "mesh_port_generation",
                "field_port_generation",
                "result_port_generation",
            )
        )
        and bool(port_id)
        and identity.get("result_port_id") == port_id
        and bool(mode_id)
        and identity.get("result_mode_id") == mode_id
        and math.isfinite(cutoff)
        and cutoff > 0.0
        and math.isclose(result_cutoff, cutoff, rel_tol=1.0e-12)
        and math.isfinite(frequency)
        and frequency > cutoff
        and math.isclose(result_frequency, frequency, rel_tol=1.0e-12)
        and identity.get("normalization") == "unit-power-wave"
        and identity.get("result_normalization") == "unit-power-wave"
        and math.isfinite(reference_plane)
        and math.isclose(
            result_reference_plane, reference_plane, rel_tol=0.0, abs_tol=1.0e-15
        )
        and _valid_sha256(identity.get("port_mesh_sha256"))
        and identity.get("result_port_mesh_sha256")
        == identity.get("port_mesh_sha256")
        and _valid_sha256(identity.get("field_eigenvector_sha256"))
        and identity.get("result_field_eigenvector_sha256")
        == identity.get("field_eigenvector_sha256")
        and _valid_sha256(identity.get("result_sha256"))
        and identity.get("accepted_result_sha256") == identity.get("result_sha256")
    )


def _wake_impedance_inputs_are_current(raw: Mapping[str, Any]) -> bool:
    identity = raw.get(
        "wake_impedance_bunch_profile_time_grid_frequency_transform_normalization_mesh_result_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    generation = str(identity.get("wake_generation", "")).strip()
    try:
        sigma = float(identity.get("bunch_sigma_s"))
        result_sigma = float(identity.get("result_bunch_sigma_s"))
        charge = float(identity.get("bunch_charge_c"))
        result_charge = float(identity.get("result_bunch_charge_c"))
        times = [float(value) for value in identity.get("time_grid_s", [])]
        result_times = [
            float(value) for value in identity.get("result_time_grid_s", [])
        ]
        wake = [float(value) for value in identity.get("wake_potential_v_c", [])]
        result_wake = [
            float(value) for value in identity.get("result_wake_potential_v_c", [])
        ]
        frequencies = [
            float(value) for value in identity.get("frequency_grid_hz", [])
        ]
        result_frequencies = [
            float(value) for value in identity.get("result_frequency_grid_hz", [])
        ]
    except (TypeError, ValueError):
        return False
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "bunch_wake_generation",
                "time_wake_generation",
                "transform_wake_generation",
                "frequency_wake_generation",
                "normalization_wake_generation",
                "mesh_wake_generation",
                "result_wake_generation",
            )
        )
        and identity.get("bunch_profile") == "gaussian"
        and identity.get("result_bunch_profile") == "gaussian"
        and math.isfinite(sigma)
        and sigma > 0.0
        and math.isclose(result_sigma, sigma, rel_tol=0.0, abs_tol=1.0e-18)
        and math.isfinite(charge)
        and charge > 0.0
        and math.isclose(result_charge, charge, rel_tol=1.0e-12)
        and len(times) >= 4
        and len(wake) == len(times)
        and all(math.isfinite(value) for value in times + wake)
        and all(left < right for left, right in zip(times, times[1:]))
        and all(
            math.isclose(
                right - left,
                times[1] - times[0],
                rel_tol=1.0e-12,
                abs_tol=1.0e-18,
            )
            for left, right in zip(times, times[1:])
        )
        and result_times == times
        and result_wake == wake
        and identity.get("fft_convention") == "exp-minus-i-omega-t"
        and identity.get("result_fft_convention") == "exp-minus-i-omega-t"
        and len(frequencies) >= 2
        and all(math.isfinite(value) and value >= 0.0 for value in frequencies)
        and all(left < right for left, right in zip(frequencies, frequencies[1:]))
        and result_frequencies == frequencies
        and identity.get("impedance_normalization")
        == "longitudinal-v-per-coulomb"
        and identity.get("result_impedance_normalization")
        == "longitudinal-v-per-coulomb"
        and _valid_sha256(identity.get("mesh_sha256"))
        and identity.get("result_mesh_sha256") == identity.get("mesh_sha256")
        and _valid_sha256(identity.get("result_sha256"))
        and identity.get("accepted_result_sha256") == identity.get("result_sha256")
    )


def _dispersive_vector_fit_inputs_are_current(raw: Mapping[str, Any]) -> bool:
    identity = raw.get(
        "dispersive_vector_fit_passivity_causality_temperature_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    generation = str(identity.get("fit_generation", "")).strip()
    try:
        temperature = float(identity.get("temperature_c"))
        result_temperature = float(identity.get("result_temperature_c"))
        frequencies = [float(item) for item in identity.get("frequency_grid_hz", [])]
        result_frequencies = [float(item) for item in identity.get("result_frequency_grid_hz", [])]
        poles = [[float(item) for item in row] for row in identity.get("poles_rad_s", [])]
        result_poles = [[float(item) for item in row] for row in identity.get("result_poles_rad_s", [])]
        residues = [[float(item) for item in row] for row in identity.get("residues", [])]
        result_residues = [[float(item) for item in row] for row in identity.get("result_residues", [])]
        minimum_dissipation = float(identity.get("minimum_dissipation"))
        result_minimum_dissipation = float(identity.get("result_minimum_dissipation"))
        causality_residual = float(identity.get("causality_residual"))
        result_causality_residual = float(identity.get("result_causality_residual"))
        residual_limit = float(identity.get("causality_residual_limit"))
    except (TypeError, ValueError):
        return False
    return (
        bool(generation)
        and all(identity.get(key) == generation for key in (
            "pole_fit_generation", "residue_fit_generation", "passivity_fit_generation",
            "causality_fit_generation", "temperature_fit_generation", "frequency_fit_generation",
            "material_fit_generation", "result_fit_generation"))
        and math.isfinite(temperature)
        and math.isclose(result_temperature, temperature, rel_tol=0.0, abs_tol=1.0e-12)
        and len(frequencies) >= 3
        and all(math.isfinite(item) and item > 0.0 for item in frequencies)
        and all(right > left for left, right in zip(frequencies, frequencies[1:]))
        and result_frequencies == frequencies
        and len(poles) == len(residues) >= 1
        and all(len(row) == 2 and all(math.isfinite(item) for item in row) for row in poles + residues)
        and all(row[0] < 0.0 for row in poles)
        and result_poles == poles and result_residues == residues
        and identity.get("passivity_enforced") is True
        and identity.get("result_passivity_enforced") is True
        and math.isfinite(minimum_dissipation) and minimum_dissipation >= 0.0
        and result_minimum_dissipation == minimum_dissipation
        and math.isfinite(causality_residual) and causality_residual >= 0.0
        and math.isfinite(residual_limit) and residual_limit >= 0.0
        and causality_residual <= residual_limit
        and result_causality_residual == causality_residual
        and _valid_sha256(identity.get("material_table_sha256"))
        and identity.get("result_material_table_sha256") == identity.get("material_table_sha256")
        and _valid_sha256(identity.get("result_sha256"))
        and identity.get("accepted_result_sha256") == identity.get("result_sha256")
    )


def _array_scan_inputs_are_current(raw: Mapping[str, Any]) -> bool:
    identity = raw.get(
        "array_embedded_pattern_feed_phase_active_reflection_scan_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    generation = str(identity.get("array_generation", "")).strip()
    element_order = identity.get("element_order")
    patterns = identity.get("embedded_pattern_sha256")
    try:
        scan_angles = [float(item) for item in identity.get("scan_angles_deg", [])]
        result_scan_angles = [float(item) for item in identity.get("result_scan_angles_deg", [])]
        feed_phase = [[float(item) for item in row] for row in identity.get("feed_phase_deg", [])]
        result_feed_phase = [[float(item) for item in row] for row in identity.get("result_feed_phase_deg", [])]
        reflection = [float(item) for item in identity.get("active_reflection_magnitude", [])]
        result_reflection = [float(item) for item in identity.get("result_active_reflection_magnitude", [])]
        accepted_power = [float(item) for item in identity.get("accepted_power_fraction", [])]
        result_accepted_power = [float(item) for item in identity.get("result_accepted_power_fraction", [])]
    except (TypeError, ValueError):
        return False
    element_count = len(element_order) if isinstance(element_order, list) else 0
    scan_count = len(scan_angles)
    return (
        bool(generation)
        and all(identity.get(key) == generation for key in (
            "pattern_array_generation", "element_array_generation", "phase_array_generation",
            "reflection_array_generation", "scan_array_generation", "power_array_generation",
            "mesh_array_generation", "result_array_generation"))
        and element_count >= 2
        and all(isinstance(item, int) and not isinstance(item, bool) and item > 0 for item in element_order)
        and len(set(element_order)) == element_count
        and identity.get("result_element_order") == element_order
        and isinstance(patterns, list) and len(patterns) == element_count
        and all(_valid_sha256(item) for item in patterns)
        and identity.get("result_embedded_pattern_sha256") == patterns
        and scan_count >= 3
        and all(math.isfinite(item) for item in scan_angles)
        and all(right > left for left, right in zip(scan_angles, scan_angles[1:]))
        and result_scan_angles == scan_angles
        and len(feed_phase) == scan_count
        and all(len(row) == element_count and all(math.isfinite(item) for item in row) for row in feed_phase)
        and result_feed_phase == feed_phase
        and len(reflection) == len(accepted_power) == scan_count
        and all(math.isfinite(item) and 0.0 <= item < 1.0 for item in reflection)
        and all(math.isfinite(item) and 0.0 <= item <= 1.0 for item in accepted_power)
        and all(math.isclose(power, 1.0 - gamma**2, rel_tol=1.0e-12, abs_tol=1.0e-12) for power, gamma in zip(accepted_power, reflection))
        and result_reflection == reflection and result_accepted_power == accepted_power
        and _valid_sha256(identity.get("array_mesh_sha256"))
        and identity.get("result_array_mesh_sha256") == identity.get("array_mesh_sha256")
        and _valid_sha256(identity.get("result_sha256"))
        and identity.get("accepted_result_sha256") == identity.get("result_sha256")
    )


def _waveguide_port_smatrix_inputs_are_current(raw: Mapping[str, Any]) -> bool:
    identity = raw.get(
        "waveguide_port_mode_power_deembed_impedance_frequency_port_smatrix_result_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    generation = str(identity.get("port_generation", "")).strip()
    modes = identity.get("mode_ids")
    ports = identity.get("port_order")
    try:
        power = [float(item) for item in identity.get("power_normalization_w", [])]
        planes = [float(item) for item in identity.get("deembed_plane_m", [])]
        impedance = [float(item) for item in identity.get("reference_impedance_ohm", [])]
        frequencies = [float(item) for item in identity.get("frequency_hz", [])]
        matrix = [
            [[float(part) for part in value] for value in row]
            for row in identity.get("smatrix_ri", [])
        ]
        result_matrix = [
            [[float(part) for part in value] for value in row]
            for row in identity.get("result_smatrix_ri", [])
        ]
    except (TypeError, ValueError):
        return False
    count = len(ports) if isinstance(ports, list) else 0
    return (
        bool(generation)
        and all(identity.get(key) == generation for key in (
            "mode_port_generation", "power_port_generation", "deembed_port_generation",
            "impedance_port_generation", "frequency_port_generation", "order_port_generation",
            "result_port_generation"))
        and count >= 2
        and all(isinstance(item, int) and item > 0 for item in ports)
        and len(set(ports)) == count and identity.get("result_port_order") == ports
        and isinstance(modes, list) and len(modes) == count and all(isinstance(item, str) and item for item in modes)
        and identity.get("result_mode_ids") == modes
        and len(power) == count and all(math.isclose(item, 1.0, rel_tol=0.0, abs_tol=1.0e-12) for item in power)
        and identity.get("result_power_normalization_w") == power
        and len(planes) == count and all(math.isfinite(item) for item in planes)
        and identity.get("result_deembed_plane_m") == planes
        and len(impedance) == count and all(math.isfinite(item) and item > 0.0 for item in impedance)
        and identity.get("result_reference_impedance_ohm") == impedance
        and len(frequencies) >= 2 and all(math.isfinite(item) and item > 0.0 for item in frequencies)
        and all(left < right for left, right in zip(frequencies, frequencies[1:]))
        and identity.get("result_frequency_hz") == frequencies
        and len(matrix) == count and all(len(row) == count for row in matrix)
        and all(len(value) == 2 and all(math.isfinite(part) for part in value) for row in matrix for value in row)
        and result_matrix == matrix
        and _valid_sha256(identity.get("mesh_sha256"))
        and identity.get("result_mesh_sha256") == identity.get("mesh_sha256")
        and _valid_sha256(identity.get("result_sha256"))
        and identity.get("accepted_result_sha256") == identity.get("result_sha256")
    )


def _sar_mass_average_inputs_are_current(raw: Mapping[str, Any]) -> bool:
    identity = raw.get("sar_mass_density_voxel_frequency_field_mesh_result_identity")
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    generation = str(identity.get("sar_generation", "")).strip()
    voxel_ids = identity.get("voxel_ids")
    try:
        mass = float(identity.get("averaging_mass_kg"))
        density = float(identity.get("tissue_density_kg_m3"))
        voxel_mass = [float(item) for item in identity.get("voxel_mass_kg", [])]
        frequency = float(identity.get("frequency_hz"))
        sar = float(identity.get("sar_w_kg"))
    except (TypeError, ValueError):
        return False
    return (
        bool(generation)
        and all(identity.get(key) == generation for key in (
            "mass_sar_generation", "density_sar_generation", "voxel_sar_generation",
            "frequency_sar_generation", "field_sar_generation", "mesh_sar_generation",
            "result_sar_generation"))
        and math.isfinite(mass) and mass > 0.0
        and identity.get("result_averaging_mass_kg") == mass
        and math.isfinite(density) and density > 0.0
        and identity.get("result_tissue_density_kg_m3") == density
        and isinstance(voxel_ids, list) and bool(voxel_ids)
        and all(isinstance(item, int) and item > 0 for item in voxel_ids)
        and len(set(voxel_ids)) == len(voxel_ids)
        and identity.get("result_voxel_ids") == voxel_ids
        and len(voxel_mass) == len(voxel_ids) and all(item > 0.0 for item in voxel_mass)
        and math.isclose(sum(voxel_mass), mass, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and identity.get("result_voxel_mass_kg") == voxel_mass
        and math.isfinite(frequency) and frequency > 0.0
        and identity.get("result_frequency_hz") == frequency
        and identity.get("field_normalization") == "accepted_power_1w"
        and identity.get("result_field_normalization") == "accepted_power_1w"
        and math.isfinite(sar) and sar >= 0.0 and identity.get("result_sar_w_kg") == sar
        and _valid_sha256(identity.get("mesh_sha256"))
        and identity.get("result_mesh_sha256") == identity.get("mesh_sha256")
        and _valid_sha256(identity.get("result_sha256"))
        and identity.get("accepted_result_sha256") == identity.get("result_sha256")
    )


def _wave_port_reference_inputs_are_current(raw: Mapping[str, Any]) -> bool:
    identity = raw.get(
        "wave_port_modal_power_impedance_deembed_phase_balance_result_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    generation = str(identity.get("port_generation", "")).strip()
    modes = identity.get("mode_ids")
    owners = identity.get("port_mode_owner_ids")
    try:
        power = [float(item) for item in identity.get("modal_power_normalization_w", [])]
        impedance = [float(item) for item in identity.get("reference_impedance_ohm", [])]
        planes = [float(item) for item in identity.get("deembed_plane_m", [])]
        phases = [float(item) for item in identity.get("phase_reference_rad", [])]
        incident = float(identity.get("incident_power_w"))
        reflected = float(identity.get("reflected_power_w"))
        transmitted = float(identity.get("transmitted_power_w"))
        dissipated = float(identity.get("dissipated_power_w"))
        result_balance = float(identity.get("result_power_balance_w"))
    except (TypeError, ValueError):
        return False
    count = len(modes) if isinstance(modes, list) else 0
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "mode_port_generation",
                "power_port_generation",
                "impedance_port_generation",
                "deembed_port_generation",
                "phase_port_generation",
                "owner_port_generation",
                "balance_port_generation",
                "result_port_generation",
            )
        )
        and count >= 2
        and all(isinstance(item, str) and item for item in modes)
        and len(set(modes)) == count
        and identity.get("result_mode_ids") == modes
        and len(power) == count
        and all(math.isfinite(item) and item > 0.0 for item in power)
        and identity.get("result_modal_power_normalization_w") == power
        and len(impedance) == count
        and all(math.isfinite(item) and item > 0.0 for item in impedance)
        and identity.get("result_reference_impedance_ohm") == impedance
        and len(planes) == count
        and all(math.isfinite(item) for item in planes)
        and identity.get("result_deembed_plane_m") == planes
        and len(phases) == count
        and all(math.isfinite(item) for item in phases)
        and identity.get("result_phase_reference_rad") == phases
        and isinstance(owners, list)
        and len(owners) == count
        and all(isinstance(item, str) and item for item in owners)
        and len(set(owners)) == count
        and identity.get("result_port_mode_owner_ids") == owners
        and math.isfinite(incident)
        and incident > 0.0
        and all(
            math.isfinite(item) and item >= 0.0
            for item in (reflected, transmitted, dissipated)
        )
        and math.isclose(
            reflected + transmitted + dissipated,
            incident,
            rel_tol=1.0e-9,
            abs_tol=1.0e-12,
        )
        and math.isclose(
            result_balance, incident, rel_tol=1.0e-9, abs_tol=1.0e-12
        )
        and _valid_sha256(identity.get("result_sha256"))
        and identity.get("accepted_result_sha256") == identity.get("result_sha256")
    )


def _farfield_basis_inputs_are_current(raw: Mapping[str, Any]) -> bool:
    identity = raw.get(
        "farfield_spherical_basis_handedness_polarization_phase_power_result_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    generation = str(identity.get("farfield_generation", "")).strip()
    try:
        theta_weights = [float(item) for item in identity.get("theta_weights", [])]
        phi_weights = [float(item) for item in identity.get("phi_weights", [])]
        radiated_power = float(identity.get("radiated_power_w"))
        integrated_power = float(identity.get("integrated_radiated_power_w"))
    except (TypeError, ValueError):
        return False
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "basis_farfield_generation",
                "handedness_farfield_generation",
                "order_farfield_generation",
                "polarization_farfield_generation",
                "phase_farfield_generation",
                "weights_farfield_generation",
                "power_farfield_generation",
                "owner_farfield_generation",
                "result_farfield_generation",
            )
        )
        and identity.get("spherical_basis") == "e_theta_e_phi"
        and identity.get("result_spherical_basis") == identity.get("spherical_basis")
        and identity.get("coordinate_handedness") == "right_handed"
        and identity.get("result_coordinate_handedness")
        == identity.get("coordinate_handedness")
        and identity.get("angular_order") == "theta_major_phi_minor"
        and identity.get("result_angular_order") == identity.get("angular_order")
        and identity.get("polarization_phase_convention") == "exp_plus_j_phase"
        and identity.get("result_polarization_phase_convention")
        == identity.get("polarization_phase_convention")
        and bool(theta_weights)
        and bool(phi_weights)
        and all(math.isfinite(item) and item >= 0.0 for item in theta_weights)
        and all(math.isfinite(item) and item >= 0.0 for item in phi_weights)
        and math.isclose(sum(theta_weights), 1.0, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(sum(phi_weights), 1.0, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and identity.get("result_theta_weights") == theta_weights
        and identity.get("result_phi_weights") == phi_weights
        and math.isfinite(radiated_power)
        and radiated_power >= 0.0
        and math.isclose(
            integrated_power,
            radiated_power,
            rel_tol=1.0e-9,
            abs_tol=1.0e-12,
        )
        and bool(identity.get("farfield_owner_id"))
        and identity.get("accepted_farfield_owner_id")
        == identity.get("farfield_owner_id")
        and _valid_sha256(identity.get("result_sha256"))
        and identity.get("accepted_result_sha256") == identity.get("result_sha256")
    )


def _dispersive_port_inputs_are_current(raw: Mapping[str, Any]) -> bool:
    identity = raw.get(
        "dispersive_port_mode_branch_cutoff_normalization_beta_phase_group_delay_mesh_result_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    generation = str(identity.get("port_generation", "")).strip()
    try:
        cutoff = float(identity.get("cutoff_frequency_hz"))
        frequencies = [float(item) for item in identity.get("frequency_hz", [])]
        beta = [
            float(item)
            for item in identity.get("propagation_constant_rad_per_m", [])
        ]
        phases = [float(item) for item in identity.get("deembedded_phase_rad", [])]
        group_delay = float(identity.get("group_delay_s"))
    except (TypeError, ValueError):
        return False
    if len(frequencies) < 3 or len(beta) != len(frequencies) or len(phases) != len(
        frequencies
    ):
        return False
    expected_delay = -(phases[-1] - phases[0]) / (
        2.0 * math.pi * (frequencies[-1] - frequencies[0])
    )
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "mode_port_generation",
                "branch_port_generation",
                "cutoff_port_generation",
                "normalization_port_generation",
                "beta_port_generation",
                "phase_port_generation",
                "delay_port_generation",
                "mesh_port_generation",
                "result_port_generation",
            )
        )
        and bool(str(identity.get("mode_id", "")).strip())
        and identity.get("result_mode_id") == identity.get("mode_id")
        and bool(str(identity.get("tracked_branch_id", "")).strip())
        and identity.get("result_tracked_branch_id")
        == identity.get("tracked_branch_id")
        and math.isfinite(cutoff)
        and cutoff > 0.0
        and identity.get("result_cutoff_frequency_hz") == cutoff
        and all(math.isfinite(item) and item > cutoff for item in frequencies)
        and all(left < right for left, right in zip(frequencies, frequencies[1:]))
        and identity.get("result_frequency_hz") == frequencies
        and identity.get("modal_normalization") == "unit_forward_power"
        and identity.get("result_modal_normalization") == "unit_forward_power"
        and identity.get("propagation_constant_sign") == "positive_forward"
        and identity.get("result_propagation_constant_sign") == "positive_forward"
        and all(math.isfinite(item) and item > 0.0 for item in beta)
        and all(left < right for left, right in zip(beta, beta[1:]))
        and identity.get("result_propagation_constant_rad_per_m") == beta
        and all(math.isfinite(item) for item in phases)
        and identity.get("result_deembedded_phase_rad") == phases
        and math.isfinite(group_delay)
        and group_delay >= 0.0
        and math.isclose(
            group_delay, expected_delay, rel_tol=1.0e-12, abs_tol=1.0e-18
        )
        and identity.get("result_group_delay_s") == group_delay
        and _valid_sha256(identity.get("mesh_sha256"))
        and identity.get("result_mesh_sha256") == identity.get("mesh_sha256")
        and bool(str(identity.get("result_owner", "")).strip())
        and identity.get("accepted_result_owner") == identity.get("result_owner")
        and _valid_sha256(identity.get("result_sha256"))
        and identity.get("accepted_result_sha256") == identity.get("result_sha256")
    )


def _transient_farfield_inputs_are_current(raw: Mapping[str, Any]) -> bool:
    identity = raw.get(
        "transient_farfield_time_gate_fft_window_phase_center_angular_energy_monitor_result_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    generation = str(identity.get("farfield_generation", "")).strip()
    try:
        gate = [float(item) for item in identity.get("time_gate_s", [])]
        phase_center = [float(item) for item in identity.get("phase_center_m", [])]
        theta = [float(item) for item in identity.get("theta_deg", [])]
        phi = [float(item) for item in identity.get("phi_deg", [])]
        accepted_energy = float(identity.get("accepted_energy_j"))
        radiated_energy = float(identity.get("radiated_energy_j"))
    except (TypeError, ValueError):
        return False
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "gate_farfield_generation",
                "fft_farfield_generation",
                "phase_center_farfield_generation",
                "angular_farfield_generation",
                "energy_farfield_generation",
                "monitor_farfield_generation",
                "owner_farfield_generation",
                "result_farfield_generation",
            )
        )
        and len(gate) == 2
        and all(math.isfinite(item) for item in gate)
        and 0.0 <= gate[0] < gate[1]
        and identity.get("result_time_gate_s") == gate
        and identity.get("fft_window") == "hann"
        and identity.get("result_fft_window") == "hann"
        and identity.get("fft_normalization") == "one_sided_energy_preserving"
        and identity.get("result_fft_normalization")
        == "one_sided_energy_preserving"
        and len(phase_center) == 3
        and all(math.isfinite(item) for item in phase_center)
        and identity.get("result_phase_center_m") == phase_center
        and len(theta) >= 2
        and len(phi) >= 2
        and all(math.isfinite(item) for item in (*theta, *phi))
        and all(left < right for left, right in zip(theta, theta[1:]))
        and all(left < right for left, right in zip(phi, phi[1:]))
        and identity.get("result_theta_deg") == theta
        and identity.get("result_phi_deg") == phi
        and math.isfinite(accepted_energy)
        and accepted_energy > 0.0
        and identity.get("result_accepted_energy_j") == accepted_energy
        and math.isfinite(radiated_energy)
        and 0.0 <= radiated_energy <= accepted_energy
        and identity.get("result_radiated_energy_j") == radiated_energy
        and bool(str(identity.get("monitor_owner", "")).strip())
        and identity.get("result_monitor_owner") == identity.get("monitor_owner")
        and bool(str(identity.get("result_owner", "")).strip())
        and identity.get("accepted_result_owner") == identity.get("result_owner")
        and _valid_sha256(identity.get("result_sha256"))
        and identity.get("accepted_result_sha256") == identity.get("result_sha256")
    )


def _eigenmode_q_closure_inputs_are_current(raw: Mapping[str, Any]) -> bool:
    identity = raw.get(
        "eigenmode_frequency_branch_energy_conductor_dielectric_radiation_q_mesh_owner_result_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    generation = str(identity.get("eigenmode_generation", "")).strip()
    try:
        frequency = float(identity.get("frequency_hz"))
        electric_energy = float(identity.get("electric_energy_j"))
        magnetic_energy = float(identity.get("magnetic_energy_j"))
        stored_energy = float(identity.get("stored_energy_j"))
        q_conductor = float(identity.get("q_conductor"))
        q_dielectric = float(identity.get("q_dielectric"))
        q_radiation = float(identity.get("q_radiation"))
        q_total = float(identity.get("q_total"))
    except (TypeError, ValueError):
        return False
    q_terms = (q_conductor, q_dielectric, q_radiation)
    if not all(math.isfinite(value) and value > 0.0 for value in q_terms):
        return False
    expected_q_total = 1.0 / sum(1.0 / value for value in q_terms)
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "frequency_generation",
                "branch_generation",
                "energy_generation",
                "conductor_q_generation",
                "dielectric_q_generation",
                "radiation_q_generation",
                "inverse_sum_generation",
                "mesh_generation",
                "owner_generation",
                "result_generation",
            )
        )
        and bool(str(identity.get("mode_id", "")).strip())
        and identity.get("result_mode_id") == identity.get("mode_id")
        and bool(str(identity.get("mode_branch", "")).strip())
        and identity.get("result_mode_branch") == identity.get("mode_branch")
        and math.isfinite(frequency)
        and frequency > 0.0
        and identity.get("result_frequency_hz") == frequency
        and all(
            math.isfinite(value) and value >= 0.0
            for value in (electric_energy, magnetic_energy)
        )
        and math.isfinite(stored_energy)
        and stored_energy > 0.0
        and math.isclose(
            stored_energy,
            electric_energy + magnetic_energy,
            rel_tol=1.0e-12,
            abs_tol=1.0e-15,
        )
        and identity.get("result_electric_energy_j") == electric_energy
        and identity.get("result_magnetic_energy_j") == magnetic_energy
        and identity.get("result_stored_energy_j") == stored_energy
        and identity.get("result_q_conductor") == q_conductor
        and identity.get("result_q_dielectric") == q_dielectric
        and identity.get("result_q_radiation") == q_radiation
        and math.isfinite(q_total)
        and q_total > 0.0
        and math.isclose(
            q_total, expected_q_total, rel_tol=1.0e-12, abs_tol=1.0e-12
        )
        and identity.get("result_q_total") == q_total
        and _valid_sha256(identity.get("mesh_sha256"))
        and identity.get("result_mesh_sha256") == identity.get("mesh_sha256")
        and bool(str(identity.get("mode_owner", "")).strip())
        and identity.get("accepted_mode_owner") == identity.get("mode_owner")
        and _valid_sha256(identity.get("result_sha256"))
        and identity.get("accepted_result_sha256") == identity.get("result_sha256")
    )


def _tdr_closure_inputs_are_current(raw: Mapping[str, Any]) -> bool:
    identity = raw.get(
        "tdr_reference_plane_velocity_time_zero_impedance_arrival_window_causality_energy_owner_result_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    generation = str(identity.get("tdr_generation", "")).strip()
    try:
        reference_plane = float(identity.get("reference_plane_m"))
        velocity = float(identity.get("propagation_velocity_m_per_s"))
        time_zero = float(identity.get("time_zero_s"))
        impedance = float(identity.get("characteristic_impedance_ohm"))
        distance = float(identity.get("reflection_distance_m"))
        arrival = float(identity.get("reflection_arrival_s"))
        window = [float(item) for item in identity.get("time_window_s", [])]
        times = [float(item) for item in identity.get("time_samples_s", [])]
        waveform = [
            float(item) for item in identity.get("reflection_waveform", [])
        ]
        pre_arrival_max = float(identity.get("pre_arrival_max_abs"))
        incident_energy = float(identity.get("incident_energy_j"))
        reflected_energy = float(identity.get("reflected_energy_j"))
        accepted_energy = float(identity.get("accepted_energy_j"))
    except (TypeError, ValueError):
        return False
    expected_arrival = (
        time_zero + 2.0 * distance / velocity if velocity > 0.0 else math.nan
    )
    pre_arrival_values = [
        abs(value) for time, value in zip(times, waveform) if time < arrival
    ]
    expected_pre_arrival_max = max(pre_arrival_values, default=0.0)
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "reference_generation",
                "velocity_generation",
                "time_zero_generation",
                "impedance_generation",
                "arrival_generation",
                "window_generation",
                "causality_generation",
                "energy_generation",
                "owner_generation",
                "result_generation",
            )
        )
        and math.isfinite(reference_plane)
        and identity.get("result_reference_plane_m") == reference_plane
        and math.isfinite(velocity)
        and velocity > 0.0
        and identity.get("result_propagation_velocity_m_per_s") == velocity
        and math.isfinite(time_zero)
        and time_zero >= 0.0
        and identity.get("result_time_zero_s") == time_zero
        and math.isfinite(impedance)
        and impedance > 0.0
        and identity.get("result_characteristic_impedance_ohm") == impedance
        and math.isfinite(distance)
        and distance > 0.0
        and identity.get("result_reflection_distance_m") == distance
        and math.isfinite(arrival)
        and math.isclose(arrival, expected_arrival, rel_tol=1.0e-12, abs_tol=1.0e-18)
        and identity.get("result_reflection_arrival_s") == arrival
        and len(window) == 2
        and 0.0 <= window[0] < window[1]
        and identity.get("result_time_window_s") == window
        and len(times) >= 5
        and len(waveform) == len(times)
        and all(math.isfinite(value) for value in (*times, *waveform))
        and all(left < right for left, right in zip(times, times[1:]))
        and window[0] <= times[0] <= arrival <= times[-1] <= window[1]
        and identity.get("result_time_samples_s") == times
        and identity.get("result_reflection_waveform") == waveform
        and math.isfinite(pre_arrival_max)
        and math.isclose(
            pre_arrival_max,
            expected_pre_arrival_max,
            rel_tol=1.0e-12,
            abs_tol=1.0e-15,
        )
        and identity.get("result_pre_arrival_max_abs") == pre_arrival_max
        and all(
            math.isfinite(value) and value >= 0.0
            for value in (incident_energy, reflected_energy, accepted_energy)
        )
        and reflected_energy <= incident_energy
        and math.isclose(
            accepted_energy,
            incident_energy - reflected_energy,
            rel_tol=1.0e-12,
            abs_tol=1.0e-15,
        )
        and identity.get("result_incident_energy_j") == incident_energy
        and identity.get("result_reflected_energy_j") == reflected_energy
        and identity.get("result_accepted_energy_j") == accepted_energy
        and bool(str(identity.get("waveform_owner", "")).strip())
        and identity.get("accepted_waveform_owner") == identity.get("waveform_owner")
        and _valid_sha256(identity.get("result_sha256"))
        and identity.get("accepted_result_sha256") == identity.get("result_sha256")
    )


def _sparameter_gated_inputs_are_current(raw: Mapping[str, Any]) -> bool:
    identity = raw.get(
        "sparameter_reference_plane_time_gate_causality_passivity_energy_port_frequency_owner_result_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    generation = str(identity.get("sparameter_generation", "")).strip()
    try:
        reference_plane = float(identity.get("reference_plane_shift_m"))
        gate = [float(item) for item in identity.get("time_gate_window_s", [])]
        times = [float(item) for item in identity.get("impulse_time_s", [])]
        response = [float(item) for item in identity.get("impulse_response", [])]
        pre_zero_max = float(identity.get("pre_zero_max_abs"))
        singular_values = [
            float(item) for item in identity.get("maximum_singular_values", [])
        ]
        incident = float(identity.get("incident_energy_j"))
        reflected = float(identity.get("reflected_energy_j"))
        transmitted = float(identity.get("transmitted_energy_j"))
        absorbed = float(identity.get("absorbed_energy_j"))
        impedances = [float(item) for item in identity.get("port_impedance_ohm", [])]
        frequencies = [float(item) for item in identity.get("frequency_grid_hz", [])]
    except (TypeError, ValueError):
        return False
    expected_pre_zero_max = max(
        (abs(value) for time, value in zip(times, response) if time < 0.0),
        default=0.0,
    )
    energy_sum = reflected + transmitted + absorbed
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "reference_generation",
                "gate_generation",
                "causality_generation",
                "passivity_generation",
                "energy_generation",
                "port_generation",
                "frequency_generation",
                "owner_generation",
                "result_generation",
            )
        )
        and math.isfinite(reference_plane)
        and identity.get("result_reference_plane_shift_m") == reference_plane
        and len(gate) == 2
        and all(math.isfinite(item) for item in gate)
        and 0.0 <= gate[0] < gate[1]
        and identity.get("result_time_gate_window_s") == gate
        and len(times) >= 3
        and len(response) == len(times)
        and all(math.isfinite(item) for item in (*times, *response))
        and all(left < right for left, right in zip(times, times[1:]))
        and identity.get("result_impulse_time_s") == times
        and identity.get("result_impulse_response") == response
        and math.isfinite(pre_zero_max)
        and math.isclose(
            pre_zero_max, expected_pre_zero_max, rel_tol=1.0e-12, abs_tol=1.0e-15
        )
        and pre_zero_max <= 1.0e-12
        and identity.get("result_pre_zero_max_abs") == pre_zero_max
        and len(singular_values) == len(frequencies) >= 2
        and all(
            math.isfinite(item) and 0.0 <= item <= 1.0 + 1.0e-12
            for item in singular_values
        )
        and identity.get("result_maximum_singular_values") == singular_values
        and all(
            math.isfinite(item) and item >= 0.0
            for item in (incident, reflected, transmitted, absorbed)
        )
        and incident > 0.0
        and math.isclose(
            energy_sum, incident, rel_tol=1.0e-12, abs_tol=1.0e-15
        )
        and identity.get("result_incident_energy_j") == incident
        and identity.get("result_reflected_energy_j") == reflected
        and identity.get("result_transmitted_energy_j") == transmitted
        and identity.get("result_absorbed_energy_j") == absorbed
        and len(impedances) >= 1
        and all(math.isfinite(item) and item > 0.0 for item in impedances)
        and identity.get("result_port_impedance_ohm") == impedances
        and all(math.isfinite(item) and item > 0.0 for item in frequencies)
        and all(left < right for left, right in zip(frequencies, frequencies[1:]))
        and identity.get("result_frequency_grid_hz") == frequencies
        and bool(str(identity.get("sparameter_owner", "")).strip())
        and identity.get("accepted_sparameter_owner")
        == identity.get("sparameter_owner")
        and _valid_sha256(identity.get("sparameter_sha256"))
        and identity.get("accepted_sparameter_sha256")
        == identity.get("sparameter_sha256")
    )


def _degenerate_eigenmode_subspace_inputs_are_current(
    raw: Mapping[str, Any],
) -> bool:
    identity = raw.get(
        "eigenmode_degenerate_subspace_principal_angle_mass_orthogonality_phase_tracking_residual_mesh_owner_result_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    generation = str(identity.get("degenerate_mode_generation", "")).strip()
    try:
        frequencies = [float(item) for item in identity.get("mode_frequencies_hz", [])]
        angles = [float(item) for item in identity.get("principal_angles_rad", [])]
        gram_real = [
            [float(item) for item in row]
            for row in identity.get("mass_gram_real", [])
        ]
        gram_imag = [
            [float(item) for item in row]
            for row in identity.get("mass_gram_imag", [])
        ]
        phase_anchors = [
            [float(item) for item in row]
            for row in identity.get("phase_anchor_complex", [])
        ]
        residuals = [float(item) for item in identity.get("residual_norms", [])]
        mesh_counts = [int(item) for item in identity.get("mesh_cell_counts", [])]
        mesh_frequencies = [
            float(item)
            for item in identity.get("mesh_converged_frequency_hz", [])
        ]
    except (TypeError, ValueError):
        return False
    tracking_ids = identity.get("tracking_subspace_ids")
    near_degenerate = (
        len(frequencies) == 2
        and all(math.isfinite(item) and item > 0.0 for item in frequencies)
        and abs(frequencies[1] - frequencies[0]) / max(frequencies) <= 1.0e-5
    )
    gram_ok = (
        len(gram_real) == len(gram_imag) == 2
        and all(len(row) == 2 for row in (*gram_real, *gram_imag))
        and all(
            math.isclose(
                gram_real[i][j], 1.0 if i == j else 0.0, abs_tol=1.0e-9
            )
            and math.isclose(gram_imag[i][j], 0.0, abs_tol=1.0e-9)
            for i in range(2)
            for j in range(2)
        )
    )
    phase_ok = (
        len(phase_anchors) == 2
        and all(len(anchor) == 2 for anchor in phase_anchors)
        and all(
            math.isfinite(real)
            and math.isfinite(imag)
            and real > 0.0
            and math.isclose(imag, 0.0, abs_tol=1.0e-9)
            for real, imag in phase_anchors
        )
    )
    mesh_errors = [
        abs(mesh_frequencies[index] - mesh_frequencies[-1])
        / mesh_frequencies[-1]
        for index in range(len(mesh_frequencies) - 1)
    ] if mesh_frequencies and mesh_frequencies[-1] > 0.0 else []
    mesh_ok = (
        len(mesh_counts) == len(mesh_frequencies) >= 3
        and all(count > 0 for count in mesh_counts)
        and all(left < right for left, right in zip(mesh_counts, mesh_counts[1:]))
        and all(math.isfinite(item) and item > 0.0 for item in mesh_frequencies)
        and len(mesh_errors) >= 2
        and mesh_errors[-1] <= mesh_errors[-2]
        and mesh_errors[-1] <= 5.0e-3
    )
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "subspace_generation",
                "principal_angle_generation",
                "mass_generation",
                "phase_generation",
                "tracking_generation",
                "residual_generation",
                "mesh_generation",
                "owner_generation",
                "result_generation",
            )
        )
        and near_degenerate
        and identity.get("result_mode_frequencies_hz") == frequencies
        and len(angles) == 2
        and all(math.isfinite(item) and 0.0 <= item <= 0.01 for item in angles)
        and identity.get("result_principal_angles_rad") == angles
        and gram_ok
        and identity.get("result_mass_gram_real") == gram_real
        and identity.get("result_mass_gram_imag") == gram_imag
        and phase_ok
        and identity.get("result_phase_anchor_complex") == phase_anchors
        and isinstance(tracking_ids, Sequence)
        and not isinstance(tracking_ids, (str, bytes))
        and len(tracking_ids) == 2
        and bool(str(tracking_ids[0]).strip())
        and tracking_ids[0] == tracking_ids[1]
        and identity.get("result_tracking_subspace_ids") == list(tracking_ids)
        and len(residuals) == 2
        and all(math.isfinite(item) and 0.0 <= item <= 1.0e-6 for item in residuals)
        and identity.get("result_residual_norms") == residuals
        and mesh_ok
        and identity.get("result_mesh_cell_counts") == mesh_counts
        and identity.get("result_mesh_converged_frequency_hz") == mesh_frequencies
        and _valid_sha256(identity.get("eigenmode_mesh_sha256"))
        and identity.get("result_eigenmode_mesh_sha256")
        == identity.get("eigenmode_mesh_sha256")
        and bool(str(identity.get("field_owner", "")).strip())
        and identity.get("accepted_field_owner") == identity.get("field_owner")
        and _valid_sha256(identity.get("field_sha256"))
        and identity.get("accepted_field_sha256") == identity.get("field_sha256")
    )


def _waveguide_port_modal_closure_inputs_are_current(raw: Mapping[str, Any]) -> bool:
    identity = raw.get(
        "waveguide_port_mode_power_orthogonality_impedance_deembed_cutoff_frequency_owner_result_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    generation = str(identity.get("waveguide_port_generation", "")).strip()
    try:
        powers = [float(item) for item in identity.get("modal_power_w", [])]
        gram_real = [[float(item) for item in row] for row in identity.get("mode_gram_real", [])]
        gram_imag = [[float(item) for item in row] for row in identity.get("mode_gram_imag", [])]
        impedance = [float(item) for item in identity.get("modal_impedance_ohm", [])]
        frequencies = [float(item) for item in identity.get("frequency_grid_hz", [])]
        cutoff = float(identity.get("cutoff_frequency_hz"))
        beta = [float(item) for item in identity.get("propagation_constant_rad_m", [])]
        reference_plane = float(identity.get("reference_plane_m"))
        deembedded_plane = float(identity.get("deembedded_reference_plane_m"))
        phase = [float(item) for item in identity.get("deembed_phase_rad", [])]
    except (TypeError, ValueError):
        return False
    count = len(frequencies)
    c0 = 299_792_458.0
    expected_impedance = [
        377.0 / math.sqrt(1.0 - (cutoff / frequency) ** 2)
        for frequency in frequencies
    ] if cutoff > 0.0 and all(frequency > cutoff for frequency in frequencies) else []
    expected_beta = [
        2.0 * math.pi / c0 * math.sqrt(frequency**2 - cutoff**2)
        for frequency in frequencies
    ] if expected_impedance else []
    distance = deembedded_plane - reference_plane
    mirrored_fields = (
        "mode_name", "normalization", "modal_power_w", "mode_gram_real",
        "mode_gram_imag", "impedance_definition", "modal_impedance_ohm",
        "frequency_grid_hz", "cutoff_frequency_hz", "propagation_constant_rad_m",
        "reference_plane_m", "deembedded_reference_plane_m", "deembed_phase_rad",
        "port_mesh_sha256",
    )
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "mode_generation", "power_generation", "orthogonality_generation",
                "impedance_generation", "deembed_generation", "cutoff_generation",
                "frequency_generation", "mesh_generation", "owner_generation",
                "result_generation",
            )
        )
        and identity.get("mode_name") == "TE10"
        and identity.get("normalization") == "accepted_power_1w"
        and count >= 3
        and len(powers) == len(impedance) == len(beta) == len(phase) == count
        and all(math.isfinite(item) and item > 0.0 for item in frequencies)
        and all(left < right for left, right in zip(frequencies, frequencies[1:]))
        and math.isfinite(cutoff)
        and cutoff > 0.0
        and all(frequency > cutoff for frequency in frequencies)
        and all(math.isclose(item, 1.0, rel_tol=0.0, abs_tol=1.0e-12) for item in powers)
        and len(gram_real) == len(gram_imag) == 2
        and all(len(row) == 2 for row in (*gram_real, *gram_imag))
        and all(
            math.isclose(gram_real[i][j], 1.0 if i == j else 0.0, abs_tol=1.0e-12)
            and math.isclose(gram_imag[i][j], 0.0, abs_tol=1.0e-12)
            for i in range(2) for j in range(2)
        )
        and identity.get("impedance_definition") == "te_wave_impedance"
        and all(
            math.isclose(observed, expected, rel_tol=1.0e-12, abs_tol=1.0e-9)
            for observed, expected in zip(impedance, expected_impedance)
        )
        and all(
            math.isclose(observed, expected, rel_tol=1.0e-12, abs_tol=1.0e-9)
            for observed, expected in zip(beta, expected_beta)
        )
        and math.isfinite(reference_plane)
        and math.isfinite(deembedded_plane)
        and distance > 0.0
        and all(
            math.isclose(observed, -expected * distance, rel_tol=1.0e-12, abs_tol=1.0e-12)
            for observed, expected in zip(phase, beta)
        )
        and all(identity.get(f"result_{field}") == identity.get(field) for field in mirrored_fields)
        and _valid_sha256(identity.get("port_mesh_sha256"))
        and bool(str(identity.get("port_owner", "")).strip())
        and identity.get("accepted_port_owner") == identity.get("port_owner")
        and _valid_sha256(identity.get("port_result_sha256"))
        and identity.get("accepted_port_result_sha256") == identity.get("port_result_sha256")
    )


def _nearfar_power_inputs_are_current(raw: Mapping[str, Any]) -> bool:
    identity = raw.get(
        "nearfar_sphere_power_directivity_gain_efficiency_polarization_quadrature_mesh_owner_result_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    generation = str(identity.get("nearfar_generation", "")).strip()
    try:
        frequency = float(identity.get("frequency_hz"))
        accepted = float(identity.get("accepted_power_w"))
        sphere = float(identity.get("enclosing_sphere_power_w"))
        radiated = float(identity.get("radiated_power_w"))
        efficiency = float(identity.get("radiation_efficiency"))
        directivity = float(identity.get("maximum_directivity_linear"))
        gain = float(identity.get("realized_gain_linear"))
        weights = [float(item) for item in identity.get("angular_quadrature_weights_sr", [])]
        intensity = [float(item) for item in identity.get("radiation_intensity_w_sr", [])]
    except (TypeError, ValueError):
        return False
    integrated_power = sum(weight * value for weight, value in zip(weights, intensity))
    expected_directivity = 4.0 * math.pi * max(intensity, default=math.nan) / radiated if radiated > 0.0 else math.nan
    mirrored_fields = (
        "frequency_hz", "accepted_power_w", "enclosing_sphere_power_w",
        "radiated_power_w", "radiation_efficiency", "maximum_directivity_linear",
        "realized_gain_linear", "polarization_basis", "copolar_definition",
        "angular_quadrature_weights_sr", "radiation_intensity_w_sr",
        "farfield_mesh_sha256",
    )
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "sphere_generation", "power_generation", "directivity_generation",
                "gain_generation", "efficiency_generation", "polarization_generation",
                "quadrature_generation", "mesh_generation", "owner_generation",
                "result_generation",
            )
        )
        and all(math.isfinite(item) and item > 0.0 for item in (frequency, accepted, sphere, radiated, directivity, gain))
        and math.isfinite(efficiency)
        and 0.0 < efficiency <= 1.0
        and math.isclose(sphere, radiated, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(efficiency, radiated / accepted, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(gain, directivity * efficiency, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and len(weights) == len(intensity) >= 4
        and all(math.isfinite(item) and item > 0.0 for item in weights)
        and all(math.isfinite(item) and item >= 0.0 for item in intensity)
        and math.isclose(sum(weights), 4.0 * math.pi, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(integrated_power, radiated, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(directivity, expected_directivity, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and identity.get("polarization_basis") == "theta_phi_right_handed"
        and identity.get("copolar_definition") == "ludwig3"
        and all(identity.get(f"result_{field}") == identity.get(field) for field in mirrored_fields)
        and _valid_sha256(identity.get("farfield_mesh_sha256"))
        and bool(str(identity.get("nearfar_owner", "")).strip())
        and identity.get("accepted_nearfar_owner") == identity.get("nearfar_owner")
        and _valid_sha256(identity.get("nearfar_result_sha256"))
        and identity.get("accepted_nearfar_result_sha256") == identity.get("nearfar_result_sha256")
    )


def _energy_history_restart_offsets_close(
    summary: Mapping[str, Any], run_count: int
) -> bool:
    segments = summary.get("energy_history_segments")
    if segments is None:
        return True
    if (
        not isinstance(segments, Sequence)
        or isinstance(segments, (str, bytes))
        or not segments
    ):
        return False
    previous_end = -1
    previous_offset_out = None
    generations = set()
    for segment in segments:
        if not isinstance(segment, Mapping):
            return False
        generation = str(segment.get("segment_generation", "")).strip()
        try:
            start = int(segment["start_run_index"])
            end = int(segment["end_run_index"])
            offset_in = float(segment["coenergy_offset_in_J"])
            offset_out = float(segment["coenergy_offset_out_J"])
        except (KeyError, TypeError, ValueError):
            return False
        if (
            not generation
            or generation in generations
            or start != previous_end + 1
            or end < start
            or end >= run_count
            or not math.isfinite(offset_in)
            or not math.isfinite(offset_out)
        ):
            return False
        if previous_offset_out is not None and not math.isclose(
            offset_in, previous_offset_out, rel_tol=1.0e-12, abs_tol=1.0e-15
        ):
            return False
        generations.add(generation)
        previous_end = end
        previous_offset_out = offset_out
    return previous_end == run_count - 1


def nonlinear_inductance_sweep_gate(
    summary: Mapping[str, Any],
    *,
    max_identity_relative_error: float = 1.0e-5,
    max_matrix_symmetry_relative_error: float = 1.0e-6,
    matrix_psd_relative_tolerance: float = 1.0e-10,
    max_replay_relative_error: float = 1.0e-9,
    maximum_residual_log10: float = -5.0,
    regime_margin: float = 0.05,
    minimum_saturation_drop: float = 0.25,
) -> dict[str, Any]:
    """Gate apparent/tangent matrices, nonlinear energy duality, and replay.

    The low-current differential inductance may exceed the apparent value while
    permeability is rising.  The gate therefore requires an observed crossover,
    rather than imposing the incorrect global rule ``L_incremental <= L_apparent``.
    """

    if not isinstance(summary, Mapping):
        raise ValueError("summary must be a mapping")
    runs = summary.get("runs")
    if (
        not isinstance(runs, Sequence)
        or isinstance(runs, (str, bytes))
        or len(runs) < 6
    ):
        raise ValueError("runs must contain at least three current levels with replay")
    tolerances = (
        max_identity_relative_error,
        max_matrix_symmetry_relative_error,
        matrix_psd_relative_tolerance,
        max_replay_relative_error,
        regime_margin,
        minimum_saturation_drop,
    )
    if any(not math.isfinite(value) or value < 0.0 for value in tolerances):
        raise ValueError("all tolerances and margins must be finite and nonnegative")
    if not math.isfinite(maximum_residual_log10):
        raise ValueError("maximum_residual_log10 must be finite")

    parsed = []
    groups: dict[float, list[dict[str, Any]]] = defaultdict(list)
    for index, raw in enumerate(runs):
        if not isinstance(raw, Mapping):
            raise ValueError(f"run {index} must be a mapping")
        requested = float(raw.get("current_A_requested", math.nan))
        replay = int(raw.get("replay", 0))
        apparent = _matrix(raw.get("apparent_inductance_H"), "apparent_inductance_H")
        incremental = _matrix(
            raw.get("incremental_inductance_H"), "incremental_inductance_H"
        )
        currents = _vector(raw.get("current_A"), "current_A")
        flux = _vector(raw.get("flux_linkage_Vs"), "flux_linkage_Vs")
        energy = float(raw.get("energy_J", math.nan))
        coenergy = float(raw.get("coenergy_J", math.nan))
        residual = float(raw.get("final_nonlinear_residual_log10", math.nan))
        if (
            not math.isfinite(requested)
            or requested <= 0.0
            or replay <= 0
            or not all(math.isfinite(value) for value in (energy, coenergy, residual))
        ):
            raise ValueError(f"run {index} has invalid current, replay, energy, or residual")

        app_metrics = _matrix_metrics(apparent)
        inc_metrics = _matrix_metrics(incremental)
        predicted_flux = [
            apparent[0][0] * currents[0] + apparent[0][1] * currents[1],
            apparent[1][0] * currents[0] + apparent[1][1] * currents[1],
        ]
        flux_error = max(
            _relative_error(actual, expected)
            for actual, expected in zip(flux, predicted_flux)
        )
        duality_target = currents[0] * flux[0] + currents[1] * flux[1]
        duality_error = _relative_error(energy + coenergy, duality_target)
        current_error = _relative_error(currents[0], requested)
        matrix_identity_ok, matrix_current_ok = _matrix_operating_point_identity_matches(
            raw, currents
        )

        def matrix_ok(metrics: Mapping[str, float]) -> bool:
            scale = max(abs(metrics["diagonal_product_H2"]), 1.0e-300)
            return (
                metrics["l11_H"] > 0.0
                and metrics["l22_H"] > 0.0
                and metrics["symmetry_relative_error"]
                <= max_matrix_symmetry_relative_error
                and metrics["determinant_H2"]
                >= -matrix_psd_relative_tolerance * scale
            )

        checks = {
            "apparent_matrix_is_symmetric_psd": matrix_ok(app_metrics),
            "incremental_matrix_is_symmetric_psd": matrix_ok(inc_metrics),
            "requested_primary_current_is_reproduced": current_error
            <= max_identity_relative_error,
            "secondary_is_open_circuit": abs(currents[1])
            <= max_identity_relative_error * max(abs(currents[0]), 1.0),
            "apparent_matrix_closes_flux_linkage": flux_error
            <= max_identity_relative_error,
            "energy_coenergy_legendre_duality_closes": duality_error
            <= max_identity_relative_error,
            "energy_and_coenergy_are_nonnegative": energy >= 0.0 and coenergy >= 0.0,
            "nonlinear_iteration_converged": residual <= maximum_residual_log10,
            "result_metadata_run_ids_are_consistent": (
                _result_metadata_run_ids_are_consistent(raw)
            ),
            "apparent_and_incremental_matrix_operating_point_ids_match": (
                matrix_identity_ok
            ),
            "matrix_operating_point_currents_match_run_current": matrix_current_ok,
            "reported_and_artifact_units_are_consistent_si": (
                _artifact_units_are_consistent(raw)
            ),
            "inductance_matrices_share_solve_sweep_generation": (
                _matrix_sweep_generations_match(raw)
            ),
            "matrix_rows_columns_and_vectors_share_port_order": (
                _matrix_port_orders_match(raw)
            ),
            "stored_energy_and_loss_series_share_si_basis": (
                _energy_loss_basis_is_si(raw)
            ),
            "sparameters_share_complex_reference_impedance_or_renormalization": (
                _sparameter_reference_impedance_is_bound(raw)
            ),
            "frequency_axis_unit_and_hz_scale_share_identity": (
                _frequency_axis_unit_is_bound(raw)
            ),
            "sparameter_port_modes_share_deembedded_reference_planes": (
                _sparameter_reference_planes_are_bound(raw)
            ),
            "energy_and_loss_share_q_frequency_sample": (
                _energy_q_frequency_sample_is_bound(raw)
            ),
            "mixed_mode_basis_matches_current_single_ended_port_order": (
                _mixed_mode_sparameter_basis_matches_port_order(raw)
            ),
            "realized_gain_and_accepted_power_share_frequency_sample": (
                _farfield_gain_power_frequency_sample_is_bound(raw)
            ),
            "field_monitor_interpolation_matches_current_mesh_pass": (
                _field_monitor_interpolation_matches_mesh_pass(raw)
            ),
            "port_deembed_reference_plane_uses_explicit_length_unit": (
                _port_deembed_reference_plane_unit_is_bound(raw)
            ),
            "sparameter_renormalization_matches_complex_reference_impedance": (
                _sparameter_renormalization_matches_reference_impedance(raw)
            ),
            "farfield_co_cross_uses_current_ludwig_polarization_basis": (
                _farfield_ludwig_polarization_basis_is_current(raw)
            ),
            "sparameters_use_one_power_wave_normalization_generation": (
                _sparameter_power_wave_normalization_is_current(raw)
            ),
            "fft_window_uses_current_coherent_gain_correction": (
                _fft_window_coherent_gain_is_current(raw)
            ),
            "sparameter_renormalization_uses_current_complex_impedance_calibration": (
                _sparameter_complex_impedance_renormalization_is_current(raw)
            ),
            "farfield_comparison_uses_explicit_current_polarization_transform": (
                _farfield_polarization_basis_transform_is_current(raw)
            ),
            "mixed_mode_sparameters_use_current_port_pair_order_and_polarity": (
                _mixed_mode_sparameter_pair_order_is_current(raw)
            ),
            "nearfield_farfield_phase_center_uses_one_global_coordinate_frame": (
                _nearfield_farfield_phase_center_frame_is_current(raw)
            ),
            "sparameter_deembed_uses_current_per_port_reference_planes": (
                _sparameter_deembed_reference_plane_map_is_current(raw)
            ),
            "time_domain_port_transform_uses_current_gate_window": (
                _time_domain_port_signal_gate_window_is_current(raw)
            ),
            "sparameter_renormalization_uses_current_port_reference_impedances": (
                _sparameter_renormalization_reference_impedance_is_current(raw)
            ),
            "realized_gain_uses_current_excitation_and_accepted_power": (
                _realized_gain_excitation_and_accepted_power_are_current(raw)
            ),
            "farfield_polarization_basis_and_phase_center_use_current_coordinates": (
                _farfield_polarization_phase_center_is_current(raw)
            ),
            "broadband_energy_q_uses_current_port_and_loss_normalization": (
                _broadband_energy_q_inputs_are_current(raw)
            ),
            "mixed_mode_uses_current_pairs_impedances_polarities_and_planes": (
                _mixed_mode_port_metadata_is_current(raw)
            ),
            "time_farfield_fft_uses_current_grid_window_scaling_and_phase_center": (
                _time_farfield_fft_inputs_are_current(raw)
            ),
            "degenerate_modes_use_current_mesh_phase_order_and_overlap": (
                _waveguide_degenerate_mode_tracking_is_current(raw)
            ),
            "dispersive_fields_use_current_causal_poles_temperature_and_units": (
                _dispersive_causal_pole_fit_is_current(raw)
            ),
            "broadband_sparameters_use_current_mesh_interpolation_modes_and_impedance": (
                _broadband_sparameter_inputs_are_current(raw)
            ),
            "transient_monitors_use_current_time_waveform_frame_and_mesh": (
                _transient_monitor_inputs_are_current(raw)
            ),
            "deembedded_network_uses_current_planes_phase_causality_passivity_and_grid": (
                _deembedded_network_inputs_are_current(raw)
            ),
            "field_circuit_cosim_uses_current_sign_impedance_and_power_balance": (
                _field_circuit_cosim_inputs_are_current(raw)
            ),
            "adaptive_results_use_current_mesh_pass_sparameter_energy_grid_and_stop_rule": (
                _adaptive_mesh_convergence_inputs_are_current(raw)
            ),
            "eigenmodes_use_current_subspace_phase_normalization_ports_and_mesh": (
                _eigenmode_tracking_inputs_are_current(raw)
            ),
            "sparameters_use_current_port_modes_planes_impedances_normalization_grid_and_result": (
                _port_network_inputs_are_current(raw)
            ),
            "farfields_use_current_angular_grid_polarization_coordinates_power_mesh_and_result": (
                _farfield_result_inputs_are_current(raw)
            ),
            "time_domain_sparameters_use_current_waveform_normalization_fft_grid_deembedding_and_result": (
                _time_domain_port_smatrix_inputs_are_current(raw)
            ),
            "near_to_far_results_use_current_huygens_orientation_phase_center_frequency_mesh_and_transform": (
                _huygens_near_far_inputs_are_current(raw)
            ),
            "waveguide_ports_use_current_mode_cutoff_normalization_plane_mesh_field_and_result": (
                _waveguide_port_mode_inputs_are_current(raw)
            ),
            "wake_impedance_uses_current_bunch_time_grid_transform_frequency_normalization_mesh_and_result": (
                _wake_impedance_inputs_are_current(raw)
            ),
            "dispersive_vector_fit_uses_current_stable_poles_residues_passivity_causality_temperature_and_result": (
                _dispersive_vector_fit_inputs_are_current(raw)
            ),
            "array_scan_uses_current_embedded_patterns_element_order_phases_reflection_power_mesh_and_result": (
                _array_scan_inputs_are_current(raw)
            ),
            "waveguide_ports_use_current_modes_power_deembed_impedance_frequency_order_smatrix_mesh_and_result": (
                _waveguide_port_smatrix_inputs_are_current(raw)
            ),
            "sar_uses_current_average_mass_density_voxels_frequency_field_mesh_and_result": (
                _sar_mass_average_inputs_are_current(raw)
            ),
            "wave_ports_use_current_modal_power_impedance_deembed_phase_owner_balance_and_result": (
                _wave_port_reference_inputs_are_current(raw)
            ),
            "farfields_use_current_spherical_basis_handedness_order_polarization_weights_power_owner_and_result": (
                _farfield_basis_inputs_are_current(raw)
            ),
            "dispersive_ports_use_current_mode_branch_cutoff_power_normalization_beta_phase_group_delay_mesh_and_result": (
                _dispersive_port_inputs_are_current(raw)
            ),
            "transient_farfields_use_current_time_gate_fft_phase_center_angles_energy_monitor_owner_and_result": (
                _transient_farfield_inputs_are_current(raw)
            ),
            "eigenmodes_use_current_frequency_branch_energy_q_inverse_sum_mesh_owner_and_result": (
                _eigenmode_q_closure_inputs_are_current(raw)
            ),
            "tdr_uses_current_reference_velocity_time_zero_impedance_arrival_causality_energy_owner_and_result": (
                _tdr_closure_inputs_are_current(raw)
            ),
            "sparameters_use_current_reference_gate_causality_passivity_energy_ports_frequency_owner_and_result": (
                _sparameter_gated_inputs_are_current(raw)
            ),
            "degenerate_eigenmodes_use_current_subspace_angles_mass_orthogonality_phase_tracking_residual_mesh_owner_and_result": (
                _degenerate_eigenmode_subspace_inputs_are_current(raw)
            ),
            "waveguide_port_modes_use_current_power_orthogonality_impedance_deembed_cutoff_frequency_mesh_owner_and_result": (
                _waveguide_port_modal_closure_inputs_are_current(raw)
            ),
            "nearfar_results_use_current_sphere_power_directivity_gain_efficiency_polarization_quadrature_mesh_owner_and_result": (
                _nearfar_power_inputs_are_current(raw)
            ),
        }
        row = {
            "current_A_requested": requested,
            "replay": replay,
            "apparent_inductance_H": apparent,
            "incremental_inductance_H": incremental,
            "current_A": currents,
            "flux_linkage_Vs": flux,
            "energy_J": energy,
            "coenergy_J": coenergy,
            "final_nonlinear_residual_log10": residual,
            "operating_point_id": raw.get("operating_point_id"),
            "reported_units": raw.get("reported_units"),
            "artifact_units": raw.get("artifact_units"),
            "differential_to_apparent_primary_ratio": incremental[0][0]
            / apparent[0][0],
            "flux_identity_relative_error": flux_error,
            "energy_coenergy_duality_relative_error": duality_error,
            "apparent_matrix_metrics": app_metrics,
            "incremental_matrix_metrics": inc_metrics,
            "checks": checks,
            "status": "ok" if all(checks.values()) else "needs_attention",
        }
        parsed.append(row)
        groups[requested].append(row)

    levels = sorted(groups)
    replay_errors = {}
    replay_checks = {}
    representatives = []
    for current in levels:
        group = sorted(groups[current], key=lambda row: row["replay"])
        replay_checks[current] = len(group) == 2 and {
            row["replay"] for row in group
        } == {1, 2}
        reference = _flatten_replay_values(group[0])
        errors = [
            _relative_error(actual, expected)
            for row in group[1:]
            for actual, expected in zip(_flatten_replay_values(row), reference)
        ]
        replay_errors[current] = max(errors, default=math.inf)
        representatives.append(group[0])

    ratios = [row["differential_to_apparent_primary_ratio"] for row in representatives]
    apparent_primary = [row["apparent_inductance_H"][0][0] for row in representatives]
    incremental_primary = [
        row["incremental_inductance_H"][0][0] for row in representatives
    ]
    peak_apparent = max(apparent_primary)
    peak_incremental = max(incremental_primary)
    family_checks = {
        "at_least_three_distinct_positive_current_levels": len(levels) >= 3,
        "every_level_has_independent_replay": all(replay_checks.values()),
        "all_run_identities_and_matrices_close": all(
            row["status"] == "ok" for row in parsed
        ),
        "replays_are_stable": max(replay_errors.values(), default=math.inf)
        <= max_replay_relative_error,
        "initial_magnetization_rise_is_observed": ratios[0] >= 1.0 + regime_margin,
        "saturated_differential_response_is_observed": ratios[-1]
        <= 1.0 - regime_margin,
        "differential_to_apparent_crossover_is_observed": any(
            left > 1.0 and right < 1.0
            for left, right in zip(ratios, ratios[1:])
        ),
        "high_current_apparent_inductance_drops_from_peak": apparent_primary[-1]
        <= (1.0 - minimum_saturation_drop) * peak_apparent,
        "high_current_incremental_inductance_drops_from_peak": incremental_primary[-1]
        <= (1.0 - minimum_saturation_drop) * peak_incremental,
        "restart_energy_history_offsets_are_continuous": (
            _energy_history_restart_offsets_close(summary, len(runs))
        ),
    }
    return {
        "policy": "nonlinear_inductance_sweep_gate_v1",
        "status": "ok" if all(family_checks.values()) else "needs_attention",
        "checks": family_checks,
        "issues": [name for name, ok in family_checks.items() if not ok],
        "current_levels_A": levels,
        "differential_to_apparent_primary_ratios": ratios,
        "maximum_replay_relative_error": max(replay_errors.values(), default=None),
        "maximum_flux_identity_relative_error": max(
            row["flux_identity_relative_error"] for row in parsed
        ),
        "maximum_energy_coenergy_duality_relative_error": max(
            row["energy_coenergy_duality_relative_error"] for row in parsed
        ),
        "maximum_matrix_symmetry_relative_error": max(
            metric["symmetry_relative_error"]
            for row in parsed
            for metric in (
                row["apparent_matrix_metrics"],
                row["incremental_matrix_metrics"],
            )
        ),
        "runs": parsed,
        "lesson": (
            "Use the apparent matrix for flux linkage at the operating point and "
            "the incremental matrix for small-signal response. Close W + W' = I dot "
            "psi, reciprocity, positive semidefiniteness, and replay. Differential "
            "inductance may exceed apparent inductance while permeability rises; the "
            "physically useful signature is a measured crossover followed by saturation."
        ),
    }
