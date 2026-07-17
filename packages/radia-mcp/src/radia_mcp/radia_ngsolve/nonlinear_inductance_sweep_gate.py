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
