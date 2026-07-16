"""Solver-neutral gate for nonlinear two-winding inductance sweeps."""

from __future__ import annotations

import math
from collections import defaultdict
from collections.abc import Mapping, Sequence
from typing import Any


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
