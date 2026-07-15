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
