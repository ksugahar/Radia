"""Mirror-symmetric three-group magnetization handoff gate."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any


def _finite(value: object, name: str) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be numeric") from exc
    if not math.isfinite(parsed):
        raise ValueError(f"{name} must be finite")
    return parsed


def _rows(value: object, name: str, minimum: int) -> list[Mapping[str, object]]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise ValueError(f"{name} must be an array")
    rows = list(value)
    if len(rows) < minimum or any(not isinstance(row, Mapping) for row in rows):
        raise ValueError(f"{name} must contain at least {minimum} objects")
    return rows


def _relative_l2(left: list[float], right: list[float]) -> float:
    if len(left) != len(right) or not left:
        return math.inf
    residual = math.sqrt(sum((a - b) ** 2 for a, b in zip(left, right)))
    scale = max(
        math.sqrt(sum(value * value for value in left)),
        math.sqrt(sum(value * value for value in right)),
        1.0e-300,
    )
    return residual / scale


def mirror_symmetric_three_magnet_handoff_gate(
    summary: Mapping[str, object],
) -> dict[str, Any]:
    """Recompute vector normalization, mirror symmetry, and repeatability."""

    if not isinstance(summary, Mapping):
        raise ValueError("summary must be an object")
    contract = summary.get("result_contract")
    if not isinstance(contract, Mapping):
        raise ValueError("result_contract must be an object")
    group_ids_value = contract.get("group_ids")
    if not isinstance(group_ids_value, Sequence) or isinstance(group_ids_value, (str, bytes)):
        raise ValueError("result_contract.group_ids must be an array")
    group_ids = [int(_finite(value, "group_id")) for value in group_ids_value]
    if group_ids != [1, 2, 3]:
        raise ValueError("group_ids must be [1, 2, 3]")
    vectors_per_group = int(
        _finite(contract.get("vectors_per_group"), "vectors_per_group")
    )
    if vectors_per_group <= 0:
        raise ValueError("vectors_per_group must be positive")
    mirror_axis = str(contract.get("mirror_axis") or "")
    if mirror_axis not in {"x", "y", "z"}:
        raise ValueError("mirror_axis must be x, y, or z")
    mirror_index = {"x": 0, "y": 1, "z": 2}[mirror_axis]
    runs = _rows(summary.get("runs"), "runs", 2)
    tolerances = summary.get("gate_tolerances", {})
    if not isinstance(tolerances, Mapping):
        raise ValueError("gate_tolerances must be an object")
    limits = {
        "direction_norm": _finite(
            tolerances.get("maximum_direction_norm_error", 1.0e-12),
            "maximum_direction_norm_error",
        ),
        "mirror": _finite(
            tolerances.get("maximum_outer_mirror_residual", 2.0e-12),
            "maximum_outer_mirror_residual",
        ),
        "outer_response": _finite(
            tolerances.get("maximum_outer_response_relative_mismatch", 1.0e-6),
            "maximum_outer_response_relative_mismatch",
        ),
        "repeat": _finite(
            tolerances.get("maximum_repeat_relative_l2", 1.0e-12),
            "maximum_repeat_relative_l2",
        ),
    }
    if any(value < 0.0 for value in limits.values()):
        raise ValueError("gate tolerances must be nonnegative")

    run_metrics = []
    repeat_vectors: list[list[float]] = []
    repeat_responses: list[list[float]] = []
    for run_index, run in enumerate(runs):
        vectors = _rows(
            run.get("magnetization_vectors"),
            f"runs[{run_index}].magnetization_vectors",
            1,
        )
        groups: dict[int, list[tuple[int, list[float], float]]] = {
            group: [] for group in group_ids
        }
        for vector_index, row in enumerate(vectors):
            group = int(_finite(row.get("group"), f"vectors[{vector_index}].group"))
            if group not in groups:
                raise ValueError(f"vectors[{vector_index}].group is not declared")
            direction_value = row.get("direction")
            if not isinstance(direction_value, Sequence) or isinstance(
                direction_value, (str, bytes)
            ) or len(direction_value) != 3:
                raise ValueError(f"vectors[{vector_index}].direction must contain three values")
            direction = [
                _finite(value, f"vectors[{vector_index}].direction[{axis}]")
                for axis, value in enumerate(direction_value)
            ]
            magnitude = _finite(
                row.get("magnitude"), f"vectors[{vector_index}].magnitude"
            )
            element_id = int(
                _finite(row.get("element_id"), f"vectors[{vector_index}].element_id")
            )
            groups[group].append((element_id, direction, magnitude))

        counts = {group: len(rows) for group, rows in groups.items()}
        maximum_norm_error = max(
            abs(math.sqrt(sum(value * value for value in direction)) - 1.0)
            for rows in groups.values()
            for _, direction, _ in rows
        )
        left_signatures = []
        for _, direction, magnitude in groups[group_ids[0]]:
            mirrored = list(direction)
            mirrored[mirror_index] *= -1.0
            left_signatures.append((*mirrored, magnitude))
        right_signatures = [
            (*direction, magnitude)
            for _, direction, magnitude in groups[group_ids[-1]]
        ]
        left_signatures.sort()
        right_signatures.sort()
        outer_mirror_residual = max(
            abs(left - right)
            for left_row, right_row in zip(left_signatures, right_signatures)
            for left, right in zip(left_row, right_row)
        ) if len(left_signatures) == len(right_signatures) else math.inf

        response_value = run.get("bmax_t_by_material")
        if not isinstance(response_value, Mapping):
            raise ValueError(f"runs[{run_index}].bmax_t_by_material must be an object")
        response = [
            _finite(
                response_value.get(str(group), response_value.get(group)),
                f"runs[{run_index}].response[{group}]",
            )
            for group in group_ids
        ]
        if any(value <= 0.0 for value in response):
            raise ValueError("all group responses must be positive")
        outer_response_mismatch = abs(response[0] - response[2]) / max(
            abs(response[0]), abs(response[2]), 1.0e-300
        )
        center_to_outer_ratio = response[1] / (0.5 * (response[0] + response[2]))
        flattened = [
            value
            for group in group_ids
            for _, direction, magnitude in sorted(groups[group], key=lambda row: row[0])
            for value in (*direction, magnitude)
        ]
        repeat_vectors.append(flattened)
        repeat_responses.append(response)
        run_metrics.append(
            {
                "group_vector_counts": counts,
                "maximum_direction_norm_error": maximum_norm_error,
                "outer_mirror_residual": outer_mirror_residual,
                "outer_response_relative_mismatch": outer_response_mismatch,
                "center_to_outer_response_ratio": center_to_outer_ratio,
                "solution_complete": run.get("end_of_moment_solution") is True
                and int(run.get("solver_returncode", -1)) == 0
                and int(run.get("error_count", -1)) == 0,
            }
        )

    maximum_vector_repeat_error = max(
        _relative_l2(repeat_vectors[0], values) for values in repeat_vectors[1:]
    )
    maximum_response_repeat_error = max(
        _relative_l2(repeat_responses[0], values) for values in repeat_responses[1:]
    )
    checks = {
        "two_or_more_complete_solver_runs": len(runs) >= 2
        and all(row["solution_complete"] for row in run_metrics),
        "three_groups_have_expected_vector_count": all(
            row["group_vector_counts"]
            == {group: vectors_per_group for group in group_ids}
            for row in run_metrics
        ),
        "all_magnetization_directions_are_unit_vectors": max(
            row["maximum_direction_norm_error"] for row in run_metrics
        )
        <= limits["direction_norm"],
        "outer_group_vectors_are_mirror_symmetric": max(
            row["outer_mirror_residual"] for row in run_metrics
        )
        <= limits["mirror"],
        "outer_group_responses_are_symmetric": max(
            row["outer_response_relative_mismatch"] for row in run_metrics
        )
        <= limits["outer_response"],
        "center_group_response_is_positive_and_lower": all(
            0.0 < row["center_to_outer_response_ratio"] < 1.0
            for row in run_metrics
        ),
        "fresh_runs_repeat_vector_package": maximum_vector_repeat_error
        <= limits["repeat"],
        "fresh_runs_repeat_group_response": maximum_response_repeat_error
        <= limits["repeat"],
    }
    return {
        "policy": "mirror_symmetric_three_magnet_handoff_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "handoff_ready": all(checks.values()),
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "limits": limits,
        "metrics": {
            "run_count": len(runs),
            "per_run": run_metrics,
            "maximum_vector_repeat_relative_l2": maximum_vector_repeat_error,
            "maximum_response_repeat_relative_l2": maximum_response_repeat_error,
        },
        "lesson": (
            "A multi-magnet handoff must preserve normalized element directions, declared group cardinality, "
            "the intended mirror transform, symmetric outer response, and repeatability across fresh solves."
        ),
    }
