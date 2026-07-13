"""Multi-conductor Maxwell-capacitance cross-formulation validation gate."""
from __future__ import annotations

import math
from collections.abc import Mapping, Sequence



def _numpy():
    import numpy as np

    return np


def _matrix_family(value: object, name: str) -> list[np.ndarray]:
    np = _numpy()
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)) or len(value) < 2:
        raise ValueError(f"{name} must contain at least two matrices")
    matrices: list[np.ndarray] = []
    dimension: int | None = None
    for index, item in enumerate(value):
        matrix = np.asarray(item, dtype=float)
        if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
            raise ValueError(f"{name}[{index}] must be square")
        if matrix.shape[0] < 3:
            raise ValueError(f"{name}[{index}] must have at least three conductors")
        if dimension is None:
            dimension = int(matrix.shape[0])
        elif matrix.shape != (dimension, dimension):
            raise ValueError(f"all {name} matrices must have the same shape")
        if not np.all(np.isfinite(matrix)):
            raise ValueError(f"{name}[{index}] must contain finite values")
        matrices.append(matrix)
    return matrices


def multiconductor_capacitance_cross_formulation_gate(
    summary: Mapping[str, object],
    *,
    max_reciprocity_relative_error: float = 2.0e-5,
    max_cross_formulation_relative_error: float = 8.0e-2,
    max_cross_formulation_error_span: float = 2.0e-2,
    minimum_position_relative_change: float = 1.0e-4,
) -> dict[str, object]:
    """Gate reciprocal passive Maxwell matrices from two formulations."""

    np = _numpy()
    if not isinstance(summary, Mapping):
        raise TypeError("summary must be a mapping")
    positions = summary.get("positions")
    if not isinstance(positions, Sequence) or isinstance(positions, (str, bytes)):
        raise ValueError("positions must be a sequence")
    positions_float = [float(value) for value in positions]
    if len(positions_float) < 2 or not all(math.isfinite(value) for value in positions_float):
        raise ValueError("positions must contain at least two finite values")
    if any(right <= left for left, right in zip(positions_float, positions_float[1:])):
        raise ValueError("positions must be strictly increasing")

    formulations = summary.get("formulations")
    if not isinstance(formulations, Mapping) or set(formulations) != {
        "volume_fem",
        "boundary_integral",
    }:
        raise ValueError("formulations must contain volume_fem and boundary_integral")
    families = {
        name: _matrix_family(value, f"formulations.{name}")
        for name, value in formulations.items()
    }
    matrix_count = len(families["volume_fem"])
    if matrix_count != len(positions_float) or any(
        len(family) != matrix_count for family in families.values()
    ):
        raise ValueError("positions and both matrix families must have equal length")
    if families["volume_fem"][0].shape != families["boundary_integral"][0].shape:
        raise ValueError("both formulations must have the same conductor count")

    family_metrics: dict[str, list[dict[str, float]]] = {}
    all_matrices: list[np.ndarray] = []
    for name, family in families.items():
        rows = []
        for matrix in family:
            scale = max(float(np.linalg.norm(matrix)), 1.0e-300)
            symmetric = 0.5 * (matrix + matrix.T)
            offdiag = matrix[~np.eye(matrix.shape[0], dtype=bool)]
            rows.append(
                {
                    "reciprocity_relative_frobenius": float(
                        np.linalg.norm(matrix - matrix.T) / scale
                    ),
                    "minimum_symmetric_eigenvalue": float(
                        np.linalg.eigvalsh(symmetric).min()
                    ),
                    "minimum_diagonal": float(np.diag(matrix).min()),
                    "maximum_offdiagonal_relative": float(offdiag.max() / scale),
                    "minimum_row_sum": float(matrix.sum(axis=1).min()),
                }
            )
            all_matrices.append(matrix)
        family_metrics[name] = rows

    cross_errors = [
        float(np.linalg.norm(boundary - volume) / np.linalg.norm(volume))
        for volume, boundary in zip(
            families["volume_fem"], families["boundary_integral"]
        )
    ]
    sensitivities = {
        name: (np.diag(family[-1]) - np.diag(family[0])) / np.diag(family[0])
        for name, family in families.items()
    }
    maximum_reciprocity_error = max(
        row["reciprocity_relative_frobenius"]
        for rows in family_metrics.values()
        for row in rows
    )
    minimum_eigenvalue = min(
        row["minimum_symmetric_eigenvalue"]
        for rows in family_metrics.values()
        for row in rows
    )
    minimum_diagonal = min(
        row["minimum_diagonal"] for rows in family_metrics.values() for row in rows
    )
    maximum_offdiagonal_relative = max(
        row["maximum_offdiagonal_relative"]
        for rows in family_metrics.values()
        for row in rows
    )
    minimum_row_sum = min(
        row["minimum_row_sum"] for rows in family_metrics.values() for row in rows
    )
    volume_sign = np.sign(sensitivities["volume_fem"])
    boundary_sign = np.sign(sensitivities["boundary_integral"])
    checks = {
        "capacitance_unit_explicit": str(summary.get("capacitance_unit") or "").strip()
        in {"F", "pF", "nF", "uF"},
        "all_matrices_reciprocal": maximum_reciprocity_error
        <= float(max_reciprocity_relative_error),
        "all_matrices_positive_definite": minimum_eigenvalue > 0.0,
        "all_diagonals_positive": minimum_diagonal > 0.0,
        "all_offdiagonals_nonpositive": maximum_offdiagonal_relative <= 0.0,
        "all_ground_row_sums_positive": minimum_row_sum > 0.0,
        "cross_formulation_endpoint_agreement": max(cross_errors)
        <= float(max_cross_formulation_relative_error),
        "cross_formulation_error_is_stable": max(cross_errors) - min(cross_errors)
        <= float(max_cross_formulation_error_span),
        "position_sensitivity_is_resolved": min(
            float(np.min(np.abs(value))) for value in sensitivities.values()
        )
        >= float(minimum_position_relative_change),
        "position_sensitivity_direction_agrees": np.array_equal(
            volume_sign, boundary_sign
        ),
    }
    return {
        "schema": "radia-ngsolve.multiconductor-capacitance-cross-formulation.v1",
        "policy": "multiconductor_capacitance_cross_formulation_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "conductor_count": int(all_matrices[0].shape[0]),
            "position_count": matrix_count,
            "maximum_reciprocity_relative_error": maximum_reciprocity_error,
            "minimum_symmetric_eigenvalue": minimum_eigenvalue,
            "minimum_diagonal": minimum_diagonal,
            "maximum_offdiagonal_relative": maximum_offdiagonal_relative,
            "minimum_row_sum": minimum_row_sum,
            "cross_formulation_relative_errors": cross_errors,
            "cross_formulation_error_span": max(cross_errors) - min(cross_errors),
            "position_sensitivity": {
                name: value.tolist() for name, value in sensitivities.items()
            },
        },
        "lesson": (
            "Validate an N-conductor sensor matrix by reciprocity, passivity, "
            "Maxwell sign structure, grounded row sums, formulation agreement, "
            "and position-sensitivity direction rather than one scalar capacitance."
        ),
    }
