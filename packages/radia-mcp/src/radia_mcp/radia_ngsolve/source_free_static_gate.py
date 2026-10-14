"""Null-solution gate for source-free static Maxwell models."""

from __future__ import annotations

import math
from collections.abc import Mapping


def source_free_static_null_solution_gate(
    summary: Mapping[str, object],
    *,
    absolute_tolerance: float = 1.0e-14,
) -> dict[str, object]:
    """Check a source-free static solve without treating it as driven validation."""
    if not isinstance(summary, Mapping):
        raise TypeError("summary must be a mapping")

    tolerance = float(absolute_tolerance)
    if not math.isfinite(tolerance) or tolerance < 0:
        raise ValueError("absolute_tolerance must be finite and nonnegative")

    tables = list(summary.get("tables") or [])
    if not tables or not all(isinstance(table, Mapping) for table in tables):
        raise ValueError("tables must contain mappings")

    values: list[float] = []
    shapes: list[tuple[int, int, int, int]] = []
    for index, table in enumerate(tables):
        matrix = table.get("values") or []
        if not isinstance(matrix, list):
            raise ValueError(f"tables[{index}].values must be a matrix")

        row_lengths: list[int] = []
        for row in matrix:
            if not isinstance(row, list):
                raise ValueError(f"tables[{index}].values must be a matrix")
            row_lengths.append(len(row))
            values.extend(float(value) for value in row)

        matrix_cols = row_lengths[0] if row_lengths else 0
        rectangular = bool(row_lengths) and all(length == matrix_cols for length in row_lengths)
        shapes.append(
            (
                int(table.get("rows", -1)),
                int(table.get("cols", -1)),
                len(matrix),
                matrix_cols if rectangular else -1,
            )
        )

    required_quantities = {"joule_loss", "hysteresis_loss", "iron_loss"}
    quantity_kinds = {str(table.get("quantity_kind") or "") for table in tables}
    checks = {
        "static_source_free_contract": (
            str(summary.get("analysis_kind")) == "magnetostatic"
            and int(summary.get("source_count", -1)) == 0
            and int(summary.get("condition_count", -1)) == 0
        ),
        "mesh_and_result_present": (
            summary.get("has_mesh") is True and summary.get("has_result") is True
        ),
        "required_loss_quantities_present": required_quantities.issubset(quantity_kinds),
        "table_shapes_consistent": all(
            rows > 0 and cols > 0 and rows == matrix_rows and cols == matrix_cols
            for rows, cols, matrix_rows, matrix_cols in shapes
        ),
        "all_observables_finite": bool(values) and all(math.isfinite(value) for value in values),
        "null_solution_observables_zero": (
            bool(values) and max(abs(value) for value in values) <= tolerance
        ),
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "source_free_static_null_solution_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "maximum_absolute_observable": max((abs(value) for value in values), default=None),
        "absolute_tolerance": tolerance,
        "lesson": (
            "With homogeneous constraints and no impressed source, the static Maxwell "
            "solution is the zero field. A successful mesh/solve plus zero loss "
            "observables is a useful end-to-end smoke, but does not validate a driven model."
        ),
    }
