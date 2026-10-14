"""Validation gates for complex vector-field component exports."""

from __future__ import annotations

import math
from typing import Any


def _relative_error(actual: float, expected: float) -> float:
    return abs(actual - expected) / max(abs(actual), abs(expected), 1.0e-300)


def complex_vector_field_maximum_gate(
    summary: dict[str, Any],
    *,
    max_component_magnitude_relative_error: float = 5.0e-5,
    max_reported_maximum_relative_error: float = 1.0e-12,
) -> dict[str, Any]:
    """Check ``|B|`` rows and per-material maxima in real/imaginary exports."""
    if not isinstance(summary, dict):
        raise ValueError("summary must be a mapping")
    cases = summary.get("cases") or []
    if not isinstance(cases, list) or not cases:
        raise ValueError("cases must contain at least one field export")

    checks: dict[str, bool] = {}
    case_metrics = []
    all_case_ids = []
    for case_index, case in enumerate(cases, start=1):
        if not isinstance(case, dict):
            raise ValueError(f"case {case_index} must be a mapping")
        case_id = str(case.get("case_id") or "").strip()
        rows = case.get("rows") or []
        maxima = case.get("maxima") or []
        if not isinstance(rows, list) or not isinstance(maxima, list):
            raise ValueError(f"case {case_index} rows/maxima must be lists")
        parsed_rows = []
        row_errors = []
        for row_index, row in enumerate(rows):
            try:
                part = str(row["part"])
                element = int(row["element"])
                material_id = int(row["material_id"])
                bx = float(row["bx_t"])
                by = float(row["by_t"])
                bz = float(row["bz_t"])
                reported = float(row["bmag_t"])
                calculated = math.sqrt(bx * bx + by * by + bz * bz)
                if not all(math.isfinite(value) for value in (bx, by, bz, reported)):
                    raise ValueError("non-finite field component")
                parsed_rows.append(
                    {
                        "part": part,
                        "element": element,
                        "material_id": material_id,
                        "reported": reported,
                        "magnitude_error": _relative_error(reported, calculated),
                    }
                )
            except (KeyError, TypeError, ValueError):
                row_errors.append(f"row {row_index} is incomplete")

        maximum_errors = []
        maximum_element_matches = []
        maximum_parse_errors = []
        for maximum_index, maximum in enumerate(maxima):
            try:
                part = str(maximum["part"])
                material_id = int(maximum["material_id"])
                reported_element = int(maximum["element"])
                reported_max = float(maximum["bmax_t"])
                candidates = [
                    row
                    for row in parsed_rows
                    if row["part"] == part and row["material_id"] == material_id
                ]
                if not candidates:
                    raise ValueError("maximum has no matching field rows")
                calculated_row = max(candidates, key=lambda row: row["reported"])
                maximum_errors.append(
                    _relative_error(reported_max, calculated_row["reported"])
                )
                reported_element_rows = [
                    row for row in candidates if row["element"] == reported_element
                ]
                maximum_element_matches.append(
                    bool(reported_element_rows)
                    and min(
                        _relative_error(reported_max, row["reported"])
                        for row in reported_element_rows
                    )
                    <= float(max_reported_maximum_relative_error)
                )
            except (KeyError, TypeError, ValueError):
                maximum_parse_errors.append(f"maximum {maximum_index} is incomplete")

        parts = {row["part"] for row in parsed_rows}
        prefix = f"case_{case_index}"
        checks[f"{prefix}_id_recorded"] = bool(case_id)
        checks[f"{prefix}_rows_parsed"] = bool(parsed_rows) and not row_errors
        checks[f"{prefix}_real_and_imaginary_parts_present"] = parts == {
            "real",
            "imaginary",
        }
        checks[f"{prefix}_tesla_frequency_excitation_recorded"] = (
            case.get("field_unit") == "T"
            and math.isfinite(float(case.get("frequency_hz", math.nan)))
            and float(case.get("frequency_hz")) > 0.0
            and math.isfinite(float(case.get("ampere_turns", math.nan)))
            and abs(float(case.get("ampere_turns"))) > 0.0
        )
        checks[f"{prefix}_component_magnitudes_close"] = bool(parsed_rows) and max(
            row["magnitude_error"] for row in parsed_rows
        ) <= float(max_component_magnitude_relative_error)
        checks[f"{prefix}_maxima_parsed"] = bool(maximum_errors) and not maximum_parse_errors
        checks[f"{prefix}_reported_maxima_close"] = bool(maximum_errors) and max(
            maximum_errors
        ) <= float(max_reported_maximum_relative_error)
        checks[f"{prefix}_maximum_elements_match"] = bool(maximum_element_matches) and all(
            maximum_element_matches
        )
        all_case_ids.append(case_id)
        case_metrics.append(
            {
                "case_id": case_id,
                "field_row_count": len(parsed_rows),
                "maximum_count": len(maximum_errors),
                "maximum_component_magnitude_relative_error": max(
                    (row["magnitude_error"] for row in parsed_rows), default=math.inf
                ),
                "maximum_reported_maximum_relative_error": max(
                    maximum_errors, default=math.inf
                ),
                "parse_errors": row_errors + maximum_parse_errors,
            }
        )
    checks["case_ids_distinct"] = len(set(all_case_ids)) == len(all_case_ids) and all(
        all_case_ids
    )
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "complex_vector_field_maximum_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "metrics": {
            "case_count": len(cases),
            "field_row_count": sum(row["field_row_count"] for row in case_metrics),
            "maximum_count": sum(row["maximum_count"] for row in case_metrics),
            "maximum_component_magnitude_relative_error": max(
                row["maximum_component_magnitude_relative_error"] for row in case_metrics
            ),
            "maximum_reported_maximum_relative_error": max(
                row["maximum_reported_maximum_relative_error"] for row in case_metrics
            ),
            "cases": case_metrics,
        },
        "tolerances": {
            "max_component_magnitude_relative_error": float(
                max_component_magnitude_relative_error
            ),
            "max_reported_maximum_relative_error": float(
                max_reported_maximum_relative_error
            ),
        },
        "notes": [
            "Validate real and imaginary component blocks independently; do not combine their maxima as if they were simultaneous scalar samples.",
            "Reported per-material maxima must identify an element that actually attains the row maximum.",
        ],
    }
