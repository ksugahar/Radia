"""Solver-neutral consistency checks for piecewise-linear nonlinear B-H data."""

from __future__ import annotations

import math
from typing import Any


MU0 = 4.0e-7 * math.pi


def _finite(value: Any, name: str) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be numeric") from exc
    if not math.isfinite(parsed):
        raise ValueError(f"{name} must be finite")
    return parsed


def _rows(
    values: Any,
    *,
    name: str,
    x_key: str,
    y_key: str,
) -> tuple[list[float], list[float], list[str]]:
    if not isinstance(values, list):
        raise ValueError(f"{name} must be a list")
    x_values: list[float] = []
    y_values: list[float] = []
    errors: list[str] = []
    for index, row in enumerate(values):
        if not isinstance(row, dict):
            errors.append(f"{name}[{index}] is not a mapping")
            continue
        try:
            x_values.append(_finite(row[x_key], f"{name}[{index}].{x_key}"))
            y_values.append(_finite(row[y_key], f"{name}[{index}].{y_key}"))
        except (KeyError, TypeError, ValueError) as exc:
            errors.append(str(exc))
    return x_values, y_values, errors


def _max_relative_error(observed: list[float], expected: list[float]) -> float:
    if len(observed) != len(expected) or not observed:
        return math.inf
    return max(
        abs(value - reference) / max(abs(value), abs(reference), 1.0e-300)
        for value, reference in zip(observed, expected, strict=True)
    )


def nonlinear_bh_piecewise_material_gate(
    summary: dict[str, Any],
    *,
    maximum_relative_identity_error: float = 1.0e-8,
    origin_atol: float = 1.0e-15,
    high_field_mu_diff_target: float = 1.0,
    high_field_mu_diff_atol: float = 2.0e-2,
    minimum_bh_point_count: int = 4,
) -> dict[str, Any]:
    """Gate secant and differential permeability derived from one B-H table.

    For a piecewise-linear B-H table, the relative secant permeability at a
    positive knot is ``B/(mu0*H)``. The differential permeability associated
    with that knot is the slope of the interval ending at the knot,
    ``delta B/(mu0*delta H)``. A central, spline, or PCHIP derivative is a
    different constitutive law and must not be silently substituted.
    """

    if not isinstance(summary, dict):
        raise ValueError("summary must be a mapping")
    contract = summary.get("contract")
    if not isinstance(contract, dict):
        raise ValueError("contract must be a mapping")
    h_b, b_t, bh_errors = _rows(
        summary.get("bh_rows"),
        name="bh_rows",
        x_key="h_a_per_m",
        y_key="b_t",
    )
    h_mu, mu_r, mu_errors = _rows(
        summary.get("secant_rows"),
        name="secant_rows",
        x_key="h_a_per_m",
        y_key="relative_mu",
    )
    h_diff, mu_diff_r, diff_errors = _rows(
        summary.get("differential_rows"),
        name="differential_rows",
        x_key="h_a_per_m",
        y_key="relative_mu_diff",
    )
    parse_errors = bh_errors + mu_errors + diff_errors

    positive_h = h_b[1:] if len(h_b) >= 2 else []
    secant_expected = [
        b / (MU0 * h) for h, b in zip(positive_h, b_t[1:], strict=True)
    ]
    differential_expected = [
        (b_right - b_left) / (MU0 * (h_right - h_left))
        for h_left, h_right, b_left, b_right in zip(
            h_b[:-1], h_b[1:], b_t[:-1], b_t[1:], strict=True
        )
        if h_right > h_left
    ]
    secant_error = _max_relative_error(mu_r, secant_expected)
    differential_error = _max_relative_error(mu_diff_r, differential_expected)
    saturation_tail_expected = contract.get("saturation_tail_expected") is True
    high_field_check = (
        bool(mu_diff_r)
        and abs(mu_diff_r[-1] - float(high_field_mu_diff_target))
        <= float(high_field_mu_diff_atol)
        and bool(mu_r)
        and mu_r[-1] > mu_diff_r[-1]
        if saturation_tail_expected
        else True
    )

    checks = {
        "rows_parsed_and_finite": not parse_errors,
        "bh_point_count_sufficient": len(h_b) >= int(minimum_bh_point_count),
        "origin_is_zero": bool(h_b)
        and abs(h_b[0]) <= float(origin_atol)
        and abs(b_t[0]) <= float(origin_atol),
        "h_grid_strictly_increasing": len(h_b) >= 2
        and all(right > left for left, right in zip(h_b, h_b[1:])),
        "b_is_monotone_nondecreasing": len(b_t) >= 2
        and all(right >= left for left, right in zip(b_t, b_t[1:])),
        "derived_grids_match_positive_bh_knots": positive_h == h_mu == h_diff
        and len(mu_r) == len(mu_diff_r) == len(h_b) - 1,
        "piecewise_linear_contract_recorded": contract.get("interpolation")
        == "piecewise_linear"
        and contract.get("differential_interval") == "left_interval_ending_at_knot",
        "permeability_definitions_recorded": contract.get("secant_definition")
        == "B/(mu0*H)"
        and contract.get("differential_definition")
        == "deltaB/(mu0*deltaH)",
        "secant_permeability_identity": secant_error
        <= float(maximum_relative_identity_error),
        "differential_permeability_identity": differential_error
        <= float(maximum_relative_identity_error),
        "differential_permeability_positive": bool(mu_diff_r)
        and all(value > 0.0 for value in mu_diff_r),
        "declared_saturation_tail_is_physical": high_field_check,
    }
    issues = parse_errors + [name for name, ok in checks.items() if not ok]
    return {
        "schema": "radia-nonlinear-bh-piecewise-material/v1",
        "policy": "piecewise_bh_secant_and_left_interval_differential_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "metrics": {
            "bh_point_count": len(h_b),
            "derived_point_count": len(mu_r),
            "secant_mu_maximum_relative_error": secant_error,
            "differential_mu_maximum_relative_error": differential_error,
            "high_field_secant_mu_r": mu_r[-1] if mu_r else math.nan,
            "high_field_differential_mu_r": mu_diff_r[-1]
            if mu_diff_r
            else math.nan,
        },
        "tolerances": {
            "maximum_relative_identity_error": float(maximum_relative_identity_error),
            "origin_atol": float(origin_atol),
            "high_field_mu_diff_target": float(high_field_mu_diff_target),
            "high_field_mu_diff_atol": float(high_field_mu_diff_atol),
            "minimum_bh_point_count": int(minimum_bh_point_count),
        },
        "notes": [
            "Use secant permeability for B/H reporting and differential permeability for a material tangent; they are not interchangeable in saturation.",
            "For a piecewise-linear table, the differential value at a positive knot belongs to the interval ending at that knot.",
            "Central or spline differentiation changes the constitutive law and can destabilize or mis-scale nonlinear iterations.",
            "A high-field vacuum-slope check is applied only when the input contract declares a saturation tail.",
        ],
    }
