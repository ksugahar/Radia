"""Canonical-unit gate for single-valued nonlinear magnetic material tables."""
from __future__ import annotations

import hashlib
import json
import math
from typing import Any


MU0 = 4.0e-7 * math.pi
SCHEMA = "radia-nonlinear-bh-table/v1"


def _rows(value: Any) -> tuple[list[list[float]], list[str]]:
    if not isinstance(value, list):
        return [], ["rows must be a list of [H, B] pairs"]
    parsed: list[list[float]] = []
    issues: list[str] = []
    for index, row in enumerate(value):
        if not isinstance(row, (list, tuple)) or len(row) != 2:
            issues.append(f"rows[{index}] must contain exactly H and B")
            continue
        try:
            h_value, b_value = (float(item) for item in row)
        except (TypeError, ValueError):
            issues.append(f"rows[{index}] must be numeric")
            continue
        if not math.isfinite(h_value) or not math.isfinite(b_value):
            issues.append(f"rows[{index}] must be finite")
            continue
        parsed.append([h_value, b_value])
    return parsed, issues


def _canonical_digest(rows: list[list[float]]) -> str:
    payload = json.dumps(rows, separators=(",", ":"), ensure_ascii=True).encode("ascii")
    return hashlib.sha256(payload).hexdigest()


def nonlinear_bh_canonical_table_gate(
    contract: dict[str, Any],
    *,
    minimum_point_count: int = 2,
    origin_atol: float = 1.0e-15,
) -> dict[str, Any]:
    """Accept only explicit SI ``[H, B]`` tables used by the nonlinear solver.

    Unit conversion belongs to the importing adapter. This gate prevents a
    plausible-looking ``[B, H]`` or CGS table from entering a solver that
    expects ``H`` in A/m and ``B`` in tesla. The material model is deliberately
    limited to a single-valued, isotropic soft-magnetic law; hysteresis and
    permanent-magnet recoil need stateful or vector constitutive contracts.
    """

    if not isinstance(contract, dict):
        raise ValueError("contract must be a mapping")
    parsed_rows, parse_issues = _rows(contract.get("rows"))
    h_values = [row[0] for row in parsed_rows]
    b_values = [row[1] for row in parsed_rows]
    digest = _canonical_digest(parsed_rows)
    declared_digest = str(contract.get("canonical_table_sha256", "")).casefold()

    checks = {
        "schema_is_canonical": contract.get("schema") == SCHEMA,
        "material_model_is_single_valued_isotropic_soft_magnetic": (
            contract.get("material_model")
            == "single_valued_isotropic_soft_magnetic"
        ),
        "column_order_is_h_then_b": contract.get("column_order") == ["H", "B"],
        "units_are_explicit_si": contract.get("units") == {"H": "A/m", "B": "T"},
        "interpolation_matches_solver": contract.get("interpolation")
        == "monotone_pchip",
        "extrapolation_matches_solver": contract.get("extrapolation")
        == "vacuum_slope",
        "rows_parse_and_are_finite": not parse_issues,
        "point_count_sufficient": len(parsed_rows) >= int(minimum_point_count),
        "origin_is_zero": bool(parsed_rows)
        and abs(h_values[0]) <= float(origin_atol)
        and abs(b_values[0]) <= float(origin_atol),
        "h_is_strictly_increasing": len(parsed_rows) >= 2
        and all(right > left for left, right in zip(h_values, h_values[1:])),
        "b_is_monotone_nondecreasing": len(parsed_rows) >= 2
        and all(right >= left for left, right in zip(b_values, b_values[1:])),
        "soft_magnetic_secant_mu_is_at_least_vacuum": len(parsed_rows) >= 2
        and all(b + 1.0e-15 >= MU0 * h for h, b in parsed_rows),
        "canonical_digest_is_present_and_matches": len(declared_digest) == 64
        and declared_digest == digest,
    }
    issues = parse_issues + [name for name, passed in checks.items() if not passed]
    return {
        "schema": "radia-nonlinear-bh-canonical-table-gate/v1",
        "policy": "explicit_si_h_b_single_valued_material_gate_v1",
        "status": "accepted" if not issues else "rejected",
        "checks": checks,
        "issues": issues,
        "canonical_table_sha256": digest,
        "metrics": {
            "point_count": len(parsed_rows),
            "h_max_a_per_m": max(h_values) if h_values else math.nan,
            "b_max_t": max(b_values) if b_values else math.nan,
        },
        "solver_contract": {
            "column_order": ["H", "B"],
            "units": {"H": "A/m", "B": "T"},
            "interpolation": "monotone_pchip",
            "extrapolation": "vacuum_slope",
        },
        "notes": [
            "Convert source tables before this gate; unit guessing from numeric magnitude is forbidden.",
            "A single-valued table is not a hysteresis model and must not be replayed as one.",
        ],
    }
