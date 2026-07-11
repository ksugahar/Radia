"""Validation gate for families of linear two-winding inductance matrices."""
from __future__ import annotations

import math


def inductance_matrix_family_gate(
    cases,
    *,
    expected_strongest_coupling_case: str | None = None,
    max_reciprocity_relative_error: float = 0.02,
    psd_relative_tolerance: float = 1.0e-12,
):
    if not isinstance(cases, list) or len(cases) < 2:
        raise ValueError("cases must contain at least two matrix rows")
    if max_reciprocity_relative_error < 0.0 or psd_relative_tolerance < 0.0:
        raise ValueError("tolerances must be nonnegative")

    rows = []
    ids = []
    for index, case in enumerate(cases):
        if not isinstance(case, dict):
            raise ValueError(f"case {index} must be an object")
        case_id = str(case.get("case_id") or "").strip()
        matrix = case.get("matrix_H")
        if not case_id or not isinstance(matrix, list) or len(matrix) != 2:
            raise ValueError(f"case {index} needs case_id and a 2x2 matrix_H")
        if any(not isinstance(row, list) or len(row) != 2 for row in matrix):
            raise ValueError(f"case {case_id} matrix_H must be 2x2")
        l11, m12 = (float(value) for value in matrix[0])
        m21, l22 = (float(value) for value in matrix[1])
        finite = all(math.isfinite(value) for value in (l11, m12, m21, l22))
        mutual = 0.5 * (m12 + m21)
        reciprocity = abs(m12 - m21) / max(abs(m12), abs(m21), 1.0e-300)
        diagonal_product = l11 * l22
        determinant = diagonal_product - mutual * mutual
        coupling = (
            abs(mutual) / math.sqrt(diagonal_product)
            if finite and l11 > 0.0 and l22 > 0.0
            else math.inf
        )
        checks = {
            "all_finite": finite,
            "positive_self_inductances": l11 > 0.0 and l22 > 0.0,
            "mutual_terms_have_consistent_sign": m12 == 0.0 or m21 == 0.0 or m12 * m21 > 0.0,
            "reciprocity_within_tolerance": reciprocity <= max_reciprocity_relative_error,
            "symmetrized_matrix_positive_semidefinite": determinant
            >= -psd_relative_tolerance * max(abs(diagonal_product), 1.0e-300),
            "coupling_coefficient_bounded": coupling <= 1.0 + psd_relative_tolerance,
        }
        ids.append(case_id)
        rows.append(
            {
                "case_id": case_id,
                "topology_class": case.get("topology_class"),
                "matrix_H": [[l11, m12], [m21, l22]],
                "mutual_mean_H": mutual,
                "reciprocity_relative_error": reciprocity,
                "determinant_H2": determinant,
                "coupling_abs": coupling,
                "checks": checks,
                "status": "ok" if all(checks.values()) else "needs_attention",
            }
        )

    strongest = max(rows, key=lambda row: row["coupling_abs"])
    unique_ids = len(set(ids)) == len(ids)
    strongest_ok = (
        expected_strongest_coupling_case is None
        or strongest["case_id"] == expected_strongest_coupling_case
    )
    family_checks = {
        "case_ids_unique": unique_ids,
        "all_case_matrices_valid": all(row["status"] == "ok" for row in rows),
        "expected_case_has_strongest_coupling": strongest_ok,
    }
    return {
        "policy": "inductance_matrix_family_gate_v1",
        "status": "ok" if all(family_checks.values()) else "needs_attention",
        "case_count": len(rows),
        "strongest_coupling_case": strongest["case_id"],
        "strongest_coupling_abs": strongest["coupling_abs"],
        "maximum_reciprocity_relative_error": max(
            row["reciprocity_relative_error"] for row in rows
        ),
        "checks": family_checks,
        "cases": rows,
        "lesson": (
            "Build a two-winding inductance matrix by exciting each winding separately. "
            "Gate reciprocal mutual terms, positive self terms, positive semidefiniteness, "
            "and |k|<=1 before comparing magnetic-circuit topologies."
        ),
    }
