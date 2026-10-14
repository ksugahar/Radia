"""Two-conductor Maxwell and mutual capacitance matrix identity gate."""
from __future__ import annotations

import math


def _matrix2(value, name):
    if not isinstance(value, list) or len(value) != 2 or any(not isinstance(row, list) or len(row) != 2 for row in value):
        raise ValueError(f"{name} must be a 2x2 array")
    matrix = [[float(item) for item in row] for row in value]
    if not all(math.isfinite(item) for row in matrix for item in row):
        raise ValueError(f"{name} must contain finite values")
    return matrix


def two_conductor_capacitance_matrix_gate(summary, *, max_symmetry_relative_error=1.0e-8, max_representation_relative_error=1.0e-8):
    """Check reciprocity, passivity, and Maxwell-to-mutual conversion."""
    if not isinstance(summary, dict):
        raise ValueError("summary must be a mapping")
    maxwell = _matrix2(summary.get("maxwell_matrix"), "maxwell_matrix")
    mutual = _matrix2(summary.get("mutual_matrix"), "mutual_matrix")
    scale = max(max(abs(item) for row in maxwell for item in row), 1.0e-300)
    symmetry_error = abs(maxwell[0][1] - maxwell[1][0]) / scale
    mutual_symmetry_error = abs(mutual[0][1] - mutual[1][0]) / scale
    expected = [[maxwell[0][0] + maxwell[0][1], -maxwell[0][1]], [-maxwell[1][0], maxwell[1][1] + maxwell[1][0]]]
    representation_error = max(abs(mutual[i][j] - expected[i][j]) for i in range(2) for j in range(2)) / scale
    determinant = maxwell[0][0] * maxwell[1][1] - maxwell[0][1] * maxwell[1][0]
    checks = {
        "unit_explicit": str(summary.get("capacitance_unit") or "").strip() in {"F", "pF", "nF", "uF"},
        "maxwell_reciprocal": symmetry_error <= float(max_symmetry_relative_error),
        "mutual_reciprocal": mutual_symmetry_error <= float(max_symmetry_relative_error),
        "maxwell_diagonal_positive": maxwell[0][0] > 0.0 and maxwell[1][1] > 0.0,
        "maxwell_offdiagonal_nonpositive": maxwell[0][1] <= 0.0 and maxwell[1][0] <= 0.0,
        "maxwell_positive_definite": determinant > 0.0,
        "ground_capacitances_positive": expected[0][0] > 0.0 and expected[1][1] > 0.0,
        "mutual_entries_nonnegative": all(item >= 0.0 for row in mutual for item in row),
        "representations_agree": representation_error <= float(max_representation_relative_error),
    }
    return {"policy":"two_conductor_capacitance_matrix_gate_v1","status":"ok" if all(checks.values()) else "needs_attention","checks":checks,"issues":[name for name,ok in checks.items() if not ok],"metrics":{"maxwell_symmetry_relative_error":symmetry_error,"mutual_symmetry_relative_error":mutual_symmetry_error,"representation_relative_error":representation_error,"maxwell_determinant":determinant,"ground_capacitance_1":expected[0][0],"ground_capacitance_2":expected[1][1],"coupling_capacitance":0.5*(expected[0][1]+expected[1][0])},"expected_mutual_matrix":expected,"lesson":"For two conductors, the reciprocal Maxwell matrix has positive diagonal and nonpositive coupling; the mutual representation uses positive coupling and Maxwell row sums as capacitances to ground."}
