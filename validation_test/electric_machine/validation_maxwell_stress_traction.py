"""Validation-class Maxwell stress tensor / traction identities.

This lightweight force example pins the local identity underneath the FEM
surface and weighted-stress extractors:

    T = (B B - 0.5 |B|^2 I) / mu0,    traction = T n.

It checks the normal-field air-gap pressure, tangential-field magnetic tension,
and an oblique-field decomposition into normal and tangential traction.  Run:

    python validation_test/electric_machine/validation_maxwell_stress_traction.py
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
SRC = REPO / "packages" / "radia-mcp" / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from radia_mcp.radia_ngsolve.force import (  # noqa: E402
    MU0,
    air_gap_maxwell_pressure,
    maxwell_stress_tensor_air,
    maxwell_traction_summary,
)


OUT_JSON = HERE / "validation_maxwell_stress_traction_summary.json"

CASES = [
    {
        "name": "normal_1T",
        "B": [0.0, 0.0, 1.0],
        "normal": [0.0, 0.0, 1.0],
        "area_m2": 2.5e-4,
    },
    {
        "name": "tangential_1T",
        "B": [1.0, 0.0, 0.0],
        "normal": [0.0, 0.0, 1.0],
        "area_m2": 2.5e-4,
    },
    {
        "name": "oblique_3_4T",
        "B": [3.0, 4.0, 0.0],
        "normal": [1.0, 0.0, 0.0],
        "area_m2": 1.0e-6,
    },
]


def _dot(a: list[float], b: list[float]) -> float:
    return sum(x * y for x, y in zip(a, b))


def _norm(a: list[float]) -> float:
    return math.sqrt(_dot(a, a))


def _assert_close(actual, expected, rtol=1.0e-12, atol=1.0e-9) -> None:
    if isinstance(actual, list):
        if len(actual) != len(expected):
            raise AssertionError(f"length mismatch: {len(actual)} != {len(expected)}")
        for a, e in zip(actual, expected):
            _assert_close(a, e, rtol=rtol, atol=atol)
        return
    if abs(actual - expected) > max(atol, rtol * max(abs(actual), abs(expected))):
        raise AssertionError(f"{actual!r} != {expected!r}")


def _matrix_symmetry_error(matrix: list[list[float]]) -> float:
    return max(
        abs(matrix[i][j] - matrix[j][i])
        for i in range(len(matrix))
        for j in range(len(matrix))
    )


def build_rows() -> list[dict]:
    rows = []
    for case in CASES:
        tensor = maxwell_stress_tensor_air(case["B"])
        summary = maxwell_traction_summary(
            case["B"],
            case["normal"],
            area_m2=case["area_m2"],
        )
        b2 = _dot(case["B"], case["B"])
        dim = len(case["B"])
        trace_expected = (1.0 - 0.5 * dim) * b2 / MU0
        trace = sum(tensor[i][i] for i in range(dim))
        tangential = summary["tangential_traction_Pa"]
        row = {
            "name": case["name"],
            "B": case["B"],
            "normal": summary["normal"],
            "area_m2": case["area_m2"],
            "stress_tensor_Pa": tensor,
            "traction": summary,
            "tensor_symmetry_abs_error": _matrix_symmetry_error(tensor),
            "trace_identity_abs_error": abs(trace - trace_expected),
            "normal_identity_abs_error": abs(
                summary["normal_traction_Pa"]
                - summary["normal_traction_identity_Pa"]
            ),
            "tangential_identity_abs_error": abs(
                summary["tangential_traction_magnitude_Pa"]
                - abs(summary["B_normal_T"] * summary["B_tangent_T"]) / MU0
            ),
            "force_magnitude_N": _norm(summary["force_N"]),
            "tangential_force_magnitude_N": case["area_m2"] * _norm(tangential),
        }
        rows.append(row)
    return rows


def validate(rows: list[dict]) -> dict:
    by_name = {row["name"]: row for row in rows}
    pressure = air_gap_maxwell_pressure(1.0)
    checks = {
        "pressure_at_1T_Pa": pressure,
        "normal_1T_normal_traction_Pa": by_name["normal_1T"]["traction"]["normal_traction_Pa"],
        "normal_1T_force_N": by_name["normal_1T"]["traction"]["force_N"],
        "tangential_1T_normal_traction_Pa": by_name["tangential_1T"]["traction"]["normal_traction_Pa"],
        "oblique_normal_traction_Pa": by_name["oblique_3_4T"]["traction"]["normal_traction_Pa"],
        "oblique_tangential_traction_magnitude_Pa": (
            by_name["oblique_3_4T"]["traction"]["tangential_traction_magnitude_Pa"]
        ),
        "max_tensor_symmetry_abs_error": max(row["tensor_symmetry_abs_error"] for row in rows),
        "max_trace_identity_abs_error": max(row["trace_identity_abs_error"] for row in rows),
        "max_normal_identity_abs_error": max(row["normal_identity_abs_error"] for row in rows),
        "max_tangential_identity_abs_error": max(row["tangential_identity_abs_error"] for row in rows),
    }

    _assert_close(checks["normal_1T_normal_traction_Pa"], pressure)
    _assert_close(checks["normal_1T_force_N"], [0.0, 0.0, pressure * CASES[0]["area_m2"]])
    _assert_close(checks["tangential_1T_normal_traction_Pa"], -pressure)
    _assert_close(checks["oblique_normal_traction_Pa"], -3.5 / MU0)
    _assert_close(checks["oblique_tangential_traction_magnitude_Pa"], 12.0 / MU0)
    _assert_close(checks["max_tensor_symmetry_abs_error"], 0.0)
    _assert_close(checks["max_trace_identity_abs_error"], 0.0)
    _assert_close(checks["max_normal_identity_abs_error"], 0.0)
    _assert_close(checks["max_tangential_identity_abs_error"], 0.0)
    return checks


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    rows = build_rows()
    checks = validate(rows)
    summary = {
        "kind": "maxwell_stress_traction_identities",
        "validation_class": True,
        "rows": rows,
        "checks": checks,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print("[Maxwell stress tensor / traction identities]")
    for row in rows:
        traction = row["traction"]
        print(
            f"  {row['name']}: normal={traction['normal_traction_Pa']:.6g} Pa, "
            f"tangent={traction['tangential_traction_magnitude_Pa']:.6g} Pa, "
            f"|F|={row['force_magnitude_N']:.6g} N"
        )
    print("[checks]")
    for key, value in checks.items():
        print(f"  {key}: {value}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
