"""Validation-class electrostatic Maxwell stress / traction identities.

For a field normal to a conductor/dielectric interface, the electrostatic
Maxwell stress gives the familiar pressure

    p = 0.5 * eps0 * E^2.

This example mirrors the magnetic Maxwell-traction validation in a form that is
easy to compare with electrostatic FEM post-processing.

Run:

    python validation_test/electrostatics/validation_electrostatic_maxwell_traction.py
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
SRC = REPO / "packages" / "radia-mcp" / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from radia_mcp.radia_ngsolve.force import (  # noqa: E402
    EPS0,
    electrostatic_stress_tensor,
    electrostatic_traction_summary,
)
from result_metadata import add_result_metadata  # noqa: E402


OUT_JSON = HERE / "validation_electrostatic_maxwell_traction_summary.json"

CASES = [
    {
        "name": "normal_1MV_per_m",
        "E": [0.0, 0.0, 1.0e6],
        "normal": [0.0, 0.0, 1.0],
        "area_m2": 0.25,
    },
    {
        "name": "tangential_1MV_per_m",
        "E": [1.0e6, 0.0, 0.0],
        "normal": [0.0, 0.0, 1.0],
        "area_m2": 0.25,
    },
    {
        "name": "oblique_3_4MV_per_m",
        "E": [3.0e6, 4.0e6, 0.0],
        "normal": [1.0, 0.0, 0.0],
        "area_m2": 1.0e-4,
    },
]


def _assert_close(value: float, expected: float, rel: float = 1.0e-12, abs_tol: float = 1.0e-12) -> None:
    if abs(value - expected) > max(abs_tol, rel * max(abs(value), abs(expected), 1.0)):
        raise AssertionError(f"{value!r} != {expected!r}")


def build_rows() -> list[dict]:
    rows = []
    for case in CASES:
        tensor = electrostatic_stress_tensor(case["E"])
        traction = electrostatic_traction_summary(
            case["E"],
            case["normal"],
            area_m2=case["area_m2"],
        )
        rows.append({
            "name": case["name"],
            "E": case["E"],
            "normal": traction["normal"],
            "area_m2": case["area_m2"],
            "stress_tensor_Pa": tensor,
            "traction": traction,
        })
    return rows


def validate(rows: list[dict]) -> dict:
    by_name = {row["name"]: row for row in rows}
    pressure_1mv = 0.5 * EPS0 * 1.0e12
    checks = {
        "pressure_1MV_per_m_Pa": pressure_1mv,
        "normal_traction_Pa": by_name["normal_1MV_per_m"]["traction"]["normal_traction_Pa"],
        "normal_force_N": by_name["normal_1MV_per_m"]["traction"]["force_N"],
        "tangential_normal_traction_Pa": (
            by_name["tangential_1MV_per_m"]["traction"]["normal_traction_Pa"]
        ),
        "oblique_normal_traction_Pa": by_name["oblique_3_4MV_per_m"]["traction"]["normal_traction_Pa"],
        "oblique_tangential_traction_magnitude_Pa": (
            by_name["oblique_3_4MV_per_m"]["traction"]["tangential_traction_magnitude_Pa"]
        ),
    }
    _assert_close(checks["normal_traction_Pa"], pressure_1mv)
    _assert_close(checks["normal_force_N"][2], 0.25 * pressure_1mv)
    _assert_close(checks["tangential_normal_traction_Pa"], -pressure_1mv)
    _assert_close(checks["oblique_normal_traction_Pa"], 0.5 * EPS0 * (9.0e12 - 16.0e12))
    _assert_close(checks["oblique_tangential_traction_magnitude_Pa"], EPS0 * 12.0e12)
    return checks


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    rows = build_rows()
    checks = validate(rows)
    summary = {
        "kind": "electrostatic_maxwell_traction_identities",
        "validation_class": True,
        "eps0": EPS0,
        "rows": rows,
        "checks": checks,
    }
    summary = add_result_metadata(summary, __file__)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print("[Electrostatic Maxwell traction identities]")
    for row in rows:
        traction = row["traction"]
        print(
            f"  {row['name']}: normal={traction['normal_traction_Pa']:.6g} Pa, "
            f"tangent={traction['tangential_traction_magnitude_Pa']:.6g} Pa"
        )
    print("[checks]")
    for key, value in checks.items():
        print(f"  {key}: {value}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
