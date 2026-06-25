"""Validation-class two-wire Lorentz/virtual-work force identity.

This example pins a common magnetostatic post-processing lesson without
depending on any commercial solver:

    Lorentz force on a current filament == d(coenergy)/d(separation)

For two long parallel wires, the separation-dependent mutual coenergy per unit
length is ``-mu0 I1 I2 log(d/d_ref)/(2*pi)``.  Its derivative is the radial
force per unit length.  Like currents therefore have a negative force in the
increasing-separation coordinate, i.e. attraction.

Run:

    python examples/force_validation/validation_parallel_wire_virtual_work_force.py
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

from radia_mcp.radia_ngsolve.force import parallel_wire_virtual_work_force_summary  # noqa: E402


OUT_JSON = HERE / "validation_parallel_wire_virtual_work_force_summary.json"

CASES = [
    {
        "label": "like_currents_horizontal_attract",
        "current1_A": 3.0,
        "current2_A": 2.0,
        "separation_xy_m": (0.035, 0.0),
        "expected_interaction": "attraction",
    },
    {
        "label": "opposite_currents_vertical_repel",
        "current1_A": 3.0,
        "current2_A": -2.0,
        "separation_xy_m": (0.0, 0.035),
        "expected_interaction": "repulsion",
    },
]


def _case_record(case: dict[str, object]) -> dict[str, object]:
    sx, sy = case["separation_xy_m"]
    separation = (sx * sx + sy * sy) ** 0.5
    row = parallel_wire_virtual_work_force_summary(
        case["current1_A"],
        case["current2_A"],
        case["separation_xy_m"],
        displacement_step_m=separation * 1.0e-4,
    )
    return {
        "label": case["label"],
        "expected_interaction": case["expected_interaction"],
        "summary": row,
    }


def build_summary() -> dict[str, object]:
    records = [_case_record(case) for case in CASES]
    max_rel_error = max(record["summary"]["force_rel_error"] for record in records)
    max_abs_error = max(record["summary"]["force_vector_abs_error_N_per_m"] for record in records)
    interactions_ok = all(
        record["summary"]["interaction"] == record["expected_interaction"]
        for record in records
    )
    directions_ok = (
        records[0]["summary"]["virtual_work_radial_force_per_length_N_per_m"] < 0.0
        and records[1]["summary"]["virtual_work_radial_force_per_length_N_per_m"] > 0.0
    )

    checks = {
        "n_cases": len(records),
        "max_force_rel_error": max_rel_error,
        "max_force_vector_abs_error_N_per_m": max_abs_error,
        "interactions_ok": interactions_ok,
        "directions_ok": directions_ok,
        "passed": max_rel_error < 1.0e-8 and interactions_ok and directions_ok,
    }
    assert checks["passed"]

    return {
        "kind": "parallel_wire_virtual_work_force_validation",
        "validation_class": True,
        "force_learning": (
            "fixed-current coenergy derivative and Lorentz force give the same "
            "parallel-wire force per unit length"
        ),
        "checks": checks,
        "cases": records,
    }


def _json_clean(value):
    if isinstance(value, float):
        return 0.0 if value == 0.0 else value
    if isinstance(value, list):
        return [_json_clean(item) for item in value]
    if isinstance(value, tuple):
        return [_json_clean(item) for item in value]
    if isinstance(value, dict):
        return {key: _json_clean(item) for key, item in value.items()}
    return value


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    summary = _json_clean(build_summary())
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    checks = summary["checks"]
    print("[parallel wire virtual work force]")
    print(f"  cases={checks['n_cases']}, max_rel_error={checks['max_force_rel_error']:.3e}")
    for record in summary["cases"]:
        row = record["summary"]
        print(
            f"  {record['label']}: interaction={row['interaction']}, "
            f"radial_virtual={row['virtual_work_radial_force_per_length_N_per_m']:.12g}, "
            f"radial_analytic={row['analytic_radial_force_per_length_N_per_m']:.12g}"
        )
    print(f"[OK] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
