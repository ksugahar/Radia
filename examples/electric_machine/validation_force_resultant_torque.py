"""Validation-class force-row resultant and torque identities.

This example checks the common post-processing step after force extraction:
collapse a table of point forces into a net force and a torque about a pivot,

    F = sum_i F_i,
    M_p = sum_i (r_i - p) x F_i.

The force rows can come from Maxwell-stress surface patches, Lorentz body-force
elements, pressure faces, or nodal loads.  The algebra here is deliberately
small enough to inspect before using a full field solve.

Run:

    python examples/electric_machine/validation_force_resultant_torque.py
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

from radia_mcp.radia_ngsolve.force import force_moment_resultant_summary  # noqa: E402


OUT_JSON = HERE / "validation_force_resultant_torque_summary.json"

RADIUS_M = 0.05
FORCE_N = 10.0


def _assert_close(actual: float, expected: float, *, rtol: float = 1.0e-12, atol: float = 1.0e-12) -> None:
    if abs(actual - expected) > max(atol, rtol * max(abs(actual), abs(expected))):
        raise AssertionError(f"{actual!r} != {expected!r}")


def build_summary() -> dict:
    couple = force_moment_resultant_summary(
        [(RADIUS_M, 0.0), (-RADIUS_M, 0.0)],
        [(0.0, FORCE_N), (0.0, -FORCE_N)],
    )
    shifted_pivot = force_moment_resultant_summary(
        [(RADIUS_M, 0.0), (-RADIUS_M, 0.0)],
        [(0.0, FORCE_N), (0.0, -FORCE_N)],
        pivot_m=(0.2, -0.1),
    )
    single_3d = force_moment_resultant_summary(
        [(0.0, 2.0 * RADIUS_M, 0.0)],
        [(FORCE_N, 0.0, 0.0)],
    )

    expected_couple = 2.0 * RADIUS_M * FORCE_N
    expected_single = -2.0 * RADIUS_M * FORCE_N
    checks = {
        "expected_couple_torque": expected_couple,
        "couple_total_force": couple["total_force"],
        "couple_total_moment": couple["total_moment"],
        "couple_torque_abs_error": abs(couple["total_moment"] - expected_couple),
        "shifted_pivot_torque_abs_error": abs(shifted_pivot["total_moment"] - expected_couple),
        "single_3d_total_moment": single_3d["total_moment"],
        "single_3d_mz_abs_error": abs(single_3d["total_moment"][2] - expected_single),
    }

    _assert_close(couple["total_force_magnitude"], 0.0)
    _assert_close(checks["couple_torque_abs_error"], 0.0)
    _assert_close(checks["shifted_pivot_torque_abs_error"], 0.0)
    _assert_close(checks["single_3d_mz_abs_error"], 0.0)

    return {
        "kind": "force_resultant_torque",
        "validation_class": True,
        "force_learning": "discrete force rows collapse to net force and pivot torque by sum(r x F)",
        "parameters": {
            "radius_m": RADIUS_M,
            "force_N": FORCE_N,
        },
        "cases": {
            "two_dimensional_force_couple": couple,
            "shifted_pivot_force_couple": shifted_pivot,
            "single_three_dimensional_force": single_3d,
        },
        "checks": checks,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    summary = build_summary()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    checks = summary["checks"]
    print("[Force resultant and torque]")
    print(f"  expected_couple_torque: {checks['expected_couple_torque']:.12g}")
    print(f"  couple_total_moment: {checks['couple_total_moment']:.12g}")
    print(f"  couple_torque_abs_error: {checks['couple_torque_abs_error']:.3e}")
    print(f"  shifted_pivot_torque_abs_error: {checks['shifted_pivot_torque_abs_error']:.3e}")
    print(f"  single_3d_mz_abs_error: {checks['single_3d_mz_abs_error']:.3e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
