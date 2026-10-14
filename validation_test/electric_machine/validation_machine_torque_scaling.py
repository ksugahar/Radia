"""Validation-class 2D sector torque to whole-machine torque scaling.

This example checks the arithmetic post-processing step used after a 2D motor
force/torque extraction:

    T_whole = T_2d * length_unit_m^2 * stack_length_m * symmetry_factor.

Run:

    python validation_test/electric_machine/validation_machine_torque_scaling.py
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

from radia_mcp.radia_ngsolve.machine_scaling import (  # noqa: E402
    MachineScaling,
    torque_scaling_summary,
)


OUT_JSON = HERE / "validation_machine_torque_scaling_summary.json"

TORQUE_2D = 2.5e6
STACK_LENGTH_M = 0.08
SYMMETRY_FACTOR = 6.0
LENGTH_UNIT_M = 1.0e-3


def _assert_close(actual: float, expected: float, *, rtol: float = 1.0e-12, atol: float = 1.0e-12) -> None:
    if abs(actual - expected) > max(atol, rtol * max(abs(actual), abs(expected))):
        raise AssertionError(f"{actual!r} != {expected!r}")


def build_summary() -> dict:
    row = torque_scaling_summary(
        TORQUE_2D,
        stack_length_m=STACK_LENGTH_M,
        symmetry_factor=SYMMETRY_FACTOR,
        length_unit_m=LENGTH_UNIT_M,
        label="one_sixth_sector",
    )
    direct = MachineScaling(
        stack_length=STACK_LENGTH_M,
        symmetry_factor=SYMMETRY_FACTOR,
        length_unit_m=LENGTH_UNIT_M,
    ).torque(TORQUE_2D)
    full_model = torque_scaling_summary(
        TORQUE_2D,
        stack_length_m=STACK_LENGTH_M,
        symmetry_factor=1.0,
        length_unit_m=LENGTH_UNIT_M,
        label="full_model_same_2d_value",
    )
    braking = torque_scaling_summary(
        -TORQUE_2D,
        stack_length_m=STACK_LENGTH_M,
        symmetry_factor=SYMMETRY_FACTOR,
        length_unit_m=LENGTH_UNIT_M,
        label="negative_torque_preserves_sign",
    )

    checks = {
        "expected_sector_torque_Nm": 0.2,
        "expected_whole_machine_torque_Nm": 1.2,
        "direct_machine_scaling_torque_Nm": direct,
        "summary_whole_machine_torque_Nm": row["whole_machine_torque_Nm"],
        "summary_minus_direct_abs_error_Nm": abs(row["whole_machine_torque_Nm"] - direct),
        "sector_to_whole_ratio": row["whole_over_sector"],
        "full_model_same_2d_value_torque_Nm": full_model["whole_machine_torque_Nm"],
        "sector_scaled_over_full_model": row["whole_machine_torque_Nm"] / full_model["whole_machine_torque_Nm"],
        "negative_torque_Nm": braking["whole_machine_torque_Nm"],
    }

    _assert_close(row["torque_one_model_sector_Nm"], checks["expected_sector_torque_Nm"])
    _assert_close(row["whole_machine_torque_Nm"], checks["expected_whole_machine_torque_Nm"])
    _assert_close(checks["summary_minus_direct_abs_error_Nm"], 0.0)
    _assert_close(checks["sector_to_whole_ratio"], SYMMETRY_FACTOR)
    _assert_close(checks["sector_scaled_over_full_model"], SYMMETRY_FACTOR)
    _assert_close(checks["negative_torque_Nm"], -checks["expected_whole_machine_torque_Nm"])

    return {
        "kind": "machine_torque_scaling",
        "validation_class": True,
        "force_learning": "2D per-depth sector torque scales by length_unit^2, stack length, and symmetry factor",
        "parameters": {
            "torque_2d_per_depth": TORQUE_2D,
            "stack_length_m": STACK_LENGTH_M,
            "symmetry_factor": SYMMETRY_FACTOR,
            "length_unit_m": LENGTH_UNIT_M,
        },
        "cases": {
            "one_sixth_sector": row,
            "full_model_same_2d_value": full_model,
            "negative_torque": braking,
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
    print("[Machine torque scaling]")
    print(f"  sector torque: {checks['expected_sector_torque_Nm']:.12g} N m")
    print(f"  whole torque:  {checks['summary_whole_machine_torque_Nm']:.12g} N m")
    print(f"  ratio:         {checks['sector_to_whole_ratio']:.12g}")
    print(f"  direct error:  {checks['summary_minus_direct_abs_error_Nm']:.3e} N m")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
