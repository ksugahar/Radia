"""Validation-class planar Lorentz block-force identity.

This lightweight example pins the 2D current-source force extraction used in
planar magnetostatics:

    F' = integral_A Jz zhat x B dA

where ``F'`` is force per out-of-plane depth [N/m].  The check uses the field of
one long round wire evaluated over a second small wire and compares the block
integral with Ampere's two-wire force law.

Run:

    python examples/electric_machine/validation_planar_lorentz_block_force.py
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

from radia_mcp.radia_ngsolve.force import MU0, planar_lorentz_force_summary  # noqa: E402
from radia_mcp.radia_ngsolve.solve import two_wire_force_per_length  # noqa: E402


OUT_JSON = HERE / "validation_planar_lorentz_block_force_summary.json"

WIRE_RADIUS_M = 0.002
WIRE_SPACING_M = 0.02
CURRENT_A = 1.0


def _norm2(values: list[float]) -> float:
    return math.hypot(values[0], values[1])


def _assert_close(actual: float, expected: float, *, rtol: float = 1.0e-12, atol: float = 1.0e-15) -> None:
    if abs(actual - expected) > max(atol, rtol * max(abs(actual), abs(expected))):
        raise AssertionError(f"{actual!r} != {expected!r}")


def build_summary() -> dict:
    area = math.pi * WIRE_RADIUS_M * WIRE_RADIUS_M
    jz = CURRENT_A / area
    by_from_left_wire = MU0 * CURRENT_A / (2.0 * math.pi * WIRE_SPACING_M)

    block = planar_lorentz_force_summary(
        jz,
        (0.0, by_from_left_wire),
        area_m2=area,
    )
    expected_magnitude = two_wire_force_per_length(CURRENT_A, CURRENT_A, WIRE_SPACING_M)
    expected_force = [-expected_magnitude, 0.0]
    force_error = [
        block["force_per_depth_N_per_m"][0] - expected_force[0],
        block["force_per_depth_N_per_m"][1] - expected_force[1],
    ]
    checks = {
        "wire_area_m2": area,
        "current_density_A_per_m2": jz,
        "external_By_T": by_from_left_wire,
        "expected_two_wire_force_magnitude_N_per_m": expected_magnitude,
        "expected_force_per_depth_N_per_m": expected_force,
        "block_force_per_depth_N_per_m": block["force_per_depth_N_per_m"],
        "force_abs_error_N_per_m": _norm2(force_error),
        "current_abs_error_A": abs(block["current_A"] - CURRENT_A),
        "mu0": MU0,
    }

    _assert_close(checks["current_abs_error_A"], 0.0)
    _assert_close(checks["force_abs_error_N_per_m"], 0.0)
    _assert_close(checks["expected_two_wire_force_magnitude_N_per_m"], 1.0e-5)

    return {
        "kind": "planar_lorentz_block_force",
        "validation_class": True,
        "force_learning": "2D Jz x B area integrals reproduce Ampere two-wire force per unit depth",
        "parameters": {
            "wire_radius_m": WIRE_RADIUS_M,
            "wire_spacing_m": WIRE_SPACING_M,
            "current_A": CURRENT_A,
        },
        "cases": {
            "right_wire_in_left_wire_field": block,
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
    print("[Planar Lorentz block force]")
    print(f"  external_By_T: {checks['external_By_T']:.12g}")
    print(f"  current_density_A_per_m2: {checks['current_density_A_per_m2']:.12g}")
    print(
        "  block_force_per_depth_N_per_m: "
        f"{checks['block_force_per_depth_N_per_m'][0]:.12g}, "
        f"{checks['block_force_per_depth_N_per_m'][1]:.12g}"
    )
    print(f"  force_abs_error_N_per_m: {checks['force_abs_error_N_per_m']:.3e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
