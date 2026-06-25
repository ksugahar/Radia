"""Validation-class 2D Maxwell-stress contour force.

This lightweight example pins the line-integral form often used in planar
magnetostatic post-processing:

    F' = integral_C T(B) n ds

where ``F'`` is force per out-of-plane depth [N/m].  It checks a single pole
face against the air-gap pressure identity and a closed rectangular contour in
a uniform field, whose net force must cancel to zero.

Run:

    python examples/electric_machine/validation_planar_maxwell_contour_force.py
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
    maxwell_contour_force_2d,
    maxwell_line_segment_force_2d,
)


OUT_JSON = HERE / "validation_planar_maxwell_contour_force_summary.json"


def _norm2(values: list[float]) -> float:
    return math.hypot(values[0], values[1])


def _assert_close(actual: float, expected: float, *, rtol: float = 1.0e-12, atol: float = 1.0e-9) -> None:
    if abs(actual - expected) > max(atol, rtol * max(abs(actual), abs(expected))):
        raise AssertionError(f"{actual!r} != {expected!r}")


def build_summary() -> dict:
    pressure = air_gap_maxwell_pressure(1.0)

    pole_face = maxwell_line_segment_force_2d(
        (0.0, -0.5),
        (0.0, 0.5),
        (1.0, 0.0),
        normal_side="right",
    )
    rectangle = [(-1.0, -0.5), (1.0, -0.5), (1.0, 0.5), (-1.0, 0.5)]
    closed_uniform = maxwell_contour_force_2d(rectangle, (1.0, 0.0), orientation="ccw")

    checks = {
        "pressure_at_1T_Pa": pressure,
        "pole_face_force_per_depth_N_per_m": pole_face["force_per_depth_N_per_m"],
        "pole_face_expected_force_per_depth_N_per_m": [pressure, 0.0],
        "pole_face_force_abs_error": _norm2([
            pole_face["force_per_depth_N_per_m"][0] - pressure,
            pole_face["force_per_depth_N_per_m"][1],
        ]),
        "closed_contour_total_force_per_depth_N_per_m": closed_uniform["total_force_per_depth_N_per_m"],
        "closed_contour_total_force_magnitude_per_depth_N_per_m": (
            closed_uniform["total_force_magnitude_per_depth_N_per_m"]
        ),
        "closed_contour_signed_area_m2": closed_uniform["polygon_signed_area_m2"],
        "closed_contour_sum_abs_normal_force_per_depth_N_per_m": (
            closed_uniform["sum_abs_normal_force_per_depth_N_per_m"]
        ),
        "mu0": MU0,
    }

    _assert_close(checks["pole_face_force_abs_error"], 0.0)
    _assert_close(checks["closed_contour_total_force_magnitude_per_depth_N_per_m"], 0.0)
    _assert_close(checks["closed_contour_signed_area_m2"], 2.0)

    return {
        "kind": "planar_maxwell_contour_force",
        "validation_class": True,
        "force_learning": "2D contour Maxwell stress integrates T n ds as force per unit depth",
        "cases": {
            "pole_face": pole_face,
            "closed_uniform_rectangle": closed_uniform,
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
    print("[Planar Maxwell contour force]")
    print(f"  pressure_at_1T_Pa: {checks['pressure_at_1T_Pa']:.12g}")
    print(f"  pole_face_force_abs_error: {checks['pole_face_force_abs_error']:.3e}")
    print(
        "  closed_contour_total_force_magnitude_per_depth_N_per_m: "
        f"{checks['closed_contour_total_force_magnitude_per_depth_N_per_m']:.3e}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
