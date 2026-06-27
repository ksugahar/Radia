"""Validation-class 2D Maxwell stress contour balance example.

Run:

    python validation_test/force_validation/validation_maxwell_contour_segment_balance.py

The example integrates a uniform magnetic field around a rectangular closed
contour.  Each segment has a nonzero stress contribution, but the closed
contour's net force cancels to zero.  This is a compact sign/orientation check
before using the same identity on FEM contour data.
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
    air_gap_maxwell_pressure,
    maxwell_contour_segment_balance_summary_2d,
)
from result_metadata import add_result_metadata  # noqa: E402


OUT_JSON = HERE / "validation_maxwell_contour_segment_balance_summary.json"
CONTOUR = [(-1.0, -0.5), (1.0, -0.5), (1.0, 0.5), (-1.0, 0.5)]
FIELD_T = (1.0, 0.0)


def build_summary() -> dict[str, object]:
    pressure = air_gap_maxwell_pressure(1.0)
    balance = maxwell_contour_segment_balance_summary_2d(
        CONTOUR,
        FIELD_T,
        expected_force_per_depth_N_per_m=(0.0, 0.0),
    )
    checks = {
        "n_segments": balance["n_segments"],
        "polygon_signed_area_m2": balance["polygon_signed_area_m2"],
        "pressure_Pa": pressure,
        "total_force_per_depth_N_per_m": balance["total_force_per_depth_N_per_m"],
        "sum_abs_normal_force_per_depth_N_per_m": balance["sum_abs_normal_force_per_depth_N_per_m"],
        "sum_abs_tangential_force_per_depth_N_per_m": balance["sum_abs_tangential_force_per_depth_N_per_m"],
        "cancellation_ratio": balance["cancellation_ratio"],
        "dominant_segment_index": balance["dominant_segment_index"],
        "status": balance["status"],
    }

    assert checks["n_segments"] == 4
    assert abs(float(checks["polygon_signed_area_m2"]) - 2.0) < 1.0e-12
    assert max(abs(value) for value in checks["total_force_per_depth_N_per_m"]) < 1.0e-9
    assert abs(float(checks["sum_abs_normal_force_per_depth_N_per_m"]) - 6.0 * pressure) < 1.0e-9
    assert abs(float(checks["sum_abs_tangential_force_per_depth_N_per_m"])) < 1.0e-12
    assert abs(float(checks["cancellation_ratio"])) < 1.0e-14
    assert checks["dominant_segment_index"] == 1
    assert checks["status"] == "ok"

    return {
        "kind": "maxwell_contour_segment_balance_validation",
        "validation_class": True,
        "learning_theme": (
            "Closed 2D Maxwell stress contours can have large local segment "
            "forces while the net force cancels by symmetry"
        ),
        "contour_vertices": CONTOUR,
        "uniform_B_T": FIELD_T,
        "checks": checks,
        "balance": balance,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    summary = add_result_metadata(build_summary(), __file__)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    checks = summary["checks"]
    print("[Maxwell contour segment balance]")
    print(
        f"  segments={checks['n_segments']} area={checks['polygon_signed_area_m2']:.12g} "
        f"status={checks['status']}"
    )
    print(
        f"  net_force={checks['total_force_per_depth_N_per_m']} "
        f"sum_abs_normal={checks['sum_abs_normal_force_per_depth_N_per_m']:.12g}"
    )
    print(f"  cancellation_ratio={checks['cancellation_ratio']:.3e}")
    print(f"[OK] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
