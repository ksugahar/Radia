"""Validation-class P1 surface-triangle constant traction load.

For one flat P1 boundary triangle with area A and constant vector traction t,

    F_e = A t,       f_i = F_e / 3

The equal nodal loads preserve both resultant force and moment because the
triangle centroid is the mean of the three vertices.  This is the small
readable boundary-load block used by FEM/BEM coupling examples before replacing
the prescribed traction with Maxwell stress.

Run:

    python validation_test/fem_readable/validation_surface_triangle_constant_traction_load.py
"""

from __future__ import annotations

import argparse
import json
from result_metadata import add_result_metadata
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
SRC = REPO / "packages" / "radia-mcp" / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from radia_mcp.radia_ngsolve.force import surface_triangle_constant_traction_load_summary  # noqa: E402


OUT_JSON = HERE / "validation_surface_triangle_constant_traction_load_summary.json"

TRIANGLE = [
    (1.0, 0.0, 0.0),
    (1.0, 2.0, 0.0),
    (1.0, 0.0, 3.0),
]
TRACTION_N_PER_M2 = (2.0, -1.0, 4.0)


def _assert_close_vector(actual, expected, *, atol: float = 1.0e-12) -> None:
    if len(actual) != len(expected):
        raise AssertionError(f"length mismatch: {len(actual)} != {len(expected)}")
    for a, e in zip(actual, expected):
        if abs(a - e) > atol:
            raise AssertionError(f"{actual!r} != {expected!r}")


def build_summary() -> dict:
    row = surface_triangle_constant_traction_load_summary(
        TRIANGLE,
        TRACTION_N_PER_M2,
    )
    shifted = surface_triangle_constant_traction_load_summary(
        TRIANGLE,
        TRACTION_N_PER_M2,
        pivot_m=row["centroid_m"],
    )

    _assert_close_vector(row["integrated_force_N"], [6.0, -3.0, 12.0])
    _assert_close_vector(row["nodal_resultant"]["total_force"], row["integrated_force_N"])
    _assert_close_vector(row["patch_resultant"]["total_moment"], [11.0, -6.0, -7.0])
    _assert_close_vector(row["nodal_resultant"]["total_moment"], row["patch_resultant"]["total_moment"])
    _assert_close_vector(shifted["patch_resultant"]["total_moment"], [0.0, 0.0, 0.0])
    _assert_close_vector(shifted["nodal_resultant"]["total_moment"], [0.0, 0.0, 0.0])

    return {
        "kind": "p1_surface_triangle_constant_traction_load",
        "validation_class": True,
        "force_learning": "constant P1 boundary traction load preserves resultant force and moment",
        "triangle_vertices_m": TRIANGLE,
        "traction_N_per_m2": TRACTION_N_PER_M2,
        "summary": row,
        "centroid_pivot_summary": shifted,
        "checks": {
            "expected_area_m2": 3.0,
            "expected_integrated_force_N": [6.0, -3.0, 12.0],
            "expected_moment_about_origin_Nm": [11.0, -6.0, -7.0],
            "max_force_preservation_abs_error_N": max(row["force_preservation_abs_errors_N"]),
            "max_moment_preservation_abs_error_Nm": max(row["moment_preservation_abs_errors_Nm"]),
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    summary = build_summary()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(add_result_metadata(summary, __file__), indent=2, sort_keys=True) + "\n", encoding="utf-8")

    checks = summary["checks"]
    print("[P1 surface-triangle constant traction load]")
    print(f"  area_m2: {summary['summary']['area']:.12g}")
    print(f"  integrated_force_N: {summary['summary']['integrated_force_N']}")
    print(f"  moment_about_origin_Nm: {summary['summary']['patch_resultant']['total_moment']}")
    print(f"  max_force_error_N: {checks['max_force_preservation_abs_error_N']:.3e}")
    print(f"  max_moment_error_Nm: {checks['max_moment_preservation_abs_error_Nm']:.3e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
