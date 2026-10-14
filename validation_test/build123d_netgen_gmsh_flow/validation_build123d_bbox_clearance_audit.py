"""Validation-class build123d bbox clearance audit.

Run:

    python validation_test/build123d_netgen_gmsh_flow/validation_build123d_bbox_clearance_audit.py

This example builds three labelled boxes before meshing.  Two pairs are
provably separated by at least one bounding-box axis.  One pair has overlapping
bounding boxes and is therefore marked for precise geometry or boolean checks.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from build123d import Box, Pos


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
SRC = REPO / "packages" / "radia-mcp" / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from radia_mcp.build123d.modeling import (  # noqa: E402
    assembly,
    shape_bbox_pair_clearance_summary,
    shape_measurement_rows,
)


OUT_JSON = HERE / "validation_build123d_bbox_clearance_audit_summary.json"


def _box(name: str, x: float):
    shape = (Pos(x, 0, 0) * Box(1, 1, 1)).solid()
    shape.label = name
    return shape


def build_summary() -> dict[str, object]:
    left = _box("left", 0.0)
    right = _box("right", 2.0)
    overlap = _box("overlap", 0.4)
    rows = shape_measurement_rows(assembly(left, right, overlap, label="bbox_audit"))
    audit = shape_bbox_pair_clearance_summary(rows)
    pairs = {row["pair"]: row for row in audit["pair_rows"]}

    checks = {
        "n_shapes": audit["n_shapes"],
        "n_pairs": audit["n_pairs"],
        "separated_pair_count": audit["separated_pair_count"],
        "bbox_overlap_pair_count": audit["bbox_overlap_pair_count"],
        "touching_pair_count": audit["touching_pair_count"],
        "min_positive_gap": audit["min_positive_gap"],
        "left_right_x_gap": pairs["left::right"]["axis_gaps"]["x"],
        "right_overlap_x_gap": pairs["right::overlap"]["axis_gaps"]["x"],
        "left_overlap_intersection_volume": pairs["left::overlap"]["bbox_intersection_volume"],
        "status": audit["status"],
    }

    assert checks["n_shapes"] == 3
    assert checks["n_pairs"] == 3
    assert checks["separated_pair_count"] == 2
    assert checks["bbox_overlap_pair_count"] == 1
    assert checks["touching_pair_count"] == 0
    assert abs(float(checks["min_positive_gap"]) - 0.6) < 1.0e-12
    assert abs(float(checks["left_right_x_gap"]) - 1.0) < 1.0e-12
    assert abs(float(checks["right_overlap_x_gap"]) - 0.6) < 1.0e-12
    assert abs(float(checks["left_overlap_intersection_volume"]) - 0.6) < 1.0e-12
    assert checks["status"] == "needs_attention"

    return {
        "kind": "build123d_bbox_clearance_audit_validation",
        "validation_class": True,
        "learning_theme": (
            "A bounding-box gap on one axis proves separation; all-axis bbox "
            "overlap only requests a precise geometry check"
        ),
        "checks": checks,
        "measurement_rows": rows,
        "audit": audit,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    summary = build_summary()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    checks = summary["checks"]
    print("[build123d bbox clearance audit]")
    print(
        f"  shapes={checks['n_shapes']} pairs={checks['n_pairs']} "
        f"separated={checks['separated_pair_count']} overlap={checks['bbox_overlap_pair_count']}"
    )
    print(
        f"  min_positive_gap={checks['min_positive_gap']:.12g} "
        f"left_overlap_intersection_volume={checks['left_overlap_intersection_volume']:.12g}"
    )
    print(f"  status={checks['status']}")
    print(f"[OK] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
