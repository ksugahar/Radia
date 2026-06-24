"""Validation-class build123d boundary normals against a Netgen .vol mesh.

The build123d side uses the analytic face normals of an axis-aligned box.  If
``--vol`` is provided, the script reads the exported triangular/tetrahedral
Netgen ``.vol`` mesh and compares each named boundary's scalar area and
oriented area vector.

This is the CAD-side companion to Maxwell-stress force integration: before
surface tractions are summed, the mesh boundary must preserve the intended
normal direction and area.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
SRC = REPO / "packages" / "radia-mcp" / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from build123d import Box  # noqa: E402

from radia_mcp.build123d.modeling import (  # noqa: E402
    box_face_vector_area_rows,
    compare_boundary_vector_area_rows,
    shape_measurement_row,
)
from radia_mcp.radia_ngsolve.netgen_vol import read_netgen_tri_tet_vol  # noqa: E402


SUMMARY_JSON = HERE / "validation_build123d_cubit_boundary_normals_summary.json"


def _compact_comparison_summary(rows: list[dict], measured_label: str) -> dict:
    return {
        "measured_label": measured_label,
        "n_cases": len(rows),
        "n_passed": sum(1 for row in rows if row["passed"]),
        "all_passed": all(row["passed"] for row in rows),
        "max_area_rel_error": max(
            (row["area_rel_error"] or 0.0 for row in rows),
            default=0.0,
        ),
        "max_vector_abs_error": max(
            (row["vector_abs_error"] or 0.0 for row in rows),
            default=0.0,
        ),
        "max_unit_normal_abs_error": max(
            (row["unit_normal_abs_error"] or 0.0 for row in rows),
            default=0.0,
        ),
    }


def build_summary(vol_path: Path | None, vector_atol: float, area_rtol: float) -> dict:
    t0 = time.perf_counter()
    size = (2.0, 3.0, 5.0)
    shape = Box(*size).solid()
    shape.label = "box_2x3x5"
    measurement = shape_measurement_row(shape)
    reference_rows = box_face_vector_area_rows(
        measurement["bounding_box"]["size"],
        center=measurement["bounding_box"]["center"],
    )

    summary: dict[str, object] = {
        "case": "axis_aligned_box_boundary_normals",
        "size": size,
        "build123d_measurement": measurement,
        "reference_boundary_rows": reference_rows,
        "reference_total_surface_area": sum(row["surface_area"] for row in reference_rows),
        "reference_closed_surface_vector_area": [
            sum(row["vector_area"][axis] for row in reference_rows)
            for axis in range(3)
        ],
        "vol": {
            "provided": vol_path is not None,
            "path": str(vol_path) if vol_path is not None else None,
            "comparison_rows": [],
            "comparison_summary": None,
        },
        "elapsed_seconds": None,
    }

    if vol_path is not None:
        mesh = read_netgen_tri_tet_vol(vol_path)
        measured_rows = list(mesh.boundary_normal_summary_rows())
        comparison_rows = compare_boundary_vector_area_rows(
            reference_rows,
            measured_rows,
            vector_atol=vector_atol,
            area_rtol=area_rtol,
            measured_label="netgen_vol",
        )
        summary["vol"] = {
            "provided": True,
            "path": str(vol_path),
            "points": len(mesh.points),
            "tetrahedra": len(mesh.tetrahedra),
            "surface_triangles": len(mesh.surface_triangles),
            "boundary_rows": measured_rows,
            "comparison_rows": comparison_rows,
            "comparison_summary": _compact_comparison_summary(comparison_rows, "netgen_vol"),
        }

    summary["elapsed_seconds"] = time.perf_counter() - t0
    return summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vol", type=Path, default=None, help="optional Netgen .vol mesh")
    parser.add_argument("--out", type=Path, default=SUMMARY_JSON)
    parser.add_argument("--vector-atol", type=float, default=1.0e-8)
    parser.add_argument("--area-rtol", type=float, default=1.0e-8)
    args = parser.parse_args(argv)

    summary = build_summary(args.vol, vector_atol=args.vector_atol, area_rtol=args.area_rtol)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    vol_summary = summary["vol"].get("comparison_summary") if isinstance(summary["vol"], dict) else None
    if vol_summary is not None and not vol_summary["all_passed"]:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
