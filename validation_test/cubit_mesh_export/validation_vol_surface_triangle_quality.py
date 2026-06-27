"""Validation-class Netgen .vol surface triangle quality example.

Coreform Cubit 2026.6 emphasizes stronger triangle/tet meshing and richer
quality assessment.  This public-safe example keeps the downstream side honest:
after a Cubit/Coreform ``export netgen`` run, radia can independently inspect
the boundary triangles that become scalar-BEM/RWG trace elements.

Run:

    python validation_test/cubit_mesh_export/validation_vol_surface_triangle_quality.py
    python validation_test/cubit_mesh_export/validation_vol_surface_triangle_quality.py --vol C:\\temp\\mesh.vol
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

from radia_mcp.radia_ngsolve.netgen_vol import parse_netgen_tri_tet_vol, read_netgen_tri_tet_vol  # noqa: E402


OUT_JSON = HERE / "validation_vol_surface_triangle_quality_summary.json"


UNIT_TET_VOL = """\
mesh3d
dimension
3
geomtype
0
facedescriptors
1
1 1 0 1 1
surfaceelements
4
1 1 1 0 3 1 2 3
1 1 1 0 3 1 4 2
1 1 1 0 3 1 3 4
1 1 1 0 3 2 4 3
volumeelements
1
1 4 1 2 3 4
points
4
0 0 0
1 0 0
0 1 0
0 0 1
pointelements
0
materials
1
1 air
bcnames
1
1 outer
endmesh
"""


SLIVER_SURFACE_VOL = """\
mesh3d
dimension
3
geomtype
0
facedescriptors
1
1 1 0 1 1
surfaceelements
1
1 1 1 0 3 1 2 3
volumeelements
0
points
3
0 0 0
1 0 0
0.999 0.000001 0
pointelements
0
materials
0
bcnames
1
1 sliver_patch
endmesh
"""


def _record(label: str, text: str) -> dict[str, object]:
    mesh = parse_netgen_tri_tet_vol(text, source=f"{label}.vol")
    return {
        "label": label,
        "mesh_summary": mesh.summary(),
        "quality_rows": list(mesh.surface_triangle_quality_rows()),
        "quality_summary": mesh.surface_triangle_quality_summary(),
        "boundary_rows": list(mesh.boundary_summary_rows()),
        "closure": mesh.surface_closure_summary(),
    }


def _external_record(path: Path) -> dict[str, object]:
    mesh = read_netgen_tri_tet_vol(path)
    return {
        "source_name": path.name,
        "mesh_summary": mesh.summary(),
        "bounding_box": mesh.bounding_box(),
        "surface_quality_summary": mesh.surface_triangle_quality_summary(),
        "tet_quality_summary": mesh.tetrahedron_quality_summary(),
        "boundary_rows": list(mesh.boundary_summary_rows()),
        "closure": mesh.surface_closure_summary(),
        "manifold": mesh.surface_edge_manifold_summary(),
    }


def build_summary(external_vol: Path | None = None) -> dict[str, object]:
    unit = _record("unit_tetrahedron_boundary", UNIT_TET_VOL)
    sliver = _record("single_sliver_patch", SLIVER_SURFACE_VOL)

    right_quality = 2.0 * (1.0 / (2.0 + math.sqrt(2.0))) / (math.sqrt(2.0) / 2.0)
    unit_quality = unit["quality_summary"]
    sliver_quality = sliver["quality_summary"]
    checks = {
        "unit_min_radius_ratio_quality": unit_quality["min_radius_ratio_quality"],
        "unit_expected_right_face_quality": right_quality,
        "unit_max_radius_ratio_quality": unit_quality["max_radius_ratio_quality"],
        "unit_max_edge_ratio": unit_quality["max_edge_ratio"],
        "unit_min_angle_deg": unit_quality["min_angle_deg"],
        "unit_max_angle_deg": unit_quality["max_angle_deg"],
        "sliver_min_radius_ratio_quality": sliver_quality["min_radius_ratio_quality"],
        "sliver_max_edge_ratio": sliver_quality["max_edge_ratio"],
        "sliver_min_angle_deg": sliver_quality["min_angle_deg"],
    }

    assert abs(checks["unit_min_radius_ratio_quality"] - right_quality) < 1.0e-15
    assert abs(checks["unit_max_radius_ratio_quality"] - 1.0) < 1.0e-15
    assert abs(checks["unit_max_edge_ratio"] - math.sqrt(2.0)) < 1.0e-15
    assert abs(checks["unit_min_angle_deg"] - 45.0) < 1.0e-12
    assert abs(checks["unit_max_angle_deg"] - 90.0) < 1.0e-12
    assert checks["sliver_min_radius_ratio_quality"] < 1.0e-5
    assert checks["sliver_max_edge_ratio"] > 900.0
    assert checks["sliver_min_angle_deg"] < 1.0e-3

    summary: dict[str, object] = {
        "kind": "netgen_vol_surface_triangle_quality_validation",
        "validation_class": True,
        "quality_definition": (
            "radius_ratio_quality = 2 * inradius / circumradius; "
            "1.0 is equilateral, near 0 marks sliver or degenerate boundary triangles"
        ),
        "release_learning": (
            "Coreform Cubit 2026.6 highlights triangle/tet meshing robustness and "
            "quality assessment; radia independently checks .vol boundary triangles "
            "before FEM/BEM trace use."
        ),
        "builtin_cases": [unit, sliver],
        "checks": checks,
    }
    if external_vol is not None:
        summary["external_vol"] = _external_record(external_vol)
    return summary


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vol", type=Path, default=None, help="Optional real Netgen .vol export to evaluate")
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    summary = build_summary(args.vol)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print("[surface triangle quality]")
    for case in summary["builtin_cases"]:
        q = case["quality_summary"]
        print(
            f"  {case['label']}: triangles={q['surface_triangles']}, "
            f"min_quality={q['min_radius_ratio_quality']:.12g}, "
            f"max_edge_ratio={q['max_edge_ratio']:.6g}, "
            f"angle_range=[{q['min_angle_deg']:.6g}, {q['max_angle_deg']:.6g}] deg"
        )
    if "external_vol" in summary:
        q = summary["external_vol"]["surface_quality_summary"]
        print(
            f"  external {summary['external_vol']['source_name']}: "
            f"triangles={q['surface_triangles']}, "
            f"min_quality={q['min_radius_ratio_quality']:.6g}, "
            f"max_edge_ratio={q['max_edge_ratio']:.6g}"
        )
    print(f"[OK] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
