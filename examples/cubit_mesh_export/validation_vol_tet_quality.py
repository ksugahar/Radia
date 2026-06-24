r"""Validation-class Netgen .vol tetrahedron quality example.

This is an example/validation run, not a pytest test.  It keeps the Cubit/Coreform
mesh-quality checks public-safe and solver-independent:

* parse only tri/tet Netgen ``.vol`` topology
* compute edge ratio, inradius, circumradius, and radius-ratio quality
* validate analytic anchors on right, equilateral, and sliver tetrahedra
* optionally evaluate a real exported ``.vol`` with ``--vol PATH``

Run:

    python examples/cubit_mesh_export/validation_vol_tet_quality.py
    python examples/cubit_mesh_export/validation_vol_tet_quality.py --vol C:\temp\mesh.vol --out C:\temp\quality.json
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


OUT_JSON = HERE / "validation_vol_tet_quality_summary.json"


VOL_TEMPLATE = """\
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
1 1 1 0 3 2 4 3
1 1 1 0 3 3 4 1
volumeelements
1
1 4 1 2 3 4
points
4
{points}
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


def _vol_from_points(points: list[tuple[float, float, float]]) -> str:
    return VOL_TEMPLATE.format(points="\n".join(f"{x:.17g} {y:.17g} {z:.17g}" for x, y, z in points))


def _case(label: str, points: list[tuple[float, float, float]]) -> dict:
    mesh = parse_netgen_tri_tet_vol(_vol_from_points(points), source=f"{label}.vol")
    row = mesh.tetrahedron_quality_rows()[0]
    return {
        "label": label,
        "mesh_summary": mesh.summary(),
        "quality_row": row,
        "quality_summary": mesh.tetrahedron_quality_summary(),
        "closure": mesh.surface_closure_summary(),
    }


def _external_record(path: Path) -> dict:
    mesh = read_netgen_tri_tet_vol(path)
    return {
        "source_name": path.name,
        "mesh_summary": mesh.summary(),
        "bounding_box": mesh.bounding_box(),
        "edge_summary": mesh.tetrahedron_edge_length_summary(),
        "quality_summary": mesh.tetrahedron_quality_summary(),
        "closure": mesh.surface_closure_summary(),
        "manifold": mesh.surface_edge_manifold_summary(),
    }


def build_summary(external_vol: Path | None = None) -> dict:
    h = math.sqrt(3.0) / 2.0
    z = math.sqrt(2.0 / 3.0)
    cases = [
        _case("right_unit_tetrahedron", [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)]),
        _case("equilateral_unit_edge", [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.5, h, 0.0), (0.5, math.sqrt(3.0) / 6.0, z)]),
        _case("thin_sliver_gate", [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (1.0e-3, 1.0e-3, 1.0e-3)]),
    ]
    by_label = {case["label"]: case for case in cases}
    right = by_label["right_unit_tetrahedron"]["quality_row"]
    equi = by_label["equilateral_unit_edge"]["quality_row"]
    sliver = by_label["thin_sliver_gate"]["quality_row"]

    expected_surface_area = 1.5 + 0.5 * math.sqrt(3.0)
    expected_inradius = 0.5 / expected_surface_area
    expected_circumradius = math.sqrt(3.0) / 2.0
    expected_quality = 3.0 * expected_inradius / expected_circumradius

    checks = {
        "right_volume": right["volume"],
        "right_surface_area": right["surface_area"],
        "right_inradius": right["inradius"],
        "right_circumradius": right["circumradius"],
        "right_radius_ratio_quality": right["radius_ratio_quality"],
        "right_expected_radius_ratio_quality": expected_quality,
        "right_edge_ratio": right["edge_ratio"],
        "equilateral_radius_ratio_quality": equi["radius_ratio_quality"],
        "equilateral_edge_ratio": equi["edge_ratio"],
        "sliver_radius_ratio_quality": sliver["radius_ratio_quality"],
        "sliver_edge_ratio": sliver["edge_ratio"],
    }

    assert abs(checks["right_volume"] - 1.0 / 6.0) < 1.0e-15
    assert abs(checks["right_surface_area"] - expected_surface_area) < 1.0e-15
    assert abs(checks["right_inradius"] - expected_inradius) < 1.0e-15
    assert abs(checks["right_circumradius"] - expected_circumradius) < 1.0e-15
    assert abs(checks["right_radius_ratio_quality"] - expected_quality) < 1.0e-15
    assert abs(checks["right_edge_ratio"] - math.sqrt(2.0)) < 1.0e-15
    assert abs(checks["equilateral_radius_ratio_quality"] - 1.0) < 1.0e-14
    assert abs(checks["equilateral_edge_ratio"] - 1.0) < 1.0e-14
    assert checks["sliver_radius_ratio_quality"] < 0.01
    assert checks["sliver_edge_ratio"] > 800.0

    summary = {
        "kind": "netgen_vol_tetrahedron_quality_validation",
        "validation_class": True,
        "quality_definition": "radius_ratio_quality = 3 * inradius / circumradius",
        "builtin_cases": cases,
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

    print("[tet quality]")
    for case in summary["builtin_cases"]:
        q = case["quality_row"]
        print(
            f"  {case['label']}: quality={q['radius_ratio_quality']:.12f}, "
            f"edge_ratio={q['edge_ratio']:.6f}, volume={q['volume']:.12g}"
        )
    if "external_vol" in summary:
        q = summary["external_vol"]["quality_summary"]
        print(
            f"  external {summary['external_vol']['source_name']}: "
            f"tets={q['tetrahedra']}, min_quality={q['min_radius_ratio_quality']:.6f}, "
            f"mean_quality={q['mean_radius_ratio_quality']:.6f}, "
            f"max_edge_ratio={q['max_edge_ratio']:.6f}"
        )
    print(f"[OK] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
