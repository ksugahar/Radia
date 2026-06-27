"""Validation-class Netgen .vol named-boundary inventory checks.

This example is intentionally outside the fast unit-test path.  It records the
boundary-condition bookkeeping that matters when a Coreform/Cubit sideset export
is used as a shared FEM/BEM mesh:

* boundary numbers and names are preserved;
* per-boundary triangle counts and areas are recoverable;
* trace node ids can be grouped by boundary for readable condition maps;
* optional external ``--vol`` input can analyze a freshly exported tri/tet mesh.

Run:

    python validation_test/cubit_mesh_export/validation_vol_boundary_inventory.py
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "packages" / "radia-mcp" / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from radia_mcp.radia_ngsolve.netgen_vol import parse_netgen_tri_tet_vol, read_netgen_tri_tet_vol  # noqa: E402


OUT_JSON = Path(__file__).with_name("validation_vol_boundary_inventory_summary.json")

BOX_SIX_BOUNDARY_VOL = """\
mesh3d
dimension
3
geomtype
0
facedescriptors
6
1 1 0 1 1
2 1 0 1 1
3 1 0 1 1
4 1 0 1 1
5 1 0 1 1
6 1 0 1 1
surfaceelements
12
5 5 1 0 3 1 3 2
5 5 1 0 3 1 4 3
6 6 1 0 3 5 6 7
6 6 1 0 3 5 7 8
3 3 1 0 3 1 2 6
3 3 1 0 3 1 6 5
4 4 1 0 3 4 7 3
4 4 1 0 3 4 8 7
1 1 1 0 3 1 5 8
1 1 1 0 3 1 8 4
2 2 1 0 3 2 3 7
2 2 1 0 3 2 7 6
volumeelements
12
1 4 1 3 2 9
1 4 1 4 3 9
1 4 5 6 7 9
1 4 5 7 8 9
1 4 1 2 6 9
1 4 1 6 5 9
1 4 4 7 3 9
1 4 4 8 7 9
1 4 1 5 8 9
1 4 1 8 4 9
1 4 2 3 7 9
1 4 2 7 6 9
points
9
0 0 0
2 0 0
2 3 0
0 3 0
0 0 5
2 0 5
2 3 5
0 3 5
1 1.5 2.5
pointelements
0
materials
1
1 air
bcnames
6
1 xmin
2 xmax
3 ymin
4 ymax
5 zmin
6 zmax
endmesh
"""


EXPECTED_AREAS = {
    "xmin": 15.0,
    "xmax": 15.0,
    "ymin": 10.0,
    "ymax": 10.0,
    "zmin": 6.0,
    "zmax": 6.0,
}


def _record(name, mesh):
    rows = list(mesh.boundary_summary_rows())
    return {
        "name": name,
        "mesh_summary": mesh.summary(),
        "bounding_box": mesh.bounding_box(),
        "boundary_rows": rows,
        "surface_area_by_boundary_number": mesh.surface_area_by_boundary_number(),
        "trace_node_ids_by_boundary_number": {
            bcnr: list(nodes) for bcnr, nodes in mesh.trace_node_ids_by_boundary_number().items()
        },
        "closure": mesh.surface_closure_summary() if mesh.tetrahedra else None,
        "manifold": mesh.surface_edge_manifold_summary(),
    }


def build_summary(vol_path=None):
    records = [_record("builtin_named_box", parse_netgen_tri_tet_vol(BOX_SIX_BOUNDARY_VOL))]
    if vol_path is not None:
        records.append(_record("external_vol", read_netgen_tri_tet_vol(vol_path)))

    builtin = records[0]
    by_name = {row["name"]: row for row in builtin["boundary_rows"]}
    checks = {
        "builtin_boundary_names": sorted(by_name),
        "builtin_boundary_count": len(builtin["boundary_rows"]),
        "builtin_total_surface_area": sum(row["surface_area"] for row in builtin["boundary_rows"]),
        "builtin_expected_total_surface_area": sum(EXPECTED_AREAS.values()),
        "builtin_total_volume": builtin["closure"]["tetrahedron_total_volume"],
        "builtin_expected_volume": 30.0,
        "builtin_named_area_errors": {
            name: abs(by_name[name]["surface_area"] - expected)
            for name, expected in EXPECTED_AREAS.items()
        },
        "builtin_all_faces_have_two_triangles": all(
            row["surface_triangles"] == 2 for row in builtin["boundary_rows"]
        ),
        "builtin_all_faces_have_four_trace_nodes": all(
            row["trace_node_count"] == 4 for row in builtin["boundary_rows"]
        ),
        "builtin_closed_manifold": builtin["manifold"]["is_closed_manifold"],
        "builtin_surface_volume_rel_error": builtin["closure"]["surface_abs_volume_rel_error"],
    }

    if checks["builtin_boundary_names"] != sorted(EXPECTED_AREAS):
        raise AssertionError("named boundary inventory drifted")
    if abs(checks["builtin_total_surface_area"] - checks["builtin_expected_total_surface_area"]) > 1.0e-14:
        raise AssertionError("boundary area total does not match analytic box area")
    if abs(checks["builtin_total_volume"] - checks["builtin_expected_volume"]) > 1.0e-14:
        raise AssertionError("box tet volume does not match analytic volume")
    if max(checks["builtin_named_area_errors"].values()) > 1.0e-14:
        raise AssertionError("per-boundary area drifted")
    if not checks["builtin_all_faces_have_two_triangles"]:
        raise AssertionError("each box face should have two boundary triangles")
    if not checks["builtin_all_faces_have_four_trace_nodes"]:
        raise AssertionError("each box face should expose four trace nodes")
    if not checks["builtin_closed_manifold"]:
        raise AssertionError("box boundary is not closed")
    if checks["builtin_surface_volume_rel_error"] > 1.0e-14:
        raise AssertionError("surface/tet volume mismatch")

    if vol_path is not None:
        external = records[-1]
        external_area_from_rows = sum(row["surface_area"] for row in external["boundary_rows"])
        checks.update({
            "external_boundary_count": len(external["boundary_rows"]),
            "external_area_from_rows": external_area_from_rows,
            "external_total_surface_area": external["closure"]["total_surface_area"]
            if external["closure"] is not None
            else external_area_from_rows,
            "external_named_boundaries": [row["name"] for row in external["boundary_rows"]],
        })
        if len(external["boundary_rows"]) == 0:
            raise AssertionError("external .vol had no boundary rows")
        if abs(checks["external_area_from_rows"] - checks["external_total_surface_area"]) > 1.0e-10:
            raise AssertionError("external boundary area rows do not sum to total area")

    return {
        "problem": "Netgen .vol named-boundary inventory validation",
        "validation_class": True,
        "checks": checks,
        "records": records,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vol", help="Optional external tri/tet Netgen .vol file")
    parser.add_argument("--output", default=str(OUT_JSON))
    args = parser.parse_args()

    summary = build_summary(args.vol)
    out = Path(args.output)
    out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(summary["checks"], indent=2, sort_keys=True))
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
