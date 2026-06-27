"""Validation-class `.vol` boundary normal/vector-area rows.

Electromagnetic force extraction needs surface orientation before it needs a
large solver: Maxwell stress integrates ``T n dS``.  This example checks that
named Coreform/Cubit sidesets exported through Netgen `.vol` preserve enough
triangle orientation information to recover each boundary's vector area and
planar unit normal.

Run:

    python validation_test/cubit_mesh_export/validation_vol_boundary_normal_vectors.py
    python validation_test/cubit_mesh_export/validation_vol_boundary_normal_vectors.py --vol C:\\temp\\box.vol
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "packages" / "radia-mcp" / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from radia_mcp.radia_ngsolve.netgen_vol import parse_netgen_tri_tet_vol, read_netgen_tri_tet_vol  # noqa: E402


OUT_JSON = Path(__file__).with_name("validation_vol_boundary_normal_vectors_summary.json")

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

EXPECTED_VECTOR_AREAS = {
    "xmin": (-15.0, 0.0, 0.0),
    "xmax": (15.0, 0.0, 0.0),
    "ymin": (0.0, -10.0, 0.0),
    "ymax": (0.0, 10.0, 0.0),
    "zmin": (0.0, 0.0, -6.0),
    "zmax": (0.0, 0.0, 6.0),
}


def _norm(values):
    return math.sqrt(sum(value * value for value in values))


def _record(label, mesh):
    rows = list(mesh.boundary_normal_summary_rows())
    return {
        "label": label,
        "mesh_summary": mesh.summary(),
        "bounding_box": mesh.bounding_box(),
        "normal_rows": rows,
        "closure": mesh.surface_closure_summary() if mesh.tetrahedra else None,
    }


def build_summary(external_vol=None):
    records = [_record("builtin_named_box", parse_netgen_tri_tet_vol(BOX_SIX_BOUNDARY_VOL))]
    if external_vol is not None:
        records.append(_record("external_vol", read_netgen_tri_tet_vol(external_vol)))

    builtin = records[0]
    by_name = {row["name"]: row for row in builtin["normal_rows"]}
    vector_errors = {
        name: _norm(
            [
                by_name[name]["vector_area"][axis] - expected[axis]
                for axis in range(3)
            ]
        )
        for name, expected in EXPECTED_VECTOR_AREAS.items()
    }
    ratio_errors = {
        row["name"]: abs(row["vector_area_norm_over_area"] - 1.0)
        for row in builtin["normal_rows"]
    }
    checks = {
        "builtin_boundary_names": sorted(by_name),
        "builtin_vector_area_errors": vector_errors,
        "builtin_max_vector_area_error": max(vector_errors.values()),
        "builtin_max_planar_ratio_error": max(ratio_errors.values()),
        "builtin_surface_vector_area_norm": builtin["closure"]["surface_vector_area_norm"],
        "builtin_surface_volume_rel_error": builtin["closure"]["surface_abs_volume_rel_error"],
        "builtin_total_area_from_normal_rows": sum(row["surface_area"] for row in builtin["normal_rows"]),
        "builtin_expected_total_area": 62.0,
    }
    if checks["builtin_boundary_names"] != sorted(EXPECTED_VECTOR_AREAS):
        raise AssertionError("boundary names drifted")
    if checks["builtin_max_vector_area_error"] > 1.0e-14:
        raise AssertionError("builtin boundary vector areas drifted")
    if checks["builtin_max_planar_ratio_error"] > 1.0e-14:
        raise AssertionError("builtin sidesets should be planar")
    if checks["builtin_surface_vector_area_norm"] > 1.0e-14:
        raise AssertionError("closed box vector area should cancel")
    if checks["builtin_surface_volume_rel_error"] > 1.0e-14:
        raise AssertionError("surface/tet volume mismatch")
    if abs(checks["builtin_total_area_from_normal_rows"] - checks["builtin_expected_total_area"]) > 1.0e-14:
        raise AssertionError("total boundary area drifted")

    if external_vol is not None:
        external = records[-1]
        ext_rows = external["normal_rows"]
        ext_area = sum(row["surface_area"] for row in ext_rows)
        ext_vector_norm = external["closure"]["surface_vector_area_norm"] if external["closure"] else None
        checks.update({
            "external_boundary_count": len(ext_rows),
            "external_total_area_from_normal_rows": ext_area,
            "external_surface_vector_area_norm": ext_vector_norm,
            "external_boundary_names": [row["name"] for row in ext_rows],
            "external_min_vector_area_norm_over_area": min(
                row["vector_area_norm_over_area"] for row in ext_rows
                if row["vector_area_norm_over_area"] is not None
            ),
        })
        if not ext_rows:
            raise AssertionError("external .vol had no boundary normal rows")
        if external["closure"] and ext_vector_norm > 1.0e-8:
            raise AssertionError("external closed surface vector area does not cancel")

    return {
        "kind": "netgen_vol_boundary_normal_vectors_validation",
        "validation_class": True,
        "force_learning": "boundary vector areas provide the n dS data for Maxwell-stress force integration",
        "checks": checks,
        "records": records,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vol", type=Path, default=None, help="Optional external tri/tet Netgen .vol file")
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    summary = build_summary(args.vol)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print("[boundary normal vectors]")
    for row in summary["records"][0]["normal_rows"]:
        print(
            f"  {row['name']}: area={row['surface_area']:.6g} "
            f"vector={row['vector_area']} normal={row['unit_normal']}"
        )
    if "external_boundary_count" in summary["checks"]:
        print(
            f"  external boundaries={summary['checks']['external_boundary_count']} "
            f"area={summary['checks']['external_total_area_from_normal_rows']:.12g}"
        )
    print(f"[OK] wrote {args.out}")


if __name__ == "__main__":
    main()
