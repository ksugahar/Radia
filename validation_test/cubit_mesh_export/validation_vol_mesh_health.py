r"""Validation-class Netgen .vol mesh health example.

This is an example/validation run, not a pytest test.  It combines the small
tri/tet-only checks that matter before a Netgen ``.vol`` export is used as a
first-order FEM/BEM trace:

* boundary triangles and volume tetrahedra are both present
* the boundary surface is a closed edge manifold
* boundary triangles match tetrahedron faces
* surface and volume element shape quality stay above a chosen gate
* the worst elements are named explicitly for debugging

Run:

    python validation_test/cubit_mesh_export/validation_vol_mesh_health.py
    python validation_test/cubit_mesh_export/validation_vol_mesh_health.py --vol C:\temp\mesh.vol --out C:\temp\mesh_health.json
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

from radia_mcp.radia_ngsolve.netgen_vol import parse_netgen_tri_tet_vol, read_netgen_tri_tet_vol  # noqa: E402


OUT_JSON = HERE / "validation_vol_mesh_health_summary.json"

HEALTH_GATE = {
    "min_surface_triangle_quality": 1.0e-6,
    "min_tetrahedron_quality": 1.0e-2,
    "closure_relative_tolerance": 1.0e-9,
    "worst_limit": 3,
}

RIGHT_TET_VOL = """\
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

SLIVER_TET_VOL = RIGHT_TET_VOL.replace("0 0 1", "0.001 0.001 0.001")

OPEN_PATCH_VOL = """\
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
0 1 0
pointelements
0
materials
0
bcnames
1
1 patch
endmesh
"""


def _record(label: str, vol_text: str) -> dict[str, object]:
    mesh = parse_netgen_tri_tet_vol(vol_text, source=f"{label}.vol")
    return {
        "label": label,
        "mesh_summary": mesh.summary(),
        "health": mesh.mesh_health_summary(**HEALTH_GATE),
    }


def _external_record(path: Path) -> dict[str, object]:
    mesh = read_netgen_tri_tet_vol(path)
    return {
        "source_name": path.name,
        "mesh_summary": mesh.summary(),
        "health": mesh.mesh_health_summary(**HEALTH_GATE),
    }


def build_summary(external_vol: Path | None = None) -> dict[str, object]:
    cases = [
        _record("clean_right_tetrahedron", RIGHT_TET_VOL),
        _record("thin_sliver_quality_gate", SLIVER_TET_VOL),
        _record("open_surface_topology_gate", OPEN_PATCH_VOL),
    ]
    by_label = {case["label"]: case for case in cases}
    clean = by_label["clean_right_tetrahedron"]["health"]
    sliver = by_label["thin_sliver_quality_gate"]["health"]
    open_patch = by_label["open_surface_topology_gate"]["health"]

    checks = {
        "clean_status": clean["status"],
        "clean_ok": clean["ok_for_first_order_fem_bem"],
        "sliver_status": sliver["status"],
        "sliver_tet_quality_ok": sliver["checks"]["tetrahedron_quality_above_threshold"],
        "open_status": open_patch["status"],
        "open_closed_manifold_ok": open_patch["checks"]["surface_is_closed_manifold"],
        "open_boundary_face_match_ok": open_patch["checks"]["boundary_faces_match_tetrahedra"],
    }

    assert checks["clean_ok"] is True
    assert checks["sliver_tet_quality_ok"] is False
    assert checks["open_closed_manifold_ok"] is False
    assert checks["open_boundary_face_match_ok"] is False

    summary: dict[str, object] = {
        "kind": "netgen_vol_mesh_health_validation",
        "validation_class": True,
        "gate": HEALTH_GATE,
        "lesson": (
            "A first-order FEM/BEM .vol mesh needs shape quality, closed boundary "
            "topology, and boundary-to-volume face consistency at the same time."
        ),
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

    print("[mesh health]")
    for case in summary["builtin_cases"]:
        health = case["health"]
        print(
            f"  {case['label']}: status={health['status']}, "
            f"ok={health['ok_for_first_order_fem_bem']}, "
            f"issues={len(health['issues'])}"
        )
    if "external_vol" in summary:
        health = summary["external_vol"]["health"]
        print(
            f"  external {summary['external_vol']['source_name']}: "
            f"status={health['status']}, ok={health['ok_for_first_order_fem_bem']}"
        )
    print(f"[OK] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
