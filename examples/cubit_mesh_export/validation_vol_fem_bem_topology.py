"""Validation-class Netgen .vol FEM/BEM topology checks.

This example is intentionally outside the fast unit-test path. It records the
topology gates needed by a first-order FEM/BEM coupling view:

* boundary triangles should form a closed 2-manifold for RWG basis functions;
* the boundary Euler characteristic catches missing or duplicated triangles;
* compact scalar-BEM nodes keep their one-based mapping back to H1 nodes;
* RWG surface edges should map onto HCurl volume edges.

Optional external ``--vol`` input can analyze a freshly exported Coreform/Cubit
tri/tet mesh.
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

from radia_mcp.radia_ngsolve.netgen_vol import parse_netgen_tri_tet_vol, read_netgen_tri_tet_vol


OUT_JSON = Path(__file__).with_name("validation_vol_fem_bem_topology_summary.json")

POINTS = [
    (0.0, 0.0, 0.0),
    (1.0, 0.0, 0.0),
    (1.0, 1.0, 0.0),
    (0.0, 1.0, 0.0),
    (0.0, 0.0, 1.0),
    (1.0, 0.0, 1.0),
    (1.0, 1.0, 1.0),
    (0.0, 1.0, 1.0),
    (0.5, 0.5, 0.5),
]

OUTWARD_SURFACE_TRIANGLES = [
    (1, 3, 2),
    (1, 4, 3),
    (5, 6, 7),
    (5, 7, 8),
    (1, 2, 6),
    (1, 6, 5),
    (4, 7, 3),
    (4, 8, 7),
    (1, 5, 8),
    (1, 8, 4),
    (2, 3, 7),
    (2, 7, 6),
]


def _star_cube_vol() -> str:
    lines = [
        "mesh3d",
        "dimension",
        "3",
        "geomtype",
        "0",
        "facedescriptors",
        "1",
        "1 1 0 1 1",
        "surfaceelements",
        str(len(OUTWARD_SURFACE_TRIANGLES)),
    ]
    for tri in OUTWARD_SURFACE_TRIANGLES:
        lines.append(f"1 1 1 0 3 {tri[0]} {tri[1]} {tri[2]}")
    lines.extend(["volumeelements", str(len(OUTWARD_SURFACE_TRIANGLES))])
    for tri in OUTWARD_SURFACE_TRIANGLES:
        lines.append(f"1 4 {tri[0]} {tri[1]} {tri[2]} 9")
    lines.extend(["points", str(len(POINTS))])
    lines.extend(f"{x} {y} {z}" for x, y, z in POINTS)
    lines.extend(["pointelements", "0", "materials", "1", "1 air", "bcnames", "1", "1 outer", "endmesh"])
    return "\n".join(lines) + "\n"


def _open_triangle_vol() -> str:
    return """\
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


def _topology_counts(topology: dict) -> dict:
    return {
        "h1_nodes": len(topology["h1"]["node_ids"]),
        "h1_trace_nodes": len(topology["h1"]["trace_node_ids"]),
        "hcurl_edges": len(topology["hcurl"]["edges"]),
        "scalar_bem_nodes": len(topology["scalar_bem"]["node_ids"]),
        "scalar_bem_triangles": len(topology["scalar_bem"]["triangles"]),
        "rwg_edges": len(topology["rwg"]["edges_local"]),
        "rwg_dof_edges": len(topology["rwg"]["dof_edge_ids"]),
        "rwg_to_hcurl_edge_ids": len(topology["trace"]["rwg_to_hcurl_edge_ids"]),
    }


def _record(name: str, mesh, include_first_order: bool) -> dict:
    record = {
        "name": name,
        "mesh_summary": mesh.summary(),
        "bounding_box": mesh.bounding_box(),
        "manifold": mesh.surface_edge_manifold_summary(),
    }
    if mesh.tetrahedra:
        record["closure"] = mesh.surface_closure_summary()
    if include_first_order:
        topology = mesh.first_order_fem_bem_topology()
        record["topology_counts"] = _topology_counts(topology)
        record["trace_policy"] = topology["policy"]
    return record


def build_summary(vol_path: str | None = None) -> dict:
    records = [
        _record("builtin_star_cube_closed", parse_netgen_tri_tet_vol(_star_cube_vol()), True),
        _record("builtin_open_triangle_patch", parse_netgen_tri_tet_vol(_open_triangle_vol()), False),
    ]
    if vol_path is not None:
        records.append(_record("external_vol", read_netgen_tri_tet_vol(vol_path), True))

    closed = records[0]
    open_patch = records[1]
    checks = {
        "builtin_closed_euler_characteristic": closed["manifold"]["euler_characteristic"],
        "builtin_closed_surface_edges": closed["manifold"]["surface_edges"],
        "builtin_closed_rwg_dof_edges": closed["topology_counts"]["rwg_dof_edges"],
        "builtin_closed_hcurl_edges": closed["topology_counts"]["hcurl_edges"],
        "builtin_open_patch_open_edges": open_patch["manifold"]["open_edges"],
        "builtin_open_patch_is_closed": open_patch["manifold"]["is_closed_manifold"],
    }

    if not closed["manifold"]["is_closed_manifold"]:
        raise AssertionError("builtin closed surface is not a 2-manifold")
    if closed["manifold"]["euler_characteristic"] != 2:
        raise AssertionError("builtin closed surface Euler characteristic drifted")
    if closed["topology_counts"]["rwg_dof_edges"] != closed["manifold"]["surface_edges"]:
        raise AssertionError("closed surface RWG dof count does not match surface edge count")
    if closed["topology_counts"]["rwg_to_hcurl_edge_ids"] != closed["topology_counts"]["rwg_dof_edges"]:
        raise AssertionError("RWG-to-HCurl trace map is incomplete")
    if open_patch["manifold"]["is_closed_manifold"] or open_patch["manifold"]["open_edges"] != 3:
        raise AssertionError("open patch was not detected")

    if vol_path is not None:
        external = records[-1]
        checks.update({
            "external_is_closed_manifold": external["manifold"]["is_closed_manifold"],
            "external_euler_characteristic": external["manifold"]["euler_characteristic"],
            "external_surface_edges": external["manifold"]["surface_edges"],
            "external_rwg_dof_edges": external["topology_counts"]["rwg_dof_edges"],
            "external_hcurl_edges": external["topology_counts"]["hcurl_edges"],
            "external_surface_abs_volume_rel_error": (
                external["closure"]["surface_abs_volume_rel_error"]
            ),
        })
        if not external["manifold"]["is_closed_manifold"]:
            raise AssertionError("external .vol boundary is not a closed 2-manifold")
        if external["manifold"]["euler_characteristic"] != 2:
            raise AssertionError("external .vol boundary is not sphere-like")
        if external["topology_counts"]["rwg_dof_edges"] != external["manifold"]["surface_edges"]:
            raise AssertionError("external .vol RWG dof count does not match surface edges")
        if external["closure"]["surface_abs_volume_rel_error"] > 1.0e-12:
            raise AssertionError("external .vol boundary volume does not match tet volume")

    return {
        "problem": "Netgen .vol first-order FEM/BEM topology validation",
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
