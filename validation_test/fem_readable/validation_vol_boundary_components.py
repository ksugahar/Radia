"""Validation-class `.vol` boundary-component trace example.

This is a MATLAB/Gypsilab-style teaching script for FEM/BEM block setup.  It
uses a tri/tet Netgen ``.vol`` mesh containing two disconnected tetrahedral
objects and checks that the boundary triangles split into independent surface
components:

* each component has its own trace node block and boundary name;
* each closed tetra boundary has Euler characteristic 2;
* the component volumes add back to the volume tetrahedra;
* RWG edge traces remain mapped to first-order HCurl volume edges.

Run:

    python validation_test/fem_readable/validation_vol_boundary_components.py
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

from radia_mcp.radia_ngsolve.netgen_vol import parse_netgen_tri_tet_vol  # noqa: E402


OUT_JSON = HERE / "validation_vol_boundary_components_summary.json"

TWO_BODY_VOL = """\
mesh3d
dimension
3
geomtype
0
facedescriptors
2
1 1 0 1 1
2 1 0 1 1
surfaceelements
8
1 1 1 0 3 1 2 3
1 1 1 0 3 1 4 2
1 1 1 0 3 2 4 3
1 1 1 0 3 3 4 1
2 2 1 0 3 5 6 7
2 2 1 0 3 5 8 6
2 2 1 0 3 6 8 7
2 2 1 0 3 7 8 5
volumeelements
2
1 4 1 2 3 4
1 4 5 6 7 8
points
8
0 0 0
1 0 0
0 1 0
0 0 1
3 0 0
4 0 0
3 1 0
3 0 1
pointelements
0
materials
1
1 air
bcnames
2
1 left_outer
2 right_outer
endmesh
"""


def _component_record(row: dict) -> dict:
    return {
        "component": row["component"],
        "boundary_names": row["boundary_names"],
        "surface_triangles": row["surface_triangles"],
        "surface_edges": row["surface_edges"],
        "trace_node_count": row["trace_node_count"],
        "trace_node_ids": row["trace_node_ids"],
        "surface_area": row["surface_area"],
        "surface_abs_volume": row["surface_abs_volume"],
        "surface_vector_area_norm": row["surface_vector_area_norm"],
        "closed_edges": row["closed_edges"],
        "open_edges": row["open_edges"],
        "is_closed_manifold": row["is_closed_manifold"],
        "euler_characteristic": row["euler_characteristic"],
    }


def validate(mesh, components: list[dict], topology: dict) -> dict:
    trace_sets = [set(row["trace_node_ids"]) for row in components]
    component_volume_sum = sum(row["surface_abs_volume"] for row in components)
    rwg_hcurl = topology["trace"]["rwg_to_hcurl_edge_ids"]
    checks = {
        "component_count": len(components),
        "component_names": [row["boundary_names"] for row in components],
        "component_trace_nodes": [row["trace_node_ids"] for row in components],
        "trace_blocks_are_disjoint": trace_sets[0].isdisjoint(trace_sets[1]),
        "all_components_closed": all(row["is_closed_manifold"] for row in components),
        "all_euler_characteristic_two": [row["euler_characteristic"] for row in components],
        "component_abs_volume_sum": component_volume_sum,
        "tetrahedron_total_volume": mesh.total_volume(),
        "component_volume_abs_error": abs(component_volume_sum - mesh.total_volume()),
        "total_surface_area": mesh.total_surface_area(),
        "hcurl_edge_count": len(topology["hcurl"]["edges"]),
        "rwg_edge_count": len(topology["rwg"]["dof_edge_ids"]),
        "rwg_to_hcurl_edge_ids": rwg_hcurl,
        "rwg_to_hcurl_is_one_to_one": len(set(rwg_hcurl)) == len(rwg_hcurl),
    }
    assert checks["component_count"] == 2
    assert checks["component_names"] == [["left_outer"], ["right_outer"]]
    assert checks["component_trace_nodes"] == [[1, 2, 3, 4], [5, 6, 7, 8]]
    assert checks["trace_blocks_are_disjoint"]
    assert checks["all_components_closed"]
    assert checks["all_euler_characteristic_two"] == [2, 2]
    assert checks["component_volume_abs_error"] < 1.0e-15
    assert checks["hcurl_edge_count"] == 12
    assert checks["rwg_edge_count"] == 12
    assert checks["rwg_to_hcurl_is_one_to_one"]
    return checks


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    mesh = parse_netgen_tri_tet_vol(TWO_BODY_VOL, source="embedded_two_body.vol")
    components = [_component_record(row) for row in mesh.surface_connected_components()]
    topology = mesh.first_order_fem_bem_topology()
    checks = validate(mesh, components, topology)
    summary = {
        "kind": "vol_boundary_component_trace_validation",
        "validation_class": True,
        "mesh_summary": mesh.summary(),
        "components": components,
        "trace": topology["trace"],
        "checks": checks,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(add_result_metadata(summary, __file__), indent=2) + "\n", encoding="utf-8")

    print("[.vol boundary components]")
    print(
        f"  components={checks['component_count']} "
        f"hcurl_edges={checks['hcurl_edge_count']} "
        f"rwg_edges={checks['rwg_edge_count']}"
    )
    for row in components:
        print(
            f"  component {row['component']}: names={row['boundary_names']} "
            f"tris={row['surface_triangles']} edges={row['surface_edges']} "
            f"nodes={row['trace_node_ids']} chi={row['euler_characteristic']} "
            f"volume={row['surface_abs_volume']:.12f}"
        )
    print("[checks]")
    for key, value in checks.items():
        print(f"  {key}: {value}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
