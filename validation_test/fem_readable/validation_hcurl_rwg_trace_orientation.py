"""Validation-class first-order HCurl/RWG trace-orientation example.

This example is intentionally small and explicit enough to translate into a
MATLAB/Gypsilab-style teaching script.  It uses a Netgen ``.vol`` mesh with one
interior node and four tetrahedra, then checks the first-order coupling views:

* volume H1 nodes keep the original one-based ``.vol`` ids
* scalar BEM compacts only boundary nodes
* volume HCurl uses oriented tetrahedron edges
* boundary RWG uses closed-surface edges and maps into HCurl edge ids
* interior HCurl spoke edges are excluded from the RWG trace

Run:

    python validation_test/fem_readable/validation_hcurl_rwg_trace_orientation.py
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


OUT_JSON = HERE / "validation_hcurl_rwg_trace_orientation_summary.json"


FOUR_TET_STAR_VOL = """\
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
4
1 4 1 2 3 5
1 4 1 4 2 5
1 4 1 3 4 5
1 4 2 4 3 5
points
5
0 0 0
1 0 0
0 1 0
0 0 1
0.25 0.25 0.25
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


def _as_edges(rows: list[list[int]]) -> set[tuple[int, int]]:
    return {tuple(row) for row in rows}


def _edge_sign_balance(topology: dict[str, object]) -> list[dict[str, object]]:
    rwg = topology["rwg"]
    tri_edges = rwg["tri_edges"]
    tri_signs = rwg["tri_edge_signs"]
    edge_signs: dict[int, list[int]] = {edge_id: [] for edge_id in rwg["dof_edge_ids"]}
    for edges, signs in zip(tri_edges, tri_signs):
        for edge_id, sign in zip(edges, signs):
            edge_signs[edge_id].append(sign)
    return [
        {
            "edge_id": edge_id,
            "signs": signs,
            "sign_sum": sum(signs),
            "has_opposite_orientations": len(signs) == 2 and sum(signs) == 0,
        }
        for edge_id, signs in sorted(edge_signs.items())
    ]


def _trace_triplets(topology: dict[str, object]) -> dict[str, list[dict[str, float | int]]]:
    trace = topology["trace"]
    return {
        "h1_to_scalar_bem": [
            {"row": row, "col": col, "value": 1.0}
            for row, col in zip(trace["h1_to_scalar_bem_rows"], trace["h1_to_scalar_bem_cols"])
        ],
        "rwg_to_hcurl": [
            {"row": row, "col": col, "value": 1.0}
            for row, col in enumerate(trace["rwg_to_hcurl_edge_ids"], start=1)
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    mesh = parse_netgen_tri_tet_vol(FOUR_TET_STAR_VOL, source="embedded_four_tet_star.vol")
    topology = mesh.first_order_fem_bem_topology()
    manifold = mesh.surface_edge_manifold_summary()
    closure = mesh.surface_closure_summary()

    hcurl_edges = _as_edges(topology["hcurl"]["edges"])
    rwg_edges = _as_edges(topology["rwg"]["dof_edges_global"])
    interior_hcurl_edges = sorted(hcurl_edges - rwg_edges)
    edge_balance = _edge_sign_balance(topology)
    tet_edge_signs = [sign for row in topology["hcurl"]["tet_edge_signs"] for sign in row]
    trace_triplets = _trace_triplets(topology)

    expected_rwg_to_hcurl = [1, 2, 3, 5, 6, 8]
    expected_interior_spokes = [(1, 5), (2, 5), (3, 5), (4, 5)]
    checks = {
        "trace_nodes_exclude_interior_node": topology["h1"]["trace_node_ids"] == [1, 2, 3, 4],
        "scalar_bem_compacts_boundary_nodes": topology["scalar_bem"]["global_node_ids"] == [1, 2, 3, 4],
        "hcurl_edge_count": len(hcurl_edges),
        "rwg_edge_count": len(rwg_edges),
        "rwg_to_hcurl_edge_ids": topology["trace"]["rwg_to_hcurl_edge_ids"],
        "rwg_to_hcurl_edge_ids_match_expected": (
            topology["trace"]["rwg_to_hcurl_edge_ids"] == expected_rwg_to_hcurl
        ),
        "interior_hcurl_edges": [list(edge) for edge in interior_hcurl_edges],
        "interior_edges_are_spokes_to_node5": interior_hcurl_edges == expected_interior_spokes,
        "all_rwg_edges_are_hcurl_edges": rwg_edges.issubset(hcurl_edges),
        "rwg_edge_signs_are_balanced": all(row["has_opposite_orientations"] for row in edge_balance),
        "tet_edge_signs_include_reversed_local_orientations": -1 in tet_edge_signs,
        "closed_boundary_manifold": bool(manifold["is_closed_manifold"]),
        "euler_characteristic_is_two": manifold["euler_characteristic"] == 2,
        "surface_volume_matches_tets": closure["surface_abs_volume_rel_error"] < 1.0e-15,
    }

    assert checks["trace_nodes_exclude_interior_node"]
    assert checks["scalar_bem_compacts_boundary_nodes"]
    assert checks["hcurl_edge_count"] == 10
    assert checks["rwg_edge_count"] == 6
    assert checks["rwg_to_hcurl_edge_ids_match_expected"]
    assert checks["interior_edges_are_spokes_to_node5"]
    assert checks["all_rwg_edges_are_hcurl_edges"]
    assert checks["rwg_edge_signs_are_balanced"]
    assert checks["tet_edge_signs_include_reversed_local_orientations"]
    assert checks["closed_boundary_manifold"]
    assert checks["euler_characteristic_is_two"]
    assert checks["surface_volume_matches_tets"]

    summary = {
        "kind": "hcurl_rwg_trace_orientation_validation",
        "validation_class": True,
        "mesh_summary": mesh.summary(),
        "topology_policy": topology["policy"],
        "hcurl_edges": topology["hcurl"]["edges"],
        "hcurl_tet_edges": topology["hcurl"]["tet_edges"],
        "hcurl_tet_edge_signs": topology["hcurl"]["tet_edge_signs"],
        "scalar_bem": topology["scalar_bem"],
        "rwg_edges_global": topology["rwg"]["dof_edges_global"],
        "rwg_edge_triangles": topology["rwg"]["edge_triangles"],
        "rwg_tri_edges": topology["rwg"]["tri_edges"],
        "rwg_tri_edge_signs": topology["rwg"]["tri_edge_signs"],
        "rwg_opposite_vertices_local": topology["rwg"]["opposite_vertices_local"],
        "trace": topology["trace"],
        "trace_triplets": trace_triplets,
        "edge_sign_balance": edge_balance,
        "manifold": manifold,
        "closure": closure,
        "checks": checks,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(add_result_metadata(summary, __file__), indent=2), encoding="utf-8")

    print("[hcurl/rwg trace orientation]")
    print(
        f"  points={mesh.summary()['points']} boundary_nodes={len(topology['h1']['trace_node_ids'])} "
        f"tets={mesh.summary()['tetrahedra']} hcurl_edges={len(hcurl_edges)} rwg_edges={len(rwg_edges)}"
    )
    print(f"  rwg -> hcurl edge ids={topology['trace']['rwg_to_hcurl_edge_ids']}")
    print(f"  interior hcurl edges={interior_hcurl_edges}")
    print(
        f"  closed={manifold['is_closed_manifold']} euler={manifold['euler_characteristic']} "
        f"balanced_edges={checks['rwg_edge_signs_are_balanced']}"
    )
    print(f"[OK] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
