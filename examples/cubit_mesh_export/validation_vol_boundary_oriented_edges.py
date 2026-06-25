"""Validation-class Netgen .vol boundary oriented-edge table.

Run:

    python examples/cubit_mesh_export/validation_vol_boundary_oriented_edges.py

This example keeps the first-order boundary trace small and readable: one
tetrahedron, four boundary triangles, six closed surface edges, and twelve
oriented triangle-edge rows.  The table is the direct pre-assembly object for a
first-order RWG-style boundary trace and the matching HCurl edge ids.
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

from radia_mcp.radia_ngsolve.netgen_vol import parse_netgen_tri_tet_vol  # noqa: E402


OUT_JSON = HERE / "validation_vol_boundary_oriented_edges_summary.json"

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


def build_summary() -> dict[str, object]:
    mesh = parse_netgen_tri_tet_vol(UNIT_TET_VOL)
    edge_summary = mesh.boundary_oriented_edge_summary()
    rows = edge_summary["rows"]

    checks = {
        "surface_triangles": edge_summary["surface_triangles"],
        "surface_edges": edge_summary["surface_edges"],
        "rwg_dof_edges": edge_summary["rwg_dof_edges"],
        "oriented_edge_rows": edge_summary["oriented_edge_rows"],
        "open_edges": edge_summary["open_edges"],
        "is_closed_manifold": edge_summary["is_closed_manifold"],
        "orientation_sign_counts": edge_summary["orientation_sign_counts"],
        "first_triangle_oriented_edges": [
            row["oriented_edge_nodes_global"] for row in rows[:3]
        ],
        "first_triangle_signs": [row["orientation_sign"] for row in rows[:3]],
        "hcurl_edge_ids": sorted({row["hcurl_edge_id"] for row in rows}),
    }

    assert checks["surface_triangles"] == 4
    assert checks["surface_edges"] == 6
    assert checks["rwg_dof_edges"] == 6
    assert checks["oriented_edge_rows"] == 12
    assert checks["open_edges"] == 0
    assert checks["is_closed_manifold"] is True
    assert checks["orientation_sign_counts"] == {"-1": 6, "1": 6}
    assert checks["first_triangle_oriented_edges"] == [[1, 2], [2, 3], [3, 1]]
    assert checks["first_triangle_signs"] == [1, 1, -1]
    assert checks["hcurl_edge_ids"] == [1, 2, 3, 4, 5, 6]

    return {
        "kind": "vol_boundary_oriented_edges_validation",
        "validation_class": True,
        "learning_theme": (
            "boundary triangle orientation should be expanded into readable "
            "RWG-style oriented edge rows before BEM assembly"
        ),
        "checks": checks,
        "edge_summary": edge_summary,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    summary = build_summary()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    checks = summary["checks"]
    print("[vol boundary oriented edges]")
    print(
        f"  triangles={checks['surface_triangles']} "
        f"surface_edges={checks['surface_edges']} rwg_dofs={checks['rwg_dof_edges']}"
    )
    print(
        f"  oriented_rows={checks['oriented_edge_rows']} "
        f"signs={checks['orientation_sign_counts']}"
    )
    print(f"  first_triangle={checks['first_triangle_oriented_edges']}")
    print(f"[OK] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
