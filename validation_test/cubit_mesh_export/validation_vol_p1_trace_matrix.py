"""Validation-class P1 FEM/BEM trace matrix from a tri/tet .vol mesh.

Run:

    python validation_test/cubit_mesh_export/validation_vol_p1_trace_matrix.py

For first-order H1 FEM and first-order scalar BEM on the same tri/tet mesh, the
boundary trace is a boolean gather matrix.  In one-based sparse COO notation:

    T = sparse(rows, cols, values, nTraceNodes, nVolumeNodes)

where each compact boundary-node row points to the original one-based volume
node id in ``cols``.
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


OUT_JSON = HERE / "validation_vol_p1_trace_matrix_summary.json"

FOUR_TET_WITH_INTERIOR_NODE_VOL = """\
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


def build_summary() -> dict[str, object]:
    mesh = parse_netgen_tri_tet_vol(
        FOUR_TET_WITH_INTERIOR_NODE_VOL,
        source="embedded_four_tet_with_interior_node.vol",
    )
    trace = mesh.p1_fem_bem_trace_matrix_summary()
    topology = mesh.first_order_fem_bem_topology()

    checks = {
        "n_volume_nodes": trace["n_volume_nodes"],
        "n_trace_nodes": trace["n_trace_nodes"],
        "matrix_shape": trace["matrix_shape"],
        "nnz": trace["nnz"],
        "rows": trace["rows"],
        "cols": trace["cols"],
        "interior_node_ids": trace["interior_node_ids"],
        "surface_triangles_local": trace["surface_triangles_local"],
        "topology_trace_cols": topology["trace"]["h1_to_scalar_bem_cols"],
        "topology_trace_rows": topology["trace"]["h1_to_scalar_bem_rows"],
    }

    assert checks["n_volume_nodes"] == 5
    assert checks["n_trace_nodes"] == 4
    assert checks["matrix_shape"] == [4, 5]
    assert checks["nnz"] == 4
    assert checks["rows"] == [1, 2, 3, 4]
    assert checks["cols"] == [1, 2, 3, 4]
    assert checks["interior_node_ids"] == [5]
    assert checks["topology_trace_rows"] == checks["rows"]
    assert checks["topology_trace_cols"] == checks["cols"]
    assert all(5 not in tri for tri in checks["surface_triangles_local"])

    return {
        "kind": "vol_p1_fem_bem_trace_matrix_validation",
        "validation_class": True,
        "learning_theme": (
            "for first-order shared tri/tet meshes, the scalar FEM-to-BEM trace "
            "is a one-based boolean gather matrix"
        ),
        "mesh_summary": mesh.summary(),
        "checks": checks,
        "trace_matrix": trace,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    summary = build_summary()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    checks = summary["checks"]
    print("[vol P1 FEM/BEM trace matrix]")
    print(
        f"  shape={checks['matrix_shape']} nnz={checks['nnz']} "
        f"volume_nodes={checks['n_volume_nodes']} trace_nodes={checks['n_trace_nodes']}"
    )
    print(f"  rows={checks['rows']}")
    print(f"  cols={checks['cols']}")
    print(f"  interior_node_ids={checks['interior_node_ids']}")
    print(f"[OK] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
