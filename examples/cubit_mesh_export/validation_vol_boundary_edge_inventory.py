"""Validation-class Netgen .vol boundary-local edge inventory.

Run:

    python examples/cubit_mesh_export/validation_vol_boundary_edge_inventory.py

For a CAD sideset exported as triangles, boundary-local perimeter edges and
within-boundary split/diagonal edges should be easy to distinguish.  This
example uses a 2 x 3 x 5 box with each rectangular face split into two
triangles: every named face has four perimeter edges and one internal diagonal.
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


OUT_JSON = HERE / "validation_vol_boundary_edge_inventory_summary.json"

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


def build_summary() -> dict[str, object]:
    mesh = parse_netgen_tri_tet_vol(BOX_SIX_BOUNDARY_VOL)
    inventory = mesh.boundary_edge_inventory_summary()
    rows = {row["name"]: row for row in inventory["rows"]}
    zmax = rows["zmax"]

    checks = {
        "boundary_count": inventory["boundary_count"],
        "surface_triangles": inventory["surface_triangles"],
        "unique_boundary_edges_total": inventory["unique_boundary_edges_total"],
        "perimeter_edges_total": inventory["perimeter_edges_total"],
        "shared_diagonal_edges_total": inventory["shared_diagonal_edges_total"],
        "overused_edges_total": inventory["overused_edges_total"],
        "zmax_perimeter_edges": zmax["perimeter_edges"],
        "zmax_shared_diagonal_edges": zmax["shared_diagonal_edges"],
        "zmax_perimeter_length": zmax["perimeter_edge_length_sum_m"],
        "zmax_diagonal_length": zmax["shared_diagonal_edge_length_sum_m"],
        "zmax_diagonal_nodes": zmax["shared_diagonal_edge_nodes"],
    }

    assert checks["boundary_count"] == 6
    assert checks["surface_triangles"] == 12
    assert checks["unique_boundary_edges_total"] == 30
    assert checks["perimeter_edges_total"] == 24
    assert checks["shared_diagonal_edges_total"] == 6
    assert checks["overused_edges_total"] == 0
    assert checks["zmax_perimeter_edges"] == 4
    assert checks["zmax_shared_diagonal_edges"] == 1
    assert abs(checks["zmax_perimeter_length"] - 10.0) < 1.0e-12
    assert abs(checks["zmax_diagonal_length"] - 13.0**0.5) < 1.0e-12
    assert checks["zmax_diagonal_nodes"] == [[5, 7]]

    return {
        "kind": "vol_boundary_edge_inventory_validation",
        "validation_class": True,
        "learning_theme": (
            "boundary-local edge inventory separates physical perimeter edges "
            "from triangulation diagonals in named .vol sidesets"
        ),
        "checks": checks,
        "inventory": inventory,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    summary = build_summary()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    checks = summary["checks"]
    print("[vol boundary edge inventory]")
    print(
        f"  boundaries={checks['boundary_count']} triangles={checks['surface_triangles']} "
        f"unique_edges={checks['unique_boundary_edges_total']}"
    )
    print(
        f"  perimeter_edges={checks['perimeter_edges_total']} "
        f"shared_diagonals={checks['shared_diagonal_edges_total']}"
    )
    print(
        f"  zmax perimeter={checks['zmax_perimeter_length']:.12g} "
        f"diagonal={checks['zmax_diagonal_length']:.12g}"
    )
    print(f"[OK] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
