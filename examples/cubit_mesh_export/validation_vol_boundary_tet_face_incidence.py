"""Validation-class `.vol` boundary-to-tetrahedron face incidence.

Run:

    python examples/cubit_mesh_export/validation_vol_boundary_tet_face_incidence.py

The check verifies that boundary triangles exported in a tri/tet Netgen `.vol`
mesh actually correspond to volume tetrahedron faces.  Exterior triangles
should match one tetrahedron face; material-interface triangles should match
two.  This catches broken sidesets before FEM/BEM trace coupling or boundary
force integration uses the mesh.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
SRC = REPO / "packages" / "radia-mcp" / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from radia_mcp.radia_ngsolve.netgen_vol import parse_netgen_tri_tet_vol  # noqa: E402


OUT_JSON = HERE / "validation_vol_boundary_tet_face_incidence_summary.json"

TWO_MATERIAL_INTERFACE_VOL = """\
mesh3d
dimension
3
geomtype
0
facedescriptors
3
1 1 0 1 1
2 1 0 2 1
3 1 2 2 1
surfaceelements
7
1 1 1 0 3 1 4 2
1 1 1 0 3 2 4 3
1 1 1 0 3 3 4 1
2 2 2 0 3 1 2 5
2 2 2 0 3 2 3 5
2 2 2 0 3 3 1 5
3 3 1 2 3 1 2 3
volumeelements
2
1 4 1 2 3 4
2 4 1 3 2 5
points
5
0 0 0
1 0 0
0 1 0
0 0 1
0 0 -1
pointelements
0
materials
2
1 air
2 core
bcnames
3
1 air_outer
2 core_outer
3 air_core_interface
endmesh
"""


def build_summary() -> dict:
    mesh = parse_netgen_tri_tet_vol(
        TWO_MATERIAL_INTERFACE_VOL,
        source="embedded_two_material_interface.vol",
    )
    summary = mesh.boundary_tet_face_incidence_summary()
    interface_rows = [
        row for row in summary["rows"]
        if row["kind"] == "interface"
    ]
    exterior_rows = [
        row for row in summary["rows"]
        if row["kind"] == "exterior"
    ]

    checks = {
        "surface_triangles": summary["surface_triangles"],
        "tetrahedra": summary["tetrahedra"],
        "exterior_surface_triangles": summary["exterior_surface_triangles"],
        "interface_surface_triangles": summary["interface_surface_triangles"],
        "orphan_surface_triangles": summary["orphan_surface_triangles"],
        "overconnected_surface_triangles": summary["overconnected_surface_triangles"],
        "domain_material_mismatch_count": summary["domain_material_mismatch_count"],
        "max_adjacent_tetrahedra": summary["max_adjacent_tetrahedra"],
        "is_volume_boundary_consistent": summary["is_volume_boundary_consistent"],
        "interface_row_count": len(interface_rows),
        "exterior_row_count": len(exterior_rows),
        "interface_adjacent_material_numbers": (
            interface_rows[0]["adjacent_material_numbers"] if interface_rows else []
        ),
        "interface_declared_domain_numbers": (
            interface_rows[0]["declared_domain_numbers"] if interface_rows else []
        ),
    }

    assert checks["surface_triangles"] == 7
    assert checks["tetrahedra"] == 2
    assert checks["exterior_surface_triangles"] == 6
    assert checks["interface_surface_triangles"] == 1
    assert checks["orphan_surface_triangles"] == 0
    assert checks["overconnected_surface_triangles"] == 0
    assert checks["domain_material_mismatch_count"] == 0
    assert checks["max_adjacent_tetrahedra"] == 2
    assert checks["is_volume_boundary_consistent"]
    assert checks["interface_adjacent_material_numbers"] == [1, 2]
    assert checks["interface_declared_domain_numbers"] == [1, 2]

    return {
        "kind": "vol_boundary_tet_face_incidence_validation",
        "validation_class": True,
        "mesh_summary": mesh.summary(),
        "checks": checks,
        "incidence_summary": summary,
    }


def main() -> int:
    summary = build_summary()
    OUT_JSON.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    checks = summary["checks"]
    print("[vol boundary tet-face incidence]")
    print(
        f"  surface={checks['surface_triangles']} tet={checks['tetrahedra']} "
        f"exterior={checks['exterior_surface_triangles']} interface={checks['interface_surface_triangles']}"
    )
    print(
        f"  orphan={checks['orphan_surface_triangles']} "
        f"overconnected={checks['overconnected_surface_triangles']} "
        f"mismatch={checks['domain_material_mismatch_count']}"
    )
    print(
        f"  interface adjacent materials={checks['interface_adjacent_material_numbers']} "
        f"declared={checks['interface_declared_domain_numbers']}"
    )
    print(f"  wrote {OUT_JSON}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
