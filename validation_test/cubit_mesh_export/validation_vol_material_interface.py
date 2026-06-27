"""Validation-class Netgen .vol material/interface inventory example.

This example keeps the mesh small enough to read by eye: two tetrahedra use
different material numbers and share one triangular material interface.  It
validates the helper tables that turn Coreform/Cubit blocks and sidesets into
readable FEM/BEM setup data:

* volume material rows: material id, name, tetra count, volume, touching
  boundaries, exterior area, and interface area;
* domain-boundary incidence rows: boundary name plus ``domin/domout`` material
  adjacency, so exterior and interface boundaries are not confused.

Run:

    python validation_test/cubit_mesh_export/validation_vol_material_interface.py
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

from radia_mcp.radia_ngsolve.netgen_vol import (  # noqa: E402
    parse_netgen_tri_tet_vol,
    read_netgen_tri_tet_vol,
)


OUT_JSON = Path(__file__).with_name("validation_vol_material_interface_summary.json")

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


def _built_in_summary() -> dict:
    mesh = parse_netgen_tri_tet_vol(TWO_MATERIAL_INTERFACE_VOL)
    material_rows = list(mesh.material_summary_rows())
    incidence_rows = list(mesh.domain_boundary_incidence_rows())
    by_material = {row["name"]: row for row in material_rows}
    by_boundary = {row["name"]: row for row in incidence_rows}
    outer_area = 1.0 + 0.5 * math.sqrt(3.0)

    checks = {
        "materials": [row["name"] for row in material_rows],
        "boundary_kinds": {row["name"]: row["kind"] for row in incidence_rows},
        "total_volume": mesh.total_volume(),
        "air_volume": by_material["air"]["volume"],
        "core_volume": by_material["core"]["volume"],
        "air_volume_fraction": by_material["air"]["volume_fraction"],
        "core_volume_fraction": by_material["core"]["volume_fraction"],
        "expected_one_tet_volume": 1.0 / 6.0,
        "expected_outer_area": outer_area,
        "air_exterior_area": by_material["air"]["exterior_surface_area"],
        "core_exterior_area": by_material["core"]["exterior_surface_area"],
        "air_interface_area": by_material["air"]["interface_surface_area"],
        "core_interface_area": by_material["core"]["interface_surface_area"],
        "interface_area": by_boundary["air_core_interface"]["surface_area"],
        "incidence_area_sum": sum(row["surface_area"] for row in incidence_rows),
        "total_surface_area": mesh.total_surface_area(),
        "interface_domin_domout": [
            by_boundary["air_core_interface"]["domin"],
            by_boundary["air_core_interface"]["domout"],
        ],
        "air_neighbor_materials": by_material["air"]["neighboring_material_numbers"],
        "core_neighbor_materials": by_material["core"]["neighboring_material_numbers"],
    }

    assert checks["materials"] == ["air", "core"]
    assert checks["boundary_kinds"] == {
        "air_outer": "exterior",
        "core_outer": "exterior",
        "air_core_interface": "interface",
    }
    assert abs(checks["total_volume"] - 1.0 / 3.0) < 1.0e-15
    assert abs(checks["air_volume"] - checks["expected_one_tet_volume"]) < 1.0e-15
    assert abs(checks["core_volume"] - checks["expected_one_tet_volume"]) < 1.0e-15
    assert abs(checks["air_volume_fraction"] - 0.5) < 1.0e-15
    assert abs(checks["core_volume_fraction"] - 0.5) < 1.0e-15
    assert abs(checks["air_exterior_area"] - outer_area) < 1.0e-15
    assert abs(checks["core_exterior_area"] - outer_area) < 1.0e-15
    assert abs(checks["air_interface_area"] - 0.5) < 1.0e-15
    assert abs(checks["core_interface_area"] - 0.5) < 1.0e-15
    assert abs(checks["interface_area"] - 0.5) < 1.0e-15
    assert abs(checks["incidence_area_sum"] - checks["total_surface_area"]) < 1.0e-15
    assert checks["interface_domin_domout"] == [1, 2]
    assert checks["air_neighbor_materials"] == [2]
    assert checks["core_neighbor_materials"] == [1]

    return {
        "kind": "netgen_vol_material_interface_inventory_validation",
        "validation_class": True,
        "built_in": {
            "checks": checks,
            "material_rows": material_rows,
            "domain_boundary_incidence_rows": incidence_rows,
        },
    }


def _external_summary(path: Path) -> dict:
    mesh = read_netgen_tri_tet_vol(path)
    material_rows = list(mesh.material_summary_rows())
    incidence_rows = list(mesh.domain_boundary_incidence_rows())
    total_material_volume = sum(row["volume"] for row in material_rows)
    incidence_area_sum = sum(row["surface_area"] for row in incidence_rows)
    checks = {
        "path": str(path),
        "summary": mesh.summary(),
        "material_count": len(material_rows),
        "incidence_row_count": len(incidence_rows),
        "total_volume": mesh.total_volume(),
        "total_material_volume": total_material_volume,
        "total_surface_area": mesh.total_surface_area(),
        "incidence_area_sum": incidence_area_sum,
        "max_material_volume_error": abs(total_material_volume - mesh.total_volume()),
        "max_incidence_area_error": abs(incidence_area_sum - mesh.total_surface_area()),
        "interface_rows": sum(1 for row in incidence_rows if row["kind"] == "interface"),
        "exterior_rows": sum(1 for row in incidence_rows if row["kind"] == "exterior"),
    }
    assert checks["summary"]["tetrahedra"] > 0
    assert checks["material_count"] > 0
    assert checks["incidence_row_count"] > 0
    assert checks["max_material_volume_error"] < 1.0e-9 * max(1.0, mesh.total_volume())
    assert checks["max_incidence_area_error"] < 1.0e-9 * max(1.0, mesh.total_surface_area())
    return {
        "checks": checks,
        "material_rows": material_rows,
        "domain_boundary_incidence_rows": incidence_rows,
    }


def build_summary(external_vol: Path | None = None) -> dict:
    summary = _built_in_summary()
    if external_vol is not None:
        summary["external"] = _external_summary(external_vol)
    return summary


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--external-vol", type=Path)
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    summary = build_summary(args.external_vol)
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    checks = summary["built_in"]["checks"]
    print("[Netgen .vol material/interface inventory]")
    print(
        f"  materials={checks['materials']}, total volume={checks['total_volume']:.15f}, "
        f"interface area={checks['interface_area']:.15f}"
    )
    print(
        f"  exterior areas air/core={checks['air_exterior_area']:.15f}/"
        f"{checks['core_exterior_area']:.15f}, incidence area error="
        f"{abs(checks['incidence_area_sum'] - checks['total_surface_area']):.3e}"
    )
    if "external" in summary:
        external = summary["external"]["checks"]
        print(
            f"  external tet/material/incidence={external['summary']['tetrahedra']}/"
            f"{external['material_count']}/{external['incidence_row_count']}, "
            f"interfaces={external['interface_rows']}"
        )
    print(f"[OK] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
