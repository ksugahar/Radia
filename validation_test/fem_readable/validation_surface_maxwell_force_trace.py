"""Validation-class `.vol` surface Maxwell traction trace.

This readable example is the FEM/BEM teaching bridge from a boundary triangle
mesh to force post-processing:

* parse the same first-order tri/tet `.vol` view used by the MATLAB prototype;
* evaluate local Maxwell traction on each oriented boundary triangle;
* distribute each constant triangle traction as P1 equivalent nodal loads;
* check that a closed surface in a uniform Maxwell stress field has zero net
  force, independent of triangle-by-triangle signs.

Run:

    python validation_test/fem_readable/validation_surface_maxwell_force_trace.py
"""

import json
from result_metadata import add_result_metadata
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "packages" / "radia-mcp" / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from radia_mcp.radia_ngsolve.force import (  # noqa: E402
    air_gap_maxwell_pressure,
    surface_triangle_maxwell_traction_summary,
)
from radia_mcp.radia_ngsolve.netgen_vol import parse_netgen_tri_tet_vol  # noqa: E402


HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "validation_surface_maxwell_force_trace_summary.json"

UNIT_TET_VOL = """\
mesh3d
dimension
3
geomtype
0
facedescriptors
4
1 1 0 1 1
2 2 0 1 1
3 3 0 1 1
4 4 0 1 1
surfaceelements
4
1 1 1 0 3 1 2 3
1 2 1 0 3 1 4 2
1 3 1 0 3 2 4 3
1 4 1 0 3 3 4 1
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
4
1 base
2 xz
3 hyp
4 yz
endmesh
"""


def _add(a, b):
    return [a[i] + b[i] for i in range(3)]


def _norm(values):
    return math.sqrt(sum(value * value for value in values))


def main():
    mesh = parse_netgen_tri_tet_vol(UNIT_TET_VOL, source="embedded_unit_tet_surface_force.vol")
    B = (0.0, 0.0, 1.0)
    face_rows = []
    nodal_loads = [[0.0, 0.0, 0.0] for _ in mesh.points]
    total_force = [0.0, 0.0, 0.0]

    for tri in mesh.surface_triangles:
        coords = [mesh.points[node - 1] for node in tri.nodes]
        row = surface_triangle_maxwell_traction_summary(coords, B)
        total_force = _add(total_force, row["integrated_force_N"])
        for local, node in enumerate(tri.nodes):
            nodal_loads[node - 1] = _add(nodal_loads[node - 1], row["nodal_force_loads_N"][local])
        face_rows.append({
            "boundary_number": tri.bcnr,
            "boundary_name": mesh.boundary_names.get(tri.bcnr, f"boundary_{tri.bcnr}"),
            "nodes": list(tri.nodes),
            "area": row["area"],
            "unit_normal": row["unit_normal"],
            "traction_Pa": row["traction_Pa"],
            "integrated_force_N": row["integrated_force_N"],
            "nodal_force_loads_N": row["nodal_force_loads_N"],
        })

    nodal_total = [sum(load[i] for load in nodal_loads) for i in range(3)]
    pressure = air_gap_maxwell_pressure(1.0)
    base = next(row for row in face_rows if row["boundary_name"] == "base")
    checks = {
        "pressure_at_1T_Pa": pressure,
        "surface_triangles": len(face_rows),
        "total_force_N": total_force,
        "total_force_norm_N": _norm(total_force),
        "nodal_total_force_N": nodal_total,
        "nodal_total_force_norm_N": _norm(nodal_total),
        "total_vs_nodal_abs_errors": [abs(total_force[i] - nodal_total[i]) for i in range(3)],
        "base_area_m2": base["area"],
        "base_force_N": base["integrated_force_N"],
        "base_force_identity_abs_error": abs(base["integrated_force_N"][2] - 0.5 * pressure),
        "surface_area_vector_norm": mesh.surface_closure_summary()["surface_vector_area_norm"],
        "surface_signed_volume": mesh.surface_closure_summary()["surface_signed_volume"],
    }

    assert checks["surface_triangles"] == 4
    assert checks["total_force_norm_N"] < 1.0e-9
    assert checks["nodal_total_force_norm_N"] < 1.0e-9
    assert max(checks["total_vs_nodal_abs_errors"]) < 1.0e-12
    assert checks["base_force_identity_abs_error"] < 1.0e-12
    assert checks["surface_area_vector_norm"] < 1.0e-15

    summary = {
        "kind": "surface_maxwell_force_trace_validation",
        "validation_class": True,
        "mesh_summary": mesh.summary(),
        "field_B_T": B,
        "face_rows": face_rows,
        "nodal_force_loads_N": nodal_loads,
        "checks": checks,
    }
    OUT_JSON.write_text(json.dumps(add_result_metadata(summary, __file__), indent=2), encoding="utf-8")

    print("[surface Maxwell traction trace]")
    print(
        f"  tris={checks['surface_triangles']} p(1T)={pressure:.12g} Pa "
        f"net_force={total_force} |net|={checks['total_force_norm_N']:.3e} N"
    )
    print(
        f"  base area={checks['base_area_m2']:.6g} "
        f"base force={checks['base_force_N']} "
        f"identity_error={checks['base_force_identity_abs_error']:.3e}"
    )
    print(f"  nodal total={nodal_total} |nodal|={checks['nodal_total_force_norm_N']:.3e} N")
    print(f"[OK] wrote {OUT_JSON}")


if __name__ == "__main__":
    main()
