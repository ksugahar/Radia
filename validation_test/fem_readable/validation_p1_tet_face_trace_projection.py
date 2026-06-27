"""Validation-class P1 tetrahedron face trace projection.

Run:

    python validation_test/fem_readable/validation_p1_tet_face_trace_projection.py

This is the local FEM/BEM coupling block used when a volume H1/P1 solution is
restricted to a boundary triangle and projected with the boundary P1 mass
matrix:

    b_j = int_face u_h N_j dS,    ||u_h||^2_face = int_face u_h^2 dS

The example uses a unit tetrahedron from a tiny Netgen `.vol` mesh, so the
volume and boundary unknowns share one-based node ids exactly as a MATLAB
teaching prototype should.
"""

from __future__ import annotations

import json
from result_metadata import add_result_metadata
import math
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
SRC = REPO / "packages" / "radia-mcp" / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from radia_mcp.radia_ngsolve.netgen_vol import parse_netgen_tri_tet_vol  # noqa: E402
from radia_mcp.radia_ngsolve.scalar_fem3d import (  # noqa: E402
    p1_tetrahedron_face_trace_summary,
)


OUT_JSON = HERE / "validation_p1_tet_face_trace_projection_summary.json"

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


def _assert_close(actual: float, expected: float, *, rtol: float = 1.0e-12, atol: float = 1.0e-15) -> None:
    if abs(actual - expected) > max(atol, rtol * max(abs(actual), abs(expected))):
        raise AssertionError(f"{actual!r} != {expected!r}")


def _surface_l2_exact(area: float, values: list[float]) -> float:
    pair_sum = sum(values[i] * values[j] for i in range(3) for j in range(i + 1, 3))
    return area / 6.0 * (sum(value * value for value in values) + pair_sum)


def _json_clean(value):
    if isinstance(value, float):
        return 0.0 if value == 0.0 else value
    if isinstance(value, list):
        return [_json_clean(item) for item in value]
    if isinstance(value, tuple):
        return [_json_clean(item) for item in value]
    if isinstance(value, dict):
        return {key: _json_clean(item) for key, item in value.items()}
    return value


def build_summary() -> dict:
    mesh = parse_netgen_tri_tet_vol(UNIT_TET_VOL, source="embedded_unit_tet_trace.vol")
    tet = mesh.tetrahedra[0]
    vertices = [mesh.points[node - 1] for node in tet.nodes]
    nodal_values = [1.0 + x + 2.0 * y + 3.0 * z for x, y, z in vertices]

    rows = []
    total_integral = 0.0
    total_l2 = 0.0
    max_integral_error = 0.0
    max_l2_error = 0.0
    max_load_sum_error = 0.0

    for tri in mesh.surface_triangles:
        face_local_nodes = tuple(tet.nodes.index(node) for node in tri.nodes)
        row = p1_tetrahedron_face_trace_summary(vertices, nodal_values, face_local_nodes)
        trace = row["trace_nodal_values"]
        expected_integral = row["area"] * sum(trace) / 3.0
        expected_l2 = _surface_l2_exact(row["area"], trace)
        load_sum = sum(row["projected_trace_load"])
        integral_error = abs(row["trace_integral"] - expected_integral)
        l2_error = abs(row["trace_l2_norm_squared"] - expected_l2)
        load_sum_error = abs(load_sum - row["trace_integral"])
        _assert_close(row["trace_integral"], expected_integral)
        _assert_close(row["trace_l2_norm_squared"], expected_l2)
        _assert_close(load_sum, row["trace_integral"])
        total_integral += row["trace_integral"]
        total_l2 += row["trace_l2_norm_squared"]
        max_integral_error = max(max_integral_error, integral_error)
        max_l2_error = max(max_l2_error, l2_error)
        max_load_sum_error = max(max_load_sum_error, load_sum_error)
        rows.append({
            "bcnr": tri.bcnr,
            "one_based_surface_nodes": list(tri.nodes),
            "face_local_nodes": list(face_local_nodes),
            "expected_integral": expected_integral,
            "expected_l2_norm_squared": expected_l2,
            "summary": row,
        })

    expected_total_integral = 3.5 + 1.5 * math.sqrt(3.0)
    _assert_close(total_integral, expected_total_integral)

    return {
        "kind": "p1_tet_face_trace_projection_validation",
        "validation_class": True,
        "mesh_summary": mesh.summary(),
        "nodal_values": nodal_values,
        "rows": rows,
        "checks": {
            "total_trace_integral": total_integral,
            "expected_total_trace_integral": expected_total_integral,
            "total_trace_l2_norm_squared": total_l2,
            "max_integral_error": max_integral_error,
            "max_l2_error": max_l2_error,
            "max_load_sum_error": max_load_sum_error,
            "row_count": len(rows),
            "surface_triangle_count": mesh.summary()["surface_triangles"],
            "policy": "shared_vol_node_ids_volume_p1_to_boundary_p1_trace",
        },
    }


def main() -> int:
    summary = build_summary()
    OUT_JSON.write_text(json.dumps(_json_clean(add_result_metadata(summary, __file__)), indent=2, sort_keys=True) + "\n", encoding="utf-8")

    checks = summary["checks"]
    print("[P1 tet face trace projection]")
    print(f"  rows={checks['row_count']} total int(u)={checks['total_trace_integral']:.12g}")
    print(f"  expected int(u)={checks['expected_total_trace_integral']:.12g}")
    print(f"  total int(u^2)={checks['total_trace_l2_norm_squared']:.12g}")
    print(
        f"  max errors: integral={checks['max_integral_error']:.3e}, "
        f"l2={checks['max_l2_error']:.3e}, load-sum={checks['max_load_sum_error']:.3e}"
    )
    print(f"  wrote {OUT_JSON}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
