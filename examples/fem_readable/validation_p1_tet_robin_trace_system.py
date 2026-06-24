"""Validation-class P1 tet Robin trace assembly example.

This is a small MATLAB/Gypsilab-style teaching script for a scalar FEM/BEM
interface view.  A Netgen ``.vol`` unit tetrahedron supplies the shared
one-based volume and boundary node ids.  The script assembles

    int_Omega grad(u).grad(v) + int_Gamma r u v
      = int_Gamma g v

with ``g = r`` on each boundary, so the constant field ``u = 1`` is the exact
solution.  The boundary term is just the P1 surface mass matrix restricted to
the `.vol` trace nodes.

Run:

    python examples/fem_readable/validation_p1_tet_robin_trace_system.py
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
SRC = REPO / "packages" / "radia-mcp" / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from radia_mcp.radia_ngsolve.netgen_vol import parse_netgen_tri_tet_vol  # noqa: E402
from radia_mcp.radia_ngsolve.scalar_fem3d import assemble_p1_tet_robin_system  # noqa: E402


OUT_JSON = HERE / "validation_p1_tet_robin_trace_system_summary.json"


UNIT_TET_ROBIN_VOL = """\
mesh3d
dimension
3
geomtype
0
facedescriptors
2
1 1 0 1 1
2 2 0 1 1
surfaceelements
4
1 1 1 0 3 1 2 3
1 2 1 0 3 1 4 2
1 2 1 0 3 2 4 3
1 2 1 0 3 3 4 1
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
2
1 base
2 sides
endmesh
"""


def _max_abs(values) -> float:
    return float(np.max(np.abs(values)))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    mesh = parse_netgen_tri_tet_vol(UNIT_TET_ROBIN_VOL, source="embedded_unit_tet_robin.vol")
    robin = {1: 1.0, 2: 2.0}
    system = assemble_p1_tet_robin_system(
        mesh.points,
        mesh.tetrahedra,
        mesh.surface_triangles,
        volume_coeff={1: 1.0},
        source={1: 0.0},
        robin_coeff=robin,
        boundary_flux=robin,
    )
    pure_neumann = assemble_p1_tet_robin_system(
        mesh.points,
        mesh.tetrahedra,
        mesh.surface_triangles,
        source=7.0,
    )
    topology = mesh.first_order_fem_bem_topology()
    closure = mesh.surface_closure_summary()
    manifold = mesh.surface_edge_manifold_summary()

    matrix = np.array(system["matrix"], dtype=float)
    rhs = np.array(system["rhs"], dtype=float)
    ones = np.ones(system["node_count"])
    residual = matrix @ ones - rhs
    energy = float(ones @ matrix @ ones)
    rhs_sum = float(rhs.sum())
    pure_matrix = np.array(pure_neumann["matrix"], dtype=float)
    pure_rhs = np.array(pure_neumann["rhs"], dtype=float)

    base_area = mesh.surface_area_by_boundary_number()[1]
    side_area = mesh.surface_area_by_boundary_number()[2]
    expected_side_area = 1.0 + 0.5 * math.sqrt(3.0)
    expected_robin_weight = base_area + 2.0 * side_area
    volume = mesh.total_volume()

    checks = {
        "base_area": base_area,
        "side_area": side_area,
        "expected_side_area": expected_side_area,
        "side_area_abs_error": abs(side_area - expected_side_area),
        "volume": volume,
        "constant_solution_residual_max_abs": _max_abs(residual),
        "energy_equals_robin_weight_abs_error": abs(energy - expected_robin_weight),
        "rhs_sum_equals_flux_weight_abs_error": abs(rhs_sum - expected_robin_weight),
        "pure_neumann_row_sum_max_abs": _max_abs(pure_matrix.sum(axis=1)),
        "pure_source_rhs_sum_abs_error": abs(float(pure_rhs.sum()) - 7.0 * volume),
        "trace_identity": (
            topology["trace"]["h1_to_scalar_bem_rows"] == [1, 2, 3, 4]
            and topology["trace"]["h1_to_scalar_bem_cols"] == [1, 2, 3, 4]
        ),
        "closed_surface": bool(manifold["is_closed_manifold"]),
        "surface_volume_abs_error": abs(abs(closure["surface_signed_volume"]) - volume),
    }

    assert checks["side_area_abs_error"] < 1.0e-15
    assert checks["constant_solution_residual_max_abs"] < 1.0e-15
    assert checks["energy_equals_robin_weight_abs_error"] < 1.0e-15
    assert checks["rhs_sum_equals_flux_weight_abs_error"] < 1.0e-15
    assert checks["pure_neumann_row_sum_max_abs"] < 1.0e-15
    assert checks["pure_source_rhs_sum_abs_error"] < 1.0e-15
    assert checks["trace_identity"]
    assert checks["closed_surface"]
    assert checks["surface_volume_abs_error"] < 1.0e-15

    summary = {
        "kind": "p1_tet_robin_trace_system_validation",
        "validation_class": True,
        "mesh_summary": mesh.summary(),
        "boundary_rows": mesh.boundary_summary_rows(),
        "robin_coeff_by_boundary": robin,
        "system_policy": system["policy"],
        "matrix": system["matrix"],
        "rhs": system["rhs"],
        "matrix_triplets": system["matrix_triplets"],
        "volume_by_material": system["volume_by_material"],
        "boundary_area_by_number": system["boundary_area_by_number"],
        "robin_area_weight": system["robin_area_weight"],
        "flux_area_weight": system["flux_area_weight"],
        "constant_solution_residual": residual.tolist(),
        "pure_neumann_row_sums": pure_matrix.sum(axis=1).tolist(),
        "trace": topology["trace"],
        "checks": checks,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    print("[P1 tet Robin trace system]")
    print(
        f"  nodes={system['node_count']} tets={mesh.summary()['tetrahedra']} "
        f"tris={mesh.summary()['surface_triangles']} volume={volume:.15f}"
    )
    print(
        f"  base area={base_area:.15f} side area={side_area:.15f} "
        f"robin weight={system['robin_area_weight']:.15f}"
    )
    print(
        f"  constant residual={checks['constant_solution_residual_max_abs']:.3e} "
        f"pure Neumann row sum={checks['pure_neumann_row_sum_max_abs']:.3e}"
    )
    print(f"  triplets={len(system['matrix_triplets'])} trace_identity={checks['trace_identity']}")
    print(f"[OK] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
