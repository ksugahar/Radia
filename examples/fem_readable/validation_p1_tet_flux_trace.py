"""Validation-class P1 tetrahedron gradient and Neumann trace example.

This is the smallest readable bridge from a volume P1 FEM solution to the
boundary flux data that a BEM/trace script consumes: a single tetrahedron, an
affine scalar field, constant element flux, and four outward face flux rows.
"""

import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "packages" / "radia-mcp" / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from radia_mcp.radia_ngsolve.netgen_vol import parse_netgen_tri_tet_vol  # noqa: E402
from radia_mcp.radia_ngsolve.scalar_fem3d import (  # noqa: E402
    p1_tetrahedron_boundary_fluxes,
    p1_tetrahedron_flux,
    p1_tetrahedron_geometry,
    p1_tetrahedron_gradient,
    p1_tetrahedron_stiffness,
)


HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "validation_p1_tet_flux_trace_summary.json"

UNIT_TET_VOL = """
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


def _quad(values, matrix):
    return sum(values[i] * sum(matrix[i][j] * values[j] for j in range(4)) for i in range(4))


def main():
    mesh = parse_netgen_tri_tet_vol(UNIT_TET_VOL, source="embedded_unit_tet_flux.vol")
    tet = mesh.tetrahedra[0]
    vertices = [mesh.points[node - 1] for node in tet.nodes]
    coeff = 2.75
    exact_grad = (1.0, -2.0, 0.5)
    nodal = [
        0.25 + exact_grad[0] * x + exact_grad[1] * y + exact_grad[2] * z
        for x, y, z in vertices
    ]
    grad = p1_tetrahedron_gradient(vertices, nodal)
    flux = p1_tetrahedron_flux(vertices, nodal, coeff=coeff)
    face_rows = p1_tetrahedron_boundary_fluxes(vertices, nodal, coeff=coeff)
    stiffness = p1_tetrahedron_stiffness(vertices, coeff=coeff)
    volume = p1_tetrahedron_geometry(vertices)["volume"]
    energy_matrix = _quad(nodal, stiffness)
    energy_direct = coeff * volume * sum(value * value for value in exact_grad)
    total_integrated_flux = sum(row["integrated_flux"] for row in face_rows)
    max_face_identity_error = max(
        abs(
            row["integrated_flux"]
            - sum(flux[i] * row["outward_area_vector"][i] for i in range(3))
        )
        for row in face_rows
    )
    checks = {
        "volume": volume,
        "gradient_abs_errors": [abs(grad[i] - exact_grad[i]) for i in range(3)],
        "flux_abs_errors": [abs(flux[i] + coeff * exact_grad[i]) for i in range(3)],
        "energy_abs_error": abs(energy_matrix - energy_direct),
        "total_integrated_flux_abs": abs(total_integrated_flux),
        "max_face_identity_error": max_face_identity_error,
        "face_count": len(face_rows),
        "surface_triangle_count": mesh.summary()["surface_triangles"],
    }
    assert max(checks["gradient_abs_errors"]) < 1.0e-15
    assert max(checks["flux_abs_errors"]) < 1.0e-15
    assert checks["energy_abs_error"] < 1.0e-15
    assert checks["total_integrated_flux_abs"] < 1.0e-15
    assert checks["max_face_identity_error"] < 1.0e-15
    assert checks["face_count"] == checks["surface_triangle_count"] == 4

    summary = {
        "kind": "p1_tet_flux_trace_validation",
        "validation_class": True,
        "mesh_summary": mesh.summary(),
        "coeff": coeff,
        "nodal_values": nodal,
        "gradient": grad,
        "flux": flux,
        "face_flux_rows": face_rows,
        "energy_matrix": energy_matrix,
        "energy_direct": energy_direct,
        "checks": checks,
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print("[P1 tet flux trace]")
    print(f"  nodes={mesh.summary()['points']} tets={mesh.summary()['tetrahedra']} volume={volume:.15f}")
    print(f"  grad={grad} flux={flux}")
    print(f"  energy matrix/direct={energy_matrix:.15f}/{energy_direct:.15f}")
    print(f"  total integrated flux={total_integrated_flux:.3e}")
    for row in face_rows:
        print(
            f"  face opposite {row['opposite_local_node']}: nodes={row['face_local_nodes']} "
            f"area={row['area']:.6f} int_flux={row['integrated_flux']:.6f}"
        )
    print(f"[OK] wrote {OUT_JSON}")


if __name__ == "__main__":
    main()
