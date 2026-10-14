"""Validation-class P1 surface-triangle trace coupling example.

This is an example/validation run, not a pytest test.  It uses a tiny Netgen
``.vol`` unit tetrahedron to show the FEM/BEM interface view that the readable
MATLAB prototype should mirror:

* volume H1 nodes and boundary P1 nodes share the same one-based node ids
* boundary triangles assemble a SurfaceL2/P1 mass matrix and load vector
* surface stiffness integrates tangential gradients on 3D triangles
* closed-surface area-vector and volume checks catch normal-orientation issues

Run:

    python validation_test/fem_readable/validation_p1_surface_trace_coupling.py
"""

from __future__ import annotations

import argparse
import json
from result_metadata import add_result_metadata
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
from radia_mcp.radia_ngsolve.scalar_fem3d import (  # noqa: E402
    p1_surface_triangle_constant_load,
    p1_surface_triangle_geometry,
    p1_surface_triangle_mass,
    p1_surface_triangle_stiffness,
)


OUT_JSON = HERE / "validation_p1_surface_trace_coupling_summary.json"


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


def _assemble_boundary(mesh):
    n = len(mesh.points)
    mass = np.zeros((n, n))
    stiffness = np.zeros((n, n))
    load = np.zeros(n)
    face_rows = []
    area_vector = np.zeros(3)

    for tri in mesh.surface_triangles:
        ids = [node - 1 for node in tri.nodes]
        coords = [mesh.points[i] for i in ids]
        geom = p1_surface_triangle_geometry(coords)
        mloc = p1_surface_triangle_mass(coords)
        kloc = p1_surface_triangle_stiffness(coords)
        floc = p1_surface_triangle_constant_load(coords)
        area_vector += np.array(geom["area_vector"])
        for a, ia in enumerate(ids):
            load[ia] += floc[a]
            for b, ib in enumerate(ids):
                mass[ia, ib] += mloc[a][b]
                stiffness[ia, ib] += kloc[a][b]
        face_rows.append({
            "nodes": list(tri.nodes),
            "area": geom["area"],
            "unit_normal": list(geom["unit_normal"]),
            "area_vector": list(geom["area_vector"]),
        })

    return mass, stiffness, load, face_rows, area_vector


def _coordinate_moments(mesh, mass):
    ones = np.ones(len(mesh.points))
    values = np.array(mesh.points)
    via_mass = ones @ mass @ values
    direct = np.zeros(3)
    for tri in mesh.surface_triangles:
        ids = [node - 1 for node in tri.nodes]
        coords = np.array([mesh.points[i] for i in ids])
        area = p1_surface_triangle_geometry(coords)["area"]
        direct += area * coords.mean(axis=0)
    return via_mass, direct


def _coordinate_surface_energies(mesh, stiffness):
    values = np.array(mesh.points)
    via_matrix = np.array([values[:, i] @ stiffness @ values[:, i] for i in range(3)])
    direct = np.zeros(3)
    basis = np.eye(3)
    for tri in mesh.surface_triangles:
        ids = [node - 1 for node in tri.nodes]
        coords = [mesh.points[i] for i in ids]
        geom = p1_surface_triangle_geometry(coords)
        normal = np.array(geom["unit_normal"])
        for i, e in enumerate(basis):
            tangent = e - np.dot(e, normal) * normal
            direct[i] += geom["area"] * np.dot(tangent, tangent)
    return via_matrix, direct


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    mesh = parse_netgen_tri_tet_vol(UNIT_TET_VOL, source="embedded_unit_tetrahedron.vol")
    mass, stiffness, load, face_rows, area_vector = _assemble_boundary(mesh)
    topology = mesh.first_order_fem_bem_topology()
    closure = mesh.surface_closure_summary()
    manifold = mesh.surface_edge_manifold_summary()

    ones = np.ones(len(mesh.points))
    total_area = mesh.total_surface_area()
    mass_area = float(ones @ mass @ ones)
    load_area = float(load.sum())
    moment_mass, moment_direct = _coordinate_moments(mesh, mass)
    energy_matrix, energy_direct = _coordinate_surface_energies(mesh, stiffness)

    expected_area = 1.5 + 0.5 * math.sqrt(3.0)
    expected_volume = 1.0 / 6.0
    expected_first_moment = (2.0 + math.sqrt(3.0)) / 6.0
    expected_coordinate_energy = 1.0 + math.sqrt(3.0) / 3.0

    checks = {
        "total_area": total_area,
        "expected_area": expected_area,
        "area_abs_error": abs(total_area - expected_area),
        "mass_area_abs_error": abs(mass_area - total_area),
        "load_area_abs_error": abs(load_area - total_area),
        "coordinate_moment_abs_errors": list(np.abs(moment_mass - moment_direct)),
        "coordinate_energy_abs_errors": list(np.abs(energy_matrix - energy_direct)),
        "stiffness_row_sum_max_abs": float(np.max(np.abs(stiffness.sum(axis=1)))),
        "surface_area_vector_norm": float(np.linalg.norm(area_vector)),
        "surface_volume_abs_error": abs(abs(closure["surface_signed_volume"]) - expected_volume),
        "expected_coordinate_first_moment": expected_first_moment,
        "expected_coordinate_energy": expected_coordinate_energy,
        "trace_is_identity_on_unit_tet": (
            topology["trace"]["h1_to_scalar_bem_rows"] == [1, 2, 3, 4]
            and topology["trace"]["h1_to_scalar_bem_cols"] == [1, 2, 3, 4]
        ),
        "rwg_hcurl_edge_id_match": topology["trace"]["rwg_to_hcurl_edge_ids"] == [1, 2, 3, 4, 5, 6],
    }

    assert checks["area_abs_error"] < 1.0e-15
    assert checks["mass_area_abs_error"] < 1.0e-15
    assert checks["load_area_abs_error"] < 1.0e-15
    assert max(checks["coordinate_moment_abs_errors"]) < 1.0e-15
    assert max(abs(v - expected_first_moment) for v in moment_mass) < 1.0e-15
    assert max(checks["coordinate_energy_abs_errors"]) < 1.0e-15
    assert max(abs(v - expected_coordinate_energy) for v in energy_matrix) < 1.0e-15
    assert checks["stiffness_row_sum_max_abs"] < 1.0e-15
    assert checks["surface_area_vector_norm"] < 1.0e-15
    assert checks["surface_volume_abs_error"] < 1.0e-15
    assert checks["trace_is_identity_on_unit_tet"]
    assert checks["rwg_hcurl_edge_id_match"]
    assert manifold["is_closed_manifold"]

    summary = {
        "kind": "p1_surface_trace_coupling_validation",
        "validation_class": True,
        "mesh_summary": mesh.summary(),
        "topology_policy": topology["policy"],
        "closure": closure,
        "manifold": manifold,
        "face_rows": face_rows,
        "surface_mass_matrix": mass.tolist(),
        "surface_stiffness_matrix": stiffness.tolist(),
        "surface_load_vector": load.tolist(),
        "coordinate_moments_via_mass": moment_mass.tolist(),
        "coordinate_moments_direct": moment_direct.tolist(),
        "coordinate_energies_via_stiffness": energy_matrix.tolist(),
        "coordinate_energies_direct": energy_direct.tolist(),
        "trace": topology["trace"],
        "checks": checks,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(add_result_metadata(summary, __file__), indent=2), encoding="utf-8")

    print("[surface trace]")
    print(
        f"  points={mesh.summary()['points']} tris={mesh.summary()['surface_triangles']} "
        f"tets={mesh.summary()['tetrahedra']} area={total_area:.15f}"
    )
    print(
        f"  mass ones={mass_area:.15f}, load sum={load_area:.15f}, "
        f"row_sum={checks['stiffness_row_sum_max_abs']:.3e}"
    )
    print(
        f"  moments={moment_mass.tolist()}, energies={energy_matrix.tolist()}"
    )
    print(
        f"  closure vector norm={checks['surface_area_vector_norm']:.3e}, "
        f"signed volume={closure['surface_signed_volume']:.15f}, "
        f"trace identity={checks['trace_is_identity_on_unit_tet']}"
    )
    print(f"[OK] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
