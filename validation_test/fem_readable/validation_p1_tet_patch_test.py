"""Validation-class P1 tetrahedron patch test.

The mesh is a unit cube split into 12 tetrahedra by connecting every boundary
triangle to one interior node.  With Dirichlet values from an affine function,
the scalar P1 stiffness assembly should recover the exact affine value at the
interior node.  This is the small 3D patch test that a readable MATLAB FEM file
should pass before growing into a full solver.

Run:

    python validation_test/fem_readable/validation_p1_tet_patch_test.py
"""

from __future__ import annotations

import json
from result_metadata import add_result_metadata
import sys
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
SRC = REPO / "packages" / "radia-mcp" / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from radia_mcp.radia_ngsolve.scalar_fem3d import (  # noqa: E402
    p1_tetrahedron_geometry,
    p1_tetrahedron_stiffness,
)


OUT_JSON = HERE / "validation_p1_tet_patch_test_summary.json"


BOUNDARY_VERTICES = [
    (0.0, 0.0, 0.0),
    (1.0, 0.0, 0.0),
    (1.0, 1.0, 0.0),
    (0.0, 1.0, 0.0),
    (0.0, 0.0, 1.0),
    (1.0, 0.0, 1.0),
    (1.0, 1.0, 1.0),
    (0.0, 1.0, 1.0),
]


FACE_TRIANGLES = [
    (0, 1, 2), (0, 2, 3),
    (4, 6, 5), (4, 7, 6),
    (0, 5, 1), (0, 4, 5),
    (3, 2, 6), (3, 6, 7),
    (0, 3, 7), (0, 7, 4),
    (1, 5, 6), (1, 6, 2),
]


AFFINE_FIELDS = [
    {"label": "x", "coeff": (0.0, 1.0, 0.0, 0.0)},
    {"label": "y", "coeff": (0.0, 0.0, 1.0, 0.0)},
    {"label": "z", "coeff": (0.0, 0.0, 0.0, 1.0)},
    {"label": "mixed", "coeff": (0.25, 1.0, 2.0, 3.0)},
]


PATCHES = [
    {"label": "centered", "interior": (0.5, 0.5, 0.5)},
    {"label": "distorted_inside", "interior": (0.43, 0.57, 0.48)},
]


def affine_value(coeff, p):
    c0, cx, cy, cz = coeff
    return c0 + cx * p[0] + cy * p[1] + cz * p[2]


def build_patch(interior):
    vertices = BOUNDARY_VERTICES + [tuple(interior)]
    center = len(vertices) - 1
    tets = [(center, a, b, c) for (a, b, c) in FACE_TRIANGLES]
    return vertices, tets


def assemble_stiffness(vertices, tets):
    k = np.zeros((len(vertices), len(vertices)))
    volumes = []
    for tet in tets:
        coords = [vertices[i] for i in tet]
        ke = p1_tetrahedron_stiffness(coords)
        volumes.append(p1_tetrahedron_geometry(coords)["volume"])
        for a, ia in enumerate(tet):
            for b, ib in enumerate(tet):
                k[ia, ib] += ke[a][b]
    return k, volumes


def solve_one_interior(k, boundary_values):
    center = k.shape[0] - 1
    rhs = -sum(k[center, j] * boundary_values[j] for j in range(center))
    return rhs / k[center, center]


def run_patch(patch):
    vertices, tets = build_patch(patch["interior"])
    k, volumes = assemble_stiffness(vertices, tets)
    rows = []
    for field in AFFINE_FIELDS:
        coeff = field["coeff"]
        boundary_values = [affine_value(coeff, p) for p in vertices[:-1]]
        computed = solve_one_interior(k, boundary_values)
        exact = affine_value(coeff, vertices[-1])
        rows.append({
            "field": field["label"],
            "computed": computed,
            "exact": exact,
            "abs_error": abs(computed - exact),
        })
    return {
        "label": patch["label"],
        "interior": list(patch["interior"]),
        "n_vertices": len(vertices),
        "n_tetrahedra": len(tets),
        "volume_sum": sum(volumes),
        "min_tet_volume": min(volumes),
        "max_tet_volume": max(volumes),
        "stiffness_row_sum_max_abs": float(np.max(np.abs(k.sum(axis=1)))),
        "fields": rows,
        "max_abs_error": max(row["abs_error"] for row in rows),
    }


def main() -> int:
    patches = [run_patch(patch) for patch in PATCHES]
    summary = {
        "kind": "p1_tetrahedron_patch_test_validation",
        "validation_class": True,
        "patches": patches,
        "max_abs_error": max(patch["max_abs_error"] for patch in patches),
    }
    OUT_JSON.write_text(json.dumps(add_result_metadata(summary, __file__), indent=2), encoding="utf-8")
    for patch in patches:
        print(
            f"[{patch['label']}] tets={patch['n_tetrahedra']} "
            f"volume={patch['volume_sum']:.12f} "
            f"row_sum={patch['stiffness_row_sum_max_abs']:.3e} "
            f"max_err={patch['max_abs_error']:.3e}"
        )
        for row in patch["fields"]:
            print(
                f"  {row['field']:5s}: computed={row['computed']:.15g} "
                f"exact={row['exact']:.15g} err={row['abs_error']:.3e}"
            )
    print(f"[OK] wrote {OUT_JSON}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
