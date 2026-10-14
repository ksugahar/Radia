"""Validation-class readable P1 tetrahedron Poisson convergence example.

This is an example/validation run rather than a pytest test.  It assembles a
small 3D scalar Poisson problem from explicit P1 tetrahedron element matrices:

    -Delta u = 3*pi^2 sin(pi x) sin(pi y) sin(pi z),  u=0 on the cube boundary.

The exact solution is ``u = sin(pi x) sin(pi y) sin(pi z)``.  The example is
written to be easy to translate into MATLAB/Gypsilab-style teaching scripts:
explicit mesh generation, element loop, Dirichlet elimination, dense solve,
and quadrature-based L2/H1 error estimates.

Run:

    python validation_test/fem_readable/validation_p1_tet_poisson_convergence.py
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

from radia_mcp.radia_ngsolve.scalar_fem3d import (  # noqa: E402
    p1_tetrahedron_geometry,
    p1_tetrahedron_stiffness,
)


OUT_JSON = HERE / "validation_p1_tet_poisson_convergence_summary.json"
MESH_DIVISIONS = (4, 6, 8)

# Six tetrahedra around the body diagonal of each cube cell.
CELL_TETS = (
    (0, 1, 3, 7),
    (0, 3, 2, 7),
    (0, 2, 6, 7),
    (0, 6, 4, 7),
    (0, 4, 5, 7),
    (0, 5, 1, 7),
)

# Symmetric four-point degree-2 tetrahedron quadrature.
_A = 0.5854101966249685
_B = 0.1381966011250105
TET_QUAD_BARY = tuple(
    tuple(_A if i == heavy else _B for i in range(4))
    for heavy in range(4)
)


def exact_u(point: tuple[float, float, float]) -> float:
    x, y, z = point
    return math.sin(math.pi * x) * math.sin(math.pi * y) * math.sin(math.pi * z)


def exact_grad(point: tuple[float, float, float]) -> np.ndarray:
    x, y, z = point
    pi = math.pi
    return np.array([
        pi * math.cos(pi * x) * math.sin(pi * y) * math.sin(pi * z),
        pi * math.sin(pi * x) * math.cos(pi * y) * math.sin(pi * z),
        pi * math.sin(pi * x) * math.sin(pi * y) * math.cos(pi * z),
    ])


def source_f(point: tuple[float, float, float]) -> float:
    return 3.0 * math.pi ** 2 * exact_u(point)


def structured_tet_mesh(n: int) -> tuple[list[tuple[float, float, float]], list[tuple[int, int, int, int]]]:
    points: list[tuple[float, float, float]] = []
    point_id: dict[tuple[int, int, int], int] = {}
    for k in range(n + 1):
        for j in range(n + 1):
            for i in range(n + 1):
                point_id[(i, j, k)] = len(points)
                points.append((i / n, j / n, k / n))

    def cube_ids(i: int, j: int, k: int) -> list[int]:
        return [
            point_id[(i, j, k)],
            point_id[(i + 1, j, k)],
            point_id[(i, j + 1, k)],
            point_id[(i + 1, j + 1, k)],
            point_id[(i, j, k + 1)],
            point_id[(i + 1, j, k + 1)],
            point_id[(i, j + 1, k + 1)],
            point_id[(i + 1, j + 1, k + 1)],
        ]

    tets: list[tuple[int, int, int, int]] = []
    for k in range(n):
        for j in range(n):
            for i in range(n):
                ids = cube_ids(i, j, k)
                for tet in CELL_TETS:
                    tets.append(tuple(ids[q] for q in tet))
    return points, tets


def _quad_point(vertices: np.ndarray, bary: tuple[float, float, float, float]) -> np.ndarray:
    out = np.zeros(3)
    for i, lam in enumerate(bary):
        out += lam * vertices[i]
    return out


def _element_load(vertices: list[tuple[float, float, float]]) -> np.ndarray:
    geom = p1_tetrahedron_geometry(vertices)
    volume = geom["volume"]
    verts = np.array(vertices, dtype=float)
    fe = np.zeros(4)
    for bary in TET_QUAD_BARY:
        p = _quad_point(verts, bary)
        val = source_f((float(p[0]), float(p[1]), float(p[2])))
        for i, lam in enumerate(bary):
            fe[i] += 0.25 * volume * val * lam
    return fe


def _is_boundary(point: tuple[float, float, float]) -> bool:
    return any(abs(v) < 1.0e-14 or abs(v - 1.0) < 1.0e-14 for v in point)


def solve_mesh(n: int) -> dict:
    points, tets = structured_tet_mesh(n)
    n_points = len(points)
    K = np.zeros((n_points, n_points))
    F = np.zeros(n_points)

    for tet in tets:
        vertices = [points[i] for i in tet]
        Ke = np.array(p1_tetrahedron_stiffness(vertices))
        fe = _element_load(vertices)
        for a, ia in enumerate(tet):
            F[ia] += fe[a]
            for b, ib in enumerate(tet):
                K[ia, ib] += Ke[a, b]

    boundary = np.array([_is_boundary(p) for p in points])
    free = np.where(~boundary)[0]
    U = np.zeros(n_points)
    U[free] = np.linalg.solve(K[np.ix_(free, free)], F[free])

    l2_sq = 0.0
    h1_sq = 0.0
    for tet in tets:
        vertices = [points[i] for i in tet]
        geom = p1_tetrahedron_geometry(vertices)
        volume = geom["volume"]
        grads = np.array(geom["gradients"], dtype=float)
        verts = np.array(vertices, dtype=float)
        u_local = np.array([U[i] for i in tet])
        grad_uh = sum(u_local[i] * grads[i] for i in range(4))
        for bary in TET_QUAD_BARY:
            p = _quad_point(verts, bary)
            uh = float(sum(bary[i] * u_local[i] for i in range(4)))
            exact = exact_u((float(p[0]), float(p[1]), float(p[2])))
            grad_err = grad_uh - exact_grad((float(p[0]), float(p[1]), float(p[2])))
            l2_sq += 0.25 * volume * (uh - exact) ** 2
            h1_sq += 0.25 * volume * float(grad_err @ grad_err)

    return {
        "n": n,
        "h": 1.0 / n,
        "n_points": n_points,
        "n_tets": len(tets),
        "n_free_dofs": int(len(free)),
        "l2_error": math.sqrt(l2_sq),
        "h1_seminorm_error": math.sqrt(h1_sq),
    }


def _rate(e0: float, e1: float, h0: float, h1: float) -> float:
    return math.log(e0 / e1) / math.log(h0 / h1)


def _validate(rows: list[dict]) -> dict:
    rates = []
    for left, right in zip(rows, rows[1:]):
        rates.append({
            "from_n": left["n"],
            "to_n": right["n"],
            "l2_rate": _rate(left["l2_error"], right["l2_error"], left["h"], right["h"]),
            "h1_rate": _rate(left["h1_seminorm_error"], right["h1_seminorm_error"], left["h"], right["h"]),
        })
    checks = {
        "l2_error_ratio_first_over_last": rows[0]["l2_error"] / rows[-1]["l2_error"],
        "h1_error_ratio_first_over_last": rows[0]["h1_seminorm_error"] / rows[-1]["h1_seminorm_error"],
        "last_l2_rate": rates[-1]["l2_rate"],
        "last_h1_rate": rates[-1]["h1_rate"],
        "rates": rates,
    }
    assert checks["l2_error_ratio_first_over_last"] > 3.0
    assert checks["h1_error_ratio_first_over_last"] > 1.8
    assert checks["last_l2_rate"] > 1.80
    assert checks["last_h1_rate"] > 0.90
    return checks


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    args = parser.parse_args()

    rows = [solve_mesh(n) for n in MESH_DIVISIONS]
    checks = _validate(rows)
    summary = {
        "kind": "p1_tet_poisson_convergence_validation",
        "validation_class": True,
        "manufactured_solution": "sin(pi*x)*sin(pi*y)*sin(pi*z)",
        "mesh_divisions": list(MESH_DIVISIONS),
        "rows": rows,
        "checks": checks,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(add_result_metadata(summary, __file__), indent=2), encoding="utf-8")

    print("[P1 tetrahedron Poisson convergence]")
    for row in rows:
        print(
            f"  n={row['n']:2d}  dofs={row['n_free_dofs']:4d}  "
            f"tets={row['n_tets']:5d}  "
            f"L2={row['l2_error']:.6e}  H1={row['h1_seminorm_error']:.6e}"
        )
    for rate in checks["rates"]:
        print(
            f"  rate n={rate['from_n']}->{rate['to_n']}: "
            f"L2={rate['l2_rate']:.3f}, H1={rate['h1_rate']:.3f}"
        )
    print(f"[OK] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
