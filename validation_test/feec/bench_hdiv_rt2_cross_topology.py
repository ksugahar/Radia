"""Compare RT2 TET/HEX/WEDGE at approximately 29k HDiv DoF.

Run one topology per process on hibino first, or on mdx only when hibino is
unavailable and its CI queue is idle. All cases use the
same unit cube and solver settings.  ``--curve-order 2`` exercises each Q2
geometry path on an affine cube, isolating implementation cost from CAD shape
and mesh-quality differences.  Field timing records the cold first call and
the steady-state median separately.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import platform
import time

import ngsolve as ng
from ngsolve.meshes import MakeStructured3DMesh
import numpy as np

import radia
from radia.vim import FieldFromSolution, Solve


_CASES = {
    "tet": {"kwargs": {"hexes": False, "nx": 7, "ny": 6, "nz": 6},
            "vertices": 4, "expected_dofs": 28656},
    "hex": {"kwargs": {"hexes": True, "nx": 7, "ny": 7, "nz": 7},
            "vertices": 8, "expected_dofs": 29106},
    "wedge": {"kwargs": {"prism": True, "nx": 6, "ny": 6, "nz": 7},
              "vertices": 6, "expected_dofs": 29160},
}


def _jsonable(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    return value


def _require_compute_host() -> str:
    hostname = platform.node()
    short = hostname.lower().split(".", 1)[0]
    if short not in {"mdx", "hibino"}:
        raise SystemExit(
            "bench_hdiv_rt2_cross_topology.py must run on mdx or hibino; "
            f"got {hostname!r}")
    return hostname


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--topology", choices=sorted(_CASES), required=True)
    parser.add_argument("--curve-order", type=int, choices=(1, 2), default=2)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    hostname = _require_compute_host()
    spec = _CASES[args.topology]
    mesh = MakeStructured3DMesh(**spec["kwargs"])
    if args.curve_order == 2:
        mesh.Curve(2)
    vertex_counts = {len(element.vertices) for element in mesh.Elements(ng.VOL)}
    if vertex_counts != {spec["vertices"]}:
        raise RuntimeError(
            f"{args.topology}: expected vertex count {spec['vertices']}, "
            f"got {sorted(vertex_counts)}")

    t0 = time.perf_counter()
    with ng.TaskManager():
        result = Solve(
            mesh,
            mu_r=100.0,
            H_ext=ng.CF((0.0, 0.0, 1000.0)),
            order=2,
            curve_order=2 if args.curve_order == 2 else 0,
            tol=1.0e-9,
        )
    external_solve_s = time.perf_counter()-t0
    if int(result["ndof"]) != spec["expected_dofs"]:
        raise RuntimeError(
            f"{args.topology}: expected {spec['expected_dofs']} HDiv DoF, "
            f"got {result['ndof']}")

    points = np.asarray([[2.5, 0.5, 0.5], [0.5, 0.5, 2.5]])
    field_samples = []
    field_reference = None
    for _ in range(7):
        t0 = time.perf_counter()
        field = np.asarray(
            FieldFromSolution(result, points, algorithm="direct"), dtype=float)
        field_samples.append(time.perf_counter()-t0)
        if field_reference is None:
            field_reference = field
        elif not np.array_equal(field, field_reference):
            raise RuntimeError("repeated direct field evaluation changed numerically")

    record = {
        "schema": "radia.hdiv_rt2_cross_topology.v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "hostname": hostname,
        "platform": platform.platform(),
        "python_version": platform.python_version(),
        "radia_version": radia.__version__,
        "ngsolve_version": ng.__version__,
        "ngsolve_threads": int(ng.ngsglobals.numthreads),
        "topology": args.topology,
        "geometry": "unit cube",
        "structured_divisions": spec["kwargs"],
        "curve_order": int(mesh.GetCurveOrder()),
        "hdiv_order": 2,
        "elements": int(mesh.ne),
        "fes_dofs": int(result["ndof"]),
        "charge_dofs": int(result["n_charge"]),
        "solve_external_s": external_solve_s,
        "field_eval_2_points_first_s": field_samples[0],
        "field_eval_2_points_steady_median_s": float(np.median(field_samples[2:])),
        "field_eval_2_points_samples_s": field_samples,
        "field_finite": bool(np.isfinite(field_reference).all()),
        "solve_settings": {
            "mu_r": 100.0,
            "H_ext_A_per_m": [0.0, 0.0, 1000.0],
            "tol": 1.0e-9,
        },
    }
    for key in (
        "iters", "demag", "M_avg", "linear_solver", "preconditioner",
        "setup_wall_s", "solve_wall_s", "post_wall_s", "total_wall_s_internal",
        "fes_wall_s", "charge_gram_wall_s", "charge_basis_wall_s",
        "charge_gram_cpp_wall_s", "hex_state_check_wall_s", "projection_wall_s",
        "demag_probe_wall_s", "cpp_solve_timings", "hmat_stats",
        "field_evaluator_build_wall_s", "field_evaluator_stats",
    ):
        if key in result and result[key] is not None:
            record[key] = _jsonable(result[key])

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(record, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(record, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
