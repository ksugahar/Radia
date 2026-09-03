"""Validate curve-order-2 HEX/WEDGE RT2 on a compute host.

The input ``.vol`` files are supplied explicitly because curved HEX/WEDGE
meshes are produced by the Cubit export workflow and are not tiny CI fixtures.
Use ``--solve`` on hibino first, or on mdx only when hibino is unavailable and
the mdx CI queue is idle. The JSON records the host, topology, dimensions,
build/solve timings, and finite field-evaluation gate.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
import platform
import socket
import time
from pathlib import Path

import numpy as np
import ngsolve as ng
import radia

from radia.vim import ChargeGram, FieldFromSolution, MagnetizationSource, Solve


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


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


def _ngsolve_surface_field(source, mesh, points, *, intorder=12):
    """Independent exterior H from the projected M.n boundary charge."""
    normal = ng.specialcf.normal(3)
    position = ng.CoefficientFunction((ng.x, ng.y, ng.z))
    sigma = ng.InnerProduct(source.magnetization, normal)
    values = []
    for point in np.asarray(points, dtype=float):
        delta = ng.CoefficientFunction(tuple(float(value) for value in point)) - position
        radius = ng.sqrt(ng.InnerProduct(delta, delta))
        integrand = sigma * delta / (4.0*np.pi*radius**3)
        values.append(ng.Integrate(integrand, mesh, ng.BND, order=intorder))
    return np.asarray(values, dtype=float)


def _run(path: Path, topology: str, solve: bool) -> dict:
    mesh = ng.Mesh(str(path))
    vertex_count = {"hex": 8, "wedge": 6}[topology]
    counts = {len(el.vertices) for el in mesh.Elements(ng.VOL)}
    if counts != {vertex_count}:
        raise ValueError(f"{path}: expected pure {topology}, got vertex counts {sorted(counts)}")
    if int(mesh.GetCurveOrder()) != 2:
        raise ValueError(f"{path}: expected curve order 2, got {mesh.GetCurveOrder()}")

    t0 = time.perf_counter()
    with ng.TaskManager():
        fes = ng.HDiv(mesh, order=2)
        B, gram, _ = ChargeGram(
            fes, _materialize_mass=False, _build_hmatrix=False)
    entries = [float(gram.entry(i, j)) for i, j in ((0, 0), (0, 1), (1, 0))]
    record = {
        "mesh": path.name,
        "mesh_sha256": _sha256(path),
        "topology": topology,
        "curve_order": 2,
        "hdiv_order": 2,
        "elements": int(mesh.ne),
        "fes_dofs": int(fes.ndof),
        "charge_dofs": int(gram.ndof()),
        "charge_shape": list(B.shape),
        "direct_entry_setup_and_sample_s": time.perf_counter() - t0,
        "sample_entries": entries,
        "sample_symmetry_error": abs(entries[1] - entries[2]),
        "sample_entries_finite": bool(np.isfinite(entries).all()),
        "state_canary": gram.hex_state_check(),
    }
    if solve:
        t0 = time.perf_counter()
        with ng.TaskManager():
            result = Solve(
                mesh, mu_r=100.0, H_ext=ng.CF((0, 0, 1000.0)),
                order=2, curve_order=2, tol=1e-9)
        record.update({
            "solve_s": time.perf_counter() - t0,
            "iterations": int(result["iters"]),
            "demag": float(result["demag"]),
            "M_avg": np.asarray(result["M_avg"], dtype=float).tolist(),
            "solve_settings": {
                "mu_r": 100.0,
                "H_ext_A_per_m": [0.0, 0.0, 1000.0],
                "tol": 1.0e-9,
            },
        })
        for key in (
            "linear_solver", "preconditioner", "setup_wall_s", "solve_wall_s",
            "post_wall_s", "total_wall_s_internal", "fes_wall_s",
            "charge_gram_wall_s", "charge_basis_wall_s", "charge_gram_cpp_wall_s",
            "hex_state_check_wall_s", "projection_wall_s", "demag_probe_wall_s",
            "cpp_solve_timings", "hmat_stats",
        ):
            if key in result and result[key] is not None:
                record[key] = _jsonable(result[key])
        vertices = np.asarray([mesh[v].point for v in mesh.vertices], dtype=float)
        center = 0.5 * (vertices.min(axis=0) + vertices.max(axis=0))
        span = np.maximum(vertices.max(axis=0) - vertices.min(axis=0), 1.0e-12)
        points = np.asarray([
            center + (2.0 * span[0], 0.0, 0.0),
            center + (0.0, 0.0, 2.0 * span[2]),
        ])
        field_t0 = time.perf_counter()
        field = np.asarray(FieldFromSolution(result, points, algorithm="direct"))
        record["field_eval_2_points_s"] = time.perf_counter() - field_t0
        record["field_finite"] = bool(np.isfinite(field).all())

        prescribed = ng.CoefficientFunction((1.0e5, 2.0e5, 3.0e5))
        with ng.TaskManager():
            source = MagnetizationSource(
                mesh, prescribed, order=2, curve_order=2, curve_gauss=8)
            field_reference = _ngsolve_surface_field(source, mesh, points)
        source_field = source.Field(points, algorithm="direct")
        reference_scale = np.maximum(np.linalg.norm(field_reference, axis=1), 1.0)
        reference_error = np.linalg.norm(source_field-field_reference, axis=1) / reference_scale
        record.update({
            "prescribed_source_kind": source.stats["field_evaluator"]["source_kind"],
            "prescribed_projection_relative_residual": source.stats["projection_relative_residual"],
            "prescribed_field_ngsolve_reference_relative_error_max": float(reference_error.max()),
            "prescribed_field_ngsolve_reference_pass": bool(reference_error.max() < 1.0e-5),
        })
        if not record["prescribed_field_ngsolve_reference_pass"]:
            raise RuntimeError(
                f"{topology} curved RT2 prescribed-source field differs from the independent "
                f"NGSolve surface integral by {reference_error.max():.3e}")
    return record


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--hex-vol", type=Path)
    parser.add_argument("--wedge-vol", type=Path)
    parser.add_argument("--solve", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not args.hex_vol and not args.wedge_vol:
        parser.error("at least one of --hex-vol or --wedge-vol is required")

    data = {
        "schema": "radia.hdiv_rt2_curved.v1",
        "timing_scope": "single-run compute-host validation; do not infer host performance from one run",
        "generated_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "hostname": socket.gethostname(),
        "platform": platform.platform(),
        "python_version": platform.python_version(),
        "radia_version": radia.__version__,
        "ngsolve_version": getattr(ng, "__version__", "unknown"),
        "cpu_count": os.cpu_count(),
        "ngsolve_threads": int(ng.ngsglobals.numthreads),
        "cases": [],
    }
    if args.hex_vol:
        data["cases"].append(_run(args.hex_vol, "hex", args.solve))
    if args.wedge_vol:
        data["cases"].append(_run(args.wedge_vol, "wedge", args.solve))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(data, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
