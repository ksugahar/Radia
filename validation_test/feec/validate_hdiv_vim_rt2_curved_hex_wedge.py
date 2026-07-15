"""Validate curve-order-2 HEX/WEDGE RT2 on a compute host.

The input ``.vol`` files are supplied explicitly because curved HEX/WEDGE
meshes are produced by the Cubit export workflow and are not tiny CI fixtures.
Use ``--solve`` only on an idle mdx or hibino host.  The JSON records the host,
topology, dimensions, build/solve timings, and finite field-evaluation gate.
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import platform
import socket
import time
from pathlib import Path

import numpy as np
import ngsolve as ng
import radia

from radia.vim import ChargeGram, FieldFromSolution, Solve


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
        })
        vertices = np.asarray([mesh[v].point for v in mesh.vertices], dtype=float)
        center = 0.5 * (vertices.min(axis=0) + vertices.max(axis=0))
        span = np.maximum(vertices.max(axis=0) - vertices.min(axis=0), 1.0e-12)
        points = np.asarray([
            center + (2.0 * span[0], 0.0, 0.0),
            center + (0.0, 0.0, 2.0 * span[2]),
        ])
        field = np.asarray(FieldFromSolution(result, points, algorithm="direct"))
        record["field_finite"] = bool(np.isfinite(field).all())
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
        "generated_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "hostname": socket.gethostname(),
        "platform": platform.platform(),
        "python_version": platform.python_version(),
        "radia_version": radia.__version__,
        "ngsolve_version": getattr(ng, "__version__", "unknown"),
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
