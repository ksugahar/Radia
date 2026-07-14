"""Compute-host scaling benchmark for solved HDiv ``rad.Fld``.

Run on an idle mdx (default) or hibino after a normal package release.  The
driver measures the public field call and checks the forced tree against the
exact direct source sum on a deterministic observation subset.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import platform
import subprocess
import time

os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

import ngsolve as ng
from ngsolve.meshes import MakeStructured3DMesh
import numpy as np
import radia as rad
from radia import vim


def _best(call, repeats=3):
    value = call()
    best = float("inf")
    for _ in range(repeats):
        started = time.perf_counter()
        value = call()
        best = min(best, time.perf_counter()-started)
    return best, np.asarray(value, float)


def _git_sha():
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], text=True,
            stderr=subprocess.DEVNULL).strip()
    except Exception:
        return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=8, help="structured hex cells per axis")
    parser.add_argument("--n-obs", type=int, default=10000)
    parser.add_argument("--reference-obs", type=int, default=128)
    parser.add_argument("--threads", type=int, default=os.cpu_count() or 1)
    parser.add_argument("--output", type=Path, default=Path(__file__).with_name(
        "hdiv_field_evaluator_scaling.json"))
    args = parser.parse_args()

    ng.SetNumThreads(args.threads)
    rad.UtiDelAll()
    from radia.vim import _radsolve
    _radsolve.clear_registry()
    mesh = MakeStructured3DMesh(
        hexes=True, nx=args.n, ny=args.n, nz=args.n,
        mapping=lambda x, y, z: (x-0.5, y-0.5, z-0.5))
    iron = vim.MeshSoftIron(mesh, mu_r=1000.0)
    source = rad.ObjBckg(lambda _point: [0.0, 0.0, 4.0e-7*np.pi*1.0e4])
    model = rad.ObjCnt([iron, source])
    solve_started = time.perf_counter()
    with ng.TaskManager():
        result = rad.Solve(model, 1.0e-6, 2000, 0)
    solve_wall = time.perf_counter()-solve_started
    evaluator = result["_field_evaluator"]

    rng = np.random.default_rng(20260714)
    pool = rng.uniform(-1.5, 1.5, (max(2*args.n_obs, 1000), 3))
    observations = np.ascontiguousarray(
        pool[np.any(np.abs(pool) > 0.55, axis=1)][:args.n_obs])
    if len(observations) != args.n_obs:
        raise RuntimeError("failed to generate the requested exterior observation count")
    reference = observations[np.linspace(
        0, len(observations)-1, min(args.reference_obs, len(observations)), dtype=int)]

    direct_ref_s, direct_ref = _best(lambda: evaluator.field(reference, "direct"))
    tree_ref_s, tree_ref = _best(lambda: evaluator.field(reference, "tree"))
    reference_scale = max(np.linalg.norm(direct_ref, axis=1).max(), 1e-300)
    tree_error = np.linalg.norm(tree_ref-direct_ref, axis=1)
    h_wall_s, h_value = _best(lambda: rad.Fld(iron, "h", observations))
    selected = evaluator.last_algorithm()
    b_wall_s, b_value = _best(lambda: rad.Fld(iron, "b", observations))

    payload = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "validation_host": platform.node(),
        "git_sha": _git_sha(),
        "radia_version": rad.__version__,
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "threads": args.threads,
        "problem": {
            "element": "structured HEX RT1",
            "n_per_axis": args.n,
            "n_elements": args.n**3,
            "ndof": result["ndof"],
            "n_observations": len(observations),
            "reference_observations": len(reference),
            "source_observation_work": int(
                result["field_evaluator_stats"]["source_count"]*len(observations)),
        },
        "solve_wall_s": solve_wall,
        "field_evaluator": result["field_evaluator_stats"],
        "field": {
            "selected_algorithm_h": selected,
            "h_public_wall_s": h_wall_s,
            "b_public_wall_s": b_wall_s,
            "direct_reference_wall_s": direct_ref_s,
            "tree_reference_wall_s": tree_ref_s,
            "tree_reference_speedup": direct_ref_s/max(tree_ref_s, 1e-300),
            "tree_error_max_over_direct_scale": float(tree_error.max()/reference_scale),
            "tree_error_p95_over_direct_scale": float(
                np.quantile(tree_error, 0.95)/reference_scale),
            "h_norm": float(np.linalg.norm(h_value)),
            "b_norm": float(np.linalg.norm(b_value)),
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    print("saved:", args.output)
    rad.UtiDelAll()
    _radsolve.clear_registry()


if __name__ == "__main__":
    main()
