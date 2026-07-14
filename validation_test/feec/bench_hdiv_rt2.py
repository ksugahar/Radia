"""Measure the public pure-TET RT1/RT2 material path on mdx or hibino.

This is a validation benchmark, not a pytest test.  It uses the same mesh for
both orders at each size and records external wall time, solver stage timings,
H-matrix statistics, iterations, and process memory.  Publication timing must
come from an idle compute host.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import platform
import time

import ngsolve as ng
from netgen.occ import Box, OCCGeometry, Pnt
import numpy as np

import radia
from radia.vim import Solve


def _require_compute_host() -> str:
    node = platform.node()
    short = node.lower().split(".", 1)[0]
    if not (
        short == "mdx"
        or short.startswith("mdx-")
        or short == "hibino"
        or short.startswith("hibino-")
    ):
        raise SystemExit("bench_hdiv_rt2.py must run on mdx or hibino; got hostname %r" % node)
    return node


def _memory() -> dict:
    try:
        import psutil

        info = psutil.Process().memory_info()
        return {
            "rss_bytes": int(info.rss),
            "peak_wset_bytes": int(getattr(info, "peak_wset", info.rss)),
        }
    except ImportError:
        return {"rss_bytes": None, "peak_wset_bytes": None}


def _mesh(maxh: float):
    return ng.Mesh(OCCGeometry(Box(Pnt(0, 0, 0), Pnt(1, 1, 1))).GenerateMesh(maxh=maxh))


def _run(mesh, order: int, mu_r: float) -> dict:
    before = _memory()
    t0 = time.perf_counter()
    with ng.TaskManager():
        result = Solve(
            mesh,
            mu_r=mu_r,
            H_ext=ng.CF((0, 0, 1000.0)),
            order=order,
            gram_eps=1e-10,
        )
    wall = time.perf_counter() - t0
    after = _memory()
    return {
        "order": order,
        "ndof": int(result["ndof"]),
        "n_charge": int(result["n_charge"]),
        "demag": float(result["demag"]),
        "M_avg": np.asarray(result["M_avg"], float).tolist(),
        "iters": int(result["iters"]),
        "wall_s_external": wall,
        "total_wall_s_internal": float(result["total_wall_s_internal"]),
        "charge_gram_wall_s": float(result["charge_gram_wall_s"]),
        "solve_wall_s": float(result["solve_wall_s"]),
        "memory_before": before,
        "memory_after": after,
        "hmat_stats": result.get("hmat_stats"),
        "cpp_solve_timings": result.get("cpp_solve_timings"),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--maxh", default="0.6,0.4,0.3")
    parser.add_argument("--mu-r", type=float, default=100.0)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).with_name("hdiv_rt2_scaling.json"),
    )
    args = parser.parse_args()

    host = _require_compute_host()
    cases = []
    for maxh in (float(value) for value in args.maxh.split(",")):
        t0 = time.perf_counter()
        mesh = _mesh(maxh)
        mesh_wall = time.perf_counter() - t0
        rows = [_run(mesh, order, args.mu_r) for order in (1, 2)]
        cases.append(
            {
                "maxh": maxh,
                "n_elements": int(mesh.ne),
                "mesh_wall_s": mesh_wall,
                "orders": rows,
                "rt2_over_rt1_wall": rows[1]["wall_s_external"] / rows[0]["wall_s_external"],
            }
        )

    payload = {
        "schema": "radia.hdiv_rt2_scaling.v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "validation_host": host,
        "platform": platform.platform(),
        "python_version": platform.python_version(),
        "radia_version": radia.__version__,
        "mu_r": args.mu_r,
        "cases": cases,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
