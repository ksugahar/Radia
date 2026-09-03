#!/usr/bin/env python3
"""Record Radia field-evaluation performance as validation evidence."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import platform
import socket
import sys
import time

import radia as rad


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = ROOT / "validation_test" / "benchmarks" / "results" / "field_parallel.json"


def create_test_magnet():
    magnet = rad.magnet_box([0, 0, 0], [0.01, 0.01, 0.01], [0, 0, 0])
    rad.ObjSetM(magnet, [0, 0, 1000])
    return magnet


def field_computation_1d(magnet, num_points: int) -> tuple[float, int]:
    started = time.perf_counter()
    for index in range(num_points):
        z = 0.1 * index / (num_points - 1)
        rad.Fld(magnet, "b", [0, 0, z])
    return time.perf_counter() - started, num_points


def field_computation_2d(magnet, grid_size: int) -> tuple[float, int]:
    points = [
        [-0.025 + 0.050 * i / (grid_size - 1),
         -0.025 + 0.050 * j / (grid_size - 1), 0.020]
        for i in range(grid_size)
        for j in range(grid_size)
    ]
    started = time.perf_counter()
    for point in points:
        rad.Fld(magnet, "b", point)
    return time.perf_counter() - started, len(points)


def run_case(name: str, function, *args) -> dict[str, object]:
    seconds, points = function(*args)
    if not math.isfinite(seconds) or seconds <= 0:
        raise RuntimeError(f"{name} produced invalid elapsed time: {seconds}")
    return {
        "name": name,
        "points": points,
        "elapsed_seconds": seconds,
        "points_per_second": points / seconds,
    }


def run_benchmark() -> dict[str, object]:
    rad.UtiDelAll()
    try:
        magnet = create_test_magnet()
        cases = [
            run_case("line_1000", field_computation_1d, magnet, 1000),
            run_case("grid_50x50", field_computation_2d, magnet, 50),
            run_case("line_5000", field_computation_1d, magnet, 5000),
            run_case("line_10000", field_computation_1d, magnet, 10000),
        ]
        return {
            "schema": "radia.validation.benchmark-field-parallel.v1",
            "generated_at_utc": datetime.now(timezone.utc).isoformat(),
            "host": socket.gethostname(),
            "platform": platform.platform(),
            "python_version": platform.python_version(),
            "radia_version": str(rad.UtiVer()),
            "thread_environment": {
                name: os.environ.get(name)
                for name in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "NGS_NUM_THREADS")
            },
            "acceptance": {
                "rule": "all timings are finite and positive",
                "status": "pass",
            },
            "cases": cases,
        }
    finally:
        rad.UtiDelAll()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)

    result = run_benchmark()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    for case in result["cases"]:
        print(
            f"{case['name']}: {case['elapsed_seconds']:.6f}s "
            f"({case['points_per_second']:.1f} points/s)"
        )
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
