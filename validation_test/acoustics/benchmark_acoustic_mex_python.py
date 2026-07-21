"""Benchmark the same acoustic C++ kernel through pybind11 and MATLAB MEX."""

from __future__ import annotations

import argparse
import importlib
import json
import os
import platform
import statistics
import sys
import time
from pathlib import Path

import numpy as np


def points(count: int) -> np.ndarray:
    index = np.arange(count, dtype=np.int64)
    return np.column_stack(
        (
            1.2 + np.mod(index, 997) / 997.0,
            (np.mod(index * 17, 991) - 495) / 4000.0,
            (np.mod(index * 31, 983) - 491) / 3500.0,
        )
    )


def timed_samples(function, repeats: int) -> list[float]:
    samples = []
    for _ in range(repeats):
        start = time.perf_counter()
        function()
        samples.append(time.perf_counter() - start)
    return samples


def run_python(args) -> dict:
    setup_start = time.perf_counter()
    sys.path.insert(0, str(args.module_root.resolve()))
    import netgen
    import radia as installed_radia

    netgen_root = Path(netgen.__file__).resolve().parent
    dll_candidates = [
        netgen_root,
        Path(sys.prefix) / "Library" / "bin",
        Path(sys.prefix),
    ]
    dll_handles = []
    if os.name == "nt":
        for directory in dll_candidates:
            if directory.is_dir():
                dll_handles.append(os.add_dll_directory(str(directory)))
    native = importlib.import_module("radia._radia_pybind")
    setup_s = time.perf_counter() - setup_start

    point_values = np.ascontiguousarray(points(args.points))
    zeta = np.ascontiguousarray(
        np.linspace(-0.8, 0.8, args.points)
        + 1j * np.linspace(0.4, -0.4, args.points)
    )
    scattering = lambda: native._AcousticSoftSphere(
        3.1, 1.0, point_values, args.terms
    )
    transfer = lambda: native._AcousticBDFDelta(zeta, "BDF2")

    first_start = time.perf_counter()
    first = scattering()
    first_s = time.perf_counter() - first_start
    for _ in range(args.warmup):
        scattering()
        transfer()
    scattering_samples = timed_samples(scattering, args.repeats)
    transfer_samples = timed_samples(transfer, args.repeats)
    scattered = np.asarray(first["scattered"])
    return {
        "schema": "radia.acoustic-backend-benchmark/v1",
        "backend": "python-pybind11",
        "host": platform.node(),
        "python_version": platform.python_version(),
        "radia_version": installed_radia.__version__,
        "platform": platform.platform(),
        "point_count": args.points,
        "terms": args.terms,
        "warmup": args.warmup,
        "repeats": args.repeats,
        "setup_s": setup_s,
        "first_scattering_s": first_s,
        "scattering_median_s": statistics.median(scattering_samples),
        "scattering_min_s": min(scattering_samples),
        "bdf_transfer_median_s": statistics.median(transfer_samples),
        "bdf_transfer_min_s": min(transfer_samples),
        "checksum_real": float(scattered.real.sum()),
        "checksum_imag": float(scattered.imag.sum()),
    }


def compare(args) -> dict:
    python = json.loads(args.python_json.read_text(encoding="utf-8"))
    matlab = json.loads(args.matlab_json.read_text(encoding="utf-8"))
    python_checksum = complex(python["checksum_real"], python["checksum_imag"])
    matlab_checksum = complex(matlab["checksum_real"], matlab["checksum_imag"])
    checksum_error = abs(python_checksum - matlab_checksum) / max(
        abs(python_checksum), 1e-300
    )
    result = {
        "schema": "radia.acoustic-python-matlab-comparison/v1",
        "generated_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "host": python["host"],
        "point_count": python["point_count"],
        "terms": python["terms"],
        "warmup": python["warmup"],
        "repeats": python["repeats"],
        "checksum_relative_error": checksum_error,
        "python": python,
        "matlab": matlab,
        "ratios": {
            "mex_over_pybind_scattering": (
                matlab["scattering_median_s"] / python["scattering_median_s"]
            ),
            "mex_over_pybind_bdf_transfer": (
                matlab["bdf_transfer_median_s"]
                / python["bdf_transfer_median_s"]
            ),
        },
    }
    if args.python_e2e is not None:
        result["python"]["process_end_to_end_s"] = args.python_e2e
    if args.matlab_e2e is not None:
        result["matlab"]["process_end_to_end_s"] = args.matlab_e2e
    return result


def parse_args():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="mode", required=True)
    python_parser = subparsers.add_parser("python")
    python_parser.add_argument("--module-root", type=Path, required=True)
    python_parser.add_argument("--output", type=Path, required=True)
    python_parser.add_argument("--points", type=int, default=20000)
    python_parser.add_argument("--terms", type=int, default=28)
    python_parser.add_argument("--warmup", type=int, default=5)
    python_parser.add_argument("--repeats", type=int, default=31)

    compare_parser = subparsers.add_parser("compare")
    compare_parser.add_argument("--python-json", type=Path, required=True)
    compare_parser.add_argument("--matlab-json", type=Path, required=True)
    compare_parser.add_argument("--output", type=Path, required=True)
    compare_parser.add_argument("--python-e2e", type=float)
    compare_parser.add_argument("--matlab-e2e", type=float)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    result = run_python(args) if args.mode == "python" else compare(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
