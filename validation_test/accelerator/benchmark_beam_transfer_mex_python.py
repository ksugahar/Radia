"""Benchmark and compare the nonlinear transfer-map pybind11/MEX backends.

The pybind lane intentionally imports only a native module built in this
checkout's ``src/radia`` directory.  It never falls back to an installed wheel.
The comparison lane is importable without a native build so that its input and
provenance contracts can be tested before an expensive paired run.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import platform
import statistics
import subprocess
import sys
import time
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CASE_SCHEMA = "radia.beam-transfer-benchmark-case/v1"
RESULT_SCHEMA = "radia.beam-transfer-backend-benchmark/v1"
COMPARISON_SCHEMA = "radia.beam-transfer-pybind-mex-comparison/v1"
DEFAULT_OUTPUT_DIR = Path(r"C:\temp\radia-validation\beam_transfer")


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _git_head(root: Path = ROOT) -> str:
    return subprocess.check_output(
        ["git", "-C", str(root), "rev-parse", "HEAD"], text=True
    ).strip()


def _require_clean_source(root: Path = ROOT) -> None:
    for arguments in (("diff", "--quiet", "HEAD", "--"),
                      ("diff", "--cached", "--quiet", "HEAD", "--")):
        if subprocess.run(["git", "-C", str(root), *arguments], check=False).returncode:
            raise RuntimeError("beam backend validation requires a clean source checkout")


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _validate_case(value: dict) -> None:
    _require(value.get("schema") == CASE_SCHEMA, f"case schema must be {CASE_SCHEMA}")
    lengths = np.asarray(value.get("lengths_m"), dtype=float)
    _require(lengths.ndim == 1 and lengths.size > 0, "lengths_m must be a nonempty vector")
    _require(np.isfinite(lengths).all() and np.all(lengths > 0), "lengths_m must be finite and positive")
    count = lengths.size
    a = np.asarray(value.get("A_per_m"), dtype=float)
    _require(a.shape == (count, 6, 6), f"A_per_m must have shape ({count}, 6, 6)")
    _require(np.isfinite(a).all(), "A_per_m must contain only finite values")
    names = value.get("names")
    _require(isinstance(names, list) and len(names) == count, "names must match lengths_m")
    _require(all(isinstance(name, str) and name for name in names), "names must be nonempty strings")
    _require(value.get("maximum_order") in (1, 2, 3), "maximum_order must be 1, 2, or 3")
    maximum_step = value.get("maximum_step_m")
    _require(isinstance(maximum_step, (int, float)) and np.isfinite(maximum_step) and maximum_step > 0,
             "maximum_step_m must be finite and positive")
    for key, width, axes in (("F2_entries", 5, (count, 6, 6, 6)),
                             ("F3_entries", 6, (count, 6, 6, 6, 6))):
        entries = value.get(key)
        _require(isinstance(entries, list), f"{key} must be a list")
        for row_number, entry in enumerate(entries):
            _require(isinstance(entry, list) and len(entry) == width,
                     f"{key}[{row_number}] must contain {width} values")
            indices = entry[:-1]
            _require(all(isinstance(index, int) and not isinstance(index, bool) for index in indices),
                     f"{key}[{row_number}] indices must be integers")
            _require(all(0 <= index < limit for index, limit in zip(indices, axes)),
                     f"{key}[{row_number}] index is out of range")
            coefficient = entry[-1]
            _require(isinstance(coefficient, (int, float)) and np.isfinite(coefficient),
                     f"{key}[{row_number}] coefficient must be finite")


def _case(path: Path):
    raw = path.read_bytes()
    value = json.loads(raw)
    _validate_case(value)
    count = len(value["lengths_m"])
    f2 = np.zeros((count, 6, 6, 6), dtype=float)
    f3 = np.zeros((count, 6, 6, 6, 6), dtype=float)
    for segment, row, column, depth, entry in value["F2_entries"]:
        f2[segment, row, column, depth] = entry
    for segment, row, column, depth, degree, entry in value["F3_entries"]:
        f3[segment, row, column, depth, degree] = entry
    return raw, value, np.asarray(value["lengths_m"], dtype=float), np.asarray(
        value["A_per_m"], dtype=float
    ), f2, f3


def _load_source_backend(source_root: Path = ROOT):
    package_dir = (source_root / "src" / "radia").resolve()
    binaries = sorted(package_dir.glob("_radia_pybind*.pyd")) + sorted(
        package_dir.glob("_radia_pybind*.so")
    )
    if not binaries:
        raise RuntimeError(
            "matching source pybind artifact is unavailable: build _radia_pybind "
            f"in {package_dir}; installed-wheel fallback is forbidden"
        )
    sys.path.insert(0, str(source_root / "src"))
    radia = importlib.import_module("radia")
    native = importlib.import_module("radia._radia_pybind")
    beam = importlib.import_module("radia.beam")
    native_path = Path(native.__file__).resolve()
    beam_path = Path(beam.__file__).resolve()
    if native_path.parent != package_dir or package_dir not in beam_path.parents:
        raise RuntimeError(
            "beam benchmark resolved outside the selected source tree; "
            "installed-wheel fallback is forbidden"
        )
    return radia, native, beam.propagate_variational_map


def _observables(result):
    return {
        "R_fro": float(np.linalg.norm(result["R"])),
        "T_fro": float(np.linalg.norm(result["T"])),
        "U_fro": float(np.linalg.norm(result["U"])),
        "R_11": float(result["R"][0, 0]),
        "T_211": float(result["T"][1, 0, 0]),
        "U_3111": float(result["U"][2, 0, 0, 0]),
        "U_4111": float(result["U"][3, 0, 0, 0]),
        "R_composition_error": float(result["diagnostics"]["R_composition_error"]),
        "T_reconstruction_error": float(result["diagnostics"]["T_reconstruction_error"]),
        "U_reconstruction_error": float(result["diagnostics"]["U_reconstruction_error"]),
    }


def benchmark(case_path: Path, repeats: int) -> dict:
    _require_clean_source()
    source_commit = _git_head()
    radia, native, propagate = _load_source_backend()
    case = _case(case_path)

    def run():
        _, value, lengths, a, f2, f3 = case
        return propagate(
            lengths, a, f2, f3, value["names"],
            maximum_order=value["maximum_order"],
            maximum_step_m=value["maximum_step_m"],
        )

    start = time.perf_counter()
    first = run()
    first_s = time.perf_counter() - start
    samples = []
    for _ in range(repeats):
        start = time.perf_counter()
        run()
        samples.append(time.perf_counter() - start)
    binary_path = Path(native.__file__).resolve()
    return {
        "schema": RESULT_SCHEMA,
        "executed_at_utc": datetime.now(UTC).isoformat(),
        "backend": "python-pybind11",
        "case_id": case[1]["case_id"],
        "case_sha256": hashlib.sha256(case[0]).hexdigest(),
        "radia_version": radia.__version__,
        "radia_git_head": source_commit,
        "python_version": platform.python_version(),
        "host": platform.node(),
        "binary": {
            "path": str(binary_path),
            "bytes": binary_path.stat().st_size,
            "sha256": _sha256(binary_path),
            "source_commit": source_commit,
        },
        "repeats": repeats,
        "first_s": first_s,
        "median_s": statistics.median(samples),
        "min_s": min(samples),
        "observables": _observables(first),
    }


def _is_hex(value: object, length: int) -> bool:
    return isinstance(value, str) and len(value) == length and all(
        character in "0123456789abcdefABCDEF" for character in value
    )


def validate_backend_report(report: dict, expected_backend: str) -> None:
    _require(report.get("schema") == RESULT_SCHEMA, "unexpected benchmark result schema")
    _require(report.get("backend") == expected_backend, f"backend must be {expected_backend}")
    _require(_is_hex(report.get("case_sha256"), 64), "case_sha256 must be SHA-256 hex")
    _require(_is_hex(report.get("radia_git_head"), 40), "radia_git_head must be a full Git SHA")
    binary = report.get("binary")
    _require(isinstance(binary, dict), "binary provenance is required")
    _require(isinstance(binary.get("path"), str) and binary["path"], "binary.path is required")
    _require(isinstance(binary.get("bytes"), int) and binary["bytes"] > 0, "binary.bytes must be positive")
    _require(_is_hex(binary.get("sha256"), 64), "binary.sha256 must be SHA-256 hex")
    _require(binary.get("source_commit") == report["radia_git_head"],
             "binary.source_commit must match radia_git_head")
    _require(isinstance(report.get("median_s"), (int, float)) and report["median_s"] > 0,
             "median_s must be positive")
    _require(isinstance(report.get("observables"), dict) and report["observables"],
             "observables are required")


def compare(python_path: Path, matlab_path: Path, *, rtol: float, atol: float) -> dict:
    python = json.loads(python_path.read_text(encoding="utf-8"))
    matlab = json.loads(matlab_path.read_text(encoding="utf-8"))
    validate_backend_report(python, "python-pybind11")
    validate_backend_report(matlab, "matlab-mex")
    required = {
        "schema": RESULT_SCHEMA,
        "case_id": python["case_id"],
        "case_sha256": python["case_sha256"],
        "radia_git_head": python["radia_git_head"],
    }
    checks = {}
    _require(set(matlab["observables"]) == set(python["observables"]),
             "MATLAB/Python observable keys differ")
    for key, expected in python["observables"].items():
        actual = float(matlab["observables"][key])
        relative = abs(actual - expected) / max(abs(expected), atol)
        checks[key] = {
            "python": float(expected),
            "matlab": actual,
            "relative_error": relative,
            "pass": bool(abs(actual - expected) <= atol + rtol * abs(expected)),
        }
    metadata = {
        key: {"python": python.get(key), "matlab": matlab.get(key),
              "pass": matlab.get(key) == expected}
        for key, expected in required.items()
    }
    result = {
        "schema": COMPARISON_SCHEMA,
        "executed_at_utc": datetime.now(UTC).isoformat(),
        "case_id": python["case_id"],
        "case_sha256": python["case_sha256"],
        "tolerance": {"rtol": rtol, "atol": atol},
        "metadata": metadata,
        "checks": checks,
        "binary_provenance": {"python": python["binary"], "matlab": matlab["binary"]},
        "ratios": {"mex_over_pybind": float(matlab["median_s"] / python["median_s"])},
        "pass": all(item["pass"] for item in metadata.values())
        and all(item["pass"] for item in checks.values()),
    }
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="mode", required=True)
    pybind = subparsers.add_parser("pybind")
    pybind.add_argument("--case", type=Path, default=HERE / "beam_transfer_benchmark_case.json")
    pybind.add_argument("--output", type=Path, default=DEFAULT_OUTPUT_DIR / "pybind.json")
    pybind.add_argument("--repeats", type=int, default=31)
    comparison = subparsers.add_parser("compare")
    comparison.add_argument("--python-json", type=Path, required=True)
    comparison.add_argument("--matlab-json", type=Path, required=True)
    comparison.add_argument("--output", type=Path, default=DEFAULT_OUTPUT_DIR / "comparison.json")
    comparison.add_argument("--rtol", type=float, default=1e-11)
    comparison.add_argument("--atol", type=float, default=1e-12)
    args = parser.parse_args()
    if args.mode == "pybind":
        if args.repeats < 1:
            parser.error("--repeats must be positive")
        result = benchmark(args.case, args.repeats)
    else:
        if args.rtol < 0 or args.atol < 0:
            parser.error("--rtol and --atol must be nonnegative")
        result = compare(args.python_json, args.matlab_json, rtol=args.rtol, atol=args.atol)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0 if result.get("pass", True) else 1


if __name__ == "__main__":
    raise SystemExit(main())
