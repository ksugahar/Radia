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
PROVENANCE_SCHEMA = "radia.native-build-provenance.v1"
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


def _finite_number(value: object) -> bool:
    return (
        isinstance(value, (int, float))
        and not isinstance(value, bool)
        and bool(np.isfinite(value))
    )


def _reject_json_constant(value: str):
    raise ValueError(f"non-finite JSON constant is forbidden: {value}")


def _load_json_bytes(raw: bytes):
    return json.loads(raw, parse_constant=_reject_json_constant)


def _load_json(path: Path):
    return _load_json_bytes(path.read_bytes())


def _json_text(value: object) -> str:
    return json.dumps(value, indent=2, allow_nan=False) + "\n"


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
            _require(_finite_number(coefficient),
                     f"{key}[{row_number}] coefficient must be finite")


def _case(path: Path):
    raw = path.read_bytes()
    value = _load_json_bytes(raw)
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


def _load_build_provenance(binary_path: Path, source_root: Path) -> dict:
    manifest_path = Path(f"{binary_path}.build.json")
    if not manifest_path.is_file():
        raise RuntimeError(
            "native build provenance manifest is missing: "
            f"{manifest_path}; rebuild the selected source checkout"
        )
    manifest = _load_json(manifest_path)
    _require(manifest.get("schema") == PROVENANCE_SCHEMA,
             f"native provenance schema must be {PROVENANCE_SCHEMA}")
    _require(manifest.get("source_dirty") is False,
             "native provenance must come from a clean source build")
    _require(_is_hex(manifest.get("source_commit"), 40),
             "native provenance source_commit must be a full Git SHA")
    _require(manifest["source_commit"] == _git_head(source_root),
             "native provenance source_commit does not match the selected checkout")
    _require(manifest.get("binary_name") == binary_path.name,
             "native provenance binary_name does not match the loaded binary")
    _require(manifest.get("binary_bytes") == binary_path.stat().st_size,
             "native provenance binary_bytes does not match the loaded binary")
    _require(manifest.get("binary_sha256") == _sha256(binary_path),
             "native provenance binary_sha256 does not match the loaded binary")
    return manifest


def _load_source_backend(source_root: Path = ROOT):
    package_dir = (source_root / "src" / "radia").resolve()
    native_path = package_dir / "_radia_pybind.pyd"
    if not native_path.is_file():
        raise RuntimeError(
            "matching source pybind artifact is unavailable: build _radia_pybind "
            f"in {package_dir}; installed-wheel fallback is forbidden"
        )
    provenance = _load_build_provenance(native_path, source_root)
    sys.path.insert(0, str(source_root / "src"))
    radia = importlib.import_module("radia")
    native = importlib.import_module("radia._radia_pybind")
    beam = importlib.import_module("radia.beam")
    loaded_native_path = Path(native.__file__).resolve()
    beam_path = Path(beam.__file__).resolve()
    if loaded_native_path != native_path or package_dir not in beam_path.parents:
        raise RuntimeError(
            "beam benchmark resolved outside the selected source tree; "
            "installed-wheel fallback is forbidden"
        )
    return radia, native, beam.propagate_variational_map, provenance


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
    if not isinstance(repeats, int) or isinstance(repeats, bool) or repeats < 1:
        raise ValueError("repeats must be a positive integer")
    _require_clean_source()
    radia, native, propagate, provenance = _load_source_backend()
    source_commit = provenance["source_commit"]
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
    _require(_finite_number(first_s) and first_s >= 0,
             "first_s must be finite and nonnegative")
    _require(all(_finite_number(sample) and sample > 0 for sample in samples),
             "benchmark samples must be finite and positive")
    measured_observables = _observables(first)
    _require(all(_finite_number(value) for value in measured_observables.values()),
             "benchmark observables must be finite")
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
            "build_manifest": str(Path(f"{binary_path}.build.json")),
        },
        "repeats": repeats,
        "first_s": first_s,
        "median_s": statistics.median(samples),
        "min_s": min(samples),
        "observables": measured_observables,
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
    _require(isinstance(binary.get("build_manifest"), str) and binary["build_manifest"],
             "binary.build_manifest is required")
    _require(_finite_number(report.get("median_s")) and report["median_s"] > 0,
             "median_s must be positive")
    _require(_finite_number(report.get("first_s")) and report["first_s"] >= 0,
             "first_s must be finite and nonnegative")
    _require(_finite_number(report.get("min_s")) and report["min_s"] > 0,
             "min_s must be finite and positive")
    repeats = report.get("repeats")
    _require(isinstance(repeats, int) and not isinstance(repeats, bool) and repeats > 0,
             "repeats must be a positive integer")
    _require(isinstance(report.get("observables"), dict) and report["observables"],
             "observables are required")
    _require(all(_finite_number(value)
                 for value in report["observables"].values()),
             "observables must contain only finite numeric values")


def compare(python_path: Path, matlab_path: Path, *, rtol: float, atol: float) -> dict:
    _require(_finite_number(rtol) and rtol >= 0,
             "rtol must be finite and nonnegative")
    _require(_finite_number(atol) and atol >= 0,
             "atol must be finite and nonnegative")
    python = _load_json(python_path)
    matlab = _load_json(matlab_path)
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
        difference = abs(actual - expected)
        denominator = max(abs(expected), atol)
        relative = 0.0 if difference == 0 else (
            None if denominator == 0 else difference / denominator
        )
        checks[key] = {
            "python": float(expected),
            "matlab": actual,
            "relative_error": relative,
            "pass": bool(difference <= atol + rtol * abs(expected)),
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
    args.output.write_text(_json_text(result), encoding="utf-8")
    print(_json_text(result), end="")
    return 0 if result.get("pass", True) else 1


if __name__ == "__main__":
    raise SystemExit(main())
