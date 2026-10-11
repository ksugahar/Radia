"""Assemble checked radia-optuna release evidence from raw benchmark JSON."""

from __future__ import annotations

import argparse
import hashlib
import json
import runpy
from pathlib import Path
from typing import Any


def _read(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise TypeError(f"{path} must contain a JSON object")
    return value


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def _require_private_same_host(python_result, matlab_result):
    # Same-host timing certification requires private runtime provenance.
    # Descriptive public OS metadata cannot establish a physical host identity.
    _require(
        bool(str(python_result.get("host", "")).strip())
        and bool(str(matlab_result.get("host", "")).strip()),
        "same-host certification requires private runtime provenance",
    )
    _require(
        str(python_result.get("host", "")).casefold()
        == str(matlab_result.get("host", "")).casefold(),
        "Python and MATLAB performance hosts differ",
    )


def _public_runtime_payload(value):
    """Keep private execution identity out of the assembled public record."""
    identity_keys = {"host", "hostname", "host_name", "validation_host", "machine_name", "computer_name"}
    if isinstance(value, dict):
        return {key: _public_runtime_payload(item) for key, item in value.items() if key not in identity_keys}
    if isinstance(value, list):
        return [_public_runtime_payload(item) for item in value]
    return value


def _validate_public_runtime_payload(value):
    # Use the same structural JSON/UNC/contextual-host detector as CI. This is
    # publication validation, independent of operational input/host admission.
    policy = runpy.run_path(str(Path(__file__).resolve().parents[2] / "tools/policy_lint.py"))
    hits = policy["public_runtime_identity_hits"]("result.json", json.dumps(value))
    _require(not hits, "public evidence contains identifying runtime metadata")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--installed", type=Path, required=True)
    parser.add_argument("--python-performance", type=Path, required=True)
    parser.add_argument("--matlab-performance", type=Path, required=True)
    parser.add_argument("--mex-cold", type=Path, required=True)
    parser.add_argument("--parallel", type=Path, required=True)
    parser.add_argument("--history", type=Path, required=True)
    parser.add_argument("--matlab-test-total", type=int, required=True)
    parser.add_argument("--matlab-test-passed", type=int, required=True)
    parser.add_argument("--environment-note", default="")
    parser.add_argument("--pre-run-cpu-percent", type=float, nargs="*", default=[])
    parser.add_argument("--free-memory-mib", type=float)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    paths = {
        "installed": args.installed,
        "python_performance": args.python_performance,
        "matlab_performance": args.matlab_performance,
        "mex_cold": args.mex_cold,
        "parallel": args.parallel,
        "history": args.history,
    }
    values = {name: _read(path) for name, path in paths.items()}
    installed = values["installed"]
    python_result = values["python_performance"]
    matlab_result = values["matlab_performance"]
    cold = values["mex_cold"]
    parallel = values["parallel"]
    history = values["history"]

    _require(installed.get("ok") is True, "installed-wheel evidence failed")
    _require(args.matlab_test_total > 0, "MATLAB test total must be positive")
    _require(
        args.matlab_test_passed == args.matlab_test_total,
        "not every MATLAB test passed",
    )
    _require(
        python_result.get("runtime") == "python-upstream",
        "Python performance runtime is invalid",
    )
    _require(
        matlab_result.get("runtime") == "matlab",
        "MATLAB performance runtime is invalid",
    )
    _require_private_same_host(python_result, matlab_result)
    _require(parallel.get("gate", {}).get("passed") is True, "parallel gate failed")
    history_rows = history.get("results", [])
    _require(bool(history_rows), "history benchmark has no rows")
    _require(
        all(row.get("converged") is True for row in history_rows),
        "history lookup results differ",
    )

    ratios: dict[str, float] = {}
    for workload in ("scalar", "grouped_conditional", "trials_dataframe"):
        python_seconds = float(python_result[workload]["median_warmed_seconds"])
        matlab_seconds = float(matlab_result[workload]["median_warmed_seconds"])
        ratios[workload] = matlab_seconds / python_seconds
    seeded_checksums_match = all(
        abs(
            float(python_result[workload]["checksum"])
            - float(matlab_result[workload]["checksum"])
        )
        <= 1.0e-12
        for workload in ("scalar", "grouped_conditional")
    )
    dataframe_shape_match = all(
        python_result["trials_dataframe"][key]
        == matlab_result["trials_dataframe"][key]
        for key in ("rows", "columns")
    )
    largest_history = max(history_rows, key=lambda row: int(row["trials"]))
    runtime_gate = (
        all(ratio <= 1.0 for ratio in ratios.values())
        and seeded_checksums_match
        and dataframe_shape_match
    )

    evidence = dict(installed)
    evidence["schema"] = "radia-optuna.release-evidence.v1"
    evidence["matlab_tests"] = {
        "schema": "radia.optuna.oracle-test-summary.v1",
        "total": args.matlab_test_total,
        "passed": args.matlab_test_passed,
        "failed": args.matlab_test_total - args.matlab_test_passed,
        "incomplete": 0,
    }
    evidence["performance"] = {
        "python": python_result,
        "matlab": matlab_result,
        "mex_cold": cold,
    }
    evidence["matlab_extensions"] = {
        "parallel_batch": parallel,
        "history_store": history,
    }
    evidence["benchmark_environment"] = {
        "note": args.environment_note,
        "pre_run_total_cpu_percent": args.pre_run_cpu_percent,
        "free_memory_mib": args.free_memory_mib,
    }
    evidence["summary"] = {
        "matlab_to_python_warmed_time_ratios": ratios,
        "seeded_checksums_match": seeded_checksums_match,
        "trials_dataframe_shape_match": dataframe_shape_match,
        "parallel_speedup": parallel["comparison"]["speedup"],
        "parallel_worker_efficiency": parallel["comparison"][
            "worker_efficiency"
        ],
        "history_freeze_exponent": history["freeze_exponent"],
        "largest_history_trials": largest_history["trials"],
        "largest_history_lookup_speedup": largest_history["lookup_speedup"],
        "mex_first_call_median_seconds": cold["median_first_call_seconds"],
        "runtime_gate_passed": runtime_gate,
        "parallel_gate_passed": parallel["gate"]["passed"],
        "history_gate_passed": all(
            row.get("converged") is True for row in history_rows
        ),
    }
    evidence["summary"]["passed"] = all(
        (
            evidence["summary"]["runtime_gate_passed"],
            evidence["summary"]["parallel_gate_passed"],
            evidence["summary"]["history_gate_passed"],
        )
    )
    evidence["source_sha256"] = {
        name: _sha256(path) for name, path in paths.items()
    }
    evidence["ok"] = evidence["summary"]["passed"]

    evidence = _public_runtime_payload(evidence)
    _validate_public_runtime_payload(evidence)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(_public_runtime_payload(evidence), indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(args.output)


if __name__ == "__main__":
    main()
