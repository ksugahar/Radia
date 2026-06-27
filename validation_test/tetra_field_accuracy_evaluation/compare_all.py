#!/usr/bin/env python
"""Summarize tetrahedral field-accuracy validation JSON files.

The historical compare script assumed an older ``evaluation_results.json``
schema.  The validation corpus now has several lanes:

* analytical_reference_results.json: uniform-M tetra MSC vs hexa reference.
* solver_comparison_results.json: Radia tetra solver vs Radia hexa solver.
* evaluation_results.json: older NGSolve H-formulation extraction attempt
  retained as a known issue / diagnostic lane.
* ngsolve_reference_results.json: independent A-formulation reference, which may
  produce NaNs in current local setups and is therefore reported, not asserted.
"""

from __future__ import annotations

import datetime as _dt
import importlib.metadata as _metadata
import json
import math
import platform
import sys
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def _version(name: str) -> str | None:
    try:
        return _metadata.version(name)
    except Exception:
        return None


def _load(name: str) -> dict | None:
    path = HERE / name
    if not path.exists():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def _finite_vector_count(values: list[list[float]]) -> int:
    return sum(
        all(math.isfinite(float(component)) for component in row)
        for row in values
    )


def _error_stats(errors: list[float]) -> dict:
    arr = np.array([float(v) for v in errors if math.isfinite(float(v))])
    if arr.size == 0:
        return {"count": 0, "avg_percent": None, "max_percent": None}
    return {
        "count": int(arr.size),
        "avg_percent": float(arr.mean()),
        "max_percent": float(arr.max()),
        "min_percent": float(arr.min()),
        "std_percent": float(arr.std()),
    }


def build_summary() -> dict:
    analytical = _load("analytical_reference_results.json")
    solver = _load("solver_comparison_results.json")
    extraction = _load("evaluation_results.json")
    ngsolve_ref = _load("ngsolve_reference_results.json")

    analytical_stats = _error_stats(analytical.get("errors", []) if analytical else [])
    solver_stats = _error_stats(solver.get("errors", []) if solver else [])
    extraction_stats = _error_stats(
        [row.get("error_percent", math.nan) for row in extraction.get("results", [])]
        if extraction else []
    )
    ngsolve_values = (
        ngsolve_ref.get("ngsolve", {}).get("B_values", [])
        if ngsolve_ref else []
    )
    finite_ngsolve = _finite_vector_count(ngsolve_values)

    checks = {
        "analytical_reference_pass": (
            analytical_stats["count"] > 0
            and analytical_stats["avg_percent"] is not None
            and analytical_stats["avg_percent"] < 0.01
        ),
        "radia_solver_comparison_pass": (
            solver_stats["count"] > 0
            and solver_stats["avg_percent"] is not None
            and solver_stats["avg_percent"] < 5.0
        ),
        "legacy_ngsolve_extraction_known_issue_recorded": (
            extraction_stats["count"] > 0
            and extraction_stats["avg_percent"] is not None
            and extraction_stats["avg_percent"] > 10.0
        ),
        "ngsolve_reference_finite_values_recorded": finite_ngsolve,
    }

    return {
        "schema": "radia.validation.tetra_field_accuracy.summary.v1",
        "generated_at_utc": _dt.datetime.now(
            _dt.timezone.utc
        ).isoformat(timespec="seconds").replace("+00:00", "Z"),
        "versions": {
            "python_version": sys.version.split()[0],
            "python_executable": sys.executable,
            "platform": platform.platform(),
            "radia_version": _version("radia"),
            "ngsolve_version": _version("ngsolve"),
            "numpy_version": _version("numpy"),
        },
        "inputs": {
            "analytical_reference_results": analytical is not None,
            "solver_comparison_results": solver is not None,
            "evaluation_results": extraction is not None,
            "ngsolve_reference_results": ngsolve_ref is not None,
        },
        "analytical_reference": analytical_stats,
        "radia_solver_comparison": solver_stats,
        "legacy_ngsolve_extraction": {
            **extraction_stats,
            "status": "known_issue" if checks["legacy_ngsolve_extraction_known_issue_recorded"] else "review",
        },
        "ngsolve_reference": {
            "vectors": len(ngsolve_values),
            "finite_vectors": finite_ngsolve,
            "status": "ok" if finite_ngsolve == len(ngsolve_values) and ngsolve_values else "review",
        },
        "checks": checks,
        "overall_status": (
            "PASS_WITH_KNOWN_ISSUES"
            if checks["analytical_reference_pass"] and checks["radia_solver_comparison_pass"]
            else "CHECK"
        ),
    }


def main() -> int:
    summary = build_summary()
    out = HERE / "comparison_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
                   encoding="utf-8")

    print("=" * 70)
    print("Tetrahedral field-accuracy validation summary")
    print("=" * 70)
    print(f"Analytical tetra-vs-hexa avg error: "
          f"{summary['analytical_reference']['avg_percent']:.6g}%")
    print(f"Radia solver tetra-vs-hexa avg error: "
          f"{summary['radia_solver_comparison']['avg_percent']:.6g}%")
    print(f"Legacy NGSolve extraction avg error: "
          f"{summary['legacy_ngsolve_extraction']['avg_percent']:.6g}% "
          f"({summary['legacy_ngsolve_extraction']['status']})")
    print(f"NGSolve reference finite vectors: "
          f"{summary['ngsolve_reference']['finite_vectors']} / "
          f"{summary['ngsolve_reference']['vectors']}")
    print(f"Overall: {summary['overall_status']}")
    print(f"Wrote: {out}")
    return 0 if summary["overall_status"].startswith("PASS") else 1


if __name__ == "__main__":
    raise SystemExit(main())
