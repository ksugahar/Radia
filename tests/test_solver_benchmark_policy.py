"""Fast placement guards for the PEEC solver benchmark evidence."""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BENCHMARKS = ROOT / "validation_test" / "solver_benchmarks"
LEGACY_METHOD = re.compile(r"\b(?:rad|radia)\.Solve\([^\n]*,\s*[12]\s*[,)]")


def test_solver_benchmarks_are_script_and_json_validation_artifacts():
    assert not list(BENCHMARKS.glob("*.ipynb"))

    manifest = json.loads(
        (BENCHMARKS / "peec_solver_benchmarks_results.json").read_text(
            encoding="utf-8"
        )
    )
    assert manifest["schema"] == "radia.validation.peec-solver-benchmarks.v2"
    assert manifest["policy"]["validation_json_required"] is True
    assert manifest["policy"]["validation_notebook_required"] is False

    for group in ("raw_benchmark_files", "source_files"):
        for artifact in manifest[group]:
            path = ROOT / artifact["path"]
            if group == "source_files":
                payload = path.read_text(encoding="utf-8").encode()
            else:
                payload = path.read_bytes()
            assert hashlib.sha256(payload).hexdigest() == artifact["sha256"], path


def test_current_solver_benchmarks_do_not_restore_retired_solve_methods():
    offenders = []
    for path in sorted(BENCHMARKS.glob("*.py")):
        if LEGACY_METHOD.search(path.read_text(encoding="utf-8")):
            offenders.append(path.name)
    assert not offenders, "retired rad.Solve method routes: " + ", ".join(offenders)
