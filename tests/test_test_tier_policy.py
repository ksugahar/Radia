"""Guard the physical split and bounded execution of Radia test evidence."""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_retired_examples_tier_is_not_restored():
    assert not (ROOT / "examples").exists()
    assert '"examples/*"' not in (ROOT / "pyproject.toml").read_text(
        encoding="utf-8"
    )
    ignore_lines = {
        line.strip()
        for line in (ROOT / ".gitignore").read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    }
    assert not any(line.lstrip("!").startswith("examples/") for line in ignore_lines)


def test_pytest_has_one_repository_configuration_source():
    assert not (ROOT / "pytest.ini").exists()
    assert "[tool.pytest.ini_options]" in (ROOT / "pyproject.toml").read_text(
        encoding="utf-8"
    )


def test_performance_scripts_live_in_the_validation_benchmark_lane():
    assert not list((ROOT / "tests").glob("benchmark_*.py"))
    assert not list((ROOT / "tests").glob("test_*performance*.py"))
    assert (
        ROOT / "validation_test" / "benchmarks" / "benchmark_field_parallel.py"
    ).is_file()


def test_fast_ci_profiles_have_explicit_paths_and_runtime_budgets():
    manifest = json.loads(
        (ROOT / "tests" / "test_tier_manifest.json").read_text(encoding="utf-8")
    )
    profiles = manifest["profiles"]
    assert profiles["fast-contracts"]["max_elapsed_seconds"] <= 60
    assert profiles["native-smoke"]["max_elapsed_seconds"] <= 120

    for name, profile in profiles.items():
        assert profile["paths"], name
        for path in profile["paths"]:
            assert (ROOT / path).is_file(), (name, path)

    fast = (ROOT / ".github" / "workflows" / "radia-fast.yml").read_text(
        encoding="utf-8"
    )
    native = (ROOT / ".github" / "workflows" / "build-test.yml").read_text(
        encoding="utf-8"
    )
    assert "tools/run_test_tier.py" in fast
    assert "tools/run_test_tier.py --profile native-smoke" in native
