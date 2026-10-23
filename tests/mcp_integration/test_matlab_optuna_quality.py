"""Public Optuna evidence, planning and release contracts.

Guards: Missing/altered installed evidence being accepted or silently replaced
by a repository, wrong release/planning boundaries and incomplete API closure.
Assumes: Stored reviewed JSON and declared package layout. The installed-copy
check uses a presence-only MEX placeholder and makes no native execution claim.
"""

import json
from copy import deepcopy
from radia_mcp.matlab import optuna_quality
from pathlib import Path

from radia_mcp.matlab.optuna_quality import (
    matlab_optuna_benchmark_plan,
    matlab_optuna_health,
    matlab_optuna_oracle_plan,
    matlab_optuna_release_gate,
)

REPO_ROOT = Path(__file__).resolve().parents[2]


def test_health_uses_distribution_and_upstream_manifests_as_truth(monkeypatch):
    # A clean source checkout intentionally has no ignored Windows MEX binary.
    # Source health therefore validates the standalone CUDA gateway command table;
    # installed-distribution health remains responsible for requiring the MEX.
    mex_path = REPO_ROOT / "packages/radia-optuna/matlab" / "optuna_mex.mexw64"
    path_is_file = Path.is_file
    monkeypatch.setattr(
        Path,
        "is_file",
        lambda path: False if path == mex_path else path_is_file(path),
    )
    health = matlab_optuna_health(str(REPO_ROOT))
    manifest = json.loads(
        (
            REPO_ROOT
            / "packages"
            / "radia-optuna"
            / "src"
            / "radia_optuna"
            / "manifest.json"
        ).read_text(encoding="utf-8")
    )
    coverage = json.loads(
        (REPO_ROOT / "packages/radia-optuna/matlab" / "optuna50_api_coverage.json").read_text(
            encoding="utf-8"
        )
    )
    test_manifest = json.loads(
        (
            REPO_ROOT
            / "packages"
            / "radia-optuna"
            / "oracles"
            / "optuna_test_manifest.json"
        ).read_text(encoding="utf-8")
    )

    assert health["ok"] is coverage["full_compatibility_complete"]
    assert health["source_kind"] == "repository"
    assert health["distribution"]["matlab_file_count"] == manifest[
        "matlab_file_count"
    ]
    assert health["native"]["command_count"] == manifest["native_command_count"]
    assert health["native"]["present"] is False
    assert health["native"]["runtime_ready"] is False
    assert health["native"]["build_required"] is True
    assert health["native"]["source_contract_present"] is True
    assert health["native"]["source_command_count"] == manifest[
        "native_command_count"
    ]
    assert health["simulink"]["entry_count"] == len(
        manifest["simulink_entry_points"]
    )
    assert health["public_api"]["entry_count"] == len(coverage["entries"])
    assert health["public_api"]["present_count"] == len(coverage["entries"])
    assert health["public_api"]["verified_count"] == coverage[
        "oracle_verified_count"
    ]
    assert health["public_api"]["asserted_count"] == coverage[
        "oracle_asserted_count"
    ]
    assert health["public_api"]["required_verified_count"] == coverage[
        "required_oracle_mapped_count"
    ]
    assert health["public_api"]["required_asserted_count"] == coverage["required_oracle_asserted_count"]
    assert health["public_api"]["missing_count"] == 0
    assert health["public_api"]["partial_count"] == 0
    assert health["public_api"]["unmapped_count"] == 0
    assert health["public_api"]["complete"] is coverage["full_compatibility_complete"]
    assert health["oracle"]["recorded_fixture_sha256"] == health["oracle"][
        "actual_fixture_sha256"
    ]
    assert health["oracle"]["recorded_inventory_sha256"] == health["oracle"][
        "actual_inventory_sha256"
    ]
    assert health["test_policy"]["entry_count"] == len(test_manifest["entries"])
    assert sum(health["test_policy"]["classifications"].values()) == len(
        test_manifest["entries"]
    )
    assert health["stewardship"]["complete"] is True


def _mock_complete_health(monkeypatch):
    health = deepcopy(matlab_optuna_health(str(REPO_ROOT)))
    health.update(ok=True, status="ready", errors=[])
    monkeypatch.setattr(optuna_quality, "matlab_optuna_health", lambda *a, **k: health)


def test_oracle_plan_targets_official_matlab_mcp_and_pinned_upstream(monkeypatch):
    plan = matlab_optuna_oracle_plan("all", str(REPO_ROOT))

    assert plan["ok"] is matlab_optuna_health(str(REPO_ROOT))["ok"]
    _mock_complete_health(monkeypatch)
    plan = matlab_optuna_oracle_plan("all", str(REPO_ROOT))
    assert plan["ok"] is True
    assert plan["execute_with"] == "evaluate_matlab_code"
    assert plan["runtime_owner"] == "MathWorks official MATLAB MCP Server"
    assert plan["ownership"] == {
        "shared_study_trial_tools": "optuna/optuna-mcp",
        "seeded_numeric_oracle": "direct optuna==5.0.0",
        "matlab_execution": "MathWorks official MATLAB MCP Server",
        "matlab_difference_contract": "radia-mcp.matlab",
    }
    assert "generate_optuna50_mcp_oracle.py" in " ".join(
        plan["fixture_regeneration"]["commands"]
    )
    assert "runtests" not in plan["matlab_code"]
    assert "testsuite" in plan["matlab_code"]
    assert plan["output"].endswith("radia-optuna-all-summary.json")
    assert "fopen" in plan["matlab_code"]
    assert plan["expected_test_count"] == matlab_optuna_health(str(REPO_ROOT))[
        "test_policy"
    ]["entry_count"]


def test_benchmark_plan_reads_checked_workload_settings_and_separates_startup(monkeypatch):
    plan = matlab_optuna_benchmark_plan(
        str(REPO_ROOT), r"C:\temp\radia-optuna-quality-test"
    )

    assert plan["ok"] is matlab_optuna_health(str(REPO_ROOT))["ok"]
    _mock_complete_health(monkeypatch)
    plan = matlab_optuna_benchmark_plan(str(REPO_ROOT), r"C:\temp\radia-optuna-quality-test")
    assert plan["ok"] is True
    assert plan["settings"] == {
        "trials": 100,
        "total_repeats": 11,
        "warmup_repeats": 3,
        "reported_repeats": 8,
        "scalar_seed": 37,
        "grouped_conditional_seed": 101,
    }
    assert plan["same_host_required"] is True
    assert plan["cold_start"]["output"].endswith("optuna_mex_cold.json")
    assert plan["python_warmed"]["runtime"] == "direct optuna==5.0.0"
    assert plan["matlab_warmed"]["execute_with"] == "evaluate_matlab_code"
    assert plan["acceptance"]["max_matlab_to_python_warmed_time_ratio"] == 1.0
    assert plan["acceptance"]["startup_measured_separately"] is True


def _passing_release_evidence() -> dict[str, object]:
    health = matlab_optuna_health(str(REPO_ROOT))
    version = health["distribution"]["version"]
    file_count = health["distribution"]["expected_matlab_file_count"]
    command_count = health["native"]["command_count"]
    simulink_count = health["simulink"]["entry_count"]
    test_count = health["test_policy"]["entry_count"]
    mex_sha256 = health["native"]["sha256"] or ("A" * 64)
    settings = {
        "trials": 100,
        "total_repeats": 11,
        "warmup_repeats": 3,
        "reported_repeats": 8,
    }
    python_result = {
        "schema": "radia.validation.optuna50-performance-runtime.v1",
        "runtime": "python-upstream",
        "host": "runtime-fixture",
        "versions": {"optuna": health["distribution"]["upstream_version"]},
        "settings": settings,
        "scalar": {"median_warmed_seconds": 2.0, "checksum": 20.0},
        "grouped_conditional": {
            "median_warmed_seconds": 4.0,
            "checksum": 104.0,
        },
        "trials_dataframe": {
            "median_warmed_seconds": 0.2,
            "rows": 1000,
            "columns": 7,
        },
    }
    matlab_result = {
        "schema": "radia.validation.optuna50-performance-runtime.v1",
        "runtime": "matlab",
        "host": "runtime-fixture",
        "versions": {"optuna_mex_command_count": command_count},
        "settings": settings,
        "scalar": {"median_warmed_seconds": 1.0, "checksum": 20.0},
        "grouped_conditional": {
            "median_warmed_seconds": 2.0,
            "checksum": 104.0,
        },
        "trials_dataframe": {
            "median_warmed_seconds": 0.1,
            "rows": 1000,
            "columns": 7,
        },
    }
    return {
        "wheel_verification": {
            "schema": "radia-optuna.wheel-verification.v1",
            "ok": True,
            "version": version,
            "matlab_file_count": file_count,
            "native_command_count": command_count,
            "simulink_standalone": True,
            "simulink_entry_count": simulink_count,
            "source_fidelity_verified": True,
            "source_fidelity_file_count": file_count + 1,
        },
        "doctor": {
            "schema": "radia-optuna.doctor.v1",
            "ok": True,
            "version": version,
            "matlab_file_count": file_count,
            "native_command_count": command_count,
            "simulink_standalone": True,
            "simulink_entry_count": simulink_count,
            "upstream_notices_complete": True,
            "mex_sha256": mex_sha256,
        },
        "matlab_tests": {
            "schema": "radia.optuna.oracle-test-summary.v1",
            "total": test_count,
            "passed": test_count,
            "failed": 0,
            "incomplete": 0,
        },
        "simulink_e2e": {
            "schema": "radia-optuna.standalone-simulink-test.v1",
            "ok": True,
            "complete_trials": 4,
            "failed_trials": 1,
            "simulink_block_trials": 4,
        },
        "table_resume": {
            "schema": "radia-optuna.table-resume-test.v1",
            "ok": True,
            "source_trial_count": 4,
            "restored_trial_count": 4,
            "table_names": [
                "TrialTable",
                "ObjectiveTable",
                "ParamTable",
                "IntermediateTable",
                "UserAttrTable",
                "ConstraintTable",
                "SamplerStateTable",
            ],
        },
        "performance": {
            "python": python_result,
            "matlab": matlab_result,
            "mex_cold": {
                "schema": "radia.validation.optuna-mex-first-call.v1",
                "host": "runtime-fixture",
                "median_first_call_seconds": 0.01,
                "binary_sha256": mex_sha256,
            },
        },
    }


def test_release_gate_accepts_complete_evidence_and_rejects_regression(monkeypatch):
    evidence = _passing_release_evidence()
    actual_health = matlab_optuna_health(str(REPO_ROOT))
    result = matlab_optuna_release_gate(evidence, str(REPO_ROOT))
    assert result["ok"] is actual_health["ok"]
    # Synthetic passing health tests gate mechanics, not current API closure.
    passing_health = deepcopy(actual_health)
    passing_health.update(ok=True, status="ready", errors=[])
    monkeypatch.setattr(optuna_quality, "matlab_optuna_health", lambda *a, **k: passing_health)
    result = matlab_optuna_release_gate(evidence, str(REPO_ROOT))

    assert result["ok"] is True
    assert result["performance"]["observed_matlab_to_python_ratios"] == {
        "scalar": 0.5,
        "grouped_conditional": 0.5,
        "trials_dataframe": 0.5,
    }

    evidence["performance"]["matlab"]["scalar"]["median_warmed_seconds"] = 2.1
    evidence["table_resume"]["table_names"].remove("SamplerStateTable")
    result = matlab_optuna_release_gate(
        json.dumps(evidence), str(REPO_ROOT), max_warmed_time_ratio=1.0
    )

    assert result["ok"] is False
    assert any("scalar MATLAB/Python" in error for error in result["errors"])
    assert any("SamplerStateTable" in error for error in result["errors"])


def test_installed_evidence_is_authoritative_without_repository(tmp_path, monkeypatch):
    import importlib.util
    import shutil
    from types import SimpleNamespace
    import pytest
    from radia_mcp.matlab.optuna_oracle import matlab_optuna_compatibility_contract

    source = REPO_ROOT / "packages/radia-optuna"
    target = tmp_path / "radia_optuna"
    shutil.copytree(source / "src/radia_optuna", target,
                    ignore=shutil.ignore_patterns("__pycache__"))
    shutil.copytree(source / "matlab", target / "matlab",
                    ignore=shutil.ignore_patterns("*.mexw64"))
    (target / "oracles").mkdir()
    for path in (source / "oracles").glob("*.json"):
        shutil.copy2(path, target / "oracles" / path.name)
    (target / "matlab/optuna_mex.mexw64").write_bytes(b"presence-only placeholder")
    shutil.copy2(source / "THIRD_PARTY_NOTICES.md", target / "matlab/THIRD_PARTY_NOTICES.md")
    monkeypatch.delenv("RADIA_REPO_ROOT", raising=False)
    original = importlib.util.find_spec
    monkeypatch.setattr(importlib.util, "find_spec", lambda name, *a, **k:
                        SimpleNamespace(origin=str(target / "__init__.py"))
                        if name == "radia_optuna" else original(name, *a, **k))
    def no_repository():
        raise AssertionError("installed evidence must not search a repository")
    monkeypatch.setattr(optuna_quality, "_repository_candidates", no_repository)
    health = matlab_optuna_health()
    assert health["ok"] is True, health
    assert health["source_kind"] == "installed-distribution", health
    assert health["oracle"]["actual_fixture_sha256"] == health["oracle"]["recorded_fixture_sha256"]
    contract = matlab_optuna_compatibility_contract()
    assert contract["ok"] is True, contract
    assert contract["root"] == str(target.resolve()), contract
    assert contract["source_kind"] == "installed-distribution", contract

    fixture = target / "oracles/optuna50_oracle.json"
    content = fixture.read_bytes()
    for damage in ("tamper", "missing"):
        if damage == "tamper":
            fixture.write_bytes(content + b" ")
        else:
            fixture.unlink()
        rejected = matlab_optuna_health()
        assert rejected["ok"] is False, rejected
        with pytest.raises((OSError, ValueError)):
            matlab_optuna_compatibility_contract()
        fixture.write_bytes(content)
    ledger = target / "evidence_manifest.json"
    original_ledger = ledger.read_bytes()
    unsafe = json.loads(original_ledger)
    unsafe["files"]["python_oracle"]["path"] = "../outside.json"
    ledger.write_text(json.dumps(unsafe), encoding="utf-8")
    assert matlab_optuna_health()["ok"] is False
    with pytest.raises(ValueError, match="Unsafe"):
        matlab_optuna_compatibility_contract()
    ledger.write_bytes(original_ledger)
    assert matlab_optuna_health()["ok"] is True
    assert matlab_optuna_compatibility_contract()["ok"] is True
