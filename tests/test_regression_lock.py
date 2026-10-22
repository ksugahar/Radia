"""Synthetic prior-output locks test plumbing, not a physical formulation."""
import json
import hashlib
import platform
import tools.regression_lock as lock
from pathlib import Path
import subprocess
import sys

import pytest
from tools.regression_lock import assert_regression, load_baseline, propose_regeneration


@pytest.fixture
def baseline(tmp_path, monkeypatch):
    source = tmp_path / "source"
    source.mkdir()
    subprocess.run(["git", "init", "-q", str(source)], check=True)
    (source / "tracked.txt").write_text("synthetic source\n", encoding="utf-8")
    subprocess.run(["git", "add", "tracked.txt"], cwd=source, check=True)
    subprocess.run(["git", "-c", "user.name=Synthetic", "-c", "user.email=synthetic@example.invalid",
                    "commit", "-qm", "Synthetic source"], cwd=source, check=True)
    monkeypatch.setattr(lock, "SOURCE_ROOT", source)
    record = {"schema": "radia.regression-lock.v1", "source_sha": lock._source_identity()["head"],
        "conditions": {"problem": "synthetic", "mesh": "none", "threads": 1}, "runtime_s": 0.25,
        "quantities": {"field": {"expected": [2.0, 0.0], "atol": 0.01, "rtol": 0.001, "unit": "T"}}}
    path = tmp_path / "approved.json"
    path.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    return path, record


def test_match_uses_declared_tolerances_and_conditions(baseline):
    path, record = baseline
    assert_regression({"field": [2.005, 0.005]}, load_baseline(path), conditions=record["conditions"])


def test_mismatch_shows_named_actual_expected_and_tolerance(baseline):
    _, record = baseline
    with pytest.raises(AssertionError, match=r"field\[0\].*actual=3.0.*expected=2.0.*atol=.*rtol="):
        assert_regression({"field": [3.0, 0.0]}, record, conditions=record["conditions"])


@pytest.mark.parametrize("actual", [{}, {"field": [2.0]}, {"field": [2.0, 0.0], "extra": 1}])
def test_quantity_names_and_shapes_cannot_silently_shrink(baseline, actual):
    _, record = baseline
    with pytest.raises(AssertionError, match="actual=.*expected="):
        assert_regression(actual, record, conditions=record["conditions"])


def test_changed_conditions_fail_with_both_versions(baseline):
    _, record = baseline
    with pytest.raises(AssertionError, match="conditions: actual=.*expected="):
        assert_regression({"field": [2.0, 0.0]}, record, conditions={**record["conditions"], "threads": 2})


@pytest.mark.parametrize("runtime", [0, 31, 999999])
def test_runtime_metadata_has_no_comparison_or_30_second_gate(baseline, runtime):
    _, record = baseline
    record["runtime_s"] = runtime
    assert_regression({"field": [2.0, 0.0]}, record, conditions=record["conditions"])


@pytest.mark.parametrize("value", [float("nan"), float("inf"), True])
def test_non_numeric_or_nonfinite_actual_is_not_a_pass(baseline, value):
    _, record = baseline
    with pytest.raises(ValueError, match="finite"):
        assert_regression({"field": [value, 0.0]}, record, conditions=record["conditions"])


@pytest.mark.parametrize("field,value", [("source_sha", "short"), ("runtime_s", -1), ("conditions", {})])
def test_invalid_provenance_rejected(baseline, field, value):
    path, record = baseline
    record[field] = value
    path.write_text(json.dumps(record), encoding="utf-8")
    with pytest.raises(ValueError):
        load_baseline(path)


def test_duplicate_keys_rejected(baseline):
    path, _ = baseline
    path.write_text('{"schema":"a","schema":"b"}', encoding="utf-8")
    with pytest.raises(ValueError, match="duplicate JSON key"):
        load_baseline(path)


def test_proposal_preserves_approved_bytes_tolerances_and_records_diff(baseline, tmp_path):
    path, record = baseline
    original = path.read_bytes()
    report = propose_regeneration(path, actual={"field": [3.0, 0.0]}, source_sha=lock._source_identity()["head"],
        conditions=record["conditions"], runtime_s=31, reason="Synthetic implementation-change proposal", out=tmp_path / "proposal")
    assert path.read_bytes() == original
    proposed = load_baseline(tmp_path / "proposal/baseline.json")
    assert proposed["quantities"]["field"]["atol"] == record["quantities"]["field"]["atol"]
    assert proposed["quantities"]["field"]["rtol"] == record["quantities"]["field"]["rtol"]
    assert proposed["quantities"]["field"]["unit"] == record["quantities"]["field"]["unit"]
    assert proposed["runtime_s"] == 31 and proposed["source_sha"] == lock._source_identity()["head"]
    assert report["quantity_changes"] == ["field"] and report["diff"]
    assert report["approval"].startswith("proposal-only") and report["reason"]


def test_unchanged_proposal_preserves_original_format_bytes(baseline, tmp_path):
    path, record = baseline
    original = path.read_bytes()
    report = propose_regeneration(path, actual={"field": [2.0, 0.0]}, source_sha=record["source_sha"],
        conditions=record["conditions"], runtime_s=record["runtime_s"], reason="Synthetic unchanged comparison", out=tmp_path / "proposal")
    assert (tmp_path / "proposal/baseline.json").read_bytes() == original
    assert report["changed"] is False and report["diff"] == []


def test_proposal_does_not_overwrite_existing_directory(baseline, tmp_path):
    path, record = baseline
    before = path.read_bytes()
    with pytest.raises(FileExistsError):
        propose_regeneration(path, actual={"field": [2.0, 0.0]}, source_sha=record["source_sha"],
            conditions=record["conditions"], runtime_s=0.5, reason="Synthetic", out=tmp_path)
    assert path.read_bytes() == before


def test_reason_required_before_proposal_output(baseline, tmp_path):
    path, record = baseline
    with pytest.raises(ValueError, match="reason"):
        propose_regeneration(path, actual={"field": [2.0, 0.0]}, source_sha=record["source_sha"],
            conditions=record["conditions"], runtime_s=0.5, reason=" ", out=tmp_path / "proposal")
    assert not (tmp_path / "proposal").exists()


def test_explicit_cli_writes_proposal_only(baseline, tmp_path):
    path, record = baseline
    actual = tmp_path / "actual.json"; actual.write_text('{"field":[2.1,0.0]}', encoding="utf-8")
    conditions = tmp_path / "conditions.json"; conditions.write_text(json.dumps(record["conditions"]), encoding="utf-8")
    command = [sys.executable, str(Path(__file__).resolve().parents[1] / "tools/regression_lock.py"), "propose",
        "--baseline", str(path), "--actual", str(actual), "--conditions", str(conditions),
        "--source-sha", subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=Path(__file__).resolve().parents[1], text=True).strip(), "--runtime-s", "0.4", "--reason", "Synthetic CLI exercise", "--out", str(tmp_path / "proposal")]
    result = subprocess.run(command, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "proposal/proposal.json").is_file()


def test_multiple_quantities_keep_independent_scalar_and_array_tolerances(baseline):
    _, record = baseline
    record["quantities"]["power"] = {"expected": 3.0, "atol": 0.0, "rtol": 0.1, "unit": "W"}
    assert_regression({"field": [2.0, 0.0], "power": 3.2}, record, conditions=record["conditions"])
    with pytest.raises(AssertionError, match="power: actual=3.4.*expected=3.0"):
        assert_regression({"field": [2.0, 0.0], "power": 3.4}, record, conditions=record["conditions"])


def test_invalid_tolerance_is_rejected_not_silently_repaired(baseline):
    _, record = baseline
    record["quantities"]["field"]["atol"] = -0.1
    with pytest.raises(ValueError, match="atol.*nonnegative"):
        assert_regression({"field": [2.0, 0.0]}, record, conditions=record["conditions"])


@pytest.mark.parametrize("failure", ["dirty-tracked", "dirty-untracked", "sha-mismatch"])
def test_proposal_rejects_dirty_source_or_claimed_sha_before_writing(baseline, tmp_path, failure):
    path, record = baseline
    claimed = record["source_sha"]
    if failure == "dirty-tracked":
        (lock.SOURCE_ROOT / "tracked.txt").write_text("changed\n", encoding="utf-8")
    elif failure == "dirty-untracked":
        (lock.SOURCE_ROOT / "untracked.txt").write_text("new\n", encoding="utf-8")
    else:
        claimed = "0" * 40
    with pytest.raises(ValueError, match="dirty|git HEAD"):
        propose_regeneration(path, actual={"field": [2.0, 0.0]}, source_sha=claimed,
            conditions=record["conditions"], runtime_s=0.1, reason="Synthetic", out=tmp_path / "proposal")
    assert not (tmp_path / "proposal").exists()


def test_proposal_records_automatic_git_and_python_identity(baseline, tmp_path):
    path, record = baseline
    report = propose_regeneration(path, actual={"field": [2.0, 0.0]}, source_sha=record["source_sha"],
        conditions=record["conditions"], runtime_s=0.1, reason="Synthetic", out=tmp_path / "proposal")
    assert report["generation"] == {"head": record["source_sha"], "dirty": False,
        "status_porcelain": "", "python_version": platform.python_version(),
        "python_implementation": platform.python_implementation()}


def test_approved_hash_is_checkout_line_ending_independent(baseline, tmp_path):
    path, record = baseline
    lf = path.read_bytes().replace(b"\r\n", b"\n")
    reports = []
    for name, raw in [("lf", lf), ("crlf", lf.replace(b"\n", b"\r\n"))]:
        path.write_bytes(raw)
        reports.append(propose_regeneration(path, actual={"field": [2.0, 0.0]},
            source_sha=record["source_sha"], conditions=record["conditions"], runtime_s=record["runtime_s"],
            reason="Synthetic line ending comparison", out=tmp_path / name))
        assert path.read_bytes() == raw
    assert reports[0]["approved_sha256"] == reports[1]["approved_sha256"] == hashlib.sha256(lf).hexdigest()
    assert reports[0]["proposed_sha256"] == reports[1]["proposed_sha256"]


def test_numpy_scalar_and_array_diagnostic(baseline):
    import numpy as np
    _, record = baseline
    assert_regression({"field": [np.float64(2.0), np.float64(0.0)]}, record, conditions=record["conditions"])
    with pytest.raises(AssertionError, match=r"NumPy arrays: use .tolist"):
        assert_regression({"field": np.array([2.0, 0.0])}, record, conditions=record["conditions"])
    assert_regression({"field": np.array([2.0, 0.0]).tolist()}, record, conditions=record["conditions"])


def test_conditions_keep_numeric_json_types_distinct(baseline):
    _, record = baseline
    with pytest.raises(AssertionError, match="conditions"):
        assert_regression({"field": [2.0, 0.0]}, record, conditions={**record["conditions"], "threads": 1.0})
