"""Historical evidence gate contracts use synthetic Git sources, never physics."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
from check_validation_records import (covered_hashes, validate_record, validate_tree,
                                      check_schema_definition, check_shape)
from ci_distribution_contract import owner_rule


def git(root, *args):
    return subprocess.check_output(["git", "-C", str(root), "-c", "user.name=Record Contract",
                                    "-c", "user.email=record@example.invalid", *args], text=True).strip()


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2), encoding="utf-8")


@pytest.fixture
def sample(tmp_path):
    git(tmp_path, "init")
    (tmp_path / "pyproject.toml").write_text('[project]\ndependencies=["ngsolve==6.2.2608"]\n')
    (tmp_path / "physics.py").write_text("VALUE = 1\n")
    (tmp_path / "original.py").write_text("print('original')\n")
    git(tmp_path, "add", ".")
    git(tmp_path, "commit", "-m", "Synthetic original source")
    source = git(tmp_path, "rev-parse", "HEAD")
    folder = tmp_path / "validation_records" / "principle" / "20261022_contract"
    folder.mkdir(parents=True)
    script = b"# immutable recorded entry\npass\n"
    result = b'{"error": 0.001}\n'
    (folder / "entry.py").write_bytes(script)
    (folder / "result.json").write_bytes(result)
    record = {
        "schema": "radia.validation-record.v2", "id": "principle/20261022_contract",
        "milestone": "Synthetic integrity milestone", "claim": "Recorded bytes match their source",
        "date": "2026-10-22", "radia_version": "5.3.3",
        "conditions": {"problem": "Synthetic contract", "mesh": "Not applicable: no solver", "hardware": "CPU contract fixture"},
        "comparison": {"reference": "Independent SHA256 of synthetic bytes", "independence": "Integrity only; no physics claim"},
        "kind": "investigation", "question": "Synthetic artifact integrity", "outcome": "observed",
        "source": {"repo_sha": source, "provenance": "exact", "dirty": False, "native_sha256": {}, "wheel_sha256": None, "source_only": True},
        "covered_paths": ["physics.py", "pyproject.toml"],
        "covered_sha256": covered_hashes(tmp_path, ["physics.py", "pyproject.toml"]),
        "entry": {"script": "entry.py", "script_sha256": hashlib.sha256(script).hexdigest(),
                  "command": ["python", "entry.py"]},
        "environment": {"host_role": "neutral-cpu", "os": "Windows", "python": "3.12",
                        "ngsolve": "6.2.2608", "threads": 1},
        "result": {"path": "result.json", "sha256": hashlib.sha256(result).hexdigest()},
        "key_values": {"error": {"value": .001, "unit": "1"}}, "limits": {"error": .01},
        "timing_s": {"elapsed": .5}, "status": "current", "stale_reason": None,
        "superseded_by": None, "private": False,
    }
    path = folder / "record.json"
    write(path, record)
    return tmp_path, path, record


def test_root_and_package_record_owner_rules():
    manifest = json.loads((ROOT / "tools/ci_distribution_manifest.json").read_text())
    assert owner_rule(manifest, "validation_records/case/record.json")["owners"] == ["radia"]
    for package in manifest["distributions"]:
        if package != "radia":
            assert owner_rule(manifest, f"packages/{package}/validation_records/case/record.json")["owners"] == [package]


def test_repository_records_pass_without_replaying_entries():
    result = validate_tree(ROOT)
    assert result["passed"] and result["numerical_work_executed"] is False


def test_valid_record_checks_bytes_and_source(sample):
    root, path, _ = sample
    assert validate_record(root, path)["status"] == "current"
    assert validate_tree(root)["records"][0]["covered_files"] == 2


@pytest.mark.parametrize("change,match", [
    ("result", "result sha256"), ("script", "entry script sha256"),
    ("private", "schema"), ("covered", "matches no tracked"), ("escape", "unsafe relative"),
    ("unknown", "schema"), ("host", "schema"),
])
def test_corrupt_or_unqualified_evidence_is_rejected(sample, change, match):
    root, path, record = sample
    if change == "result":
        (path.parent / "result.json").write_text('{}')
    elif change == "script":
        (path.parent / "entry.py").write_text('pass')
    elif change == "private":
        record["private"] = True
    elif change == "covered":
        record["covered_paths"] = ["missing/*.py"]
    elif change == "escape":
        record["result"]["path"] = "../outside.json"
    elif change == "unknown":
        record["unrecognized"] = 1
    elif change == "host":
        record["environment"]["host_role"] = "machine-name"
    write(path, record)
    with pytest.raises(ValueError, match=match):
        validate_record(root, path)


def test_current_source_change_fails_but_stale_retains_original_evidence(sample):
    root, path, record = sample
    before = (path.parent / "result.json").read_bytes()
    (root / "physics.py").write_text('VALUE = 2\n')
    git(root, "add", "physics.py"); git(root, "commit", "-m", "Synthetic changed source")
    with pytest.raises(ValueError, match="claims current"):
        validate_record(root, path)
    record.update(status="stale", stale_reason="covered source changed")
    write(path, record)
    assert validate_record(root, path)["changed_covered_files"] == ["physics.py"]
    assert (path.parent / "result.json").read_bytes() == before



@pytest.mark.parametrize("question", ["C:/private/input", "192.0.2.1", "commercial reference data"])
def test_private_metadata_is_rejected(sample, question):
    root, path, record = sample
    record["question"] = question
    write(path, record)
    with pytest.raises(ValueError, match="private"):
        validate_record(root, path)



def test_orphan_result_cannot_bypass_schema(sample):
    root, _, _ = sample
    write(root / "validation_records" / "orphan" / "result.json", {})
    with pytest.raises(ValueError, match="no record.json"):
        validate_tree(root)


def test_record_identity_is_its_directory(sample):
    root, path, record = sample
    record["id"] = "different/20261022_record"
    write(path, record)
    with pytest.raises(ValueError, match="id must match"):
        validate_record(root, path)


def test_uncovered_source_change_does_not_make_a_record_stale(sample):
    root, path, _ = sample
    (root / "notes.txt").write_text('unrelated source note\n')
    git(root, "add", "notes.txt"); git(root, "commit", "-m", "Synthetic unrelated note")
    assert validate_record(root, path)["status"] == "current"


def test_nan_result_cannot_pass_hash_validation(sample):
    root, path, record = sample
    payload = b'{"error": NaN}'
    (path.parent / "result.json").write_bytes(payload)
    record["result"]["sha256"] = hashlib.sha256(payload).hexdigest()
    write(path, record)
    with pytest.raises(ValueError, match="non-finite"):
        validate_record(root, path)


def test_fixed_contract_rejects_unknown_keywords():
    with pytest.raises(ValueError, match="unsupported keywords"):
        check_schema_definition({"type": "object", "$ref": "not-supported"})


@pytest.mark.parametrize("value", [0, 1, "false", None])
def test_fixed_contract_boolean_is_not_integer(value):
    with pytest.raises(ValueError, match="schema"):
        check_shape(value, {"const": False})
    check_shape(False, {"const": False})


def test_fixed_contract_rejects_nonfinite_number():
    with pytest.raises(ValueError, match="schema"):
        check_shape(float('inf'), {"type": "number"})


def test_record_gate_works_without_site_packages(sample):
    root, _, _ = sample
    code = ("import sys; from pathlib import Path; sys.path.insert(0, sys.argv[1]); "
            "from check_validation_records import validate_tree; "
            "assert validate_tree(Path(sys.argv[2]))['passed']")
    result = subprocess.run([sys.executable, "-I", "-S", "-c", code,
                             str(ROOT / "tools"), str(root)], capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr


@pytest.mark.parametrize("field", ["milestone", "claim", "date", "radia_version", "conditions", "comparison"])
def test_milestone_metadata_is_required(sample, field):
    root, path, record = sample
    del record[field]
    write(path, record)
    with pytest.raises(ValueError, match="required"):
        validate_record(root, path)


@pytest.mark.parametrize("field,value", [("claim", "first\nsecond"), ("claim", "  "),
    ("date", "2026-02-30"), ("radia_version", "unknown"),
    ("conditions", {"problem": "case"}), ("comparison", {"reference": "analytic"})])
def test_invalid_milestone_metadata_is_rejected(sample, field, value):
    root, path, record = sample
    record[field] = value
    write(path, record)
    with pytest.raises(ValueError):
        validate_record(root, path)


@pytest.mark.parametrize("text", ["https://doi.org/10.1000/example", r"\mu_0", "CUDA 12.8.0.61"])
def test_public_text_is_not_private(sample, text):
    root, path, record = sample
    record["question"] = text
    write(path, record)
    assert validate_record(root, path)["status"] == "current"


def test_source_only_or_binary_identity_required(sample):
    root, path, record = sample
    record["source"]["source_only"] = False
    write(path, record)
    with pytest.raises(ValueError, match="requires native"):
        validate_record(root, path)


def test_superseding_record_must_exist(sample):
    root, path, record = sample
    record["superseded_by"] = "principle/20261022_missing"
    write(path, record)
    with pytest.raises(ValueError, match="target must exist"):
        validate_record(root, path)


def test_unavailable_source_history_is_provenance_only(sample):
    root, path, record = sample
    record["source"]["repo_sha"] = "b" * 40
    record["environment"].pop("ngsolve")
    write(path, record)
    assert validate_record(root, path)["status"] == "current"


def test_normalized_covered_hashes_accept_crlf(sample):
    root, path, record = sample
    (root / "physics.py").write_bytes(b"VALUE = 1\r\n")
    assert validate_record(root, path)["status"] == "current"


def test_schema_error_includes_actual_expected_and_record_path(sample):
    root, path, record = sample
    record["environment"]["threads"] = "many"
    write(path, record)
    with pytest.raises(ValueError) as caught:
        validate_record(root, path)
    assert str(path) in str(caught.value)
    assert "actual='many'" in str(caught.value) and "expected type" in str(caught.value)


@pytest.mark.parametrize("change", ["added", "removed", "dependency_patch"])
def test_content_inventory_detects_changes_without_history(sample, change):
    root, path, record = sample
    if change == "added":
        record["covered_paths"].append("*.py")
        (root / "added.py").write_text("NEW = 1\n")
        git(root, "add", "added.py"); git(root, "commit", "-m", "Add covered file")
    elif change == "removed":
        (root / "physics.py").unlink()
    else:
        (root / "pyproject.toml").write_text('[project]\ndependencies=["ngsolve==6.2.2609"]\n')
    write(path, record)
    with pytest.raises(ValueError, match="claims current.*actual.*expected"):
        validate_record(root, path)


def test_pre_pin_trial_version_is_recordable(sample):
    root, path, record = sample
    record["environment"]["ngsolve"] = "6.2.2607"
    write(path, record)
    assert validate_record(root, path)["status"] == "current"


def test_invalid_id_calendar_date_is_rejected(sample):
    root, path, record = sample
    record["id"] = "principle/20260230_contract"
    write(path, record)
    with pytest.raises(ValueError, match="day"):
        validate_record(root, path)


def test_covered_metadata_is_author_selected_not_forced_to_root(sample):
    root, path, record = sample
    record["covered_paths"] = ["physics.py"]
    record["covered_sha256"] = covered_hashes(root, record["covered_paths"])
    (root / "pyproject.toml").write_text('[project]\nversion="5.3.4"\ndependencies=["ngsolve==6.2.2608"]\n')
    write(path, record)
    assert validate_record(root, path)["status"] == "current"


def pre_rewrite(record):
    record["source"].update(repo_sha=None, provenance="pre-rewrite")
    record.update(status="stale", stale_reason="pre-rewrite: execution identity lost in confirmed rewrite")


def test_pre_rewrite_keeps_null_source_and_stale_without_changing_artifacts(sample):
    root, path, record = sample
    before = {name: (path.parent / name).read_bytes() for name in ("entry.py", "result.json")}
    pre_rewrite(record)
    write(path, record)
    result = validate_record(root, path)
    assert result["source"] is None and result["provenance"] == "pre-rewrite"
    assert result["status"] == "stale" and result["changed_covered_files"] == []
    assert result["numerical_work_executed"] is False
    assert all((path.parent / name).read_bytes() == raw for name, raw in before.items())


@pytest.mark.parametrize("change,match", [
    ("exact-null", "exact requires a full"),
    ("pre-rewrite-sha", "repo_sha=null"),
    ("current", "requires status=stale"),
    ("wrong-reason", "stale_reason requires"),
    ("empty-explanation", "stale_reason requires"),
])
def test_provenance_cross_field_contract_is_explicit(sample, change, match):
    root, path, record = sample
    old = record["source"]["repo_sha"]
    pre_rewrite(record)
    if change == "exact-null":
        record["source"]["provenance"] = "exact"
    elif change == "pre-rewrite-sha":
        record["source"]["repo_sha"] = old
    elif change == "current":
        record.update(status="current", stale_reason=None)
    elif change == "wrong-reason":
        record["stale_reason"] = "covered files changed"
    else:
        record["stale_reason"] = "pre-rewrite:   "
    write(path, record)
    with pytest.raises(ValueError, match=match) as error:
        validate_record(root, path)
    if change == "pre-rewrite-sha":
        assert old not in str(error.value)


@pytest.mark.parametrize("change,match", [
    ("artifact", "requires native"), ("result", "result sha256"),
    ("script", "entry script sha256"), ("dirty", "schema"),
    ("covered", "covered_sha256 must contain"), ("private", "private"),
])
def test_pre_rewrite_does_not_bypass_existing_integrity_guards(sample, change, match):
    root, path, record = sample
    pre_rewrite(record)
    if change == "artifact":
        record["source"]["source_only"] = False
    elif change == "result":
        (path.parent / "result.json").write_bytes(b"{}")
    elif change == "script":
        (path.parent / "entry.py").write_bytes(b"pass")
    elif change == "dirty":
        record["source"]["dirty"] = True
    elif change == "covered":
        record["covered_sha256"] = {}
    else:
        record["question"] = "C:/private/original"
    write(path, record)
    with pytest.raises(ValueError, match=match):
        validate_record(root, path)


@pytest.mark.parametrize("change", ["v1", "missing-provenance", "public-old-sha"])
def test_latest_schema_has_no_legacy_or_old_identity_field(sample, change):
    root, path, record = sample
    if change == "v1":
        record["schema"] = "radia.validation-record.v1"
    elif change == "missing-provenance":
        del record["source"]["provenance"]
    else:
        pre_rewrite(record)
        record["source"]["old_sha"] = "b" * 40
    write(path, record)
    with pytest.raises(ValueError, match="schema"):
        validate_record(root, path)
