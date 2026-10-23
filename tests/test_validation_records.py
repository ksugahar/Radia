"""Historical evidence contracts use synthetic Git sources, never physics.

Guards: Broken record identity, unapproved/misbound adoption declarations and
incorrect convergence arithmetic. Assumes: Synthetic approval strings test
provenance structure, not human authentication or reference/theory validity.
"""
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
from validation_adoption import criterion_hash, scope_hash


def seal_adoption(path, record):
    """Only synthetic fixtures: deliberately construct consistent byte identities."""
    record["adoption"]["approval"]["scope_sha256"] = scope_hash(record)
    record["adoption"]["approval"]["rationale_sha256"] = criterion_hash(record["adoption"]["rationale"])
    result = record["adoption"]["result"]
    payload = {"adoption_result": {k: v for k, v in result.items() if k != "artifact_sha256"}}
    data = (json.dumps(payload, sort_keys=True) + "\n").encode()
    (path.parent / "result.json").write_bytes(data)
    digest = hashlib.sha256(data).hexdigest()
    record["result"]["sha256"] = result["artifact_sha256"] = digest
    record["adoption"]["approval"]["criterion_sha256"] = criterion_hash(record["adoption"]["criterion"])
    write(path, record)


@pytest.fixture
def adoption_sample(sample):
    root, path, record = sample
    source = record["source"]["repo_sha"]
    levels = [{"id": f"L{i}", "parameters": {"h": h}}
              for i, h in enumerate((1., .5, .25), 1)]
    criterion = {
        "reference": {"kind": "independent-route", "citation": "Synthetic alternate method",
                      "independence": "Synthetic schema fixture; no solver or physical claim",
                      "shared_components": []},
        "quantity": {"name": "Synthetic relative error", "definition": "Fixture value", "unit": "1"},
        "variables": {"h": "Synthetic mesh size"}, "levels": levels,
        "series": [{"id": "h-series", "level_ids": ["L1", "L2", "L3"], "axis": "h",
                    "direction": "decreasing", "rate_model": "log-ratio", "minimum_rate": 1.5,
                    "monotonicity": "strict-decrease", "finest_threshold": {"operator": "lt", "value": .01},
                    "expected_theory": {"model": "h^2", "norm": "synthetic", "source": "Synthetic rationale",
                                        "assumptions": "Fixture only", "limitations": "No mathematical certification"}}],
        "failure_definition": "Any failed required gate", "saturation_policy": "blocks-pass",
        "residual_threshold": 1e-6,
        "threshold_verification": {"mode": "reviewer", "review_ref": "synthetic-review-20261023"},
    }
    result = {"source_sha": source, "artifact_sha256": "a" * 64, "verdict": "PASS", "levels": [
        {**level, "error": error, "relative_residual": 1e-9, "source_sha": source,
         "mesh_sha256": "b" * 64, "runtime_s": .1}
        for level, error in zip(levels, (.04, .01, .0025))]}
    record["outcome"] = "pass"
    record["category"] = "adoption"
    record["adoption"] = {"feature": "synthetic.production_entry", "criterion": criterion,
                           "approval": {"approved_by": "user", "approved_on": "2026-10-22",
                                        "approval_ref": "synthetic-consent-20261023", "criterion_sha256": "a" * 64},
                           "result": result, "protected_by": ["tests/test_source_lock.py"]}
    record["adoption"]["rationale"] = {
        "inference": "Synthetic evidence supports this declared fixture only", "scope": "No physics certification",
        "assumptions": [{"id": "reference.checked", "kind": "reference-verification",
                         "statement": "Synthetic independent reference is correct", "evidence_refs": ["artifact:result"],
                         "invalidated_by": ["reference.invalid"], "recheck": "Recheck reference and all dependent verdicts"}],
        "depends_on": ["reference.checked"], "break_conditions": [{"id": "reference.invalid", "condition": "Reference expression found incorrect"}]}
    lock = root / "tests/test_source_lock.py"
    lock.parent.mkdir(exist_ok=True)
    lock.write_text('"""Guards: Synthetic declared route.\nAssumes: surrogate test of adoption record '
                    + record["id"] + ' (convergence level 1)\n"""\n', encoding="utf-8")
    seal_adoption(path, record)
    return root, path, record


def test_convergence_adoption_reports_rates_without_executing(adoption_sample):
    root, path, _ = adoption_sample
    result = validate_record(root, path)
    assert result["adoption_claim"] is True
    assert result["observed_rates"]["h-series"] == pytest.approx([2., 2.])
    assert result["approval_authentication"] == "independent-review-required"
    assert result["numerical_work_executed"] is False


@pytest.mark.parametrize("change", ["missing-approval", "agent", "date", "private-ref", "criterion-hash",
                                 "old-schema", "missing-category"])
def test_invalid_adoption_approval_and_latest_only(adoption_sample, change):
    root, path, record = adoption_sample
    approval = record["adoption"]["approval"]
    if change == "missing-approval": del record["adoption"]["approval"]
    elif change == "agent": approval["approved_by"] = "agent"
    elif change == "date": approval["approved_on"] = "2026-02-30"
    elif change == "private-ref": approval["approval_ref"] = "C:/private/consent"
    elif change == "criterion-hash": approval["criterion_sha256"] = "0" * 64
    elif change == "old-schema": record["schema"] = "radia.validation-record.v2"
    else: del record["category"]
    write(path, record)
    with pytest.raises(ValueError): validate_record(root, path)


@pytest.mark.parametrize("change", ["source", "level-source", "mesh-hash", "parameters", "missing-level",
                                 "duplicate-level", "result-hash", "outcome", "unsafe-lock", "missing-lock"])
def test_adoption_identity_mismatches_fail(adoption_sample, change):
    root, path, record = adoption_sample
    a = record["adoption"]; result = a["result"]
    if change == "source": result["source_sha"] = "0" * 40
    elif change == "level-source": result["levels"][0]["source_sha"] = "0" * 40
    elif change == "mesh-hash": result["levels"][0]["mesh_sha256"] = "not-a-hash"
    elif change == "parameters": result["levels"][0]["parameters"] = {"h": 2.}
    elif change == "missing-level": result["levels"].pop()
    elif change == "duplicate-level": result["levels"][0]["id"] = "L2"
    elif change == "outcome": record["outcome"] = "observed"
    elif change == "unsafe-lock": a["protected_by"] = ["tests/../../escape.py"]
    elif change == "missing-lock": a["protected_by"] = ["tests/missing.py"]
    seal_adoption(path, record)
    if change == "result-hash":
        result["artifact_sha256"] = "0" * 64
        write(path, record)
    with pytest.raises(ValueError): validate_record(root, path)


def test_hashed_artifact_values_cannot_differ_from_declaration(adoption_sample):
    root, path, record = adoption_sample
    payload = json.loads((path.parent / "result.json").read_text())
    payload["adoption_result"]["levels"][0]["error"] = .08
    write(path.parent / "result.json", payload)
    record["result"]["sha256"] = record["adoption"]["result"]["artifact_sha256"] = hashlib.sha256((path.parent / "result.json").read_bytes()).hexdigest()
    write(path, record)
    with pytest.raises(ValueError, match="artifact values"): validate_record(root, path)


@pytest.mark.parametrize("change", ["not-monotone", "below-rate", "strict-equality", "residual", "nonfinite"])
def test_convergence_cannot_claim_pass_after_failed_gate(adoption_sample, change):
    root, path, record = adoption_sample
    a = record["adoption"]; rows = a["result"]["levels"]
    if change == "not-monotone": rows[1]["error"] = .05
    elif change == "below-rate": a["criterion"]["series"][0]["minimum_rate"] = 2.5
    elif change == "strict-equality": a["criterion"]["series"][0]["finest_threshold"]["value"] = .0025
    elif change == "residual": rows[0]["relative_residual"] = 1.1e-6
    else: rows[0]["error"] = float("nan")
    seal_adoption(path, record)
    with pytest.raises(ValueError): validate_record(root, path)


def test_inclusive_threshold_is_distinct_from_strict(adoption_sample):
    root, path, record = adoption_sample
    record["adoption"]["criterion"]["series"][0]["finest_threshold"] = {"operator": "le", "value": .0025}
    seal_adoption(path, record)
    assert validate_record(root, path)["adoption_claim"] is True


def test_preserve_failed_convergence_as_fail_not_adopted(adoption_sample):
    root, path, record = adoption_sample
    record["adoption"]["result"]["levels"][1]["error"] = .05
    record["adoption"]["result"]["verdict"] = "FAIL"
    record["outcome"] = "fail"
    seal_adoption(path, record)
    result = validate_record(root, path)
    assert result["adoption_claim"] is False and result["adoption_verdict"] == "FAIL"


def test_exact_zero_error_does_not_create_infinite_rate_pass(adoption_sample):
    root, path, record = adoption_sample
    record["adoption"]["result"]["levels"][-1]["error"] = 0.
    record["adoption"]["result"]["verdict"] = "UNRESOLVED"
    record["outcome"] = "observed"
    seal_adoption(path, record)
    assert validate_record(root, path)["adoption_claim"] is False


@pytest.mark.parametrize("reference", [
    {"kind": "analytic", "citation": "Missing DOI/symbolic evidence"},
    {"kind": "independent-route", "citation": "Missing independence/shared-parts"},
    {"kind": "measurement", "citation": "Missing uncertainty"},
    {"kind": "benchmark", "citation": "Missing uncertainty"},
])
def test_reference_requires_type_specific_evidence(adoption_sample, reference):
    root, path, record = adoption_sample
    record["adoption"]["criterion"]["reference"] = reference
    seal_adoption(path, record)
    with pytest.raises(ValueError): validate_record(root, path)


def test_safe_literal_binding_never_imports_validation_script(adoption_sample):
    root, path, record = adoption_sample
    text = b"LIMIT = 0.01\nraise RuntimeError('must not execute')\n"
    (path.parent / "validation.py").write_bytes(text)
    record["adoption"]["criterion"]["threshold_verification"] = {
        "mode": "literal-constants", "review_ref": "synthetic-wiring-review",
        "script": {"path": "validation.py", "sha256": hashlib.sha256(text).hexdigest()},
        "bindings": [{"constant": "LIMIT", "criterion_pointer": "/series/0/finest_threshold/value"}],
    }
    seal_adoption(path, record)
    assert validate_record(root, path)["threshold_verification"] == "literal-constants"
    record["adoption"]["criterion"]["series"][0]["finest_threshold"]["value"] = .02
    seal_adoption(path, record)
    with pytest.raises(ValueError, match="literal binding"): validate_record(root, path)


def test_pre_rewrite_evidence_cannot_be_exact_adoption(adoption_sample):
    root, path, record = adoption_sample
    record["source"].update(repo_sha=None, provenance="pre-rewrite")
    record.update(status="stale", stale_reason="pre-rewrite: synthetic historical evidence")
    write(path, record)
    with pytest.raises(ValueError, match="source provenance"): validate_record(root, path)


@pytest.mark.parametrize("cycle", [False, True])
def test_supersedes_requires_existing_acyclic_record(sample, cycle):
    root, path, record = sample
    record["supersedes"] = [{"id": record["id"] if cycle else "principle/20261021_missing",
                              "relation": "corrects", "reason": "Synthetic relation"}]
    write(path, record)
    with pytest.raises(ValueError): validate_record(root, path)


def test_record_without_adoption_has_no_adoption_claim(sample):
    root, path, _ = sample
    assert validate_record(root, path)["adoption_claim"] is False


def test_stale_pass_does_not_establish_current_adoption(adoption_sample):
    root, path, record = adoption_sample
    record.update(status="stale", stale_reason="Synthetic result needs revalidation")
    write(path, record)
    result = validate_record(root, path)
    assert result["adoption_verdict"] == "PASS" and result["adoption_claim"] is False


def test_same_byte_approval_does_not_approve_modified_criterion(adoption_sample):
    root, path, record = adoption_sample
    record["adoption"]["criterion"]["series"][0]["minimum_rate"] = 1.4
    write(path, record)
    with pytest.raises(ValueError, match="approved criterion hash"): validate_record(root, path)


def test_log_per_degree_p_series_is_distinct_from_h_slope(adoption_sample):
    root, path, record = adoption_sample
    criterion = record["adoption"]["criterion"]
    criterion["variables"] = {"p": "Synthetic polynomial degree"}
    series = criterion["series"][0]
    series.update(axis="p", direction="increasing", rate_model="log-per-unit", minimum_rate=1.)
    for i, (declared, actual) in enumerate(zip(criterion["levels"], record["adoption"]["result"]["levels"]), 1):
        declared["parameters"] = actual["parameters"] = {"p": i}
    seal_adoption(path, record)
    assert validate_record(root, path)["observed_rates"]["h-series"] == pytest.approx([1.38629436112] * 2)


def test_varying_geometry_inside_p_series_is_rejected(adoption_sample):
    root, path, record = adoption_sample
    c = record["adoption"]["criterion"]
    c["variables"]["geometry"] = "Synthetic geometry order"
    for i, (declared, actual) in enumerate(zip(c["levels"], record["adoption"]["result"]["levels"]), 1):
        declared["parameters"]["geometry"] = actual["parameters"]["geometry"] = i
    seal_adoption(path, record)
    with pytest.raises(ValueError, match="held-fixed"): validate_record(root, path)


@pytest.fixture
def analytic_sample(adoption_sample):
    root, path, record = adoption_sample
    # Synthetic canonical entry, never a scientific approval or execution claim.
    bibliography = root / "packages/radia-mcp/src/radia_mcp/bibliography/data/references.bib"
    bibliography.parent.mkdir(parents=True)
    bibliography.write_text('@article{SyntheticReference, doi={10.0000/synthetic-schema-fixture}}\n', encoding="utf-8")
    receipt = b'{"synthetic": true, "execution_claim": false}\n'
    (path.parent / "symbolic.json").write_bytes(receipt)
    record["adoption"]["criterion"]["reference"] = {
        "kind": "analytic", "citation": "Synthetic canonical reference/receipt fixture only",
        "bib_key": "SyntheticReference", "symbolic_verification": {
            "path": "symbolic.json", "sha256": hashlib.sha256(receipt).hexdigest()},
    }
    seal_adoption(path, record)
    return root, path, record, bibliography


@pytest.mark.parametrize("identifier", ["doi={10.0000/synthetic-schema-fixture}", "isbn={978-0-13-249020-7}"])
def test_valid_analytic_receipt_is_hashed_not_executed(analytic_sample, identifier):
    root, path, record, bibliography = analytic_sample
    bibliography.write_text('@book{SyntheticReference, '+identifier+'}\n', encoding="utf-8")
    assert validate_record(root, path)["theory_applicability"] == "independent-review-required"
    (path.parent / "symbolic.json").write_bytes(b"changed")
    with pytest.raises(ValueError, match="artifact hash"): validate_record(root, path)


def test_analytic_reference_rejects_missing_canonical_key(analytic_sample):
    root, path, record, bibliography = analytic_sample
    record["adoption"]["criterion"]["reference"]["bib_key"] = "MissingReference"
    seal_adoption(path, record)
    with pytest.raises(ValueError, match="analytic bibliography key"): validate_record(root, path)


def test_analytic_reference_rejects_url_only_entry(analytic_sample):
    root, path, record, bibliography = analytic_sample
    bibliography.write_text('@book{SyntheticReference, url={https://example.org/reference}}\n', encoding="utf-8")
    with pytest.raises(ValueError, match="analytic bibliography identifier"): validate_record(root, path)


def test_analytic_reference_rejects_invalid_isbn_checksum(analytic_sample):
    root, path, record, bibliography = analytic_sample
    bibliography.write_text('@book{SyntheticReference, isbn={9780132490208}}\n', encoding="utf-8")
    with pytest.raises(ValueError, match="checksum-valid ISBN-13"): validate_record(root, path)


@pytest.mark.parametrize("extra", ["LIMIT = 2 * 0.01\n", "if True:\n    LIMIT = 0.02\n"])
def test_reassigned_literal_is_not_safe_machine_verification(adoption_sample, extra):
    root, path, record = adoption_sample
    data = ("LIMIT = 0.01\n" + extra).encode()
    (path.parent / "binding.py").write_bytes(data)
    record["adoption"]["criterion"]["threshold_verification"] = {
        "mode": "literal-constants", "review_ref": "synthetic-review", "script": {
            "path": "binding.py", "sha256": hashlib.sha256(data).hexdigest()},
        "bindings": [{"constant": "LIMIT", "criterion_pointer": "/series/0/finest_threshold/value"}],
    }
    seal_adoption(path, record)
    with pytest.raises(ValueError, match="assignment count"): validate_record(root, path)


def test_p_then_h_chain_uses_separate_rates_not_a_mixed_axis(adoption_sample):
    root, path, record = adoption_sample
    c = record["adoption"]["criterion"]
    c["variables"] = {"p": "polynomial degree", "h": "mesh size"}
    parameters = ({"p": 1, "h": 1.}, {"p": 2, "h": 1.}, {"p": 2, "h": .5})
    for declared, actual, params in zip(c["levels"], record["adoption"]["result"]["levels"], parameters):
        declared["parameters"] = actual["parameters"] = params
    original = c["series"][0]
    p_series = json.loads(json.dumps(original))
    p_series.update(id="p", level_ids=["L1", "L2"], axis="p", direction="increasing",
                    rate_model="log-per-unit", minimum_rate=1., finest_threshold={"operator": "lt", "value": .02})
    h_series = json.loads(json.dumps(original))
    h_series.update(id="h", level_ids=["L2", "L3"])
    c["series"] = [p_series, h_series]
    seal_adoption(path, record)
    actual = validate_record(root, path)["observed_rates"]
    assert actual["p"] == pytest.approx([1.38629436112])
    assert actual["h"] == pytest.approx([2.])


@pytest.mark.parametrize("count", [1, 2])
def test_series_needs_at_least_two_levels_without_imposing_a_third(adoption_sample, count):
    root, path, record = adoption_sample
    c = record["adoption"]["criterion"]
    c["levels"] = c["levels"][:count]
    record["adoption"]["result"]["levels"] = record["adoption"]["result"]["levels"][:count]
    c["series"][0]["level_ids"] = [row["id"] for row in c["levels"]]
    c["series"][0]["finest_threshold"]["value"] = .02
    seal_adoption(path, record)
    if count == 1:
        with pytest.raises(ValueError): validate_record(root, path)
    else:
        assert validate_record(root, path)["adoption_claim"] is True


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
        "schema": "radia.validation-record.v3", "category": "milestone", "supersedes": [], "id": "principle/20261022_contract",
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


def test_rationale_returns_reverse_review_dependencies(adoption_sample):
    root, path, record = adoption_sample
    result = validate_record(root, path)
    edge = result["assumption_dependencies"]["reference.checked"]
    assert edge["record_id"] == record["id"]
    assert edge["recheck"] == record["adoption"]["rationale"]["assumptions"][0]["recheck"]
    assert result["rationale_validity"] == "independent-review-required"


@pytest.mark.parametrize("damage", ["missing", "duplicate", "unknown", "omitted", "unapproved-change"])
def test_rationale_invalid_or_unbound_dependencies_fail(adoption_sample, damage):
    root, path, record = adoption_sample
    a = record["adoption"]; rationale = a["rationale"]
    if damage == "missing":
        del a["rationale"]
        write(path, record)
    else:
        if damage == "duplicate": rationale["assumptions"].append(dict(rationale["assumptions"][0]))
        elif damage == "unknown": rationale["depends_on"] = ["missing.id"]
        elif damage == "omitted": rationale["depends_on"] = []
        else: rationale["scope"] = "Expanded without user approval"
        if damage == "unapproved-change": write(path, record)
        else: seal_adoption(path, record)
    with pytest.raises(ValueError): validate_record(root, path)


@pytest.fixture
def champion_sample(adoption_sample):
    import copy
    root, origin_path, original = adoption_sample
    record = copy.deepcopy(original)
    record.pop("adoption")
    record["id"] = "principle/20261023_champion"
    record["date"] = "2026-10-23"
    record["category"] = "accuracy-champion"
    folder = root / "validation_records" / record["id"]
    folder.mkdir(parents=True)
    (folder / "entry.py").write_bytes((origin_path.parent / "entry.py").read_bytes())
    levels = []
    for number, h, error in [(1, .25, .0025), (2, .125, .000625)]:
        fields = {"id": f"C{number}", "parameters": {"h": h}, "error": error,
                  "source_sha": record["source"]["repo_sha"], "mesh_sha256": "b" * 64,
                  "runtime_s": number * 2100.}
        payload = (json.dumps({"champion_level": fields}, sort_keys=True) + "\n").encode()
        name = f"level{number}.json"
        (folder / name).write_bytes(payload)
        levels.append({**fields, "artifact": {"path": name, "sha256": hashlib.sha256(payload).hexdigest()}})
    record["champion"] = {
        "adoption_record": {"id": original["id"], "sha256": hashlib.sha256(origin_path.read_bytes()).hexdigest()},
        "feature": original["adoption"]["feature"], "quantity": {key: original["adoption"]["criterion"]["quantity"][key] for key in ("name", "unit")},
        "selected_level_id": "C2", "selection_reason": "Synthetic finest level, not a real champion",
        "levels": levels}
    path = folder / "record.json"
    seal_champion(path, record)
    return root, path, record


def seal_champion(path, record):
    champion = record["champion"]
    data = (json.dumps({"champion_result": {k: v for k, v in champion.items() if k != "adoption_record"}},
                       sort_keys=True) + "\n").encode()
    (path.parent / "result.json").write_bytes(data)
    record["result"]["sha256"] = hashlib.sha256(data).hexdigest()
    write(path, record)


@pytest.mark.parametrize("category", ["accuracy-champion", "speed", "scale"])
def test_champion_records_table_and_origin_without_replay(champion_sample, category):
    root, path, record = champion_sample
    record["category"] = category
    write(path, record)
    result = validate_record(root, path)
    assert result["champion_category"] == category
    assert result["adoption_origin"] == record["champion"]["adoption_record"]["id"]
    assert result["champion_levels"] == 2
    assert result["adoption_claim"] is False
    assert result["timing_used_as_gate"] is False and result["numerical_replay"] is False


@pytest.mark.parametrize("damage", ["missing-block", "category", "origin-hash", "missing-origin", "feature",
                                    "level-hash", "level-value", "selected", "endpoint", "source"])
def test_champion_identity_mismatches_fail(champion_sample, damage):
    root, path, record = champion_sample
    c = record["champion"]
    if damage == "missing-block": del record["champion"]
    elif damage == "category": record["category"] = "unrelated"
    elif damage == "origin-hash": c["adoption_record"]["sha256"] = "0" * 64
    elif damage == "missing-origin": c["adoption_record"]["id"] = "missing/20261023_origin"
    elif damage == "feature": c["feature"] = "different.production_entry"
    elif damage == "level-hash": c["levels"][0]["artifact"]["sha256"] = "0" * 64
    elif damage == "level-value": c["levels"][0]["error"] = .003
    elif damage == "selected": c["selected_level_id"] = "missing"
    elif damage == "endpoint": c["selected_level_id"] = "C1"
    else: c["levels"][-1]["source_sha"] = "0" * 40
    if damage == "missing-block": write(path, record)
    else: seal_champion(path, record)
    with pytest.raises(ValueError): validate_record(root, path)


def test_changed_code_requires_stale_but_preserves_champion(champion_sample):
    root, path, record = champion_sample
    (root / "physics.py").write_text("VALUE = 2\n")
    git(root, "add", "physics.py"); git(root, "commit", "-m", "Change covered source")
    with pytest.raises(ValueError, match="current after covered content changed"):
        validate_record(root, path)
    origin_path = root / "validation_records" / record["champion"]["adoption_record"]["id"] / "record.json"
    origin = json.loads(origin_path.read_text())
    origin["status"] = "stale"; origin["stale_reason"] = "Covered source changed"
    write(origin_path, origin)
    record["champion"]["adoption_record"]["sha256"] = hashlib.sha256(origin_path.read_bytes()).hexdigest()
    record["status"] = "stale"; record["stale_reason"] = "Covered source changed"
    seal_champion(path, record)
    result = validate_record(root, path)
    assert result["status"] == result["adoption_origin_status"] == "stale"
    assert origin_path.is_file() and path.is_file()


def test_one_off_extension_can_have_a_single_new_level(champion_sample):
    root, path, record = champion_sample
    record["category"] = "scale"
    record["champion"]["levels"] = record["champion"]["levels"][-1:]
    seal_champion(path, record)
    result = validate_record(root, path)
    assert result["champion_levels"] == 1
    assert result["adoption_origin"] == record["champion"]["adoption_record"]["id"]


def test_champion_receipt_cannot_replace_a_number_with_boolean(champion_sample):
    root, path, record = champion_sample
    level = record["champion"]["levels"][0]
    level["runtime_s"] = 0.0
    fields = {k: v for k, v in level.items() if k != "artifact"}
    fields["runtime_s"] = False
    data = (json.dumps({"champion_level": fields}) + "\n").encode()
    (path.parent / level["artifact"]["path"]).write_bytes(data)
    level["artifact"]["sha256"] = hashlib.sha256(data).hexdigest()
    seal_champion(path, record)
    with pytest.raises(ValueError, match="champion level artifact values"):
        validate_record(root, path)


@pytest.mark.parametrize("change", ["feature", "record-id"])
def test_sealed_decision_cannot_be_reassigned_to_another_scope(adoption_sample, change):
    root, path, record = adoption_sample
    if change == "feature": record["adoption"]["feature"] = "different.production_entry"
    else:
        old = path.parent
        record["id"] = "principle/20261023_reassigned"
        path = root / "validation_records" / record["id"] / "record.json"
        path.parent.mkdir(parents=True)
        for name in ("entry.py", "result.json"):
            (path.parent / name).write_bytes((old / name).read_bytes())
    write(path, record)
    with pytest.raises(ValueError, match="approved feature/record scope"):
        validate_record(root, path)


@pytest.mark.parametrize("change", ["artifact", "review", "condition", "approval-date"])
def test_adoption_rejects_unresolved_refs_or_postdated_consent(adoption_sample, change):
    root, path, record = adoption_sample
    premise = record["adoption"]["rationale"]["assumptions"][0]
    if change == "artifact": premise["evidence_refs"] = ["artifact:missing"]
    elif change == "review": premise["evidence_refs"] = ["review:undeclared"]
    elif change == "condition": premise["invalidated_by"] = ["undeclared.failure"]
    else: record["adoption"]["approval"]["approved_on"] = "2099-01-01"
    seal_adoption(path, record)
    with pytest.raises(ValueError): validate_record(root, path)


def test_resolved_evidence_reports_hashed_artifacts_and_review_limits(adoption_sample):
    root, path, record = adoption_sample
    record["adoption"]["rationale"]["assumptions"][0]["evidence_refs"].append("review:synthetic-review-20261023")
    seal_adoption(path, record)
    refs = validate_record(root, path)["assumption_dependencies"]["reference.checked"]["resolved_evidence"]
    assert refs["artifact:result"]["sha256"] == record["result"]["sha256"]
    assert refs["review:synthetic-review-20261023"]["verification"] == "independent-review-required"


@pytest.mark.parametrize("cited", [True, False])
def test_surrogate_assumes_must_cite_this_decision_without_import(adoption_sample, cited):
    root, path, record = adoption_sample
    test = root / "tests/test_surrogate.py"
    test.parent.mkdir(exist_ok=True)
    identity = record["id"] if cited else "different/20261022_decision"
    test.write_text('"""Guards: Declared route changes.\nAssumes: surrogate test of adoption record '
                    + identity + ' (convergence level 1)\n"""\nraise RuntimeError("must not import")\n')
    record["adoption"]["protected_by"] = ["tests/test_surrogate.py"]
    seal_adoption(path, record)
    if cited: assert validate_record(root, path)["adoption_claim"] is True
    else:
        with pytest.raises(ValueError, match="lock Assumes citation"): validate_record(root, path)


@pytest.mark.parametrize("direction", ["back", "forward"])
def test_supersession_must_be_reciprocal(sample, direction):
    import copy
    root, path, record = sample
    previous = copy.deepcopy(record)
    previous.update(id="principle/20261021_previous", date="2026-10-21", superseded_by=record["id"])
    target = path.parent.parent / "20261021_previous/record.json"
    write(target, previous)
    for name in ("entry.py", "result.json"):
        (target.parent / name).write_bytes((path.parent / name).read_bytes())
    record["supersedes"] = [{"id": previous["id"], "relation": "corrects", "reason": "Synthetic relation"}]
    if direction == "back": previous["superseded_by"] = None; write(target, previous)
    else: record["supersedes"] = []
    write(path, record)
    with pytest.raises(ValueError, match="supersedes backlink|superseded_by forward link"):
        validate_record(root, path if direction == "back" else target)


def test_category_is_a_declared_record_role(sample):
    root, path, record = sample
    record["category"] = "fastest-ever"
    write(path, record)
    with pytest.raises(ValueError, match="enum"): validate_record(root, path)


def test_declared_artifact_approval_and_record_references_resolve(adoption_sample):
    import copy
    root, path, record = adoption_sample
    other = copy.deepcopy(record)
    other.pop("adoption")
    other["category"] = "milestone"; other["id"] = "principle/20261021_evidence"
    target = root / "validation_records" / other["id"] / "record.json"
    write(target, other)
    refs = ["artifact:entry", "approval:" + record["adoption"]["approval"]["approval_ref"], "record:" + other["id"]]
    record["adoption"]["rationale"]["assumptions"][0]["evidence_refs"] = refs
    seal_adoption(path, record)
    actual = validate_record(root, path)["assumption_dependencies"]["reference.checked"]["resolved_evidence"]
    assert actual["artifact:entry"]["sha256"] == record["entry"]["script_sha256"]
    assert actual[refs[1]]["verification"] == "independent-review-required"
    assert actual[refs[2]]["sha256"] == hashlib.sha256(target.read_bytes()).hexdigest()


def reseal_champion_levels(path, record):
    """Synthetic receipts only: isolate metadata guards from hash-mismatch guards."""
    for level in record["champion"]["levels"]:
        fields = {k: v for k, v in level.items() if k != "artifact"}
        data = (json.dumps({"champion_level": fields}, sort_keys=True) + "\n").encode()
        (path.parent / level["artifact"]["path"]).write_bytes(data)
        level["artifact"]["sha256"] = hashlib.sha256(data).hexdigest()
    seal_champion(path, record)


def test_current_champion_with_stale_origin_fails_repository_gate(champion_sample):
    root, path, record = champion_sample
    target = root / "validation_records" / record["champion"]["adoption_record"]["id"] / "record.json"
    origin = json.loads(target.read_text())
    origin["status"] = "stale"; origin["stale_reason"] = "Reviewed premise invalidated"
    write(target, origin)
    record["champion"]["adoption_record"]["sha256"] = hashlib.sha256(target.read_bytes()).hexdigest()
    seal_champion(path, record)
    with pytest.raises(ValueError, match="champion origin status"): validate_tree(root)


@pytest.mark.parametrize("change", ["quantity", "parameter-keys"])
def test_champion_cannot_change_the_adopted_quantity_or_variables(champion_sample, change):
    root, path, record = champion_sample
    if change == "quantity": record["champion"]["quantity"] = {"name": "Torque", "unit": "N m"}
    else:
        for level in record["champion"]["levels"]: level["parameters"] = {"zz": 1.0}
    reseal_champion_levels(path, record)
    with pytest.raises(ValueError, match="champion quantity|champion parameter keys"): validate_record(root, path)


def test_failed_outcome_cannot_claim_champion(champion_sample):
    root, path, record = champion_sample
    record["outcome"] = "fail"
    write(path, record)
    with pytest.raises(ValueError, match="champion outcome"): validate_record(root, path)


def test_accuracy_champion_cannot_record_increasing_errors(champion_sample):
    root, path, record = champion_sample
    # Both are below the origin's finest error; only the increasing trend is wrong.
    for level, error in zip(record["champion"]["levels"], (.0005, .0009)): level["error"] = error
    reseal_champion_levels(path, record)
    with pytest.raises(ValueError, match="champion error trend"): validate_record(root, path)


def test_adoption_cannot_cite_itself_as_evidence(adoption_sample):
    root, path, record = adoption_sample
    record["adoption"]["rationale"]["assumptions"][0]["evidence_refs"] = ["record:" + record["id"]]
    seal_adoption(path, record)
    with pytest.raises(ValueError, match="self evidence reference"): validate_record(root, path)


def test_adoption_requires_a_surrogate_lock(adoption_sample):
    root, path, record = adoption_sample
    record["adoption"]["protected_by"] = []
    write(path, record)
    with pytest.raises(ValueError, match="protected_by"): validate_record(root, path)


def test_surrogate_lock_must_be_a_python_module(adoption_sample):
    root, path, record = adoption_sample
    lock = root / "tests/lock.txt"
    lock.write_bytes((root / record["adoption"]["protected_by"][0]).read_bytes())
    record["adoption"]["protected_by"] = ["tests/lock.txt"]
    write(path, record)
    with pytest.raises(ValueError, match="lock extension"): validate_record(root, path)


def test_accuracy_champion_cannot_be_worse_than_origin(champion_sample):
    root, path, record = champion_sample
    for level, error in zip(record["champion"]["levels"], (.9, .5)): level["error"] = error
    reseal_champion_levels(path, record)
    with pytest.raises(ValueError, match="champion origin error bound"): validate_record(root, path)


def test_recorded_finest_still_checks_rate_and_preserves_large_error(adoption_sample):
    root, path, record = adoption_sample
    a = record["adoption"]
    a["criterion"]["series"][0]["finest_threshold"] = {
        "operator": "recorded-only", "reason": "Synthetic user approval excludes finest accuracy only"}
    for row in a["result"]["levels"]:
        row["error"] *= 1000
    seal_adoption(path, record)
    assert validate_record(root, path)["adoption_claim"] is True
    a["criterion"]["series"][0]["minimum_rate"] = 2.5
    seal_adoption(path, record)
    with pytest.raises(ValueError, match="series verdict"):
        validate_record(root, path)


@pytest.mark.parametrize("change", ["no-reason", "hidden-threshold", "numeric-no-value"])
def test_ambiguous_finest_gate_is_rejected(adoption_sample, change):
    root, path, record = adoption_sample
    gate = {"operator": "recorded-only", "reason": "Synthetic approval"}
    if change == "no-reason": del gate["reason"]
    elif change == "hidden-threshold": gate["value"] = 1e20
    else: gate = {"operator": "le"}
    record["adoption"]["criterion"]["series"][0]["finest_threshold"] = gate
    seal_adoption(path, record)
    with pytest.raises(ValueError):
        validate_record(root, path)


def test_missing_level_timing_is_explicit_not_zero(adoption_sample):
    root, path, record = adoption_sample
    row = record["adoption"]["result"]["levels"][0]
    row["runtime_s"] = None
    row["runtime_note"] = "Historical campaign recorded only its aggregate elapsed time"
    seal_adoption(path, record)
    assert validate_record(root, path)["adoption_claim"] is True
    del row["runtime_note"]
    seal_adoption(path, record)
    with pytest.raises(ValueError, match="runtime_note"):
        validate_record(root, path)
