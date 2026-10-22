"""Adapt immutable nested records to common zero-tolerance regression locks."""
import hashlib
import json
from pathlib import Path

from tools.regression_lock import assert_regression, flatten_json_record
from validation_test.radia_mcp.netgen_vol_examples import ROOT

FIXTURES = ROOT / "tests/mcp_integration/fixtures/netgen_vol_analytic"


def approved_case(name):
    manifest = json.loads((FIXTURES / "manifest.json").read_text(encoding="utf-8"))
    reference = manifest["cases"][name]
    filename = reference["file"]
    if Path(filename).name != filename:
        raise ValueError("Expected a fixture filename, not a path")
    path = FIXTURES / filename
    raw = path.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    if digest != reference["sha256"]:
        raise AssertionError(f"{name}: actual fixture SHA256={digest}; expected={reference['sha256']}")
    provenance = {"case": name, "reference_source_sha": manifest["reference_source_sha"],
                  "fixture_sha256": digest, "git_blob": reference["git_blob"],
                  "original_path": reference["original_path"],
                  "reference_conditions": manifest["conditions"]}
    return path, json.loads(raw), provenance


def assert_summary_regression(name, summary):
    _, reference, provenance = approved_case(name)
    # Preserve the existing post-build JSON normalization (including tuples).
    # Do not move this before an example's own _json_clean logic.
    summary = json.loads(json.dumps(summary, allow_nan=False))
    values, structure = flatten_json_record(reference)
    actual, actual_structure = flatten_json_record(summary)
    baseline = {"schema": "radia.regression-lock.v1",
                "source_sha": provenance["reference_source_sha"],
                "conditions": {"provenance": provenance, "structure": structure,
                               "historical_runtime": "unknown; virtual wrapper runtime_s=0 is not a measurement"},
                # In-memory adapter, never written over the frozen record.
                "runtime_s": 0,
                "quantities": {key: {"expected": value, "atol": 0, "rtol": 0,
                                     "unit": "original record field units (unchanged)"}
                               for key, value in values.items()}}
    assert_regression(actual, baseline,
                      conditions={**baseline["conditions"], "structure": actual_structure})
