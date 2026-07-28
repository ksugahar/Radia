from __future__ import annotations

import hashlib
import json

import pytest

from radia_mcp.radia_ngsolve.validation_evidence import validate_evidence_bundle


def _artifact(*, passed: bool = True, source: bool = True) -> str:
    lane = {
        "owner": "owner",
        "capability_gain": "new decision",
        "positive_probe": "positive",
        "negative_probe": "negative",
        "protocol_probe": "initialize/list/call",
        "verification": "tests passed",
        "commit": "a" * 40,
    }
    payload = {
        "schema": "example.learning.v1",
        "created_at_utc": "2026-07-28T00:00:00Z",
        "versions": {"solver": "1.0"},
        "pass": passed,
        "retirement_capabilities": ["planar_dc", "axisymmetric_dc"],
        "mcp_balance": {
            "policy": "equal_capability_gain_v1",
            "status": "verified",
            "capability_id": "example",
            "public": lane,
            "source_tool": lane if source else {},
        },
    }
    return json.dumps(payload, sort_keys=True)


def _record(name: str, capabilities: list[str], raw: str | None = None) -> dict:
    if raw is None:
        payload = json.loads(_artifact())
        payload["retirement_capabilities"] = capabilities
        raw = json.dumps(payload, sort_keys=True)
    return {
        "artifact_id": name,
        "artifact_json": raw,
        "artifact_sha256": hashlib.sha256(raw.encode()).hexdigest(),
        "capabilities": capabilities,
    }


def _packet(records: list[dict], required: list[str] | None = None) -> dict:
    return {
        "inventory_sha256": "b" * 64,
        "required_capabilities": required or ["planar_dc", "axisymmetric_dc"],
        "artifacts": records,
    }


def test_accepts_balanced_artifact_and_explicit_multiple_capabilities():
    result = validate_evidence_bundle(
        _packet([_record("a", ["planar_dc", "axisymmetric_dc"])])
    )
    assert result["retirement_ready"] is True
    assert result["accepted_capabilities"] == ["axisymmetric_dc", "planar_dc"]


def test_canonical_digest_is_independent_of_required_and_artifact_order():
    a = _record("a", ["planar_dc"])
    b = _record("b", ["axisymmetric_dc"])
    first = validate_evidence_bundle(_packet([a, b]))
    second = validate_evidence_bundle(
        _packet([b, a], ["axisymmetric_dc", "planar_dc"])
    )
    assert first["evidence_bundle_sha256"] == second["evidence_bundle_sha256"]


def test_reports_exact_missing_coverage():
    result = validate_evidence_bundle(_packet([_record("a", ["planar_dc"])]))
    assert result["contract_valid"] is True
    assert result["retirement_ready"] is False
    assert result["missing_capabilities"] == ["axisymmetric_dc"]


def test_rejects_payload_digest_mismatch():
    row = _record("a", ["planar_dc"])
    row["artifact_sha256"] = "0" * 64
    with pytest.raises(ValueError, match="does not match"):
        validate_evidence_bundle(_packet([row]))


@pytest.mark.parametrize(
    "raw, message",
    [
        (_artifact(passed=False), "pass=true"),
        (_artifact(source=False), "source lane"),
    ],
)
def test_rejects_failed_or_single_lane_artifacts(raw: str, message: str):
    with pytest.raises(ValueError, match=message):
        validate_evidence_bundle(_packet([_record("a", ["planar_dc"], raw)]))


def test_rejects_conflicting_duplicate_capability():
    with pytest.raises(ValueError, match="claimed by both"):
        validate_evidence_bundle(
            _packet(
                [
                    _record("a", ["planar_dc"]),
                    _record("b", ["planar_dc", "axisymmetric_dc"]),
                ]
            )
        )
