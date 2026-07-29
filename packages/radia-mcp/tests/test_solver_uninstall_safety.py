from __future__ import annotations

import copy

from radia_mcp.fem.uninstall_safety import (
    REQUIRED_SCAN_SCOPES,
    _sha,
    validate_solver_uninstall_safety_evidence,
)


def _evidence() -> dict:
    digest = "a" * 64
    payload = {
        "schema": "radia.solver-uninstall-safety-evidence.v2",
        "executed_at_utc": "2026-07-29T00:00:00Z",
        "execution_version": {
            "producer": "test",
            "producer_version": "1",
            "python": "3.12",
        },
        "archive_snapshots": [
            {
                "role": "runtime",
                "file_count": 2,
                "root_manifest_sha256": digest,
                "archive_sha256": digest,
                "verified_archive_sha256": digest,
                "archive_size_bytes": 100,
                "stable_storage": True,
                "archive_verified": True,
            }
        ],
        "rollback_installer": {
            "version": "example",
            "size_bytes": 100,
            "sha256": digest,
            "verified_sha256": digest,
            "stable_storage": True,
            "hash_verified": True,
            "temporary_storage": False,
        },
        "dependency_scan": {
            "executed_at_utc": "2026-07-29T00:00:00Z",
            "scopes": [
                {
                    "category": category,
                    "scanned": True,
                    "result_sha256": digest,
                    "verified_result_sha256": digest,
                }
                for category in sorted(REQUIRED_SCAN_SCOPES)
            ],
            "unresolved_dependencies": [],
            "active_processes": [],
            "live_compatibility_default_enabled": False,
            "scan_complete": True,
            "scan_sha256": digest,
        },
        "installation_roots": {
            "discovered_root_count": 2,
            "planned_removal_root_count": 2,
            "mismatch_detected": True,
            "mismatch_plan_complete": True,
        },
        "replacement_evidence_bundle": {
            "schema": "radia.validation-evidence-bundle.v1",
            "status": "complete",
            "contract_valid": True,
            "retirement_ready": True,
            "solver_uninstall_performed": False,
            "required_capabilities": ["planar_dc", "axisymmetric_dc"],
            "accepted_capabilities": ["axisymmetric_dc", "planar_dc"],
            "missing_capabilities": [],
            "evidence_bundle_sha256": "d" * 64,
        },
        "solver_uninstall_performed": False,
    }
    payload["evidence_payload_sha256"] = _sha(payload)
    return payload


def test_accepts_complete_reversible_uninstall_evidence() -> None:
    result = validate_solver_uninstall_safety_evidence(_evidence())
    assert result["status"] == "accepted"
    assert result["ready_for_explicit_uninstall_approval"] is True
    assert result["solver_uninstall_performed"] is False
    assert result["retirement_capabilities"] == [
        "archive_snapshot",
        "replacement_capability_bundle",
        "live_dependency_scan_clean",
        "rollback_installer_preserved",
    ]


def test_rejects_stale_archive_and_temporary_rollback() -> None:
    payload = _evidence()
    payload["archive_snapshots"][0]["archive_verified"] = False
    payload["archive_snapshots"][0]["archive_sha256"] = "b" * 64
    payload["rollback_installer"]["stable_storage"] = False
    payload["rollback_installer"]["temporary_storage"] = True
    payload["evidence_payload_sha256"] = _sha(payload)
    result = validate_solver_uninstall_safety_evidence(payload)
    assert result["status"] == "rejected"
    assert result["checks"]["archive_snapshots_verified"] is False
    assert result["checks"]["rollback_installer_verified"] is False


def test_rejects_unresolved_dependency_active_process_and_missing_scope() -> None:
    payload = _evidence()
    payload["dependency_scan"]["unresolved_dependencies"] = ["automation/caller.py"]
    payload["dependency_scan"]["active_processes"] = [{"pid": 7}]
    payload["dependency_scan"]["scopes"] = payload["dependency_scan"]["scopes"][:-1]
    payload["evidence_payload_sha256"] = _sha(payload)
    result = validate_solver_uninstall_safety_evidence(payload)
    assert result["checks"]["dependency_scopes_complete"] is False
    assert result["checks"]["no_unresolved_runtime_dependencies"] is False
    assert result["checks"]["no_active_solver_processes"] is False


def test_rejects_unplanned_root_mismatch_and_stale_payload() -> None:
    payload = _evidence()
    payload["installation_roots"]["mismatch_plan_complete"] = False
    stale = payload["evidence_payload_sha256"]
    result = validate_solver_uninstall_safety_evidence(payload)
    assert result["checks"]["root_mismatch_has_explicit_plan"] is False
    assert result["checks"]["evidence_payload_sha256"] is False
    assert stale == payload["evidence_payload_sha256"]


def test_rejects_dependency_scan_from_another_execution() -> None:
    payload = _evidence()
    payload["dependency_scan"]["executed_at_utc"] = "2026-07-28T00:00:00Z"
    payload["evidence_payload_sha256"] = _sha(payload)

    result = validate_solver_uninstall_safety_evidence(payload)

    assert result["status"] == "rejected"
    assert result["checks"]["dependency_scan_timestamp_matches_evidence"] is False


def test_rejects_claim_that_uninstall_already_ran() -> None:
    payload = copy.deepcopy(_evidence())
    payload["solver_uninstall_performed"] = True
    payload["evidence_payload_sha256"] = _sha(payload)
    result = validate_solver_uninstall_safety_evidence(payload)
    assert result["checks"]["uninstall_not_performed"] is False
    assert result["solver_uninstall_performed"] is False


def test_rejects_archive_only_evidence_without_complete_replacement_bundle() -> None:
    payload = _evidence()
    payload["replacement_evidence_bundle"]["status"] = "incomplete"
    payload["replacement_evidence_bundle"]["retirement_ready"] = False
    payload["replacement_evidence_bundle"]["accepted_capabilities"] = [
        "planar_dc"
    ]
    payload["replacement_evidence_bundle"]["missing_capabilities"] = [
        "axisymmetric_dc"
    ]
    payload["evidence_payload_sha256"] = _sha(payload)

    result = validate_solver_uninstall_safety_evidence(payload)

    assert result["status"] == "rejected"
    assert result["ready_for_explicit_uninstall_approval"] is False
    assert result["checks"]["replacement_capabilities_complete"] is False
