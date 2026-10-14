"""Solver-neutral evidence gate for a reversible local uninstall."""

from __future__ import annotations

import hashlib
import json
import re
from collections.abc import Mapping, Sequence
from typing import Any


SCHEMA = "radia.solver-uninstall-safety-evidence.v3"
GATE_SCHEMA = "radia.solver-uninstall-safety-gate.v3"
REPLACEMENT_SCHEMA = "radia.validation-evidence-bundle.v1"
PRODUCTION_PROOF_SCHEMA = "radia.validation.production-replacement-proof.v1"
REQUIRED_PRODUCTION_PROOFS = {
    "motor_dual_lane",
    "native_motor_angle_family",
}
REQUIRED_PRODUCTION_CHECKS = {
    "motor_dual_lane": {
        "both_lanes_solved",
        "shared_identity_matches",
        "angle_grids_complete",
        "torque_waveforms_nonconstant",
        "correlation_gate_pass",
        "dual_mcp_gate_pass",
    },
    "native_motor_angle_family": {
        "standalone_batch",
        "all_tests_passed",
        "mex_and_sources_hashed",
        "periodic_motor_handle_tested",
        "simulink_compile_tested",
        "foreign_openmp_runtime_isolated",
    },
}
REQUIRED_SCAN_SCOPES = {
    "active_processes",
    "automation_scripts",
    "converter",
    "environment",
    "file_associations",
    "mcp_config",
    "public_solver_repo",
    "registry",
    "scheduled_tasks",
    "shortcuts",
    "source_mcp",
}
_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")


def _sha(value: object) -> str:
    payload = dict(value) if isinstance(value, Mapping) else value
    if isinstance(payload, dict):
        payload.pop("evidence_payload_sha256", None)
    raw = json.dumps(
        payload, ensure_ascii=True, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def _rows(value: object) -> list[Mapping[str, Any]]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        return []
    return [row for row in value if isinstance(row, Mapping)]


def _validate_production_replacement_proof(record: object) -> bool:
    if not isinstance(record, Mapping):
        return False
    raw = record.get("artifact_json")
    declared_sha = str(record.get("artifact_sha256", "")).lower()
    if not isinstance(raw, str) or not raw.strip():
        return False
    actual_sha = hashlib.sha256(raw.encode("utf-8")).hexdigest()
    if not _SHA256_RE.fullmatch(declared_sha) or declared_sha != actual_sha:
        return False
    try:
        artifact = json.loads(raw)
    except json.JSONDecodeError:
        return False
    if not isinstance(artifact, Mapping):
        return False
    versions = artifact.get("execution_version")
    proofs = artifact.get("proofs")
    if not isinstance(versions, Mapping) or not isinstance(proofs, Mapping):
        return False
    if set(proofs) != REQUIRED_PRODUCTION_PROOFS:
        return False
    if not all(str(versions.get(key, "")).strip() for key in ("radia", "matlab", "simulink")):
        return False

    for proof_id, required_checks in REQUIRED_PRODUCTION_CHECKS.items():
        proof = proofs.get(proof_id)
        if not isinstance(proof, Mapping) or proof.get("status") != "pass":
            return False
        checks = proof.get("checks")
        hashes = proof.get("artifact_sha256_by_role")
        if not isinstance(checks, Mapping) or not isinstance(hashes, Mapping):
            return False
        if not required_checks.issubset(checks):
            return False
        if not all(checks.get(name) is True for name in required_checks):
            return False
        if not hashes or not all(
            _SHA256_RE.fullmatch(str(value).lower()) for value in hashes.values()
        ):
            return False

    payload = dict(artifact)
    declared_payload_sha = str(payload.pop("proof_payload_sha256", "")).lower()
    return bool(
        artifact.get("schema") == PRODUCTION_PROOF_SCHEMA
        and artifact.get("status") == "pass"
        and artifact.get("pass") is True
        and str(artifact.get("executed_at_utc", "")).strip()
        and _SHA256_RE.fullmatch(declared_payload_sha)
        and declared_payload_sha == _sha(payload)
    )


def validate_solver_uninstall_safety_evidence(
    evidence: Mapping[str, Any],
) -> dict[str, Any]:
    """Validate content-addressed uninstall evidence without reading local paths."""

    if not isinstance(evidence, Mapping):
        raise ValueError("evidence must be an object")
    archives = _rows(evidence.get("archive_snapshots"))
    dependencies = evidence.get("dependency_scan")
    rollback = evidence.get("rollback_installer")
    install_roots = evidence.get("installation_roots")
    versions = evidence.get("execution_version")
    replacement = evidence.get("replacement_evidence_bundle")
    production_proof = evidence.get("production_replacement_proof")
    if not isinstance(dependencies, Mapping):
        dependencies = {}
    if not isinstance(rollback, Mapping):
        rollback = {}
    if not isinstance(install_roots, Mapping):
        install_roots = {}
    if not isinstance(versions, Mapping):
        versions = {}
    if not isinstance(replacement, Mapping):
        replacement = {}

    archive_roles = [str(row.get("role", "")).strip() for row in archives]
    archive_checks = [
        bool(
            str(row.get("role", "")).strip()
            and int(row.get("file_count", 0) or 0) > 0
            and int(row.get("archive_size_bytes", 0) or 0) > 0
            and _SHA256_RE.fullmatch(str(row.get("root_manifest_sha256", "")).lower())
            and _SHA256_RE.fullmatch(str(row.get("archive_sha256", "")).lower())
            and str(row.get("verified_archive_sha256", "")).lower()
            == str(row.get("archive_sha256", "")).lower()
            and row.get("stable_storage") is True
            and row.get("archive_verified") is True
        )
        for row in archives
    ]
    scope_rows = _rows(dependencies.get("scopes"))
    scope_names = {str(row.get("category", "")).strip() for row in scope_rows}
    scope_checks = [
        bool(
            str(row.get("category", "")).strip()
            and row.get("scanned") is True
            and _SHA256_RE.fullmatch(str(row.get("result_sha256", "")).lower())
            and str(row.get("verified_result_sha256", "")).lower()
            == str(row.get("result_sha256", "")).lower()
        )
        for row in scope_rows
    ]

    checks = {
        "schema": evidence.get("schema") == SCHEMA,
        "execution_timestamp_recorded": bool(str(evidence.get("executed_at_utc", "")).strip()),
        "execution_version_complete": all(
            str(versions.get(key, "")).strip()
            for key in ("producer", "producer_version", "python")
        ),
        "archive_snapshots_present": bool(archives),
        "archive_roles_unique": bool(archive_roles)
        and len(archive_roles) == len(set(archive_roles)),
        "archive_snapshots_verified": bool(archive_checks) and all(archive_checks),
        "rollback_installer_verified": bool(
            int(rollback.get("size_bytes", 0) or 0) > 0
            and _SHA256_RE.fullmatch(str(rollback.get("sha256", "")).lower())
            and str(rollback.get("verified_sha256", "")).lower()
            == str(rollback.get("sha256", "")).lower()
            and str(rollback.get("version", "")).strip()
            and rollback.get("stable_storage") is True
            and rollback.get("hash_verified") is True
            and rollback.get("temporary_storage") is False
        ),
        "dependency_scopes_complete": REQUIRED_SCAN_SCOPES.issubset(scope_names),
        "dependency_scope_results_hashed": bool(scope_checks) and all(scope_checks),
        "no_unresolved_runtime_dependencies": dependencies.get("unresolved_dependencies") == [],
        "no_active_solver_processes": dependencies.get("active_processes") == [],
        "live_compatibility_disabled_by_default": dependencies.get(
            "live_compatibility_default_enabled"
        )
        is False,
        "dependency_scan_complete": dependencies.get("scan_complete") is True
        and bool(_SHA256_RE.fullmatch(str(dependencies.get("scan_sha256", "")).lower())),
        "dependency_scan_timestamp_matches_evidence": bool(
            str(evidence.get("executed_at_utc", "")).strip()
        )
        and str(dependencies.get("executed_at_utc", "")).strip()
        == str(evidence.get("executed_at_utc", "")).strip(),
        "installation_roots_accounted_for": int(
            install_roots.get("discovered_root_count", 0) or 0
        )
        > 0
        and int(install_roots.get("planned_removal_root_count", 0) or 0)
        == int(install_roots.get("discovered_root_count", 0) or 0),
        "root_mismatch_has_explicit_plan": (
            install_roots.get("mismatch_detected") is False
            or install_roots.get("mismatch_plan_complete") is True
        ),
        "replacement_capabilities_complete": bool(
            replacement.get("schema") == REPLACEMENT_SCHEMA
            and replacement.get("status") == "complete"
            and replacement.get("contract_valid") is True
            and replacement.get("retirement_ready") is True
            and replacement.get("solver_uninstall_performed") is False
            and isinstance(replacement.get("required_capabilities"), Sequence)
            and not isinstance(replacement.get("required_capabilities"), (str, bytes))
            and bool(replacement.get("required_capabilities"))
            and sorted(replacement.get("accepted_capabilities", []))
            == sorted(replacement.get("required_capabilities", []))
            and replacement.get("missing_capabilities") == []
            and _SHA256_RE.fullmatch(
                str(replacement.get("evidence_bundle_sha256", "")).lower()
            )
        ),
        "production_replacement_proof_complete": (
            _validate_production_replacement_proof(production_proof)
        ),
        "uninstall_not_performed": evidence.get("solver_uninstall_performed") is False,
        "evidence_payload_sha256": str(
            evidence.get("evidence_payload_sha256", "")
        ).lower()
        == _sha(evidence),
    }
    passed = all(checks.values())
    return {
        "schema": GATE_SCHEMA,
        "status": "accepted" if passed else "rejected",
        "pass": passed,
        "ready_for_explicit_uninstall_approval": passed,
        "solver_uninstall_performed": False,
        "checks": checks,
        "archive_roles": sorted(archive_roles),
        "scan_scope_categories": sorted(scope_names),
        "evidence_payload_sha256": evidence.get("evidence_payload_sha256"),
        "retirement_capabilities": [
            "archive_snapshot",
            "replacement_capability_bundle",
            "production_replacement_proof",
            "live_dependency_scan_clean",
            "rollback_installer_preserved",
        ]
        if passed
        else [],
    }
