"""Solver-neutral evidence gate for a reversible local uninstall."""

from __future__ import annotations

import hashlib
import json
import re
from collections.abc import Mapping, Sequence
from typing import Any


SCHEMA = "radia.solver-uninstall-safety-evidence.v2"
GATE_SCHEMA = "radia.solver-uninstall-safety-gate.v2"
REPLACEMENT_SCHEMA = "radia.validation-evidence-bundle.v1"
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
            "live_dependency_scan_clean",
            "rollback_installer_preserved",
        ]
        if passed
        else [],
    }
