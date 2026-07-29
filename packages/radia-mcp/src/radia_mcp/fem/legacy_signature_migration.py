"""Closed-world gate for solver-neutral signature-migration evidence."""

from __future__ import annotations

import hashlib
import json
import re
from collections.abc import Mapping
from datetime import datetime
from typing import Any


SCHEMA = "radia.legacy-signature-migration-evidence.v1"
GATE_SCHEMA = "radia.legacy-signature-migration-gate.v1"
_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
_TOKEN_RE = re.compile(r"^[A-Za-z0-9._-]+$")
_TOP_LEVEL_KEYS = {
    "schema",
    "executed_at_utc",
    "execution_version",
    "solver_launched",
    "live_dependency_required",
    "signature_count",
    "assigned_signature_count",
    "lane_counts",
    "requirement_counts",
    "native_emitter_candidate_count",
    "false_solver_ready_count",
    "identity_mismatch_count",
    "conversion_error_count",
    "script_compile_failure_count",
    "signature_catalog_sha256",
    "evidence_payload_sha256",
}
_VERSION_KEYS = {"producer", "producer_version", "solver"}
_LANES = {
    "planar_magnetostatic_age",
    "planar_hcurl_eddy_bubble",
    "axifem_henrotte",
    "axifem_eddy",
    "scalar_h1_electrostatic",
}
_REQUIREMENTS = {
    "age_motion_interface",
    "boundary_operator_identity",
    "circuit_excitation_identity",
    "curved_profile_geometry",
    "kelvin_open_boundary",
    "nonlinear_material_curve",
    "point_source_or_probe_identity",
}


def _sha(value: Mapping[str, Any]) -> str:
    payload = dict(value)
    payload.pop("evidence_payload_sha256", None)
    raw = json.dumps(
        payload, ensure_ascii=True, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def _mapping(value: object) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


def _nonnegative_int(value: object) -> int | None:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        return None
    return value


def _counter_map(value: object, allowed: set[str]) -> dict[str, int] | None:
    mapping = _mapping(value)
    if not mapping or not set(mapping) <= allowed:
        return None
    result: dict[str, int] = {}
    for key, raw in mapping.items():
        count = _nonnegative_int(raw)
        if count is None:
            return None
        result[str(key)] = count
    return result


def _timestamp_is_parseable(value: object) -> bool:
    try:
        parsed = datetime.fromisoformat(str(value).replace("Z", "+00:00"))
    except ValueError:
        return False
    return parsed.tzinfo is not None


def validate_legacy_signature_migration_evidence(
    evidence: Mapping[str, Any],
) -> dict[str, Any]:
    """Validate aggregate routing evidence without local paths or provenance."""

    if not isinstance(evidence, Mapping):
        raise ValueError("evidence must be an object")
    versions = _mapping(evidence.get("execution_version"))
    lanes = _counter_map(evidence.get("lane_counts"), _LANES)
    requirements = _counter_map(evidence.get("requirement_counts"), _REQUIREMENTS)
    signature_count = _nonnegative_int(evidence.get("signature_count"))
    assigned = _nonnegative_int(evidence.get("assigned_signature_count"))
    candidates = _nonnegative_int(evidence.get("native_emitter_candidate_count"))
    false_ready = _nonnegative_int(evidence.get("false_solver_ready_count"))
    mismatches = _nonnegative_int(evidence.get("identity_mismatch_count"))
    conversion_errors = _nonnegative_int(evidence.get("conversion_error_count"))
    compile_failures = _nonnegative_int(evidence.get("script_compile_failure_count"))
    checks = {
        "closed_world_top_level": set(evidence) == _TOP_LEVEL_KEYS,
        "schema": evidence.get("schema") == SCHEMA,
        "execution_timestamp_recorded": _timestamp_is_parseable(
            evidence.get("executed_at_utc")
        ),
        "execution_version_closed_world": set(versions) == _VERSION_KEYS
        and all(
            bool(_TOKEN_RE.fullmatch(str(versions.get(key, ""))))
            for key in _VERSION_KEYS
        )
        and versions.get("solver") == "not-launched",
        "solver_not_launched": evidence.get("solver_launched") is False,
        "live_dependency_not_required": evidence.get("live_dependency_required") is False,
        "lane_vocabulary_closed": lanes is not None,
        "requirement_vocabulary_closed": requirements is not None,
        "all_signatures_assigned": None not in (signature_count, assigned)
        and signature_count > 0
        and assigned == signature_count
        and lanes is not None
        and sum(lanes.values()) == signature_count,
        "no_false_ready_or_conversion_failure": None not in (
            candidates,
            false_ready,
            mismatches,
            conversion_errors,
            compile_failures,
        )
        and candidates <= signature_count
        and false_ready == 0
        and mismatches == 0
        and conversion_errors == 0
        and compile_failures == 0,
        "signature_catalog_digest": bool(
            _SHA256_RE.fullmatch(
                str(evidence.get("signature_catalog_sha256", "")).lower()
            )
        ),
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
        "checks": checks,
        "signature_count": signature_count,
        "lane_counts": lanes or {},
        "requirement_counts": requirements or {},
        "native_emitter_candidate_count": candidates,
        "live_dependency_required": False,
    }
