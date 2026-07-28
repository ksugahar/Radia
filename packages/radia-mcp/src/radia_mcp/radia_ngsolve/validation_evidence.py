"""Content-addressed evidence bundles for retiring a solver dependency."""

from __future__ import annotations

import hashlib
import json
import re
from collections.abc import Mapping, Sequence
from typing import Any


SCHEMA = "radia.validation-evidence-bundle.v1"
_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
_COMMIT_RE = re.compile(r"^[0-9a-f]{40}$")


def _canonical_sha256(value: object) -> str:
    payload = json.dumps(
        value, ensure_ascii=True, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _string_list(value: object, label: str) -> list[str]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise ValueError(f"{label} must be a list")
    rows = [str(item).strip() for item in value]
    if not rows or any(not item for item in rows):
        raise ValueError(f"{label} must contain nonempty strings")
    if len(rows) != len(set(rows)):
        raise ValueError(f"{label} must not contain duplicates")
    return rows


def _lane_is_verified(lane: object) -> bool:
    if not isinstance(lane, Mapping):
        return False
    commit = str(lane.get("commit", "")).lower()
    return all(
        str(lane.get(field, "")).strip()
        for field in (
            "owner",
            "capability_gain",
            "positive_probe",
            "negative_probe",
            "protocol_probe",
            "verification",
        )
    ) and bool(_COMMIT_RE.fullmatch(commit))


def _validate_artifact(record: Mapping[str, Any]) -> dict[str, Any]:
    artifact_id = str(record.get("artifact_id", "")).strip()
    raw = record.get("artifact_json")
    declared_sha = str(record.get("artifact_sha256", "")).lower()
    if not artifact_id:
        raise ValueError("artifact_id is required")
    if not isinstance(raw, str) or not raw.strip():
        raise ValueError(f"{artifact_id}: artifact_json is required")
    actual_sha = hashlib.sha256(raw.encode("utf-8")).hexdigest()
    if not _SHA256_RE.fullmatch(declared_sha) or declared_sha != actual_sha:
        raise ValueError(f"{artifact_id}: artifact_sha256 does not match content")
    try:
        artifact = json.loads(raw)
    except json.JSONDecodeError as exc:
        raise ValueError(f"{artifact_id}: artifact_json is invalid: {exc}") from exc
    if not isinstance(artifact, Mapping):
        raise ValueError(f"{artifact_id}: artifact_json must encode an object")
    if not str(artifact.get("schema", "")).strip():
        raise ValueError(f"{artifact_id}: artifact schema is required")
    if artifact.get("pass") is not True:
        raise ValueError(f"{artifact_id}: artifact must have pass=true")
    if not str(
        artifact.get("created_at_utc") or artifact.get("executed_at_utc") or ""
    ).strip():
        raise ValueError(f"{artifact_id}: execution timestamp is required")
    versions = artifact.get("versions")
    execution_version = artifact.get("execution_version")
    if not isinstance(versions, Mapping) and not isinstance(execution_version, Mapping):
        raise ValueError(f"{artifact_id}: version identity is required")
    balance = artifact.get("mcp_balance")
    if not isinstance(balance, Mapping):
        raise ValueError(f"{artifact_id}: mcp_balance is required")
    if balance.get("policy") != "equal_capability_gain_v1":
        raise ValueError(f"{artifact_id}: equal public/source learning is required")
    if balance.get("status") != "verified":
        raise ValueError(f"{artifact_id}: mcp_balance must be verified")
    if not _lane_is_verified(balance.get("public")):
        raise ValueError(f"{artifact_id}: public lane evidence is incomplete")
    if not _lane_is_verified(balance.get("source_tool")):
        raise ValueError(f"{artifact_id}: source lane evidence is incomplete")
    capabilities = _string_list(record.get("capabilities"), "capabilities")
    declared = _string_list(
        artifact.get("retirement_capabilities"), "retirement_capabilities"
    )
    if set(declared) != set(capabilities):
        raise ValueError(
            f"{artifact_id}: record capabilities do not match artifact attestation"
        )
    return {
        "artifact_id": artifact_id,
        "artifact_sha256": actual_sha,
        "artifact_schema": str(artifact["schema"]),
        "capabilities": sorted(capabilities),
        "capability_id": str(balance.get("capability_id", "")),
        "public_commit": str(balance["public"]["commit"]).lower(),
        "source_commit": str(balance["source_tool"]["commit"]).lower(),
    }


def validate_evidence_bundle(packet: Mapping[str, Any]) -> dict[str, Any]:
    """Validate in-memory artifacts without granting local filesystem access."""

    if not isinstance(packet, Mapping):
        raise ValueError("packet must be an object")
    inventory_sha = str(packet.get("inventory_sha256", "")).lower()
    if not _SHA256_RE.fullmatch(inventory_sha):
        raise ValueError("inventory_sha256 must be a lowercase SHA-256 digest")
    required = sorted(_string_list(packet.get("required_capabilities"), "required_capabilities"))
    records = packet.get("artifacts")
    if not isinstance(records, Sequence) or isinstance(records, (str, bytes)):
        raise ValueError("artifacts must be a list")

    validated: list[dict[str, Any]] = []
    owners: dict[str, str] = {}
    for raw_record in records:
        if not isinstance(raw_record, Mapping):
            raise ValueError("each artifact record must be an object")
        record = _validate_artifact(raw_record)
        for capability in record["capabilities"]:
            previous = owners.get(capability)
            if previous is not None:
                raise ValueError(
                    f"capability {capability!r} is claimed by both {previous!r} "
                    f"and {record['artifact_id']!r}"
                )
            owners[capability] = record["artifact_id"]
        validated.append(record)

    accepted = sorted(set(required) & set(owners))
    missing = sorted(set(required) - set(owners))
    unexpected = sorted(set(owners) - set(required))
    canonical = {
        "schema": SCHEMA,
        "inventory_sha256": inventory_sha,
        "required_capabilities": required,
        "artifacts": sorted(validated, key=lambda row: row["artifact_id"]),
    }
    return {
        **canonical,
        "status": "complete" if not missing else "incomplete",
        "contract_valid": True,
        "retirement_ready": not missing,
        "solver_uninstall_performed": False,
        "accepted_capabilities": accepted,
        "missing_capabilities": missing,
        "unexpected_capabilities": unexpected,
        "capability_owners": dict(sorted(owners.items())),
        "evidence_bundle_sha256": _canonical_sha256(canonical),
    }
