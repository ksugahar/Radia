"""Closed-world gate for solver-neutral historical corpus evidence."""

from __future__ import annotations

import hashlib
import json
import re
from collections.abc import Mapping
from datetime import datetime
from typing import Any


SCHEMA = "radia.legacy-solver-corpus-evidence.v1"
GATE_SCHEMA = "radia.legacy-solver-corpus-gate.v1"
_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
_TOKEN_RE = re.compile(r"^[A-Za-z0-9._-]+$")
_TOP_LEVEL_KEYS = {
    "schema",
    "executed_at_utc",
    "execution_version",
    "solver_launched",
    "live_dependency_required",
    "model_surface",
    "automation_surface",
    "document_surface",
    "topic_surface",
    "evidence_payload_sha256",
}
_VERSION_KEYS = {"producer", "producer_version", "solver"}
_SURFACE_KEYS = {
    "model_surface": {
        "input_count",
        "unique_content_count",
        "semantic_signature_count",
        "parse_error_count",
        "catalog_sha256",
    },
    "automation_surface": {
        "script_count",
        "unique_command_count",
        "classified_command_count",
        "unknown_command_count",
        "catalog_sha256",
    },
    "document_surface": {
        "document_count",
        "unique_content_count",
        "covered_document_count",
        "failure_count",
        "catalog_sha256",
    },
    "topic_surface": {
        "discovered_topic_count",
        "catalogued_topic_count",
        "missing_topic_count",
        "catalog_sha256",
    },
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


def _surface_shape(evidence: Mapping[str, Any], name: str) -> bool:
    surface = _mapping(evidence.get(name))
    return set(surface) == _SURFACE_KEYS[name]


def _surface_hashes(evidence: Mapping[str, Any]) -> bool:
    return all(
        bool(_SHA256_RE.fullmatch(str(_mapping(evidence.get(name)).get("catalog_sha256", "")).lower()))
        for name in _SURFACE_KEYS
    )


def _timestamp_is_parseable(value: object) -> bool:
    try:
        parsed = datetime.fromisoformat(str(value).replace("Z", "+00:00"))
    except ValueError:
        return False
    return parsed.tzinfo is not None


def validate_legacy_corpus_evidence(evidence: Mapping[str, Any]) -> dict[str, Any]:
    """Validate aggregate corpus coverage without reading a local path."""

    if not isinstance(evidence, Mapping):
        raise ValueError("evidence must be an object")

    versions = _mapping(evidence.get("execution_version"))
    models = _mapping(evidence.get("model_surface"))
    automation = _mapping(evidence.get("automation_surface"))
    documents = _mapping(evidence.get("document_surface"))
    topics = _mapping(evidence.get("topic_surface"))

    model_inputs = _nonnegative_int(models.get("input_count"))
    model_unique = _nonnegative_int(models.get("unique_content_count"))
    model_signatures = _nonnegative_int(models.get("semantic_signature_count"))
    parse_errors = _nonnegative_int(models.get("parse_error_count"))
    scripts = _nonnegative_int(automation.get("script_count"))
    commands = _nonnegative_int(automation.get("unique_command_count"))
    classified = _nonnegative_int(automation.get("classified_command_count"))
    unknown = _nonnegative_int(automation.get("unknown_command_count"))
    document_count = _nonnegative_int(documents.get("document_count"))
    document_unique = _nonnegative_int(documents.get("unique_content_count"))
    covered_documents = _nonnegative_int(documents.get("covered_document_count"))
    document_failures = _nonnegative_int(documents.get("failure_count"))
    discovered_topics = _nonnegative_int(topics.get("discovered_topic_count"))
    catalogued_topics = _nonnegative_int(topics.get("catalogued_topic_count"))
    missing_topics = _nonnegative_int(topics.get("missing_topic_count"))

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
        "surface_shapes_closed_world": all(
            _surface_shape(evidence, name) for name in _SURFACE_KEYS
        ),
        "model_surface_complete": None not in (
            model_inputs,
            model_unique,
            model_signatures,
            parse_errors,
        )
        and model_inputs > 0
        and 0 < model_signatures <= model_unique <= model_inputs
        and parse_errors == 0,
        "automation_surface_complete": None not in (
            scripts,
            commands,
            classified,
            unknown,
        )
        and scripts > 0
        and commands > 0
        and classified == commands
        and unknown == 0,
        "document_surface_complete": None not in (
            document_count,
            document_unique,
            covered_documents,
            document_failures,
        )
        and document_count > 0
        and 0 < document_unique <= document_count
        and covered_documents == document_count
        and document_failures == 0,
        "topic_surface_complete": None not in (
            discovered_topics,
            catalogued_topics,
            missing_topics,
        )
        and discovered_topics > 0
        and catalogued_topics == discovered_topics
        and missing_topics == 0,
        "surface_hashes_recorded": _surface_hashes(evidence),
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
        "coverage": {
            "model_inputs": model_inputs,
            "semantic_signatures": model_signatures,
            "automation_scripts": scripts,
            "automation_commands": commands,
            "documents": document_count,
            "topics": discovered_topics,
        },
        "evidence_payload_sha256": evidence.get("evidence_payload_sha256"),
        "live_dependency_required": False,
    }
