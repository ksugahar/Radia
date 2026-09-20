"""Explicit MCP metadata shared by radia-mcp capability servers."""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterable, Mapping
from datetime import datetime
import json
import math
import re
import sys
from typing import Any

from radia_mcp.common.server_hardening import annotation_preset_name, infer_tool_annotations


SCHEMA = "cae-ai-lab.mcp-server-contract.v3"
ARTIFACT_IDENTITY_SCHEMA = "cae-ai-lab.solver-artifact-identity.v1"
CAPABILITY_MATURITY_LEVELS = (
    "knowledge_only",
    "artifact_gate",
    "solver_ready",
    "live_verified",
    "cross_validated",
)
_CLIENT_HISTORY_LIMIT = 8
_CLIENT_FIELD_LIMIT = 256
_SHA256 = re.compile(r"^[0-9a-fA-F]{64}$")
_ARTIFACT_GATE = re.compile(
    r"(?:^|_)(?:audit|check|compare|diagnos|gate|inspect|lint|preflight|"
    r"review|validate)(?:_|$)"
)
_SOLVER_READY = re.compile(
    r"(?:^|_)(?:execute|launch|run|solve|start|submit)(?:_|$)"
)


def _bounded_client_field(value: Any) -> str | None:
    if value is None:
        return None
    return str(value)[:_CLIENT_FIELD_LIMIT]


def _capture_current_client(low_level: Any, server_name: str, version: str) -> None:
    """Capture initialize clientInfo on the first request after handshake."""
    try:
        session = low_level.request_context.session
    except LookupError:
        return
    if getattr(session, "_caeai_client_info_captured", False):
        return
    params = getattr(session, "client_params", None)
    client_info = getattr(params, "clientInfo", None)
    if client_info is None:
        return

    snapshot = {
        "name": _bounded_client_field(getattr(client_info, "name", None)),
        "title": _bounded_client_field(getattr(client_info, "title", None)),
        "version": _bounded_client_field(getattr(client_info, "version", None)),
    }
    session._caeai_client_info_captured = True
    state = low_level._caeai_client_connection_state
    state["connection_count"] += 1
    state["latest"] = snapshot
    state["recent"] = (state["recent"] + [snapshot])[-_CLIENT_HISTORY_LIMIT:]
    print(
        f"[{server_name}] {version} client connected: "
        f"{json.dumps(snapshot, ensure_ascii=False)}; "
        "solver/application sessions were not modified",
        file=sys.stderr,
        flush=True,
    )


def _install_client_info_capture(mcp: Any, server_name: str, version: str) -> bool:
    """Install bounded, process-local clientInfo capture without side effects."""
    low_level = getattr(mcp, "_mcp_server", None)
    if low_level is None or getattr(low_level, "_caeai_client_info_capture", False):
        return False
    low_level._caeai_client_connection_state = {
        "connection_count": 0,
        "latest": None,
        "recent": [],
    }
    for request_type, previous in tuple(low_level.request_handlers.items()):
        async def _capture_then_handle(
            request: Any, _previous: Any = previous
        ) -> Any:
            _capture_current_client(low_level, server_name, version)
            return await _previous(request)

        low_level.request_handlers[request_type] = _capture_then_handle
    low_level._caeai_client_info_capture = True
    return True


def get_client_connection_state(mcp: Any) -> dict[str, Any]:
    """Return a detached copy of bounded client connection metadata."""
    low_level = getattr(mcp, "_mcp_server", None)
    state = getattr(low_level, "_caeai_client_connection_state", None) or {}
    latest = state.get("latest")
    return {
        "connection_count": int(state.get("connection_count", 0)),
        "latest": dict(latest) if isinstance(latest, dict) else None,
        "recent": [
            dict(item)
            for item in state.get("recent", [])
            if isinstance(item, dict)
        ],
    }


def _infer_maturity(tool_name: str) -> str:
    name = tool_name.lower()
    if _SOLVER_READY.search(name):
        return "solver_ready"
    if _ARTIFACT_GATE.search(name):
        return "artifact_gate"
    return "knowledge_only"


def build_capability_maturity_contract() -> dict[str, Any]:
    """Return the fleet-wide, evidence-ordered capability vocabulary."""
    return {
        "schema": "cae-ai-lab.mcp-capability-maturity.v1",
        "levels": list(CAPABILITY_MATURITY_LEVELS),
        "ordering": "evidence_strength",
        "automatic_ceiling": "solver_ready",
        "evidence_required_for": ["live_verified", "cross_validated"],
    }


def build_artifact_identity_contract() -> dict[str, Any]:
    """Return the neutral solver-artifact identity contract."""
    return {
        "schema": ARTIFACT_IDENTITY_SCHEMA,
        "canonical_container": ".hdf5",
        "required_fields": [
            "schema",
            "artifact_id",
            "created_at_utc",
            "producer",
            "solver",
            "status",
            "coordinate_system",
            "unit_system",
            "input_sha256",
            "result_sha256",
            "total_compute_seconds",
            "timing_breakdown_s",
        ],
        "conditional_fields": {
            "mesh_based": ["mesh_sha256"],
            "transient": ["time_axis"],
        },
        "digest": "sha256",
        "timing_policy": "one to four largest nonnegative stages",
    }


def _utc_timestamp(value: Any) -> bool:
    if not isinstance(value, str) or not value:
        return False
    try:
        parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError:
        return False
    return parsed.tzinfo is not None and parsed.utcoffset().total_seconds() == 0


def validate_solver_artifact_identity(
    artifact: Mapping[str, Any],
) -> dict[str, Any]:
    """Validate portable solver-result identity before physics-specific gates."""
    contract = build_artifact_identity_contract()
    errors: list[str] = []
    missing = [key for key in contract["required_fields"] if key not in artifact]
    errors.extend(f"missing:{key}" for key in missing)

    if artifact.get("schema") != ARTIFACT_IDENTITY_SCHEMA:
        errors.append("schema:unexpected")
    for key in ("artifact_id", "coordinate_system"):
        if not isinstance(artifact.get(key), str) or not artifact.get(key):
            errors.append(f"{key}:nonempty_string_required")
    if artifact.get("unit_system") != "SI":
        errors.append("unit_system:SI_required")
    if artifact.get("status") not in {"complete", "ok", "passed"}:
        errors.append("status:completed_result_required")
    if not _utc_timestamp(artifact.get("created_at_utc")):
        errors.append("created_at_utc:UTC_ISO8601_required")

    for key in ("producer", "solver"):
        value = artifact.get(key)
        if not isinstance(value, Mapping):
            errors.append(f"{key}:mapping_required")
            continue
        for field in ("name", "version"):
            if not isinstance(value.get(field), str) or not value.get(field):
                errors.append(f"{key}.{field}:nonempty_string_required")

    for key in ("input_sha256", "result_sha256"):
        if not _SHA256.fullmatch(str(artifact.get(key, ""))):
            errors.append(f"{key}:sha256_required")

    elapsed = artifact.get("total_compute_seconds")
    if (
        isinstance(elapsed, bool)
        or not isinstance(elapsed, (int, float))
        or not math.isfinite(float(elapsed))
        or float(elapsed) < 0.0
    ):
        errors.append("total_compute_seconds:finite_nonnegative_required")

    timing = artifact.get("timing_breakdown_s")
    if not isinstance(timing, Mapping) or not 1 <= len(timing) <= 4:
        errors.append("timing_breakdown_s:one_to_four_stages_required")
    elif any(
        not isinstance(name, str)
        or not name
        or isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(float(value))
        or float(value) < 0.0
        for name, value in timing.items()
    ):
        errors.append("timing_breakdown_s:finite_nonnegative_values_required")

    mesh_based = bool(artifact.get("mesh_based", False))
    if mesh_based and not _SHA256.fullmatch(str(artifact.get("mesh_sha256", ""))):
        errors.append("mesh_sha256:required_for_mesh_based_result")
    transient = artifact.get("study_type") == "transient"
    if transient:
        time_axis = artifact.get("time_axis")
        if not isinstance(time_axis, Mapping):
            errors.append("time_axis:required_for_transient_result")
        elif any(key not in time_axis for key in ("unit", "count", "start", "end")):
            errors.append("time_axis:unit_count_start_end_required")

    unique_errors = sorted(set(errors))
    return {
        "schema": "cae-ai-lab.solver-artifact-identity-validation.v1",
        "status": "accepted" if not unique_errors else "rejected",
        "accepted": not unique_errors,
        "errors": unique_errors,
        "contract": contract,
    }


def audit_tool_contract(mcp: Any) -> dict[str, Any]:
    """Summarize the live FastMCP registry without invoking any tool."""
    tools = getattr(getattr(mcp, "_tool_manager", None), "_tools", {})
    return audit_tool_definitions(tools)


def audit_tool_definitions(tools: dict[str, Any]) -> dict[str, Any]:
    """Audit FastMCP definitions or public SDK Tool objects identically."""
    missing_titles: list[str] = []
    missing_annotations: list[str] = []
    missing_contract_meta: list[str] = []
    invalid_maturity: list[str] = []
    presets: Counter[str] = Counter()
    sources: Counter[str] = Counter()
    maturity: Counter[str] = Counter()
    structured: list[str] = []
    unstructured: list[str] = []
    for name, tool in tools.items():
        if not getattr(tool, "title", None):
            missing_titles.append(name)
        annotations = getattr(tool, "annotations", None)
        if annotations is None or any(
            value is None
            for value in (
                getattr(annotations, "readOnlyHint", None),
                getattr(annotations, "destructiveHint", None),
                getattr(annotations, "idempotentHint", None),
                getattr(annotations, "openWorldHint", None),
            )
        ):
            missing_annotations.append(name)
        presets[annotation_preset_name(annotations)] += 1
        meta = dict(getattr(tool, "meta", None) or {})
        if meta.get("caeai.contract") != SCHEMA:
            missing_contract_meta.append(name)
        maturity_level = str(meta.get("caeai.maturity", ""))
        maturity[maturity_level or "unspecified"] += 1
        if maturity_level not in CAPABILITY_MATURITY_LEVELS:
            invalid_maturity.append(name)
        sources[str(meta.get("caeai.annotation_source", "unspecified"))] += 1
        metadata = getattr(tool, "fn_metadata", None)
        if (getattr(metadata, "output_schema", None) is not None
                or getattr(tool, "outputSchema", None) is not None):
            structured.append(name)
        else:
            unstructured.append(name)
    return {
        "schema": SCHEMA,
        "complete": not (
            missing_titles
            or missing_annotations
            or missing_contract_meta
            or invalid_maturity
        ),
        "n_tools": len(tools),
        "annotation_presets": dict(sorted(presets.items())),
        "annotation_sources": dict(sorted(sources.items())),
        "capability_maturity": dict(sorted(maturity.items())),
        "capability_maturity_contract": build_capability_maturity_contract(),
        "artifact_identity_contract": build_artifact_identity_contract(),
        "missing_titles": sorted(missing_titles),
        "missing_annotations": sorted(missing_annotations),
        "missing_contract_meta": sorted(missing_contract_meta),
        "invalid_maturity": sorted(invalid_maturity),
        "structured_output_tools": sorted(structured),
        "unstructured_output_tools": sorted(unstructured),
    }


def apply_tool_contract(
    mcp: Any,
    *,
    server_name: str,
    version: str,
    tool_prefix: str = "",
    set_server_metadata: bool = True,
) -> dict[str, Any]:
    """Attach titles, annotations, instructions, and version to one server."""

    tools = getattr(getattr(mcp, "_tool_manager", None), "_tools", {})
    selected = {
        name: tool
        for name, tool in tools.items()
        if not tool_prefix or name.startswith(tool_prefix)
    }
    for name, tool in selected.items():
        if not getattr(tool, "title", None):
            tool.title = " ".join(part.capitalize() for part in name.split("_") if part)
        meta = dict(getattr(tool, "meta", None) or {})
        if getattr(tool, "annotations", None) is None:
            annotations, source = infer_tool_annotations(
                name, str(getattr(tool, "description", "") or "")
            )
            tool.annotations = annotations
            meta["caeai.annotation_source"] = source
        else:
            meta.setdefault("caeai.annotation_source", "explicit-decorator")
        meta["caeai.annotation_preset"] = annotation_preset_name(
            tool.annotations
        )
        meta["caeai.contract"] = SCHEMA
        meta.setdefault("caeai.maturity", _infer_maturity(name))
        metadata = getattr(tool, "fn_metadata", None)
        meta["caeai.output_mode"] = (
            "structured"
            if getattr(metadata, "output_schema", None) is not None
            else "unstructured"
        )
        tool.meta = meta

    low_level = getattr(mcp, "_mcp_server", None)
    if set_server_metadata and low_level is not None:
        low_level.version = version
        if not getattr(low_level, "instructions", None):
            low_level.instructions = (
                f"Call the status/profile tool before routing {server_name}; "
                "validate solver ownership and artifacts before side effects."
            )
    _install_client_info_capture(mcp, server_name, version)
    audit = audit_tool_contract(mcp)
    mcp._radia_runtime_contract = audit
    return audit


def build_runtime_contract(
    capability_packs: Iterable[str], session_policy: str
) -> dict[str, Any]:
    return {
        "schema": SCHEMA,
        "runtime_core": "FastMCP over stdio",
        "capability_packs": list(capability_packs),
        "skill_layer": "separate workflow skills",
        "complete_tool_annotations": True,
        "annotation_policy": (
            "explicit preferred; conservative audited inference otherwise"
        ),
        "capability_maturity": build_capability_maturity_contract(),
        "artifact_identity": build_artifact_identity_contract(),
        "client_connection_visibility": {
            "capture_point": "first post-initialize request",
            "captured_fields": ["name", "title", "version"],
            "process_local_history_limit": _CLIENT_HISTORY_LIMIT,
            "source_session_side_effect": False,
        },
        "session_policy": session_policy,
        "protocol_smoke": "initialize + tools/list + representative tools/call",
    }
