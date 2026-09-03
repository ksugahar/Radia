"""Explicit MCP metadata shared by radia-mcp capability servers."""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterable
from typing import Any

from .server_hardening import (
    annotation_preset_name,
    infer_tool_annotations,
)


SCHEMA = "cae-ai-lab.mcp-server-contract.v2"


def audit_tool_contract(mcp: Any) -> dict[str, Any]:
    """Summarize the live FastMCP registry without invoking any tool."""
    tools = getattr(getattr(mcp, "_tool_manager", None), "_tools", {})
    missing_titles: list[str] = []
    missing_annotations: list[str] = []
    missing_contract_meta: list[str] = []
    presets: Counter[str] = Counter()
    sources: Counter[str] = Counter()
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
        sources[str(meta.get("caeai.annotation_source", "unspecified"))] += 1
        metadata = getattr(tool, "fn_metadata", None)
        if getattr(metadata, "output_schema", None) is not None:
            structured.append(name)
        else:
            unstructured.append(name)
    return {
        "schema": SCHEMA,
        "complete": not (
            missing_titles or missing_annotations or missing_contract_meta
        ),
        "n_tools": len(tools),
        "annotation_presets": dict(sorted(presets.items())),
        "annotation_sources": dict(sorted(sources.items())),
        "missing_titles": sorted(missing_titles),
        "missing_annotations": sorted(missing_annotations),
        "missing_contract_meta": sorted(missing_contract_meta),
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
        "session_policy": session_policy,
        "protocol_smoke": "initialize + tools/list + representative tools/call",
    }
