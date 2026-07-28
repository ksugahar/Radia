"""Explicit MCP metadata shared by radia-mcp capability servers."""

from __future__ import annotations

import re
from collections.abc import Iterable
from typing import Any

from mcp.types import ToolAnnotations


SCHEMA = "cae-ai-lab.mcp-server-contract.v1"
_WRITE = re.compile(
    r"(?:^|_)(?:add|convert|create|delete|edit|exec|export|generate|import|install|launch|mesh|remove|run|save|set|solve|start|stop|submit|update|write)(?:_|$)"
)
_READ = re.compile(
    r"(?:^|_)(?:analysis|analyze|audit|catalog|check|compare|diagnos|estimate|explain|gate|get|health|index|inspect|inventory|lesson|lint|list|overview|plan|policy|preflight|profile|query|read|readiness|reference|report|roadmap|route|search|status|summary|taxonomy|topic|topology|validate)(?:_|$)"
)
_OPEN = re.compile(
    r"(?:^|_)(?:download|fetch|github|http|live|online|remote|session|web)(?:_|$)"
)


def _annotations(name: str, description: str) -> ToolAnnotations:
    mutating = bool(
        _WRITE.search(name.lower())
        or re.search(
            r"\b(?:create|delete|edit|execute|export|launch|modify|run|save|solve|update|write)\b",
            description,
            re.IGNORECASE,
        )
    )
    read_only = not mutating and bool(
        _READ.search(name.lower())
        or re.search(
            r"\b(?:analyze|check|compare|compile|compute|inspect|list|plan|read|report|return|route|search|validate)\b",
            description,
            re.IGNORECASE,
        )
    )
    if read_only:
        return ToolAnnotations(
            readOnlyHint=True,
            destructiveHint=False,
            idempotentHint=True,
            openWorldHint=bool(_OPEN.search(name.lower())),
        )
    return ToolAnnotations(
        readOnlyHint=False,
        destructiveHint=True,
        idempotentHint=False,
        openWorldHint=True,
    )


def apply_tool_contract(
    mcp: Any,
    *,
    server_name: str,
    version: str,
    tool_prefix: str = "",
    set_server_metadata: bool = True,
) -> dict[str, int]:
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
        if getattr(tool, "annotations", None) is None:
            tool.annotations = _annotations(
                name, str(getattr(tool, "description", "") or "")
            )
        meta = dict(getattr(tool, "meta", None) or {})
        meta.setdefault("caeai.contract", SCHEMA)
        tool.meta = meta

    low_level = getattr(mcp, "_mcp_server", None)
    if set_server_metadata and low_level is not None:
        low_level.version = version
        if not getattr(low_level, "instructions", None):
            low_level.instructions = (
                f"Call the status/profile tool before routing {server_name}; "
                "validate solver ownership and artifacts before side effects."
            )
    return {"tools": len(selected)}


def build_runtime_contract(
    capability_packs: Iterable[str], session_policy: str
) -> dict[str, Any]:
    return {
        "schema": SCHEMA,
        "runtime_core": "FastMCP over stdio",
        "capability_packs": list(capability_packs),
        "skill_layer": "separate workflow skills",
        "explicit_tool_annotations": True,
        "session_policy": session_policy,
        "protocol_smoke": "initialize + tools/list + representative tools/call",
    }
