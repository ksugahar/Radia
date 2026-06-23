"""Generic per-server status / introspection tool factory.

Pattern adopted 2026-05-24 from wjc9011/COMSOL_Multiphysics_MCP
(comsol_status, pdf_search_status, pdf_list_modules — each server
gives a clear "what can I do, am I healthy" snapshot).

Why this helps LLMs (and humans):
  - 37 radia_mcp servers (incl. meta + panel_review) — picking the
    right one is a discovery problem. mcp-server-elf is a separate
    PyPI package (`pip install mcp-server-elf`) and not in this
    catalog.
  - Each server's `<name>_status()` tool returns:
      - server name + module path
      - tool list (auto-introspected from FastMCP)
      - optional-dependency probe (e.g. chromadb installed?)
      - one-line "what am I for" description
      - cross-link to companion servers
  - LLM can chain `<server>_status()` -> understand server -> call
    the right knowledge tool without trial-and-error

Usage in a server's server.py:

    from mcp.server.fastmcp import FastMCP
    from radia_mcp.common.status import register_status_tool

    mcp = FastMCP("mcp-server-bayesian-opt")

    # ... register your @mcp.tool() decorators here ...

    register_status_tool(
        mcp,
        server_name="mcp-server-bayesian-opt",
        description="Bayesian optimization for EM engineering",
        subpackage="radia_mcp.bayesian_opt",
        related_servers=["topology-optimization", "evolutionary"],
        optional_deps=["pymc", "emcee", "numpyro"],
    )
"""

from __future__ import annotations
import importlib
import sys
from typing import Optional


def _probe_dep(module_name: str) -> dict:
    """Check if an optional dep is importable. No version check (lazy)."""
    spec = importlib.util.find_spec(module_name)
    if spec is None:
        return {"installed": False, "version": None}
    try:
        m = importlib.import_module(module_name)
        return {"installed": True,
                "version": getattr(m, "__version__", "unknown")}
    except Exception as e:
        return {"installed": False, "version": None, "error": str(e)}


def build_status_payload(
    server_name: str,
    description: str,
    subpackage: str,
    related_servers: Optional[list[str]] = None,
    optional_deps: Optional[list[str]] = None,
    mcp_tools: Optional[list[str]] = None,
) -> dict:
    """Build the status dict shape that every radia_mcp.* server returns."""
    payload = {
        "server": server_name,
        "subpackage": subpackage,
        "description": description,
        "python_version": (f"{sys.version_info.major}."
                            f"{sys.version_info.minor}."
                            f"{sys.version_info.micro}"),
    }
    if mcp_tools is not None:
        payload["tools"] = sorted(mcp_tools)
        payload["n_tools"] = len(mcp_tools)
    if related_servers:
        payload["related_servers"] = related_servers
    if optional_deps:
        payload["optional_deps"] = {
            d: _probe_dep(d) for d in optional_deps
        }
        payload["all_optional_deps_installed"] = all(
            v["installed"] for v in payload["optional_deps"].values()
        )
    return payload


def _introspect_fastmcp_tools(mcp) -> list[str]:
    """Best-effort: list tool names registered on a FastMCP instance.

    FastMCP keeps tools in an internal registry; the exact attribute
    name varies by version. Falls back to None if introspection fails.
    """
    # Try the common attribute names across mcp 1.0+ versions
    for attr in ("_tools", "tools", "_tool_registry"):
        if hasattr(mcp, attr):
            obj = getattr(mcp, attr)
            if isinstance(obj, dict):
                return list(obj.keys())
            if hasattr(obj, "_tools"):   # nested
                inner = obj._tools
                if isinstance(inner, dict):
                    return list(inner.keys())
    # Try Pydantic-style or list iteration
    try:
        return [t.name for t in mcp.list_tools()]
    except Exception:
        return []


def register_status_tool(
    mcp,
    server_name: str,
    description: str,
    subpackage: str,
    related_servers: Optional[list[str]] = None,
    optional_deps: Optional[list[str]] = None,
    tool_name: Optional[str] = None,
) -> None:
    """Register a `<server>_status` MCP tool that returns the dict above.

    Args:
        mcp: FastMCP instance (already created in the server module)
        server_name: e.g. "mcp-server-bayesian-opt"
        description: one-line "what am I for"
        subpackage: e.g. "radia_mcp.bayesian_opt"
        related_servers: list of MCP server short names that pair well
                          (e.g. ["topology-optimization", "evolutionary"])
        optional_deps: pip package names to probe (chromadb, pymc, etc.)
        tool_name: override the auto-generated tool name (default:
                    derives from server_name by stripping "mcp-server-"
                    and appending "_status")
    """
    if tool_name is None:
        short = server_name.removeprefix("mcp-server-").replace("-", "_")
        tool_name = f"{short}_status"

    @mcp.tool(name=tool_name)
    def _status() -> dict:
        return build_status_payload(
            server_name=server_name,
            description=description,
            subpackage=subpackage,
            related_servers=related_servers,
            optional_deps=optional_deps,
            mcp_tools=_introspect_fastmcp_tools(mcp),
        )
    # docstring set after definition so it appears in tool description
    _status.__doc__ = (
        f"Status / introspection for {server_name}.\n\n"
        f"Returns a dict with: server name, subpackage path, tool list,\n"
        f"optional dependency probe, Python version, and related servers.\n"
        f"Call this first if you're unsure whether the server is healthy\n"
        f"or what tools it exposes."
    )
