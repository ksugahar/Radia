"""Generic per-server status / introspection tool factory.

Pattern adopted 2026-05-24 from wjc9011/COMSOL_Multiphysics_MCP
(comsol_status, pdf_search_status, pdf_list_modules — each server
gives a clear "what can I do, am I healthy" snapshot).

Why this helps LLMs (and humans):
  - 49 radia_mcp servers (incl. meta + panel_review) — picking the
    right one is a discovery problem. mcp-server-elf is a separate
    PyPI package (`pip install mcp-server-elf`) and not in this
    catalog.
  - Each server's `<name>_status()` tool returns:
      - server name + module path
      - tool list (auto-introspected from FastMCP)
      - optional-dependency probe (e.g. chromadb installed?)
      - lightweight `--selftest` command and optional heavier audit command
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
from functools import lru_cache
import importlib
import importlib.metadata
import inspect
import sys
from typing import Optional


@lru_cache(maxsize=None)
def _probe_dep(module_name: str) -> dict:
    """Check if an optional dep is importable without importing it.

    Status tools are called frequently by agents. Importing optional packages
    such as chromadb, matplotlib, or ngsolve just to read ``__version__`` makes
    status calls slow and can trigger side effects. ``find_spec`` plus package
    metadata is enough for an MCP health hint.
    """
    spec = importlib.util.find_spec(module_name)
    if spec is None:
        return {"installed": False, "version": None}
    try:
        version = importlib.metadata.version(module_name)
    except importlib.metadata.PackageNotFoundError:
        version = "unknown"
    except Exception as e:
        version = "unknown"
        return {"installed": True, "version": version,
                "version_probe_error": str(e)}
    return {"installed": True, "version": version}


def build_status_payload(
    server_name: str,
    description: str,
    subpackage: str,
    related_servers: Optional[list[str]] = None,
    optional_deps: Optional[list[str]] = None,
    mcp_tools: Optional[list[str]] = None,
    audit_command: Optional[str] = None,
) -> dict:
    """Build the status dict shape that every radia_mcp.* server returns."""
    payload = {
        "server": server_name,
        "subpackage": subpackage,
        "description": description,
        "selftest_command": f"{server_name} --selftest",
        "python_version": (f"{sys.version_info.major}."
                            f"{sys.version_info.minor}."
                            f"{sys.version_info.micro}"),
    }
    if audit_command:
        payload["audit_command"] = audit_command
    if mcp_tools is not None:
        payload["tools"] = sorted(mcp_tools)
        payload["n_tools"] = len(mcp_tools)
    if related_servers:
        payload["related_servers"] = related_servers
    if optional_deps:
        payload["optional_deps"] = {
            d: dict(_probe_dep(d)) for d in optional_deps
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
    manager = getattr(mcp, "_tool_manager", None)
    if manager is not None:
        tools = getattr(manager, "_tools", None)
        if isinstance(tools, dict):
            return list(tools.keys())

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
        tools = mcp.list_tools()
        if inspect.isawaitable(tools):
            tools.close()
            return []
        return [t.name for t in tools]
    except Exception:
        return []


def register_status_tool(
    mcp,
    server_name: str,
    description: str,
    subpackage: str,
    related_servers: Optional[list[str]] = None,
    optional_deps: Optional[list[str]] = None,
    audit_command: Optional[str] = None,
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
        audit_command: optional heavier repo-wide validation/audit command,
                       separate from the lightweight `--selftest` health check.
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
            audit_command=audit_command,
        )
    # docstring set after definition so it appears in tool description
    _status.__doc__ = (
        f"Status / introspection for {server_name}.\n\n"
        f"Returns a dict with: server name, subpackage path, tool list,\n"
        f"optional dependency probe, Python version, selftest command,\n"
        f"optional audit command, and related servers.\n"
        f"Call this first if you're unsure whether the server is healthy\n"
        f"or what tools it exposes."
    )
