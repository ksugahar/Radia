"""Factory for the uniform `<server>_topics()` MCP tool.

Pattern background
------------------
Many radia_mcp.* servers expose a single dispatcher tool that takes a
`topic: str` argument (e.g. `fusion(topic="tokamak")`,
`bayesian_opt(topic="acquisition_functions")`). For those servers the
list of valid topic values is buried inside the dispatcher's `if
topic == "x"` chain or in a docstring — neither is machine-readable
from the MCP `tools/list` response, so an LLM has to either:

  - call the dispatcher with an invalid topic and parse the error
    message (fragile, wastes a turn), or
  - call radia_mcp.meta.radia_mcp_get(name) to read the natural-
    language description (no enum guarantee).

`register_topics_tool` adds an explicit `<server>_topics()` tool that
returns the authoritative topic enumeration in a stable shape:

    {
        "server": "mcp-server-fusion-reactor",
        "n_topics": 4,
        "topics": [
            {"name": "overview",    "description": "..."},
            {"name": "tokamak",     "description": "..."},
            {"name": "stellarator", "description": "..."},
            {"name": "all",         "description": "Concatenated dump"},
        ],
    }

Not needed for multi-tool servers (cubit, radia-ngsolve, ...): each
of their tools is already its own "topic" and is exposed via the MCP
standard `tools/list` method.

Usage
-----
At the bottom of a dispatcher-style server.py, after the dispatcher
`@mcp.tool()` and before `def main()`:

    from ..common import register_topics_tool
    from .knowledge import TOPICS    # or define inline

    register_topics_tool(
        mcp,
        server_name="mcp-server-fusion-reactor",
        topics=TOPICS,                    # dict[str, str]
    )

The companion `TOPICS: dict[str, str]` lives in the knowledge module
as the single source of truth — the dispatcher chain reads from it
or developers keep both in sync manually.
"""
from __future__ import annotations

import importlib.metadata
from typing import Mapping


def register_topics_tool(
    mcp,
    server_name: str,
    topics: Mapping[str, str],
    tool_name: str | None = None,
) -> str:
    """Register a `<server>_topics` MCP tool exposing the enum.

    Args:
        mcp:          FastMCP instance.
        server_name:  Full server name (e.g. "mcp-server-fusion-reactor") — only
                      used in the payload; the registered tool name
                      derives from it via the rule
                          mcp-server-<short>  ->  <short>_topics
                      unless `tool_name` overrides.
        topics:       Mapping topic_name -> one-line description.
        tool_name:    Override the derived tool name (rarely needed).

    Returns:
        The actual tool name registered.
    """
    if tool_name is None:
        short = server_name.removeprefix("mcp-server-").replace("-", "_")
        tool_name = f"{short}_topics"

    # Stable, snapshot copy — server should declare topics in module
    # scope and never mutate.
    topics_snapshot = [
        {"name": k, "description": str(v)} for k, v in topics.items()
    ]

    def _topics_impl() -> dict[str, object]:
        return {
            "server": server_name,
            "n_topics": len(topics_snapshot),
            "topics": topics_snapshot,
            "hint": (
                f"Call the dispatcher tool with topic=<name> from this "
                f"list. Topic 'all' (if present) returns everything "
                f"concatenated."
            ),
        }

    # Attach a meaningful name so MCP tool listing shows it correctly.
    _topics_impl.__name__ = tool_name
    _topics_impl.__doc__ = (
        f"Authoritative list of topics accepted by this server's "
        f"dispatcher tool. Returns {len(topics_snapshot)} topics."
    )

    mcp.tool()(_topics_impl)

    # Dispatcher servers historically register topics after status. Apply the
    # same fleet contract immediately so registration order cannot leave this
    # late tool without annotations, title, metadata, or an output schema.
    from .mcp_contract import apply_tool_contract

    try:
        version = importlib.metadata.version("radia-mcp")
    except importlib.metadata.PackageNotFoundError:
        version = "unknown"
    apply_tool_contract(
        mcp,
        server_name=server_name,
        version=version,
        tool_prefix=tool_name,
        set_server_metadata=False,
    )
    return tool_name
