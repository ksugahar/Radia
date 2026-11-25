"""MCP server for NGSolve weak-form variational formulations (ngsolve-mcp.weak-form)."""

from __future__ import annotations

import argparse
import asyncio

from mcp.server.fastmcp import FastMCP

from radia_mcp.fem.fem_reference_topics import WEAK_FORM_TOPICS as _TOPICS

from radia_mcp.common import register_status_tool
from radia_mcp.meta.catalog import CATALOG

mcp = FastMCP("ngsolve-mcp.weak-form")


@mcp.tool()
def weak_form_usage(topic: str = "overview") -> str:
    """Return NGSolve weak-form formulation guidance for the given topic."""
    key = topic.strip().lower()
    if key in _TOPICS:
        return _TOPICS[key]
    return f"Unknown topic '{topic}'. Available: {', '.join(sorted(_TOPICS))}"


# Registered after every tool: the status registration applies the
# annotation contract to the tools that exist at this point.
register_status_tool(
    mcp,
    "ngsolve-mcp.weak-form",
    "NGSolve weak-form variational formulations (HCurl/HDiv/H1, BilinearForm, LinearForm, SymbolicBFI).",
    "ngsolve_mcp.weak_form",
    related_servers=CATALOG["ngsolve-mcp.weak-form"]["related"],
)


def _selftest() -> None:
    tools = asyncio.run(mcp.list_tools())
    if "weak_form_usage" not in {t.name for t in tools}:
        raise SystemExit("weak_form_usage tool missing")
    print("ngsolve-mcp.weak-form selftest passed")


def main() -> None:
    parser = argparse.ArgumentParser(description="ngsolve-mcp.weak-form MCP server")
    parser.add_argument("--selftest", action="store_true", help="Verify tool registration")
    if parser.parse_args().selftest:
        from radia_mcp.common.utf8_stdout import use_utf8_stdout

        use_utf8_stdout()
        _selftest()
        return
    from radia_mcp.common.windows_stdio import run_stdio

    run_stdio(mcp)


if __name__ == "__main__":
    main()
