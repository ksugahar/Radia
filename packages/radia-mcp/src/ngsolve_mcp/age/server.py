"""MCP server for Air-Gap Element (AGE) rotating-machine coupling (ngsolve-mcp.age)."""

from __future__ import annotations

import argparse
import asyncio

from mcp.server.fastmcp import FastMCP

from radia_mcp.common import register_status_tool
from radia_mcp.motor.age import AGE_TOPICS as _MOVED

mcp = FastMCP("ngsolve-mcp.age")

_TOPICS = {
    "overview": (
        'ngsolve-mcp.age is a legacy route. Its topics are now served by radia-mcp.motor: motor_age_quality(topic=<topic>), with the same text.'
        + "\nTopics: " + ", ".join(_MOVED)
    ),
    **_MOVED,
}


@mcp.tool()
def age_usage(topic: str = "overview") -> str:
    """Return Air-Gap Element (AGE) guidance for the given topic."""
    key = topic.strip().lower()
    if key in _TOPICS:
        return _TOPICS[key]
    return f"Unknown topic '{topic}'. Available: {', '.join(sorted(_TOPICS))}"


# Registered after every tool: the status registration applies the
# annotation contract to the tools that exist at this point.
register_status_tool(
    mcp,
    "ngsolve-mcp.age",
    "Air-Gap Element (AGE) analytical/harmonic rotor-stator coupling for NGSolve.",
    "ngsolve_mcp.age",
    related_servers=["radia-mcp.ngsolve", "fem", "motor"],
)


def _selftest() -> None:
    tools = asyncio.run(mcp.list_tools())
    if "age_usage" not in {t.name for t in tools}:
        raise SystemExit("age_usage tool missing")
    print("ngsolve-mcp.age selftest passed")


def main() -> None:
    parser = argparse.ArgumentParser(description="ngsolve-mcp.age MCP server")
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
