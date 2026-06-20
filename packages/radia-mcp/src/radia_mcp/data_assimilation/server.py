"""
Data Assimilation MCP Server (radia_mcp.data_assimilation)

Kalman / EnKF / 4D-Var for EM state estimation + sensor fusion.

Usage:
    mcp-server-data-assimilation              # Start MCP server
    mcp-server-data-assimilation --selftest   # Self-test
"""

import sys

from mcp.server.fastmcp import FastMCP
from ..common import register_status_tool, register_topics_tool

from .knowledge import get_data_assimilation_documentation, TOPICS

mcp = FastMCP("mcp-server-data-assimilation")


@mcp.tool()
def data_assimilation(topic: str = "all") -> str:
    """Data assimilation for EM state estimation + sensor fusion.

    Topics:
        "all"            - Everything
        "overview"       - DA family tree + decision table per EM use
        "kalman_family"  - KF / EKF / UKF / EnKF; magnet tracking,
                            eddy current state, hysteresis state
        "four_dvar"      - 4D-Var with TLM + adjoint; window-based
                            smoothing for time-evolving EM systems
    """
    return get_data_assimilation_documentation(topic)




register_status_tool(
    mcp,
    server_name='mcp-server-data-assimilation',
    description='Kalman / EnKF / 4D-Var for EM state estimation + sensor fusion',
    subpackage='radia_mcp.data_assimilation',
    related_servers=["mor", "fusion-reactor"],
    optional_deps=["filterpy"],
)


register_topics_tool(
    mcp,
    server_name='mcp-server-data-assimilation',
    topics=TOPICS,
)


def main():
    if "--selftest" in sys.argv:
        print("Data Assimilation MCP server self-test:")
        for t in ["overview", "kalman_family", "four_dvar", "all"]:
            doc = data_assimilation(t)
            if len(doc) < 200:
                print(f"  FAIL data_assimilation('{t}'): only {len(doc)} chars")
                sys.exit(1)
            print(f"  data_assimilation('{t}'): {len(doc)} chars OK")
        print("  PASSED")
        return
    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
