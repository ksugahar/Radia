"""MCP Server: radia_mcp.fusion_reactor

Fusion reactor magnet knowledge.

Distilled from W:/.../99_アプリケーション/08_核融合/ (52 files, 1.6 GB).

Cross-references:
- `radia_mcp.accelerator` — adjacent SC magnet expertise
- `radia_mcp.electromagnet` — DC magnet Hantila
- `radia_mcp.fem.potential_formulations.h_formulation` — HCurl H for SC

Usage:
    mcp-server-fusion-reactor              # stdio
    mcp-server-fusion-reactor --selftest   # self-test
"""
import sys
from mcp.server.fastmcp import FastMCP
from ..common import register_status_tool, register_topics_tool
from .knowledge import get_knowledge, TOPICS

mcp = FastMCP("mcp-server-fusion-reactor")


@mcp.tool()
def fusion_reactor(topic: str = "overview") -> str:
    """
    Fusion reactor magnet knowledge.

    Args:
        topic: One of:
            "overview"     - Confinement landscape (DEFAULT)
            "tokamak"      - Tokamak (ITER) TF/PF/CS coil system
            "stellarator"  - Stellarator (LHD, W7-X, Mitsubishi lineage)
            "all"          - Everything
    """
    return get_knowledge(topic)


register_status_tool(
    mcp,
    server_name='mcp-server-fusion-reactor',
    description='Fusion reactor magnets: tokamak ITER + stellarator LHD/W7-X/heliotron lineage',
    subpackage='radia_mcp.fusion_reactor',
    related_servers=["accelerator", "electromagnet"],
)

register_topics_tool(
    mcp,
    server_name='mcp-server-fusion-reactor',
    topics=TOPICS,
)


def main():
    if "--selftest" in sys.argv:
        print("fusion-reactor MCP server self-test:")
        print(f"  knowledge: {len(get_knowledge('all'))} chars")
        print("OK")
        return
    mcp.run()


if __name__ == "__main__":
    main()
