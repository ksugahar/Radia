"""MCP Server: radia_mcp.litz_transmission

Litz wire + transmission line knowledge for high-frequency conductor analysis.

Distilled from W:/.../39_部品材料/02_リッツ線/ (44 files) + 03_伝送線路/ (22 files).

Cross-references:
- `radia_mcp.peec.carstensen_ac_copper_loss` — Carstensen AC loss
- `radia_mcp.motor.hollaus_eddy` — Hollaus MSFEM
- `radia_mcp.pcb.coil_compensation` — Litz wire in WPT coils

Usage:
    mcp-server-litz-transmission              # stdio
    mcp-server-litz-transmission --selftest   # self-test
"""
import sys
from mcp.server.fastmcp import FastMCP
from ..common import register_status_tool, register_topics_tool
from .knowledge import get_knowledge, TOPICS

mcp = FastMCP("mcp-server-litz-transmission")


@mcp.tool()
def litz_transmission(topic: str = "overview") -> str:
    """
    Litz wire + transmission line knowledge.

    Args:
        topic: One of:
            "overview"      - Skin depth + AC loss summary (DEFAULT)
            "litz"          - Litz wire (Dowell, homogenization, magnetic-plated)
            "transmission"  - Multiconductor transmission lines (MTL)
            "all"           - Everything
    """
    return get_knowledge(topic)




register_status_tool(
    mcp,
    server_name='mcp-server-litz-transmission',
    description='Litz wire AC loss (Dowell, homogenization, magnetic-plated wire) + multiconductor transmission line theory',
    subpackage='radia_mcp.litz_transmission',
    related_servers=["peec", "ih", "pcb"],
)


register_topics_tool(
    mcp,
    server_name='mcp-server-litz-transmission',
    topics=TOPICS,
)


def main():
    if "--selftest" in sys.argv:
        print("litz_transmission MCP server self-test:")
        print(f"  knowledge: {len(get_knowledge('all'))} chars")
        print("OK")
        return
    mcp.run()


if __name__ == "__main__":
    main()
