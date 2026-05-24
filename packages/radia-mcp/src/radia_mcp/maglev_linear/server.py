"""MCP Server: radia_mcp.maglev_linear

Magnetic levitation + linear drive knowledge.

Distilled from W:/.../99_アプリケーション/07_磁気浮上/ + 09_リニアドライブ/.

Cross-references:
- `radia_mcp.wpt.applications.robot_bearingless` — bearingless motor + WPT
- `radia_mcp.motor` — analogous rotary motor knowledge
- `radia_mcp.team_benchmark.force_motion.problem_28` — levitation benchmark

Usage:
    mcp-server-maglev-linear              # stdio
    mcp-server-maglev-linear --selftest   # self-test
"""
import sys
from mcp.server.fastmcp import FastMCP
from ..common import register_status_tool, register_topics_tool
from .knowledge import get_knowledge, TOPICS

mcp = FastMCP("mcp-server-maglev-linear")


@mcp.tool()
def maglev_linear(topic: str = "overview") -> str:
    """
    Magnetic levitation + linear drive knowledge.

    Args:
        topic: One of:
            "overview" - Levitation + linear drive landscape (DEFAULT)
            "maglev"   - Levitation (EMS/EDS/SCMaglev/Halbach/bearingless ★)
            "linear"   - Linear induction (LIM) + linear synchronous (LSM)
            "all"      - Everything
    """
    return get_knowledge(topic)




register_status_tool(
    mcp,
    server_name='mcp-server-maglev-linear',
    description='Maglev (EMS/EDS/SCMaglev/Halbach/bearingless ★) + linear drives (LIM/LSM). Lab specialty: bearingless + WPT',
    subpackage='radia_mcp.maglev_linear',
    related_servers=["motor", "wpt"],
)


register_topics_tool(
    mcp,
    server_name='mcp-server-maglev-linear',
    topics=TOPICS,
)


def main():
    if "--selftest" in sys.argv:
        print("maglev_linear MCP server self-test:")
        print(f"  knowledge: {len(get_knowledge('all'))} chars")
        print("OK")
        return
    mcp.run()


if __name__ == "__main__":
    main()
