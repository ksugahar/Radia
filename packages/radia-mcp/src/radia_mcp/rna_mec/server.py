"""MCP Server: radia_mcp.rna_mec

Reluctance Network Analysis (RNA) / Magnetic Equivalent Circuit (MEC) knowledge.

Distilled from W:/.../13_RNA/ (17 files, lab Sugahara 田中/羽根 lineage).

Cross-references:
- `radia_mcp.bem.mmm_msc.rna_mmm` — RNA + MMM coupling
- `radia_mcp.magnetic_materials.hysteresis_models.lab_core` — Play model
- `radia_mcp.mor.systematic.cln` — Cauer ladder for eddy

Usage:
    mcp-server-rna-mec              # stdio
    mcp-server-rna-mec --selftest   # self-test
"""
import sys
from mcp.server.fastmcp import FastMCP
from ..common import register_status_tool, register_topics_tool
from .knowledge import get_knowledge, TOPICS

mcp = FastMCP("mcp-server-rna-mec")


@mcp.tool()
def rna_mec(topic: str = "overview") -> str:
    """
    Reluctance Network Analysis / Magnetic Equivalent Circuit.

    Args:
        topic: One of:
            "overview"            - RNA / MEC landscape + lab lineage (DEFAULT)
            "rna_mmm_coupling"    - RNA + MMM for transformer
            "dynamic_hysteresis"  - ★ Dynamic hysteresis MEC (Play + Cauer)
            "all"                 - Everything
    """
    return get_knowledge(topic)




register_status_tool(
    mcp,
    server_name='mcp-server-rna-mec',
    description='RNA / Magnetic Equivalent Circuit. ★ Lab specialty: dynamic hysteresis MEC (Play + Cauer)',
    subpackage='radia_mcp.rna_mec',
    related_servers=["magnetic-materials", "mor", "ih"],
)


register_topics_tool(
    mcp,
    server_name='mcp-server-rna-mec',
    topics=TOPICS,
)


def main():
    if "--selftest" in sys.argv:
        print("rna_mec MCP server self-test:")
        print(f"  knowledge: {len(get_knowledge('all'))} chars")
        print("OK")
        return
    mcp.run()


if __name__ == "__main__":
    main()
