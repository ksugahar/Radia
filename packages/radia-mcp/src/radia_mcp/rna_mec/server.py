"""MCP Server: radia_mcp.rna_mec

Reluctance Network Analysis (RNA) / Magnetic Equivalent Circuit (MEC) knowledge.

Distilled from public-safe curated corpus (17 files, lab Sugahara 田中/羽根 lineage).

Cross-references:
- `radia_mcp.radia_ngsolve.hdiv_vim` — magnetic-material field coupling
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
            "overview"                       - RNA / MEC landscape + lab lineage (DEFAULT)
            "mec_basics"                     - Reluctance, permeance, MMF, Hopkinson law
            "nodal_vs_mesh_analysis"         - Derbas 2009: KCL vs KVL, Jacobian
            "reluctance_network_construction" - flux tube, claw-pole example
            "lumped_extraction_fea"          - L, M from FEA (Lee 2005 TEAM-28)
            "cauer_ladder_rna"               - CLN for eddy-current MOR (Kameari 2018)
            "rna_magnetic_coupling"               - Janet 2004-2005 mixed method, CT
            "electromechanical_coupling"     - State-space RNA + ODE (RK4)
            "team28_reduced_model"           - TEAM-28 in depth, 85h -> 1h speedup
            "topology_optimization"          - Yin 2023 grid RNA + AVM
            "dynamic_hysteresis"             - Lab Play + Cauer dynamic MEC
            "vs_pec_peec"                    - Acronym map: RNA/MEC/PEEC/FEM/HDiv-VIM
            "all"                            - Everything
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
