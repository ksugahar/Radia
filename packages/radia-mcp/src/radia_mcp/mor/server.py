"""
MOR (Model Order Reduction) MCP Server (radia_mcp.mor)

Knowledge layer for model order reduction: the systematic MOR taxonomy
(projection / POD / Krylov-PRIMA / system-theoretic / data-driven,
hyperreduction) and a literature-corpus search.

Usage:
    mcp-server-mor              # Start MCP server (stdio)
    mcp-server-mor --selftest   # Self-test
"""

import sys

from mcp.server.fastmcp import FastMCP
from radia_mcp.common import register_status_tool

from .systematic_knowledge import get_systematic_mor_knowledge

mcp = FastMCP("mcp-server-mor")


@mcp.tool()
def mor_systematic(topic: str = "mor_taxonomy") -> str:
    """
    Systematic MOR knowledge -- distilled from the deGruyter 3-volume
    MOR Handbook (1221 pages) + Kiss-Orosz 2024 Energies rotating-
    machines review + Ioan Vol3 Ch5 EM-specific MOR.

    Args:
        topic: One of:
            "mor_taxonomy"          - Three families: projection / system / data
            "projection_pod_rb"     - POD + Reduced Basis (Vol 2 Ch 2, 4)
            "projection_krylov"     - Krylov / Arnoldi / Lanczos (Vol 1 Ch 3)
            "pgd"                   - Proper Generalized Decomposition (Vol 2 Ch 3)
            "hyperreduction"        - DEIM / EIM for nonlinear (Vol 2 Ch 5)
            "system_theoretic_bt"   - Balanced Truncation, H2/H_inf (Vol 1 Ch 2)
            "data_driven_dmd_oi"    - DMD, OI, Loewner (Vol 1 Ch 6 + Vol 2 Ch 7)
            "parametric_pmor"       - Parametric MOR (Vol 2 Ch 1)
            "em_specific_ioan"      - Ioan Vol 3 Ch 5 (56p EM-specific)
            "rotating_machines_kiss_orosz_2024" - 2024 review for motors
            "software_lab"          - pyMOR / MORLab / lab tools
            "lab_recommendation"    - Decision guide for radia + NGSolve
            "all"                   - Everything
    """
    return get_systematic_mor_knowledge(topic)


@mcp.tool()
def mor_bibliography(query: str = "") -> str:
    """Search the local literature corpus (RADIA_LIT_ROOT) by filename keywords."""
    from radia_mcp.literature_index.index_tools import lit_search
    return lit_search(query)


register_status_tool(
    mcp,
    server_name='mcp-server-mor',
    description='Model Order Reduction: POD, PRIMA/Krylov, hyperreduction (DEIM)',
    subpackage='radia_mcp.mor',
    related_servers=["radia-ngsolve", "rna-mec"],
)


def main():
    if "--selftest" in sys.argv:
        print("MOR MCP server self-test:")
        from .systematic_knowledge import SECTIONS as S
        for k in S:
            r = mor_systematic(k)
            assert len(r) > 500, f"{k} short ({len(r)})"
            print(f"  systematic({k!r}): {len(r)} chars")
        print("  PASSED")
        return
    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
