"""
NDT MCP Server (radia_mcp.ndt) -- Non-Destructive Evaluation knowledge.

`ndt(topic=...)` dispatches the curated knowledge.py. The accepted topic
enum is published via `ndt_topics()` (auto-registered by
register_topics_tool). Search a local literature corpus with the
`ndt_bibliography` (backed by the literature index).
"""

import sys

from mcp.server.fastmcp import FastMCP
from radia_mcp.common import register_status_tool, register_topics_tool
from .knowledge import get_knowledge, TOPICS

mcp = FastMCP("mcp-server-ndt")


@mcp.tool()
def ndt(topic: str = "overview") -> str:
    """
    Electromagnetic non-destructive testing (NDT / NDE) knowledge.

    Covers eddy current testing (ECT), magnetic flux leakage (MFL),
    microwave / EMI metal detection, and the underlying FEM / BEM /
    hysteresis / open-boundary infrastructure.

    Args:
        topic: One of:
            "overview"                  - NDT landscape (DEFAULT)
            "eddy_current_basics"       - Skin depth, impedance plane
            "probe_types"               - Pancake/bobbin/split-D/array
            "defect_models"             - Crack / void / back-side
            "signal_processing"         - Phase, multi-freq, POD
            "fem_formulation"           - A-V / T-Omega / h / edge elements
            "integral_equation"         - BEM, VIE, LIE for cracks
            "multiply_connected"        - Cut surfaces, tree-cotree
            "open_boundary"             - PML, Kelvin, ABC, FEM-BEM
            "hysteresis_in_ndt"         - Vector Play / Stop / Energy
            "mfl_magnetic_flux_leakage" - DC + Hall arrays for pipes / rails
            "microwave_emi_landmine"    - GPR, metal detector, NMR MOUSE
            "ml_for_ndt"                - CNN / U-Net / PINN for NDT
            "all"                       - Everything concatenated
    """
    return get_knowledge(topic)


@mcp.tool()
def ndt_bibliography(query: str = "") -> str:
    """Search the local literature corpus (RADIA_LIT_ROOT) by filename keywords."""
    from radia_mcp.literature_index.index_tools import lit_search
    return lit_search(query)


register_status_tool(
    mcp,
    server_name='mcp-server-ndt',
    description='Non-destructive testing: eddy current testing, '
                'magnetic flux leakage, MFL signal analysis',
    subpackage='radia_mcp.ndt',
    related_servers=["ih", "magnetic-materials"],
)

register_topics_tool(
    mcp,
    server_name='mcp-server-ndt',
    topics=TOPICS,
)


def main():
    if "--selftest" in sys.argv:
        n_topics = len(TOPICS)
        n_chars = len(get_knowledge("all"))
        print(f"NDT MCP server self-test:")
        print(f"  topics:       {n_topics}")
        print(f"  knowledge:    {n_chars} chars")
        print("OK")
        return
    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
