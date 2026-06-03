"""MCP Server: radia_mcp.maglev  (magnetic levitation)

Magnetic levitation knowledge -- transport / suspension / control.  The
headline content is the lab's own Radia-based maglev research line
(CAE-AI Lab, Yano + Sugahara): Radia IEM <-> reduced-potential FEM weak
coupling for moving-magnet eddy-current force (topic radia_iem_fem) and
Cauer Ladder Network model-order reduction for real-time control-coupled
maglev (topic cln_mor_control).  The former linear-drive material
(LIM/LSM, end effects) was removed.

For the levitation-FORCE physics (induction lift, EML melting, magnetic
bearings, superconducting, diamagnetic, Earnshaw, force computation) see
the sibling server `radia_mcp.levitation`.

Cross-references:
- `radia_mcp.levitation` — levitation-force physics + bearings + EML
- `radia_mcp.mor` (mor_cln) — Cauer Ladder Network MOR theory
- `radia_mcp.fem` (potential_formulations) — A-phi / T-Omega / A-T gauges
- `radia_mcp.team_benchmark.force_motion.problem_28` — TEAM 28 levitation benchmark
- `radia_mcp.motor` — analogous rotary motor knowledge

Usage:
    mcp-server-maglev              # stdio
    mcp-server-maglev --selftest   # self-test
"""
import sys
from mcp.server.fastmcp import FastMCP
from ..common import register_status_tool, register_topics_tool
from .knowledge import get_knowledge, TOPICS

mcp = FastMCP("mcp-server-maglev")


@mcp.tool()
def maglev(topic: str = "overview") -> str:
    """
    Magnetic levitation knowledge (transport / suspension / control).

    Args:
        topic: One of:
            "overview"            - Maglev landscape + the lab's research (DEFAULT)
            "radia_iem_fem"       - Radia IEM <-> reduced-potential FEM weak coupling (Yano)
            "cln_mor_control"     - Cauer Ladder Network MOR for control-coupled maglev (Yano)
            "pm_maglev_zero_power"- Passive PM levitation, Maxwell-Earnshaw
            "eddy_current_maglev" - Eddy-current EDS, Kansai 2D model, Arago
            "sumitomo_heavy_industrial" - JP 7-327337 PM bearing + JP 2007-215264 mover
            "kansai_research"     - Saiki/Fujii magnetic-wheel lineage
            "scmaglev_eds"        - SCMaglev (Chuo Shinkansen) -- SC-EDS levitation
            "halbach_arrays"      - Halbach + Inductrack
            "all"                 - Everything
    """
    return get_knowledge(topic)




register_status_tool(
    mcp,
    server_name='mcp-server-maglev',
    description='Magnetic levitation (EMS/EDS/PM/SC/Halbach). Lab research line: Radia IEM<->FEM weak coupling for moving-magnet eddy-current force + Cauer Ladder Network MOR for control-coupled maglev (Yano, CAE-AI). Sibling: levitation (force physics).',
    subpackage='radia_mcp.maglev',
    related_servers=["levitation", "mor", "motor"],
)


register_topics_tool(
    mcp,
    server_name='mcp-server-maglev',
    topics=TOPICS,
)


def main():
    if "--selftest" in sys.argv:
        print("maglev MCP server self-test:")
        print(f"  knowledge: {len(get_knowledge('all'))} chars")
        print("OK")
        return
    mcp.run()


if __name__ == "__main__":
    main()
