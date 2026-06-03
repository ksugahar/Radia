"""MCP Server: radia_mcp.levitation

Magnetic levitation FORCE-physics + stationary / industrial / laboratory
levitation knowledge.

The force-physics counterpart to `radia_mcp.maglev` (transport):
pick by intent -- move a vehicle along a track -> maglev; suspend,
melt, bear, or COMPUTE a levitation force -> levitation.

Cross-references:
- `radia_mcp.maglev` -- maglev systems (EMS/EDS trains, SCMaglev, wheels) + Radia-IEM/CLN research
- `radia_mcp.ih` -- induction heating (EML = levitation + IH together)
- `radia_mcp.peec` -- eddy-current / SIBC engine for induction lift
- `radia_mcp.electromagnet` -- DC magnet / pole-face B (EMS & AMB)
- `radia_mcp.differential_forms` (em_force_recipe) -- rigorous EM force
- `radia_mcp.team_benchmark` (force_motion.problem_28) -- EDS lift benchmark

Usage:
    mcp-server-levitation              # stdio
    mcp-server-levitation --selftest   # self-test
"""
import sys
from mcp.server.fastmcp import FastMCP
from ..common import register_status_tool, register_topics_tool
from .knowledge import get_knowledge, TOPICS

mcp = FastMCP("mcp-server-levitation")


@mcp.tool()
def levitation(topic: str = "overview") -> str:
    """
    Magnetic levitation FORCE-physics + stationary levitation knowledge.

    Args:
        topic: One of:
            "overview"             - Levitation-force taxonomy + Radia fit (DEFAULT)
            "induction_levitation" - Eddy-current lift: jumping ring, circuit model
            "eml_melting"          - Electromagnetic levitation MELTING (ties to IH)
            "magnetic_bearings"    - AMB force-current-displacement, neg. stiffness
            "superconducting"      - Meissner vs flux pinning, HTS bulk, frozen image
            "diamagnetic"          - grad(B^2) levitation, graphite / water-frog
            "earnshaw_stability"   - Earnshaw's theorem + its 5 loopholes
            "force_computation"    - Maxwell stress / virtual work / time-avg J x B
            "benchmarks"           - TEAM 28, jumping-ring analytic, EML lift coeff
            "all"                  - Everything
    """
    return get_knowledge(topic)


register_status_tool(
    mcp,
    server_name='mcp-server-levitation',
    description='Magnetic levitation FORCE-physics: induction/eddy-current lift, '
                'EML melting (ties to IH), magnetic bearings (AMB), superconducting '
                '(Meissner/pinning), diamagnetic, Earnshaw + loopholes, force '
                'computation. Force-physics counterpart to maglev (transport).',
    subpackage='radia_mcp.levitation',
    related_servers=["maglev", "ih", "peec", "electromagnet"],
)


register_topics_tool(
    mcp,
    server_name='mcp-server-levitation',
    topics=TOPICS,
)


def main():
    if "--selftest" in sys.argv:
        print("levitation MCP server self-test:")
        n_all = len(get_knowledge('all'))
        print(f"  knowledge: {n_all} chars")
        # Every TOPICS key must dispatch to real content (no unknown-topic).
        bad = []
        for tk in TOPICS:
            head = get_knowledge(tk)[:80].lower()
            if "unknown topic" in head or "not found" in head:
                bad.append(tk)
        if bad:
            print(f"  FAIL: topics returned unknown-topic: {bad}")
            sys.exit(1)
        print(f"  topics: {len(TOPICS)} dispatch OK")
        print("PASSED")
        return
    mcp.run()


if __name__ == "__main__":
    main()
