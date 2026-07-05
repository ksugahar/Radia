"""MCP Server: radia_mcp.metamaterial

Electromagnetic metamaterial knowledge:
- Left-handed / negative-index / double-negative (DNG) materials
- Split-ring resonator (SRR), wire-medium homogenization
- Composite right/left-handed (CRLH) transmission lines
- Magneto-optical metamaterials (Faraday rotation enhancement)
- Transformation optics and the Kelvin-transform connection
- Topology optimization for EM material design
- Perfect lens / evanescent amplification (Pendry vs Bergamin)
- Gyrator-based active LH lines (lab thesis)

Distilled from PDFs in
  public-safe curated corpus
  public-safe curated corpus
  public-safe curated corpus

Cross-references:
- `radia_mcp.electromagnet` -- Kelvin transform (same coordinate-
  transformation family as transformation optics cloaks/lenses)
- `radia_mcp.peec` -- CRLH unit cell as PEEC ladder
- `radia_mcp.ih` -- ESIM and effective surface impedance
- `radia_mcp.radia_ngsolve` -- periodic FES, BEM, analytical formulas

Usage:
    mcp-server-metamaterial              # stdio
    mcp-server-metamaterial --selftest   # self-test
"""

import sys

from mcp.server.fastmcp import FastMCP
from ..common import register_status_tool, register_topics_tool

from .knowledge import get_knowledge, TOPICS

try:
    from .bibliography_index_knowledge import get_bibliography_index
except ImportError:
    def get_bibliography_index(query: str = "") -> str:
        return ("Metamaterial bibliography index not yet generated. "
                "Run the auto-indexer.")


mcp = FastMCP("mcp-server-metamaterial")


@mcp.tool()
def metamaterial(topic: str = "overview") -> str:
    """
    Electromagnetic metamaterial knowledge.

    Args:
        topic: One of (see TOPICS for the full enum):
            "overview"                     - Field map (DEFAULT)
            "lh_materials_neg_index"       - DNG / LH / Vesalago
            "srr_split_ring"               - SRR equivalent circuit
            "effective_medium"             - Homogenization, MOM tensor
            "transformation_optics_kelvin" - Pendry/Leonhardt + Kelvin
            "magneto_optical"              - Faraday enhancement (MOM)
            "topology_optimization_design" - Density-method TO for EM
            "crlh_transmission_line"       - Caloz-Itoh CRLH cell
            "gyrator_lh_line"              - Active LH TL (lab thesis)
            "perfect_lens_evanescent"      - Pendry vs Bergamin caveat
            "applications_automotive"      - Toyota CRDL roadmap
            "acoustic_sonic_crystal"       - Sonic crystals / band gaps
            "all"                          - Everything

        Use the `metamaterial_topics` tool to list authoritative
        topic names.
    """
    return get_knowledge(topic)


@mcp.tool()
def metamaterial_bibliography(query: str = "") -> str:
    """Search the metamaterial bibliography catalog of cited PDFs."""
    return get_bibliography_index(query)


register_status_tool(
    mcp,
    server_name='mcp-server-metamaterial',
    description=(
        'Metamaterials: DNG/LH, SRR, CRLH, effective medium, '
        'magneto-optical, transformation optics (Kelvin link), '
        'topology-optimization design, perfect-lens evanescent '
        'subtlety, automotive applications'
    ),
    subpackage='radia_mcp.metamaterial',
    related_servers=[
        "electromagnet",        # Kelvin transform
        "peec",                 # CRLH as PEEC ladder
        "ih",                   # ESIM surface impedance
        "radia_ngsolve",        # periodic FES + BEM
        "pcb",                  # WPT-side ladder networks
        "litz-transmission",    # adjacent TL topic
    ],
)

register_topics_tool(
    mcp,
    server_name='mcp-server-metamaterial',
    topics=TOPICS,
)


def main():
    if "--selftest" in sys.argv:
        print("metamaterial MCP server self-test:")
        print(f"  topics:    {len(TOPICS)}")
        print(f"  knowledge: {len(get_knowledge('all'))} chars")
        n_bib = len(get_bibliography_index(""))
        print(f"  bibliography: {n_bib} chars")
        # Round-trip every topic to make sure the dispatcher is wired
        for name in TOPICS:
            body = get_knowledge(name)
            assert body and len(body) > 100, f"topic '{name}' too short"
        print("  all topics dispatch: OK")
        print("PASSED")
        return
    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
