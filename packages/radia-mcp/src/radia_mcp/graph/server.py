"""
Graph MCP Server (radia_mcp.graph)

Surfaces the Sugahara Lab publication-figure style guide via two MCP
tools:

    graph_style_guide(target='all'|'paper_single_column'|...)
        Returns the lab-standard rules (Times New Roman, units in
        parentheses, no in-figure title, font sizes per embed width,
        IEEE/IEEJ conventions) plus the numeric profile for the
        selected target.

    graph_size_for_target(target, embed_width_cm=None)
        Computes the recommended figure size + font sizes + a
        ready-to-paste MATLAB or Matplotlib snippet for the chosen
        embedding width.

Promoted on 2026-05-26 from
  s:/mcp-server/src/mcp_server_document/graph/
into radia-mcp following the standard `radia_mcp.<topic>.server`
pattern (statusable, --selftest-able, discoverable via the meta
catalog).

The companion Python helpers (apply_lab_style, lab_figsize,
lab_savefig, tighten_margins, label_curve_endpoints, add_slope_guide,
check_legend_overlap, find_best_legend_loc, plot_asymptote_ratio_sweep,
plot_basis_size_convergence) are NOT registered as MCP tools — they
have richer signatures (matplotlib objects in/out) than MCP-JSON can
carry.  They live in `radia_mcp.graph.tools` for direct import by
analysis scripts.

Usage:
    mcp-server-graph              # start MCP server (stdio transport)
    mcp-server-graph --selftest   # smoke check every profile in both tools
"""

import sys

from mcp.server.fastmcp import FastMCP

from ..common import register_status_tool

from . import tools as _tools

mcp = FastMCP("mcp-server-graph")


# ============================================================
# Tool registration
# ============================================================
# Auto-register the MCP-callable subset (graph_*).  The Python helpers
# (apply_lab_style etc.) have matplotlib-typed signatures that don't
# round-trip through MCP-JSON, so we skip them deliberately.
_REGISTERED: list[str] = []
for _name in dir(_tools):
    if _name.startswith("graph_"):
        _fn = getattr(_tools, _name)
        if callable(_fn):
            mcp.tool()(_fn)
            _REGISTERED.append(_name)


# ============================================================
# Self-introspection (uniform with other radia_mcp servers)
# ============================================================

register_status_tool(
    mcp,
    server_name="mcp-server-graph",
    description="Sugahara Lab publication-figure style guide: "
                "IEEE / IEEJ font/size profiles, MATLAB + Matplotlib "
                "snippets, lab style rules (units in parentheses, no "
                "in-figure title, Times New Roman serif).",
    subpackage="radia_mcp.graph",
    related_servers=["mathematica", "literature-index"],
    # matplotlib is needed for the python helper functions
    # (apply_lab_style / lab_figsize / etc.), NOT for the MCP tools.
    # Server loads without it.
    optional_deps=["matplotlib"],
)


# ============================================================
# Entry point
# ============================================================

def main():
    if "--selftest" in sys.argv:
        print("graph MCP server self-test:")
        print(f"  Registered tools ({len(_REGISTERED)}):")
        for name in sorted(_REGISTERED):
            print(f"    - {name}")
        # Smoke-test both tools across every profile so a typo in the
        # _PROFILES dict immediately surfaces in CI.
        profiles = sorted(_tools._PROFILES.keys())
        print(f"  Profiles ({len(profiles)}):")
        for prof in profiles:
            guide = _tools.graph_style_guide(prof)
            size = _tools.graph_size_for_target(prof)
            assert len(guide) > 200, (
                f"graph_style_guide({prof!r}) returned only "
                f"{len(guide)} chars (suspiciously short)"
            )
            assert "Unknown target" not in guide, (
                f"graph_style_guide({prof!r}) reported unknown target"
            )
            assert "Unknown target" not in size, (
                f"graph_size_for_target({prof!r}) reported unknown target"
            )
            print(f"    {prof:42s} guide={len(guide):4d} ch, "
                  f"size={len(size):4d} ch")
        # Sanity: 'all' returns everything
        full = _tools.graph_style_guide("all")
        assert all(p in full for p in profiles), \
            "graph_style_guide('all') missing some profiles"
        # Sanity: unknown target returns help text (not crash)
        unk = _tools.graph_style_guide("not-a-real-profile")
        assert "Unknown target" in unk and "Valid:" in unk
        print(f"  graph_style_guide('all')               -> "
              f"{len(full):5d} chars")
        print(f"  graph_style_guide('not-a-real-profile') -> "
              f"unknown-target help text emitted")
        print("  PASSED")
        return

    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
