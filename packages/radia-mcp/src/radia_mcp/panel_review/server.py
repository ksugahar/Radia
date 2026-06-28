"""
Panel Review MCP Server (radia_mcp.panel_review)

Surfaces the current Radia notebook panel review contract.  The old PySide6
panel review chain is retired; this server keeps the historical topic names as
compatibility redirects to the Jupyter notebook workbench route.

Why an MCP server (rather than just SKILL.md files):

- Future AI sessions / other AI agents working on Radia panels can
  query the notebook contract without first invoking a Skill or
  reading a markdown file -- the knowledge arrives in their context
  via a single `panel_review(topic="...")` call.
- Cross-team usability: someone editing a panel in a different
  codebase, helping a Radia user, or reviewing a Radia PR can pull
  the bug catalogue / 13-check list without local repo access.
- Versioning: bumping radia-mcp on PyPI ships the updated text
  everywhere (the SKILL.md files only update for cloners).

Usage:
    mcp-server-panel-review            # start MCP server (stdio transport)
    mcp-server-panel-review --selftest # run smoke check on all topics
"""

import sys

from mcp.server.fastmcp import FastMCP

from ..common import register_status_tool, register_topics_tool
from .panel_review_knowledge import (
    get_panel_review_documentation,
    TOPICS,
)

mcp = FastMCP("mcp-server-panel-review")


# ============================================================
# MCP Tool
# ============================================================

@mcp.tool()
def panel_review(topic: str = "overview") -> str:
    """
    Get Radia notebook panel review documentation.

    Historical topic names remain accepted, but all topics now point at the
    notebook workbench contract: DesignSpec, *_notebook Workbench wiring,
    result artifacts, validation_test, and no-PySide regression checks.

    Args:
        topic: Documentation topic.  Options:
            "overview"            - Notebook panel review contract
            "5_skills_chain"      - Compatibility alias to notebook contract
            "13_checks"           - Compatibility alias to notebook contract
            "bug_catalogue"       - Compatibility alias to notebook contract
            "val_checkbox_trap"   - Compatibility alias to notebook contract
            "map_value_reject"    - Compatibility alias to notebook contract
            "widget_calc_gap"     - Compatibility alias to notebook contract
            "smoke_scenarios"     - Compatibility alias to notebook contract
            "red_flags"           - Compatibility alias to notebook contract
            "workflow"            - Compatibility alias to notebook contract
            "all"                 - All topics concatenated
    """
    return get_panel_review_documentation(topic)


# ============================================================
# MCP Prompts
# ============================================================

@mcp.prompt()
def review_a_panel(panel_path: str = "src/radia/radia_ih.py") -> str:
    """Run a thorough review of a Radia notebook panel."""
    return (
        f"Please run a thorough notebook panel review for: {panel_path}\n\n"
        "Steps (use the panel_review MCP tool for reference text):\n"
        "1. panel_review(topic='overview') for the notebook contract.\n"
        "2. Inspect DesignSpec, Workbench.build_command(), and calc argparse.\n"
        "3. Run `python -m pytest validation_test/panels/test_notebook_workbench.py -q`.\n"
        "4. Confirm result-bearing notebook + sidecar JSON policy where docs are involved.\n"
        "5. Report BUG / RISK / NIT / OK with file/line references.\n\n"
        "Common gotchas:\n"
        " - JSON used as preset storage instead of run artifact\n"
        " - notebook imports PySide6/PyQt\n"
        " - Workbench argv drifts from calc_*.py argparse\n"
    )


# ============================================================
# Entry point
# ============================================================

# ============================================================
# Self-introspection (uniform with other radia_mcp servers)
# ============================================================

register_status_tool(
    mcp,
    server_name="mcp-server-panel-review",
    description="Radia notebook panel review contract (DesignSpec, "
                "Workbench, result artifacts, no-PySide gate)",
    subpackage="radia_mcp.panel_review",
    related_servers=["meta"],
    optional_deps=[],
)

register_topics_tool(
    mcp,
    server_name="mcp-server-panel-review",
    topics=TOPICS,
)


def main():
    if "--selftest" in sys.argv:
        print("panel-review MCP server self-test:")
        for topic in TOPICS:
            doc = panel_review(topic)
            assert len(doc) > 50, f"topic {topic!r} returned <50 chars"
            print(f"  panel_review({topic!r:22}) -> {len(doc):5d} chars")
        # Also verify the 'all' synonym
        doc_all = panel_review("all")
        assert len(doc_all) > sum(len(v) for v in TOPICS.values()) - 100
        print(f"  panel_review('all'                  ) -> "
              f"{len(doc_all):5d} chars")
        # Verify unknown topic returns help text
        unk = panel_review("not-a-real-topic")
        assert "Unknown topic" in unk
        assert "Available topics" in unk
        print(f"  panel_review('not-a-real-topic')     -> "
              f"unknown-topic help text emitted")
        print("  PASSED")
        return

    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
