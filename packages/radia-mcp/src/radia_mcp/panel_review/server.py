"""
Panel Review MCP Server (radia_mcp.panel_review)

Surfaces Radia's GUI panel review skill chain (panel-cli-diff,
panel-review, panel-qt-test, panel-preview, panel-smoke) and the
catalogue of known panel bug classes through one queryable MCP tool.

Why an MCP server (rather than just SKILL.md files):

- Future AI sessions / other AI agents working on Radia panels can
  query the right check / fix without first invoking a Skill or
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
    Get Radia GUI panel review skill-chain documentation and bug catalogue.

    Surfaces the 5-skill chain (panel-cli-diff / panel-review /
    panel-qt-test / panel-preview / panel-smoke) and the known bug
    patterns (val()-truthy trap, AMG-class map-value typos, widget-vs-
    calc Stage-1 gap, widget-default-x-calc-limitation, etc.).

    Args:
        topic: Documentation topic.  Options:
            "overview"            - High-level: 5-skill chain + glossary
            "5_skills_chain"      - When-to-use of each skill in detail
            "13_checks"           - Full checklist used by panel-review
                                     (checks 11-13 added 2026-05-24)
            "bug_catalogue"       - Known bug patterns (A..H), each
                                     linked to its detection skill
            "val_checkbox_trap"   - Deep dive on ModePanel.val() returning
                                     string "1"/"0" for QCheckBox
                                     (both truthy)
            "map_value_reject"    - Deep dive on _XXX_MAP × argparse
                                     choices mismatch (AMG-to-AMS class)
            "widget_calc_gap"     - Deep dive on Stage-1 (solver-switch
                                     pinning) policy + how today's IGTE
                                     paper widget gap happened
            "smoke_scenarios"     - panel-smoke skill format + workflow
            "red_flags"           - Triage table for user-reported panel
                                     issues
            "workflow"            - Recommended skill ordering before
                                     shipping a panel mode
            "all"                 - All topics concatenated
    """
    return get_panel_review_documentation(topic)


# ============================================================
# MCP Prompts
# ============================================================

@mcp.prompt()
def review_a_panel(panel_path: str = "src/radia/radia_ih.py") -> str:
    """Run a thorough review of a Radia Layer-3 panel."""
    return (
        f"Please run a thorough panel review for: {panel_path}\n\n"
        "Steps (use the panel_review MCP tool for reference text):\n"
        "1. panel_review(topic='5_skills_chain') for the workflow order.\n"
        "2. panel_review(topic='13_checks') for the static review list.\n"
        "3. Run `python tests/panels/check_panel_cli.py --panel "
        f"{panel_path.split('/')[-1]}` for static argv ↔ argparse\n"
        "   matching (Check 4 catches map-value rejects).\n"
        "4. If the panel has paper-citable or release-noted features,\n"
        "   add a panel-smoke scenario to tests/panels/smoke_scenarios.py\n"
        "   then run `python .claude/skills/panel-smoke/smoke.py --panel "
        f"{panel_path.split('/')[-1].replace('.py', '')}`.\n"
        "5. Report BUG / RISK / NIT / OK per panel-review SKILL.md\n"
        "   output format.\n\n"
        "Common gotchas (see panel_review(topic='bug_catalogue')):\n"
        " - val() returns string for QCheckBox (both truthy)\n"
        " - _SOLVER_MAP values must match calc argparse choices\n"
        " - widget defaults must not silently break the calc layer\n"
    )


# ============================================================
# Entry point
# ============================================================

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
