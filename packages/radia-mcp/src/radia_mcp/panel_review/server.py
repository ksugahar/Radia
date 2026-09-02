"""
Panel Review MCP Server (radia_mcp.panel_review)

Surfaces the current Radia Simulink application-block review contract. The old
PySide6 and notebook-workbench review chains are retired.

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
    Get Radia Simulink application-block review / construction documentation.

    The topics cover the Simulink block contract, DesignSpec/headless wiring,
    masks, typed ports, result artifacts, validation, and the Cubit boundary.

    Args:
        topic: Documentation topic.  Options:
            "overview"                - application-block review contract
            "build_application_block" - block construction workflow
            "cubit_boundary"          - CAD/mesh/Simulink ownership boundary
            "workflow"                - block construction workflow
            "all"                     - all topics concatenated
    """
    return get_panel_review_documentation(topic)


# ============================================================
# MCP Prompts
# ============================================================

@mcp.prompt()
def review_a_panel(panel_path: str = "matlab/radia_simulink_library.slx") -> str:
    """Run a thorough review of a Radia application block."""
    return (
        f"Please run a thorough Simulink application-block review for: {panel_path}\n\n"
        "Steps (use the panel_review MCP tool for reference text):\n"
        "1. panel_review(topic='overview') for the block contract.\n"
        "2. Inspect DesignSpec, calc argparse, mask, ports, and artifact runner.\n"
        "3. Run Python and MATLAB Simulink workflow tests.\n"
        "4. Confirm the block matches the headless golden and failure provenance.\n"
        "5. Report BUG / RISK / NIT / OK with file/line references.\n\n"
        "Common gotchas:\n"
        " - solver logic duplicated in a mask callback\n"
        " - Python launched on every time step\n"
        " - MEX promoted without parity/lifecycle/long-run evidence\n"
    )


@mcp.prompt()
def build_application_block(
    app_name: str = "ih",
    source_material: str = "validation_test/induction_heating",
) -> str:
    """Build or upgrade a Radia Simulink application block."""
    return (
        f"Please build or upgrade the Radia Simulink application block for app={app_name!r} "
        f"from source material under {source_material}.\n\n"
        "Use `panel_review(topic='build_application_block')` and "
        "`panel_review(topic='cubit_boundary')` first.\n\n"
        "Required shape:\n"
        "1. Move reusable computation into `src/` or a headless `calc_*.py`.\n"
        "2. Keep heavy checks in `validation_test/`.\n"
        "3. Create/update `<App>DesignSpec` and a masked block in the single "
        "Radia Simulink library.\n"
        "4. Save `command.txt`, `run.log`, `solver_result.json`, and versioned "
        "`result.json`.\n"
        "5. Run the Python application runner and MATLAB Simulink tests.\n"
        "6. Do not add PySide6/PyQt to normal Radia Python; Cubit's own "
        "toolbar runtime is the only PySide exception. Do not add a new "
        "notebook workbench.\n"
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
    description="Radia Simulink application-block review/construction contract "
                "(DesignSpec, masks, ports, result artifacts, no-PySide gate)",
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
