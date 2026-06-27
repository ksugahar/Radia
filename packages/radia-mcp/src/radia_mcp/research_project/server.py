"""MCP Server: radia_mcp.research_project

Research-project dashboard / scanner that runs health reports across all
document types in a project tree:
  * research_project_consistency_check
  * research_project_deadline_gantt
  * research_project_health_dashboard
  * research_project_scan

Promoted 2026-06-02 into radia-mcp's public PyPI package.  Its per-domain
health handlers import the radia-mcp document subpackages (grant_writing /
paper_writing / poster / doc_convert / bibliography / pdf).

Usage:
    mcp-server-research-project              # stdio
    mcp-server-research-project --selftest   # self-test
"""

import sys

from mcp.server.fastmcp import FastMCP

from ..common import register_status_tool
from . import register

mcp = FastMCP("mcp-server-research-project")

# Register all research_project_* tools from tools.py.
_n_tools = register(mcp)

register_status_tool(
    mcp,
    server_name="mcp-server-research-project",
    description=(
        "Research-project dashboard: consistency check, deadline gantt, and a "
        "health-dashboard / scan that aggregates per-document-type health "
        "reports (grant / paper / poster / pptx / bibliography / pdf via "
        "radia-mcp)."
    ),
    subpackage="radia_mcp.research_project",
    related_servers=["document-meta", "grant-writing", "paper-writing", "poster", "bibliography"],
    optional_deps=["pymupdf", "python-pptx"],
)


def main():
    """Entry point for mcp-server-research-project."""
    if "--selftest" in sys.argv:
        print(f"mcp-server-research-project self-test: registered {_n_tools} tools")
        return
    mcp.run()


if __name__ == "__main__":
    main()
