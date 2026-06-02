"""MCP Server: radia_mcp.research_project

Research-project dashboard / scanner that runs health reports across all
document types in a project tree:
  * research_project_consistency_check
  * research_project_deadline_gantt
  * research_project_health_dashboard
  * research_project_scan

Promoted 2026-06-02 from mcp-server-document.research_project (LAB-private)
to radia-mcp (public PyPI).  Its per-domain health handlers import the
radia-mcp document subpackages (paper_writing / poster / doc_convert /
bibliography / pdf).  The grant_writing health handler is graceful: it uses
the LAB-private mcp-server-document.grant_writing opportunistically if that
package is installed, otherwise reports grant health as unavailable -- so
radia-mcp has no hard dependency on the private package.

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
        "reports (paper / poster / pptx / bibliography / pdf via radia-mcp; "
        "grant health is opportunistic via the LAB-private mcp-server-document)."
    ),
    subpackage="radia_mcp.research_project",
    related_servers=["document-meta", "paper-writing", "poster", "bibliography"],
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
