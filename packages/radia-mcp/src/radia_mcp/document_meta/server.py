"""MCP Server: radia_mcp.document_meta

Cross-cutting document helpers that do not belong to one document type:
  * document_meta_deadline_countdown -- days-to-deadline + phase/action hint
  * document_meta_diff_versions      -- unified diff of two manuscript versions
  * document_meta_template_loader    -- LaTeX skeletons (IEEE/IEEJ/APS/Beamer/JSPS)
  * document_meta_lint_all           -- run every applicable radia-mcp lint over a file
  * document_meta_write_notebook_result_json -- write versioned JSON sidecar from saved outputs
  * document_meta_write_docs_notebook_result_jsons -- batch sidecars for docs notebooks
  * document_meta_notebook_result_audit -- check ipynb saved outputs + JSON sidecars
  * document_meta_examples_notebook_audit -- examples -> docs or validation_test audit
  * document_meta_panel_layout_audit -- root-level panels migration impact audit

Promoted + REDESIGNED 2026-06-02 into radia-mcp's public PyPI package.
The lint registry now covers radia-mcp's public document subpackages
(grant_writing / presentation / paper_writing / poster).

Usage:
    mcp-server-document-meta              # stdio
    mcp-server-document-meta --selftest   # self-test
"""

import sys

from mcp.server.fastmcp import FastMCP

from ..common import register_status_tool
from . import register

mcp = FastMCP("mcp-server-document-meta")

# Register all document_meta_* tools from tools.py.
_n_tools = register(mcp)

register_status_tool(
    mcp,
    server_name="mcp-server-document-meta",
    description=(
        "Cross-cutting document helpers: deadline countdown, manuscript "
        "version diff, LaTeX template scaffolding, and an all-domain lint "
        "orchestrator (registry over radia-mcp grant_writing / presentation "
        "/ paper_writing / poster lints), plus Radia repo "
        "notebook-result and panel-layout audits."
    ),
    subpackage="radia_mcp.document_meta",
    related_servers=["grant-writing", "presentation", "paper-writing", "poster", "bibliography"],
    optional_deps=["pymupdf"],
)


def main():
    """Entry point for mcp-server-document-meta."""
    if "--selftest" in sys.argv:
        print(f"mcp-server-document-meta self-test: registered {_n_tools} tools")
        return
    mcp.run()


if __name__ == "__main__":
    main()
