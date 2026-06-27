"""MCP Server: radia_mcp.grant_writing

Grant proposal writing support: Japanese technical-prose lint, section
coverage checks, KDDI Digital Innovation / social-implementation axes,
budget alignment checks, and integrated health reports.

Promoted to radia-mcp so the document-writing servers are registered in
parallel: paper-writing / figure / grant-writing / presentation.

Usage:
    mcp-server-grant-writing              # stdio
    mcp-server-grant-writing --selftest   # self-test
"""
from __future__ import annotations

import sys

from mcp.server.fastmcp import FastMCP

from ..common import register_status_tool
from . import register
from .tools import grant_writing_health_report

mcp = FastMCP("mcp-server-grant-writing")

_n_tools = register(mcp)

register_status_tool(
    mcp,
    server_name="mcp-server-grant-writing",
    description=(
        "Grant proposal lint and review helpers: Japanese technical prose, "
        "section coverage, KDDI Digital Innovation social-implementation "
        "axes, budget alignment, and integrated health reports."
    ),
    subpackage="radia_mcp.grant_writing",
    related_servers=["paper-writing", "figure", "presentation", "document-meta"],
    optional_deps=[],
)


def main():
    """Entry point for mcp-server-grant-writing."""
    if "--selftest" in sys.argv:
        sample = (
            "社会的課題は地域製造業のパワーエレクトロニクス設計である。"
            "本研究は生成AIとMCPによりLTspice、Radia、NGSolveを接続し、"
            "PoC基板で実証する。成果はOSSレポジトリ、技術プレゼン、"
            "予算内訳、年度スケジュールとして公開する。"
        )
        report = grant_writing_health_report(sample, program="kddi_digital")
        assert report["overall_score"] >= 1.0
        print(f"mcp-server-grant-writing self-test: registered {_n_tools} tools")
        print(f"  sample health score: {report['overall_score']}/10")
        return
    mcp.run()


if __name__ == "__main__":
    main()
