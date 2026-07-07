"""MCP Server: radia_mcp.grant_writing

Grant proposal writing support: Japanese technical-prose lint, section
coverage checks, KDDI Digital Innovation / social-implementation axes,
domain focus checks, budget alignment checks, and integrated health reports.

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
        "axes, domain focus, budget alignment, and integrated health reports."
    ),
    subpackage="radia_mcp.grant_writing",
    related_servers=["paper-writing", "figure", "presentation", "document-meta"],
    optional_deps=[],
)


def main():
    """Entry point for mcp-server-grant-writing."""
    if "--selftest" in sys.argv:
        sample = (
            "本研究の主題はパワーエレクトロニクス基板CAE-AI環境の社会実装である。"
            "社会的課題は地域製造業のパワーエレクトロニクス設計である。"
            "本研究は生成AIとMCPによりLTspice、SPICE、Radia、NGSolve、PEEC、"
            "熱解析を接続する。商用CAEを直ちに置換するものではなく、"
            "Python-nativeで現代的なAPIを持つ入口としてPoC基板で実証する。"
            "成果はOSSレポジトリ、技術プレゼン、"
            "予算内訳、年度スケジュールとして公開する。"
            "予算は助成上限額に近い申請額とし、単価、数量、月数、年度配分、"
            "見積根拠を積算して検証ループに対応づける。"
        )
        report = grant_writing_health_report(sample, program="kddi_digital")
        assert report["overall_score"] >= 1.0
        print(
            "mcp-server-grant-writing self-test: "
            f"registered {_n_tools} domain tools (+ status tool)"
        )
        print(f"  sample health score: {report['overall_score']}/10")
        return
    mcp.run()


if __name__ == "__main__":
    main()
