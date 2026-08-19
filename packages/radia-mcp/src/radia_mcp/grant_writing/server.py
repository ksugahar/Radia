"""MCP Server: radia_mcp.grant_writing

Grant proposal writing support: Japanese technical-prose lint, section
coverage checks, KDDI Digital Innovation / social-implementation axes,
KAKENHI AI-era OSS research-platform axes, domain focus checks, budget
alignment checks, internal-evidence-to-external-scale checks, and integrated
health reports. Collaborative-integration checks cover lifecycle cost,
ecosystem boundaries, scope, negative results, ethics, and asset provenance.
Tool-to-domain checks require infrastructure proposals to end in a measurable
field outcome and falsifiable knowledge product. Derived-metric checks require
separate calibration and held-out validation. Cross-organization pilot checks
distinguish scientific re-execution from internal use, links, and user counts.
Literature-gap checks prevent non-detection in a bounded corpus from being
promoted into field-wide adoption claims or the proposal's academic gap.
Abstraction checks keep named software out of the research concept while
retaining names where implementation or feasibility must be reproducible.
Reviewer-vocabulary checks explain OSS/AI terms, prefer readable field terms,
and keep named benchmarks in a verification rather than significance role.
Persuasion-quality checks catch self-negating evidence, abrupt equations,
undefined symbols, defensive paragraphs, optional branches, acronym piles, and
internal memo shorthand. Review-format checks encode the in-house KAKENHI
briefing realities: three review criteria, ~100 proposals per reviewer-month,
monochrome printing, researchmap-referenced records, the human-rights/legal
box, and the funding-overlap box format.

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
        "axes, KAKENHI AI-era OSS platform axes, domain focus, budget "
        "alignment, internal-to-external scale, collaborative-integration "
        "risks, tool-to-domain outcomes, derived-metric validation, "
        "cross-organization pilots, literature-gap evidence scope, "
        "named-software abstraction, "
        "reviewer vocabulary and benchmark role, "
        "persuasion hierarchy and equation introductions, "
        "KAKENHI review-format realities, "
        "and integrated health reports."
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
            "見積根拠を積算して検証ループに対応づける。公式料金表、料金年度、"
            "参照日、税込区分、最低購入単位、有効期限、為替、端数処理を記録する。"
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
