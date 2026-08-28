"""MCP Server: radia_mcp.grant_writing

Grant proposal writing support: Japanese technical-prose lint, section
coverage checks, KDDI Digital Innovation / social-implementation axes,
official KAKENHI B/C review axes, KAKENHI AI-era OSS research-platform axes,
domain focus checks, budget alignment checks,
internal-evidence-to-external-scale checks, and integrated health reports.
Collaborative-integration checks cover lifecycle cost,
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
internal memo shorthand. Adjacent-reviewer readability checks catch short but
concept-dense sentences, method-name piles, and paragraphs that mix scientific,
decision, and platform layers; they also expose takeaways delayed until after
method names and numerical evidence, vague decision objects, and required scope
with no deliverable. Reviewer-vocabulary checks distinguish an
MCP access interface from the repository or database that stores its sources.
Reviewer-momentum checks distinguish empty excitement language from a useful
opening arc: reviewer-visible stakes, a concrete bottleneck or unused
opportunity, the research move, and an observable payoff. They also catch
method inventories that appear before the problem is understood.
The official KAKENHI reference map separates
the three research-plan elements, internationality, and budget validity.
Review-format checks separately encode in-house briefing realities:
~100 proposals per reviewer-month, monochrome printing,
researchmap-referenced records, the human-rights/legal box, and the
funding-overlap box format. They also enforce the form-specific rule that a
non-applicant keeps the final-year-early-application page completely blank
instead of writing "not applicable." Central-claim checks catch one
question restated as two: keyword coverage cannot see it, because every
required word is present and the defect is that the words disagree.
Argument-evidence mapping indexes the question, gap, operations, decision
rules, knowledge output, preliminary evidence, responsibilities, and negative
results for an LLM/human close read without scoring their scientific validity.
It also locates whether preparation evidence is tied to a research item the
team can start or execute.
Checks that a section cannot answer -- budget itemization in a research plan
-- report themselves inapplicable instead of scoring it low.
Japanese readability is scored on a separate 100-point diagnostic using
Japanese sentence structure, logical order, subject-predicate proximity,
concept load, notation consistency, and committed proposal wording. English
is outside that score and is never averaged into it.

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
from .tools import (
    grant_writing_health_report,
    grant_writing_japanese_readability_score,
)

mcp = FastMCP("mcp-server-grant-writing")

_n_tools = register(mcp)

register_status_tool(
    mcp,
    server_name="mcp-server-grant-writing",
    description=(
        "Grant proposal lint and review helpers: Japanese technical prose, "
        "section coverage, KDDI Digital Innovation social-implementation "
        "axes, official KAKENHI B/C axes, KAKENHI AI-era OSS platform axes, "
        "domain focus, budget "
        "alignment, internal-to-external scale, collaborative-integration "
        "risks, tool-to-domain outcomes, derived-metric validation, "
        "cross-organization pilots, literature-gap evidence scope, "
        "named-software abstraction, "
        "reviewer vocabulary and benchmark role, "
        "persuasion hierarchy and equation introductions, "
        "adjacent-domain reviewer readability and concept density, "
        "Japanese-only 100-point readability scoring, "
        "reviewer momentum from concrete tension to observable payoff, "
        "MCP role accuracy and preparation-to-plan traceability, "
        "KAKENHI official review structure and review-format realities, "
        "form-specific blank-field rules for conditional KAKENHI sections, "
        "central-claim consistency across summary and body, "
        "non-scoring argument-evidence maps for close reading, "
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
        assert report["defect_counts"]["total"] >= 0
        readable = grant_writing_japanese_readability_score(
            "設計条件の選択には時間を要する。"
            "本研究では候補順位が一致する条件を明らかにする。"
            "二つの解析法を比較し、適用範囲を判定する。"
        )
        assert readable["status"] == "pass"
        assert readable["score"] >= 85
        english = grant_writing_japanese_readability_score(
            "This English proposal is outside the Japanese score."
        )
        assert english["status"] == "not_applicable"
        print(
            "mcp-server-grant-writing self-test: "
            f"registered {_n_tools} domain tools (+ status tool)"
        )
        print(f"  sample located defects: {report['defect_counts']['total']}"
              f" (defect_score {report['defect_score']}/10)")
        print(
            "  japanese readability: "
            f"{readable['score']}/100; English excluded as designed"
        )
        return
    mcp.run()


if __name__ == "__main__":
    main()
