"""Regression tests for the curated IEEJ technical-report knowledge."""

from radia_mcp.radia_ngsolve.knowledge.technical_reports import (
    get_technical_reports_documentation,
)


def test_report_topics_expose_actionable_guidance():
    open_boundary = get_technical_reports_documentation("open_boundary")
    mor = get_technical_reports_documentation("cln")
    optimization = get_technical_reports_documentation("ai")

    assert "Strategic Dual Image" in open_boundary
    assert "Kelvin" in open_boundary
    assert "IABC" in open_boundary
    assert "expansion points" in mor
    assert "Lanczos" in mor
    assert "current continuity" in mor
    assert "uncertainty" in optimization
    assert "field model" in optimization


def test_all_report_knowledge_contains_sources_and_cross_links():
    text = get_technical_reports_documentation("all")

    assert len(text) > 20_000
    assert "2005-12" in text
    assert "2025-04" in text
    assert "NGSolve" in text
    assert "force_validation" in text
    assert "radia" in text.lower()


def test_aliases_and_unknown_topic_are_stable():
    assert get_technical_reports_documentation("kelvin") == (
        get_technical_reports_documentation("open-boundary")
    )
    unknown = get_technical_reports_documentation("not-a-report-topic")
    assert unknown.startswith("Unknown topic")
    assert "overview" in unknown

