"""Regression contracts for the September draft-check review."""

import asyncio
import math

import pytest
from radia_mcp.grant_writing import register
from radia_mcp.grant_writing import tools as t


@pytest.mark.parametrize("ending", ["。", "．", "?", "？"])
def test_singularity_punctuation_and_short_questions(ending):
    result = t.grant_writing_central_question_singularity_check(
        f"本研究の問いは一つである。目的は何か{ending}さらに何を保つか{ending}"
    )
    assert result["applicable"]
    assert result["score"] is None
    assert result["risks"][0]["type"] == "announced_singular_but_multiple"
    assert result["risks"][0]["question_count"] == 2


def test_enumerations_do_not_make_an_inapplicable_finding():
    result = t.grant_writing_central_question_singularity_check("成果は二つある。計画は三つある。")
    assert not result["applicable"]
    assert not result["risks"]
    assert len(result["enumeration_declarations"]) == 2


def test_separate_announcements_do_not_borrow_each_others_questions():
    result = t.grant_writing_central_question_singularity_check(
        "中心の問いは一つである。何が変わるか？\n" "本研究の問いは一つである。何が変わるか？"
    )
    assert result["announcement_count"] == 2
    assert not result["risks"]


def test_health_report_preserves_unscored_question_candidates():
    text = "本研究の問いは一つである。何が変わるか？さらに何を保つか？"
    result = t.grant_writing_health_report(text)
    assert result["detailed_results"]["central_question_singularity"]["risk_count"]
    assert any(q["tool"] == "singularity" for q in result["questions"])
    assert "central_question_singularity" not in result["detailed_scores"]
    assert (
        "central_question_singularity"
        not in t.grant_writing_health_report(text, skip="singularity")["tools_run"]
    )


def test_length_capacity_is_not_a_quality_score():
    result = t.grant_writing_draft_length_budget_check(
        "本文", page_limit=2, reserved_pages=1, figure_count=2
    )
    assert result["status"] == "no_prose_capacity"
    assert result["score"] is None
    assert result["risks"][0]["severity"] == "HIGH"
    assert result["estimated_allowance"] == 0


@pytest.mark.parametrize(
    "count,status", [(99, "within_estimate"), (100, "within_estimate"), (101, "exceeds_estimate")]
)
def test_length_estimate_boundary(count, status):
    result = t.grant_writing_draft_length_budget_check(
        "字" * count, page_limit=1, chars_per_page=100
    )
    assert result["counted_characters"] == count
    assert result["status"] == status
    assert result["score"] is None


@pytest.mark.parametrize(
    "kwargs",
    [
        {"page_limit": -1},
        {"page_limit": math.nan},
        {"reserved_pages": math.inf},
        {"chars_per_page": 0},
        {"chars_per_page": 1.5},
        {"figure_count": -1},
        {"figure_count": 1.5},
        {"figure_page_cost": -1},
        {"page_limit": True},
    ],
)
def test_invalid_length_inputs_fail_loudly(kwargs):
    with pytest.raises(ValueError):
        t.grant_writing_draft_length_budget_check("本文", **kwargs)


def test_unspecified_page_limit_is_not_pass():
    result = t.grant_writing_draft_length_budget_check("本文")
    assert not result["applicable"]
    assert result["status"] == "not_applicable"
    assert result["score"] is None


@pytest.mark.parametrize(
    "text",
    [
        "人間と機械の調和。社会の価値。目標と従来。独創。マイルストーンと月。共同研究者。研究歴。論文。",
        "導体を回転させ磁束密度を測定する。",
    ],
)
def test_keywords_never_certify_presence_or_absence(text):
    result = t.grant_writing_form_field_coverage_check(text, preset="tateisi_research")
    assert result["score"] is None
    assert result["candidate_count"] + result["unmatched_count"] == 8
    assert "covered_count" not in result
    assert "missing_count" not in result
    assert not result["risks"]
    assert result["status"] == "manual_review_required"


@pytest.mark.parametrize(
    "kwargs",
    [{"preset": "unknown"}, {"fields": ", ,"}, {"preset": "tateisi_research", "fields": "目的"}],
)
def test_invalid_field_configuration_is_visible(kwargs):
    with pytest.raises(ValueError):
        t.grant_writing_form_field_coverage_check("本文", **kwargs)


def test_draft_path_and_field_deduplication(tmp_path):
    path = tmp_path / "draft.md"
    path.write_text("目標を定める。", encoding="utf-8")
    result = t.grant_writing_form_field_coverage_check(str(path), fields="目標,目標,根拠")
    assert result["field_count"] == 2
    assert result["candidate_count"] == 1
    assert result["unmatched_count"] == 1


def test_new_tools_are_registered_with_real_mcp_schemas():
    from mcp.server.fastmcp import FastMCP

    server = FastMCP("draft-regression")
    register(server)
    listed = {tool.name: tool for tool in asyncio.run(server.list_tools())}
    for suffix in (
        "central_question_singularity_check",
        "draft_length_budget_check",
        "form_field_coverage_check",
    ):
        assert "text" in listed["grant_writing_" + suffix].inputSchema["properties"]


@pytest.mark.parametrize("ending", ["。", "?", "？"])
def test_claim_heading_does_not_swallow_following_statement(ending):
    statements = t._claim_statements(
        "中心の問い\n本研究の問いは明快である。" f"磁場測定の境界条件を検証できるか{ending}"
    )
    assert len(statements) == 1
    assert "中心の問い\n" not in statements[0]["text"]
    assert "境界条件" in statements[0]["text"]


def test_kanji_advice_measures_cause_without_renaming_domain_terms():
    from radia_mcp.paper_writing._ja_lint import grant_writing_check_kanji_ratio

    result = grant_writing_check_kanji_ratio("設定電流と動的口径。" * 6)
    assert result["cause"]["formal_noun_in_kanji_count"] == 0
    assert result["cause"]["repeated_domain_terms"]
    assert result["status"] == "kanji_heavy"


def test_wrapped_claim_is_not_a_heading():
    statements = t._claim_statements(
        "本研究の問いは低周波磁場において\n磁場測定の境界条件を検証できるか？"
    )
    assert len(statements) == 1
    assert "低周波磁場" in statements[0]["text"]


def test_formal_nouns_and_compounds_are_distinguished():
    from radia_mcp.paper_writing._ja_lint import grant_writing_check_kanji_ratio

    result = grant_writing_check_kanji_ratio("行う事は測る時に決める。時間と物理。" * 6)
    assert result["cause"]["formal_noun_in_kanji_count"] == 12
