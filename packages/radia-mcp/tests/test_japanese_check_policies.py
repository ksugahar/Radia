"""Shared writing algorithms must retain caller-specific rules and adapters."""

import inspect

import pytest
from radia_mcp._shared import japanese_checks
from radia_mcp.paper_writing import _ja_lint as paper
from radia_mcp.presentation import _text_checks as slides
from radia_mcp.presentation import cross_lint


@pytest.mark.parametrize("module", [paper, slides])
@pytest.mark.parametrize("name", [
    "grant_writing_lint_bedrock",
    "grant_writing_suggest_redundancy_fixes",
    "grant_writing_check_misuse_japanese",
])
def test_compatibility_signatures_remain_text_only(module, name):
    function = getattr(module, name)
    assert list(inspect.signature(function).parameters) == ["text"]
    assert function.__doc__


def test_misuse_rules_remain_distinct_between_papers_and_slides():
    text = "新規性を解析手法に置く。"
    result = paper.grant_writing_check_misuse_japanese(text)
    assert result["total_matches"] == 1
    assert result["issues"][0]["match"] == "新規性を解析手法に置く"
    assert cross_lint.presentation_check_misuse_japanese(text)["total_matches"] == 0


@pytest.mark.parametrize("module", [paper, slides])
def test_misuse_wrapper_uses_its_current_policy(module, monkeypatch):
    monkeypatch.setattr(module, "_MISUSE", [("only-here", "local policy")])
    result = module.grant_writing_check_misuse_japanese("only-here")
    assert result["issues"] == [{
        "pattern": "only-here", "match": "only-here", "position": 0,
        "fix": "local policy",
    }]


@pytest.mark.parametrize("module", [paper, slides])
def test_redundancy_wrapper_uses_its_current_policy(module, monkeypatch):
    monkeypatch.setattr(module, "_REDUNDANCY", [("xx", "")])
    result = module.grant_writing_suggest_redundancy_fixes("xx xx")
    assert result["total_matches"] == 2
    assert result["unique_patterns_hit"] == 1
    assert result["top_fixes"][0]["action"] == "delete"
    assert result["top_fixes"][0]["suggested"] is None


@pytest.mark.parametrize("module", [paper, slides])
def test_bedrock_wrapper_uses_its_current_scanner(module, monkeypatch):
    monkeypatch.setattr(module, "_scan_hedges", lambda text: (7, {"caller rule": 7}))
    result = module.grant_writing_lint_bedrock("Clean input.")
    assert result["issues"][0]["count"] == 7
    assert result["issues"][0]["examples"] == ["caller rule"]


def test_shared_misuse_retains_total_and_truncation_contract():
    result = japanese_checks._check_misuse_japanese("x" * 40, [("x", "hint")])
    assert result["total_matches"] == 40
    assert len(result["issues"]) == 30
    assert result["issues"][-1]["position"] == 29
