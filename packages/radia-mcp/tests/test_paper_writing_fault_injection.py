"""Detector failures must not become clean composite diagnostics."""
import importlib

import pytest

from radia_mcp.paper_writing.plans.T9 import paper_writing_reviewer_2_trigger_summary
from radia_mcp.paper_writing.plans.T18 import paper_writing_run_full_workflow
from radia_mcp.paper_writing._em_paper_style import paper_writing_em_submission_gate


DETECTORS = [
    ("plans.T3", "paper_writing_claim_quantification", "unquantified_hype",
     {"observations": {"unquantified_count": 0}}),
    ("plans.T5", "paper_writing_related_work_density", "self_cite_excess",
     {"observations": {"self_cite_ratio": 0}}),
    ("plans.T4", "paper_writing_limitation_statement_presence", "no_limitation",
     {"observations": {"limitation_paragraphs": ["limited"], "discussion_found": True}}),
    ("plans.T6", "paper_writing_figure_referencing_coverage", "unreferenced_figures",
     {"observations": {"unreferenced_labels": []}}),
    ("cross_lint", "paper_writing_find_undefined_acronyms", "undefined_acronyms",
     {"undefined_count": 0}),
    ("tools", "paper_writing_count_weak_expressions", "weak_expressions",
     {"total_weak_expressions": 0}),
    ("tools", "paper_writing_check_english_redflags", "english_redflags",
     {"total_issues": 0}),
]


def patch_detector(monkeypatch, detector, payload, raises=False):
    module, name, _, _ = detector
    def run(*args, **kwargs):
        if raises:
            raise RuntimeError("injected detector crash")
        return payload
    monkeypatch.setattr(importlib.import_module("radia_mcp.paper_writing." + module), name, run)


@pytest.fixture
def clean_detectors(monkeypatch):
    for detector in DETECTORS:
        patch_detector(monkeypatch, detector, detector[3])


@pytest.mark.parametrize("detector", DETECTORS, ids=[d[2] for d in DETECTORS])
@pytest.mark.parametrize("mode", ["exception", "error", "ok_false", "skip", "empty", "none"])
def test_each_failed_detector_prevents_clean_score(monkeypatch, clean_detectors, detector, mode):
    payload = {"error": {"error": "injected"}, "ok_false": {"ok": False},
               "skip": {"status": "skip"}, "empty": {}, "none": None}.get(mode)
    patch_detector(monkeypatch, detector, payload, raises=mode == "exception")
    result = paper_writing_reviewer_2_trigger_summary("text", bib="bib")
    assert result["status"] == "partial"
    assert result["score"] is None
    assert result["observations"]["risk_percent"] is None
    assert result["observations"]["unknown_triggers"] == [detector[2]]
    assert "すべて clear" not in " ".join(result["comments"])


def test_all_detectors_failed_is_unavailable(monkeypatch):
    for detector in DETECTORS:
        patch_detector(monkeypatch, detector, {"error": "injected"})
    result = paper_writing_reviewer_2_trigger_summary("text", bib="bib")
    assert result["status"] == "unavailable"
    assert result["score"] is None
    assert len(result["observations"]["unknown_triggers"]) == 7


def test_missing_bib_is_partial_but_findings_survive(monkeypatch, clean_detectors):
    patch_detector(monkeypatch, DETECTORS[0], {"observations": {"unquantified_count": 3}})
    result = paper_writing_reviewer_2_trigger_summary("text")
    assert result["status"] == "partial"
    assert result["score"] is None
    assert result["observations"]["top_triggers"] == ["unquantified_hype"]
    assert len(result["comments"]) == 2


def test_complete_clean_detectors_can_score_ten(clean_detectors):
    result = paper_writing_reviewer_2_trigger_summary("text", bib="bib")
    assert result["status"] == "complete"
    assert result["score"] == 10
    assert result["observations"]["risk_percent"] == 0


@pytest.mark.parametrize("value", [-1, float("nan"), float("inf"), True, "0"])
def test_invalid_counts_are_unknown(monkeypatch, clean_detectors, value):
    patch_detector(monkeypatch, DETECTORS[0], {"observations": {"unquantified_count": value}})
    result = paper_writing_reviewer_2_trigger_summary("text", bib="bib")
    assert result["score"] is None
    assert result["observations"]["unknown_triggers"] == ["unquantified_hype"]


def test_workflow_keeps_active_triggers_and_forwards_bib(monkeypatch, clean_detectors):
    def related(text, bib, author_last_names):
        assert bib == "supplied bib"
        assert author_last_names == "Smith"
        return {"observations": {"self_cite_ratio": 0.8}}
    monkeypatch.setattr(importlib.import_module("radia_mcp.paper_writing.plans.T5"),
                        DETECTORS[1][1], related)
    result = paper_writing_run_full_workflow("text", bib="supplied bib", author_last_names="Smith",
                                           skip_phases="phase1,phase2,phase3,phase5,phase6")
    phase4 = result["phases"]["phase4_reviewer_2_simulate"]
    assert phase4["total_triggers"] == 1
    assert phase4["top_tier_count"] == 1
    assert phase4["risk_percent"] > 0
    assert result["status"] == "complete"


def test_workflow_reports_incomplete_review(clean_detectors):
    result = paper_writing_run_full_workflow("text", skip_phases="phase1,phase2,phase3,phase5,phase6")
    assert result["status"] == "partial"
    assert result["execution_summary"]["n_errors"] == 1
    assert "判定は保留" in result["overall_summary"]
    assert "致命傷無し" not in result["overall_summary"]


@pytest.mark.parametrize("payload", [None, {}, {"error": "injected"}])
def test_workflow_invalid_reviewer_result_is_not_success(monkeypatch, payload):
    monkeypatch.setattr(importlib.import_module("radia_mcp.paper_writing.plans.T9"),
                        "paper_writing_reviewer_2_trigger_summary", lambda *a, **kw: payload)
    result = paper_writing_run_full_workflow("text", skip_phases="phase1,phase2,phase3,phase5,phase6")
    assert result["status"] == "partial"
    assert result["execution_summary"]["n_errors"] == 1
    assert "Workflow 完了" not in result["overall_summary"]


@pytest.mark.parametrize("mode,expected", [("exception", "fail"), ("error", "fail"),
    ("ok_false", "fail"), ("skip", "skip"), ("not_applicable", "skip")])
@pytest.mark.parametrize("tool,check", [
    ("paper_writing_check_equation_numbering", "equation_numbering"),
    ("paper_writing_count_underlines", "count_underlines"),
    ("paper_writing_check_figure_forward_reference", "figure_forward_reference"),
])
def test_gate_does_not_promote_detector_failure(monkeypatch, tmp_path, mode, expected, tool, check):
    from radia_mcp.paper_writing import tools
    def failed(*args, **kwargs):
        if mode == "exception":
            raise RuntimeError("injected")
        return {"error": {"error": "injected"}, "ok_false": {"ok": False},
                "skip": {"status": "skip"}, "not_applicable": {"applicable": False}}[mode]
    monkeypatch.setattr(tools, tool, failed)
    tex = tmp_path / "paper.tex"
    bib = tmp_path / "references.bib"
    tex.write_text("Text.", encoding="utf-8")
    bib.write_text("", encoding="utf-8")
    result = paper_writing_em_submission_gate(tex_path=str(tex), bib_path=str(bib),
                                             auto_resolve_inputs=False)
    row = next(row for row in result["checks"] if row["name"] == check)
    assert row["status"] == expected
    assert result["verdict"] != "pass"
    assert "SUBMISSION-READY" not in result["advice"]
