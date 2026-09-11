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


# module, callable, gate row, minimal successful payload
GATE_DETECTORS = [
    ("tools", "paper_writing_target_venue_policy", "target_venue_policy",
     {"status": "pass", "target_category": "ieee"}),
    ("tools", "paper_writing_check_equation_numbering", "equation_numbering",
     {"dangling_refs": [], "unused_labels": []}),
    ("tools", "paper_writing_count_underlines", "count_underlines", {"total_underlines": 0}),
    ("tools", "paper_writing_check_figure_forward_reference", "figure_forward_reference",
     {"dangling_refs": [], "orphan_labels": []}),
    ("_undefined_variables", "paper_writing_check_undefined_variables", "undefined_variables", {"n_undefined": 0}),
    ("_undefined_acronyms", "paper_writing_check_undefined_acronyms", "undefined_acronyms", {"n_undefined": 0}),
    ("_digest_lints", "paper_writing_check_ref_label_consistency", "ref_label_consistency", {"status": "pass"}),
    ("_digest_lints", "paper_writing_check_ieee_keywords", "ieee_keywords", {"status": "pass"}),
    ("tools", "paper_writing_check_digest_human_review_triggers", "digest_human_review_triggers", {"status": "pass"}),
    ("tools", "paper_writing_bilingual_readability_check", "bilingual_adjacent_reviewer_readability", {"status": "pass"}),
    ("tools", "paper_writing_lint_reference_format", "lint_reference_format", {"problems": []}),
    ("tools", "paper_writing_check_citation_usage", "check_citation_usage",
     {"missing_in_bib_count": 0, "unused_in_tex_count": 0}),
    ("tools", "paper_writing_check_self_citation_ratio", "self_citation_ratio",
     {"self_citation_ratio": 0, "verdict": "OK"}),
    ("_citation_verify", "paper_writing_check_citation_keys_exist", "citation_keys_exist", {"status": "pass"}),
    ("tools", "paper_writing_validate_abstract_length", "validate_abstract_length", {"within_limit": True}),
    ("tools", "paper_writing_check_abstract_background_ratio", "abstract_background_ratio", {"background_ratio": 0}),
    ("tools", "paper_writing_count_weak_expressions", "abstract_weak_expressions", {"total_weak_expressions": 0}),
    ("tools", "paper_writing_check_abstract_no_math_no_citation", "abstract_no_math_no_citation", {"status": "pass"}),
    ("tools", "paper_writing_validate_pdf_pages", "validate_pdf_pages", {"page_count": 1}),
    ("_pdf_layout_visual", "paper_writing_detect_page_whitespace_anomalies", "page_whitespace_anomalies", {"flagged_count": 0}),
    ("_pdf_layout_visual", "paper_writing_check_floats_far_from_reference", "floats_far_from_reference", {"flagged_count": 0}),
    ("_pdf_overlap_detection", "paper_writing_detect_text_image_overlap", "text_image_overlap", {"n_overlaps": 0}),
    ("_pdf_overlap_detection", "paper_writing_detect_text_overflow_page", "text_overflow_page", {"n_overflows": 0}),
    ("_digest_lints", "paper_writing_check_pdf_unresolved_markers", "pdf_unresolved_markers", {"status": "pass"}),
]


@pytest.fixture
def gate_inputs(monkeypatch, tmp_path):
    for detector in GATE_DETECTORS:
        patch_detector(monkeypatch, detector, detector[3])
    tex = tmp_path / "paper.tex"
    bib = tmp_path / "references.bib"
    pdf = tmp_path / "paper.pdf"
    tex.write_text("Text.", encoding="utf-8")
    bib.write_text("", encoding="utf-8")
    # All PDF readers are mocked: this suite tests orchestration, not rendering.
    pdf.write_bytes(b"mock PDF input")
    return dict(tex_path=str(tex), bib_path=str(bib), pdf_path=str(pdf),
                abstract_text="An abstract.", author_last_names="Smith",
                target_venue="IEEE", page_limit=2, auto_resolve_inputs=False)


@pytest.mark.parametrize("mode,expected", [("exception", "fail"), ("error", "fail"),
    ("ok_false", "fail"), ("skip", "skip"), ("not_applicable", "skip"),
    ("empty", "fail"), ("none", "fail"), ("malformed", "fail"),
    ("status_error", "fail")])
@pytest.mark.parametrize("detector", GATE_DETECTORS, ids=[d[2] for d in GATE_DETECTORS])
def test_gate_does_not_promote_detector_failure(monkeypatch, gate_inputs, mode, expected, detector):
    def failed(*args, **kwargs):
        if mode == "exception":
            raise RuntimeError("injected")
        if mode == "malformed":
            return {key: None for key in detector[3]}
        return {"error": {"error": "injected"}, "ok_false": {"ok": False},
                "skip": {"status": "skip"}, "not_applicable": {"applicable": False},
                "empty": {}, "none": None, "status_error": {"status": "error"}}[mode]
    module, tool, check, _ = detector
    monkeypatch.setattr(importlib.import_module("radia_mcp.paper_writing." + module), tool, failed)
    result = paper_writing_em_submission_gate(**gate_inputs)
    row = next(row for row in result["checks"] if row["name"] == check)
    assert row["status"] == expected
    assert result["verdict"] != "pass"
    assert "SUBMISSION-READY" not in result["advice"]
    assert any(c["name"] == "pdf_unresolved_markers" for c in result["checks"])


def test_gate_complete_success_is_still_possible(gate_inputs):
    result = paper_writing_em_submission_gate(**gate_inputs)
    assert result["verdict"] == "pass"
    assert result["n_passed"] == len(GATE_DETECTORS) + 1  # bibliography policy
    assert result["n_skipped"] == 0


@pytest.mark.parametrize("ratio,verdict,expected", [(0.2, "OK", "pass"), (0.21, "WARN", "warn"), (0.3, "HIGH", "warn")])
def test_gate_self_citation_verdict_mapping(monkeypatch, gate_inputs, ratio, verdict, expected):
    detector = next(d for d in GATE_DETECTORS if d[2] == "self_citation_ratio")
    patch_detector(monkeypatch, detector, {"self_citation_ratio": ratio, "verdict": verdict})
    result = paper_writing_em_submission_gate(**gate_inputs)
    assert next(c for c in result["checks"] if c["name"] == detector[2])["status"] == expected
    assert result["verdict"] == expected


@pytest.mark.parametrize("ratio,expected", [(0.25, "pass"), (0.26, "warn"), (0.4, "warn"), (0.41, "fail")])
def test_gate_background_ratio_boundaries(monkeypatch, gate_inputs, ratio, expected):
    detector = next(d for d in GATE_DETECTORS if d[2] == "abstract_background_ratio")
    patch_detector(monkeypatch, detector, {"background_ratio": ratio})
    result = paper_writing_em_submission_gate(**gate_inputs)
    assert next(c for c in result["checks"] if c["name"] == detector[2])["status"] == expected
    assert result["verdict"] == expected


@pytest.mark.parametrize("value", [float("nan"), float("inf"), -1, True, "0"])
@pytest.mark.parametrize("check,key", [
    ("abstract_background_ratio", "background_ratio"),
    ("self_citation_ratio", "self_citation_ratio"),
    ("validate_pdf_pages", "page_count"),
    ("page_whitespace_anomalies", "flagged_count"),
])
def test_gate_invalid_measurements_fail(monkeypatch, gate_inputs, check, key, value):
    detector = next(d for d in GATE_DETECTORS if d[2] == check)
    patch_detector(monkeypatch, detector, {**detector[3], key: value})
    result = paper_writing_em_submission_gate(**gate_inputs)
    assert next(c for c in result["checks"] if c["name"] == check)["status"] == "fail"
    assert result["verdict"] == "fail"


@pytest.mark.parametrize("payload", [None, {}, {"error": "injected"}, {"ok": False},
    {"ok": True, "files_resolved": []},
    {"ok": True, "files_resolved": ["main.tex"], "files_missing": ["child.tex"]}])
def test_gate_resolver_invalid_results_fail(monkeypatch, gate_inputs, payload):
    resolver = importlib.import_module("radia_mcp.paper_writing._tex_resolver")
    monkeypatch.setattr(resolver, "resolve_input_chain", lambda *a, **kw: payload)
    result = paper_writing_em_submission_gate(**{**gate_inputs, "auto_resolve_inputs": True})
    assert next(c for c in result["checks"] if c["name"] == "multifile_resolved")["status"] == "fail"
    assert result["verdict"] == "fail"


@pytest.mark.parametrize("tool,check,options", [
    ("resolve_input_chain", "multifile_resolved", {"auto_resolve_inputs": True}),
    ("extract_abstract_from_tex", "abstract_extracted", {"abstract_text": ""}),
])
def test_gate_preprocessing_exception_is_reported(monkeypatch, gate_inputs, tool, check, options):
    def crashed(*args, **kwargs):
        raise RuntimeError("injected preprocessing failure")
    monkeypatch.setattr(importlib.import_module("radia_mcp.paper_writing._tex_resolver"), tool, crashed)
    result = paper_writing_em_submission_gate(**{**gate_inputs, **options})
    assert next(c for c in result["checks"] if c["name"] == check)["status"] == "fail"
    assert result["verdict"] == "fail"


@pytest.mark.parametrize("payload", [{"error": "injected"}, {}, ["text"], 1])
def test_gate_abstract_extractor_requires_text(monkeypatch, gate_inputs, payload):
    monkeypatch.setattr(importlib.import_module("radia_mcp.paper_writing._tex_resolver"),
                        "extract_abstract_from_tex", lambda *a, **kw: payload)
    result = paper_writing_em_submission_gate(**{**gate_inputs, "abstract_text": ""})
    assert next(c for c in result["checks"] if c["name"] == "abstract_extracted")["status"] == "fail"
    assert result["verdict"] == "fail"
