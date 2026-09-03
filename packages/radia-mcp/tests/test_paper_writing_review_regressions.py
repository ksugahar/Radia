from __future__ import annotations

import inspect

import pytest

from radia_mcp.paper_writing import tools as pw
from radia_mcp.paper_writing._citation_verify import (
    _parse_bib_lightweight,
    _title_already_cited,
    paper_writing_verify_citation,
)
from radia_mcp.paper_writing._em_paper_style import (
    paper_writing_em_submission_gate,
)
from radia_mcp.paper_writing._pdf_layout_visual import (
    paper_writing_check_floats_far_from_reference,
    paper_writing_detect_page_whitespace_anomalies,
)
from radia_mcp.paper_writing._pdf_overlap_detection import (
    paper_writing_detect_text_overflow_page,
)
from radia_mcp.paper_writing._tex_resolver import resolve_input_chain


def _checks_by_name(result: dict) -> dict[str, dict]:
    return {check["name"]: check for check in result["checks"]}


def test_submission_gate_rejects_supplied_missing_files(tmp_path):
    result = paper_writing_em_submission_gate(
        tex_path=str(tmp_path / "missing.tex"),
        bib_path=str(tmp_path / "missing.bib"),
        pdf_path=str(tmp_path / "missing.pdf"),
        page_limit=2,
    )
    checks = _checks_by_name(result)
    assert result["verdict"] == "fail"
    assert checks["tex_input"]["status"] == "fail"
    assert checks["bib_input"]["status"] == "fail"
    assert checks["pdf_input"]["status"] == "fail"
    assert not any("0 pages" in check["summary"] for check in result["checks"])


def test_submission_gate_maps_detector_failures_to_status(tmp_path):
    tex = tmp_path / "paper.tex"
    bib = tmp_path / "references.bib"
    abstract = " ".join(["word"] * 251)
    tex.write_text(
        "\\documentclass{article}\n"
        "\\begin{document}\n"
        f"\\begin{{abstract}}{abstract}\\end{{abstract}}\n"
        "\\section{Introduction}\n"
        "\\underline{bad} See Fig.~\\ref{fig:missing}, "
        "Eq.~\\eqref{eq:missing}, and \\cite{known,missing}.\n"
        "\\end{document}\n",
        encoding="utf-8",
    )
    bib.write_text("@article{known, title={Nested {MRI} title}}\n", encoding="utf-8")
    result = paper_writing_em_submission_gate(
        tex_path=str(tex),
        bib_path=str(bib),
        abstract_text=abstract,
        target_venue="IEEE Transactions on Magnetics",
    )
    checks = _checks_by_name(result)
    assert result["verdict"] == "fail"
    assert checks["figure_forward_reference"]["status"] == "fail"
    assert checks["equation_numbering"]["status"] == "fail"
    assert checks["count_underlines"]["status"] == "warn"
    assert checks["lint_reference_format"]["status"] == "fail"
    assert checks["check_citation_usage"]["status"] == "fail"
    assert checks["validate_abstract_length"]["status"] == "fail"


def test_display_math_scanner_does_not_treat_line_spacing_as_math():
    source = (
        r"\begin{equation} a=b\\[2pt] c=d \label{eq:one}\end{equation} "
        r"copper cylinder is discretised. \[x=y\]"
    )
    plain = pw._paper_writing_plain_text(source)
    assert "copper cylinder is discretised" in plain


def test_japanese_abstract_sentences_split_without_spaces():
    result = pw.paper_writing_check_abstract_background_ratio(
        "背景を述べる。課題を述べる。本研究では3例を解析した。"
    )
    assert result["total_sentences"] == 3


def test_conclusion_plain_text_preserves_percent_suffix():
    source = (
        "Body introduces MMPM and 12 GHz.\n\n"
        "## Conclusions\nThe error was 5% lower with MMPM at 12 GHz."
    )
    result = pw.paper_writing_check_conclusion_first_use(source)
    assert result["conclusion_found"] is True
    assert "5%" not in result["new_numeric_claims"]


def test_cref_comma_list_is_split_into_individual_keys(tmp_path):
    tex = tmp_path / "refs.tex"
    tex.write_text(
        r"\label{fig:a}\label{fig:b}\cref{fig:a,fig:b}",
        encoding="utf-8",
    )
    figure = pw.paper_writing_check_figure_forward_reference(str(tex))
    assert figure["orphan_label_count"] == 0
    assert figure["dangling_ref_count"] == 0


def test_balanced_bib_parser_and_title_match_are_exact(tmp_path):
    bib = tmp_path / "references.bib"
    bib.write_text(
        "@string{ieee = {IEEE}}\n"
        "@article{mri, title={A {MRI} Method for Nested {Fields}}, "
        "author={A, B}, year={2025}}\n",
        encoding="utf-8",
    )
    entries = _parse_bib_lightweight(str(bib))
    assert set(entries) == {"mri"}
    assert entries["mri"]["title"] == "A {MRI} Method for Nested {Fields}"
    assert _title_already_cited(entries, "Introduction") is None
    assert _title_already_cited(entries, "A MRI Method for Nested Fields") == "mri"


def test_input_chain_resolves_from_main_compile_directory(tmp_path):
    section_dir = tmp_path / "sec"
    section_dir.mkdir()
    main = tmp_path / "main.tex"
    first = section_dir / "a.tex"
    second = section_dir / "b.tex"
    main.write_text(r"\input{sec/a}", encoding="utf-8")
    first.write_text(r"A \input{sec/b}", encoding="utf-8")
    second.write_text(r"B \cite{nested}", encoding="utf-8")
    result = resolve_input_chain(str(main))
    assert result["ok"] is True
    assert not result["files_missing"]
    assert r"\cite{nested}" in result["merged_tex"]


def test_computational_electromagnetics_venues_are_supported():
    for venue in ("IGTE", "COMPUMAG", "CEFC"):
        result = pw.paper_writing_target_venue_policy(venue)
        assert result["status"] == "pass"
        assert result["target_category"] == "electromagnetics"


def test_float_distance_without_aux_is_not_reported_as_pass(tmp_path):
    tex = tmp_path / "paper.tex"
    pdf = tmp_path / "paper.pdf"
    tex.write_text(
        r"\begin{figure}\label{fig:a}\end{figure}", encoding="utf-8"
    )
    pdf.write_bytes(b"not opened because aux is absent")
    result = paper_writing_check_floats_far_from_reference(str(tex), str(pdf))
    assert result["applicable"] is False
    assert result["status"] == "not_applicable"


def test_float_distance_uses_aux_page_and_pdf_prose_reference(tmp_path):
    fitz = pytest.importorskip("fitz")
    tex = tmp_path / "paper.tex"
    aux = tmp_path / "paper.aux"
    pdf = tmp_path / "paper.pdf"
    tex.write_text(
        r"See Fig.~\ref{fig:a}.\begin{figure}\caption{Result}\label{fig:a}"
        r"\end{figure}",
        encoding="utf-8",
    )
    aux.write_text(
        r"\newlabel{fig:a}{{1}{3}{Result}{figure.1}{}}",
        encoding="utf-8",
    )
    with fitz.open() as document:
        document.new_page().insert_text((72, 72), "As shown in Fig. 1, the result agrees.")
        document.new_page().insert_text((72, 72), "Intermediate text")
        document.new_page().insert_text((72, 72), "Fig. 1. Result")
        document.save(pdf)

    result = paper_writing_check_floats_far_from_reference(
        str(tex), str(pdf), max_pages_apart=1, aux_path=str(aux)
    )

    assert result["applicable"] is True
    assert result["flagged_count"] == 1
    assert result["flagged"][0]["float_page"] == 3
    assert result["flagged"][0]["first_ref_page"] == 1


def test_cropbox_with_nonzero_pdf_origin_does_not_create_false_overflow(tmp_path):
    fitz = pytest.importorskip("fitz")
    pdf = tmp_path / "cropped.pdf"
    with fitz.open() as document:
        page = document.new_page(width=600, height=800)
        page.set_cropbox(fitz.Rect(40, 50, 560, 750))
        page.insert_text((20, 30), "Text inside the visible crop area")
        document.save(pdf)

    result = paper_writing_detect_text_overflow_page(str(pdf), use_cropbox=True)

    assert result["n_overflows"] == 0


def test_temporary_doi_failure_is_not_called_fabrication(tmp_path, monkeypatch):
    bib = tmp_path / "references.bib"
    bib.write_text("@article{known, title={Known}, year={2025}}", encoding="utf-8")
    monkeypatch.setattr(
        "radia_mcp.paper_writing.paper_download.paper_writing_resolve_doi",
        lambda _doi: {
            "ok": False,
            "error": "crossref HTTP 429",
            "temporary_failure": True,
        },
    )

    result = paper_writing_verify_citation(
        "A claim", str(bib), candidate_doi="10.1234/example"
    )

    assert result["verdict"] == "error"
    assert "temporarily unavailable" in result["advice"]
    assert "Do not infer" in result["advice"]
    assert "likely fabricated" not in result["advice"]


@pytest.mark.parametrize(
    ("status_code", "temporary"),
    [(403, True), (429, True), (503, True), (404, False)],
)
def test_crossref_http_status_classification(monkeypatch, status_code, temporary):
    from radia_mcp.paper_writing import paper_download

    class Response:
        pass

    response = Response()
    response.status_code = status_code

    class Requests:
        @staticmethod
        def get(*_args, **_kwargs):
            return response

    monkeypatch.setattr(paper_download, "_require_requests", lambda: Requests)
    result = paper_download.paper_writing_resolve_doi("DOI:10.1234/test")

    assert result["temporary_failure"] is temporary
    assert result["error_kind"] == ("temporary" if temporary else "not_found")


def test_whitespace_detector_default_does_not_use_old_75_percent_cutoff():
    default = inspect.signature(
        paper_writing_detect_page_whitespace_anomalies
    ).parameters["whitespace_threshold"].default
    assert default == 0.95


def test_triangle_uses_lowercase_english_content_words():
    text = r"""
\title{Autonomous Electromagnetic Motion Analysis}
\begin{abstract}Autonomous electromagnetic motion analysis is demonstrated.\end{abstract}
\section{Conclusions}The autonomous electromagnetic motion analysis was validated.
"""
    result = pw.paper_writing_title_abstract_conclusion_triangle(text)
    assert result["score"] >= 8.0


def test_citation_health_parses_single_line_nested_entry():
    bib = "@article{x, title={Nested {MRI}}, author={A, B}, year={2025}}"
    result = pw.paper_writing_citation_health_4_axes(bib, current_year=2026)
    assert result["score"] is not None


def test_reproducibility_does_not_read_coil_or_table_as_metadata():
    result = pw.paper_writing_reproducibility_open_science_check(
        "The coil result is listed in Table 2.1."
    )
    assert result["per_axis"]["competing_interests"]["covered"] is False
    assert result["per_axis"]["methods_replicability"]["covered"] is False


def test_physical_unit_is_not_misread_as_p_value():
    result = pw.paper_writing_statistical_reporting_compliance(
        "The loss was p = 0.05 W."
    )
    assert result["total_p_values"] == 0


def test_terminology_normalizer_preserves_line_endings(tmp_path):
    path = tmp_path / "paper.tex"
    path.write_bytes("立方体ではなくcubeを用いる。\nSecond line.\n".encode("utf-8"))
    result = pw.paper_writing_normalize_terminology_file(
        str(path), dry_run=False
    )
    assert result["changed"] is True
    assert b"\r\n" not in path.read_bytes()


def test_terminology_normalizer_never_writes_replacement_char(tmp_path):
    path = tmp_path / "paper.tex"
    path.write_bytes(b"\x81")
    result = pw.paper_writing_normalize_terminology_file(
        str(path), dry_run=False
    )
    assert "error" in result
    assert path.read_bytes() == b"\x81"
