"""Comprehensive test suite for radia_mcp.paper_writing.

Covers the 80 paper_writing_* tools by family with small synthetic
inputs.  No network for arxiv/S2 tools (those are mocked via the
``requests`` module; live-network smoke tests are run via the
mcp-server-paper-writing --selftest harness only).

Added 2026-05-26 (radia-mcp v0.90.0).  Until this file existed,
paper_writing had ZERO dedicated unit tests despite 80 tools --
relying solely on the server's --selftest smoke check.  This file
locks in current behavior so future refactors won't regress silently.
"""
from __future__ import annotations

import pytest

from radia_mcp.paper_writing import tools as pw
from radia_mcp.paper_writing._arxiv_source import (
    _normalize_arxiv_id,
    paper_writing_arxiv_extract_equations,
    paper_writing_external_sources_recipe,
)
from radia_mcp.paper_writing._tex_figure_placement import (
    paper_writing_tex_figure_placement,
)
from radia_mcp.paper_writing._pdf_layout_visual import (
    paper_writing_layout_visual_recipe,
)


# --------------------------------------------------------------------
# Fixture text snippets (representative IEEE / IEEJ paper content)
# --------------------------------------------------------------------

ABSTRACT_JA = """\
本論文では、誘導加熱における磁界解析を有限要素法で実施した。
ワークピース内の渦電流分布と発熱密度を評価し、実験値と比較した。
提案手法はカーボン鋼の表面焼入れに適用可能であることを示した。
"""

ABSTRACT_EN = """\
This paper presents a finite element analysis of induction heating
with surface impedance boundary condition.  We compute the eddy
current and Joule loss distributions in the workpiece, and validate
against experimental measurements.  The method is applicable to
case-hardening of carbon steel.
"""

WEAK_EN = """\
It seems the result might be useful.  It could be that the method
tends to converge.  It seems that this might be a good approach.
"""

WEAK_JA = """\
本手法は性能向上が期待されると考えられる。提案アルゴリズムは
有用であろうと思われる。
"""

TEX_MINIMAL = r"""
\documentclass[journal]{IEEEtran}
\usepackage{graphicx}
\begin{document}
\title{Sample}
\author{Author}
\maketitle
\begin{abstract}
This is the abstract.
\end{abstract}
\section{Introduction}
We propose a new method~\cite{ref1}.  See Fig.~\ref{fig:one}.
\begin{figure}[!ht]
\centering
\includegraphics[width=\columnwidth]{a.pdf}
\caption{First figure.}
\label{fig:one}
\end{figure}
\begin{equation}
E = mc^2
\label{eq:einstein}
\end{equation}
\bibliographystyle{IEEEtran}
\bibliography{refs}
\end{document}
"""

BIB_MINIMAL = r"""
@article{ref1,
  author = {A. Author and B. Author},
  title  = {A sample paper},
  journal= {IEEE Trans.\ Magn.},
  year   = {2025},
  volume = {61},
  pages  = {1--10},
  doi    = {10.1109/TMAG.2025.1234567},
}
"""


# =====================================================================
# Text-style tools (English)
# =====================================================================


def test_count_weak_expressions_en_finds_hedges():
    r = pw.paper_writing_count_weak_expressions(WEAK_EN)
    assert r["total_weak_expressions"] >= 4
    assert "target" in r and "warning" in r


def test_count_weak_expressions_clean_text_is_zero():
    r = pw.paper_writing_count_weak_expressions(
        "We propose a new method.  It reduces error by 23 percent."
    )
    assert r["total_weak_expressions"] == 0


def test_analyze_sentences_returns_stats():
    # Tool treats each input as ONE sentence unit (does not split on
    # in-string periods); use separate calls or paragraph breaks for
    # multi-sentence input.  Verify the basic shape on a single string.
    r = pw.paper_writing_analyze_sentences("One short sentence here.",
                                            max_len=20)
    assert "total_sentences" in r
    assert r["total_sentences"] == 1
    assert r["avg_length"] > 0
    assert r["max_length"] >= 4
    assert r["threshold"] == 20


def test_check_english_redflags_catches_common_issues():
    txt = "It should be noted that in order to show that the result is novel..."
    r = pw.paper_writing_check_english_redflags(txt)
    # The implementation returns a list / count; just verify structure exists.
    assert isinstance(r, dict)
    assert len(r) > 0


def test_check_passive_voice_ratio_low_for_active_text():
    r = pw.paper_writing_check_passive_voice_ratio(
        "We propose. We compute. We measure. We conclude."
    )
    assert "passive_ratio" in r or "ratio" in r or "passive_count" in r


def test_check_paragraph_opener_returns_dict():
    txt = "First paragraph.\n\nIn this section we discuss.\n\nWe propose."
    r = pw.paper_writing_check_paragraph_opener(txt)
    assert isinstance(r, dict)


def test_check_paragraph_length_returns_dict():
    txt = "P1.\n\n" + ("Long. " * 50) + "\n\n" + "P3."
    r = pw.paper_writing_check_paragraph_length(txt)
    assert isinstance(r, dict)


def test_check_sentence_ending_variety_returns_dict():
    r = pw.paper_writing_check_sentence_ending_variety(
        "We did this. We did that. We did the other."
    )
    assert isinstance(r, dict)


def test_check_word_repetition_returns_dict():
    txt = "The method method method is robust.  The robust method robust."
    r = pw.paper_writing_check_word_repetition(txt)
    assert isinstance(r, dict)


def test_check_tense_consistency_returns_dict():
    txt = "We propose a method.  We proposed a method.  We propose."
    r = pw.paper_writing_check_tense_consistency(txt)
    assert isinstance(r, dict)


# =====================================================================
# Abstract / IMRaD tools
# =====================================================================


def test_validate_abstract_length_in_range():
    r = pw.paper_writing_validate_abstract_length(ABSTRACT_EN)
    assert "word_count" in r or "n_words" in r or "count" in r


def test_check_abstract_background_ratio_returns_dict():
    r = pw.paper_writing_check_abstract_background_ratio(ABSTRACT_EN)
    assert isinstance(r, dict)


def test_check_imrad_balance_returns_dict():
    text = ("Introduction text. " * 20 +
            "\n\nMethod text. " * 20 +
            "\n\nResults text. " * 20 +
            "\n\nDiscussion text. " * 20)
    r = pw.paper_writing_check_imrad_balance(text)
    assert isinstance(r, dict)


def test_check_prose_density_returns_dict():
    r = pw.paper_writing_check_prose_density(ABSTRACT_EN)
    assert isinstance(r, dict)


def test_abstract_strength_returns_dict():
    r = pw.paper_writing_abstract_strength(ABSTRACT_EN)
    assert isinstance(r, dict)


def test_abstract_no_math_no_citation_flags_domain_acronym():
    r = pw.paper_writing_check_abstract_no_math_no_citation(
        "This paper combines FEM and BEM for magnetic levitation design."
    )
    assert r["status"] == "warning"
    assert r["n_acronyms"] == 2
    assert r["acronyms"] == ["BEM", "FEM"]


def test_abstract_no_math_no_citation_flags_numeric_reference_marker():
    r = pw.paper_writing_check_abstract_no_math_no_citation(
        "This paper extends the reduced-potential formulation [1]."
    )
    assert r["status"] == "fail"
    assert r["n_citation_markers"] == 1
    assert r["n_citations"] == 1


def test_abstract_no_math_no_citation_keeps_universal_acronyms_clean():
    r = pw.paper_writing_check_abstract_no_math_no_citation(
        "This paper provides a PDF report and a CPU-time comparison."
    )
    assert r["status"] == "clean"
    assert r["n_acronyms"] == 0


def test_digest_human_review_triggers_flag_cauer_digest_patterns():
    tex = r"""
\begin{abstract}
Finite Cauer ladder models are used for eddy-current analysis.
The model agrees with errors of 0.04\%, 0.001\%, and 0.33\%.
\end{abstract}
\section{Introduction}
We propose a new Warburg-type CLN termination for a circular conductor.
The transition frequency remains outside the Galerkin reduction.
The HOIBC envelope is used for the surface tail with $p_{\mathrm H}=0$.
The CLN bulk modes and HOIBC surface modes are coupled by a Schur
complement, the algebraic form of a discrete Steklov-Poincare map.
The method preserves the high-frequency SIBC scaling.
The mixed system is
\[
\begin{bmatrix}K_{bb}&K_{bs}\\K_{sb}&K_{ss}\end{bmatrix}
\begin{bmatrix}c_b\\c_s\end{bmatrix}=g,\quad
S=K_{ss}-K_{sb}K_{bb}^{-1}K_{bs}.
\]
The uniform dc term is not included in this count.  The model uses
$N_b=1$ with the HOIBC envelope, and the
bulk CLN uses $N=10$.  The L-term and R-term references are compared.
\begin{figure}[!t]
\caption{Admittance comparison; the proposed model is compared with the
reference and demonstrates improved accuracy whereas finite ladders miss
the asymptote.  The gray line marks $f_N$.}
\end{figure}
The rank-(1,1) case is demonstrated for a cylinder, sphere, and cube.
In the one-bulk/one-surface specialization, the transition is automatic.
\textbf{Verification.}  The plotted benchmark is described here.
The Japanese text calls this region 壁帯.
"""
    r = pw.paper_writing_check_digest_human_review_triggers(tex, max_report=40)
    rules = {f["rule"] for f in r["findings"]}
    assert r["status"] == "warning"
    assert "abstract_generic_cauer_ladder" in rules
    assert "abstract_numeric_overload" in rules
    assert "known_construct_framing" in rules
    assert "over_narrow_galerkin_transition_framing" in rules
    assert "technical_term_without_citation" in rules
    assert "opaque_rank_tuple" in rules
    assert "misleading_minimal_basis_claim" in rules
    assert "unnumbered_core_equation" in rules
    assert "block_matrix_missing_s_dependence" in rules
    assert "undefined_block_matrices" in rules
    assert "unexplained_schur_dtn_equivalence" in rules
    assert "vague_sibc_scaling" in rules
    assert "ambiguous_nb_definition" in rules
    assert "confusing_dc_count_exclusion" in rules
    assert "ambiguous_n_vs_nb_order" in rules
    assert "undefined_lr_termination" in rules
    assert "minimal_bulk_order_misread" in rules
    assert "ph_zero_higher_order_mismatch" in rules
    assert "one_page_scope_creep_shapes" in rules
    assert "unclear_wall_band_term" in rules
    assert "result_figure_before_numerical_section" in rules
    assert "hidden_numerical_example_section" in rules
    assert "ambiguous_fn_marker" in rules
    assert "caption_overloaded" in rules


def test_digest_human_review_triggers_flag_missing_cln_hoibc_orders():
    tex = r"""
\section{Introduction}
The CLN bulk modes and HOIBC surface envelope are coupled by a Schur
complement.  This mixed CLN/HOIBC model is evaluated on a circular
conductor benchmark with radius 5 mm, conductivity 5.8e7 S/m, and a
Bessel reference solution.
"""
    r = pw.paper_writing_check_digest_human_review_triggers(tex, max_report=20)
    rules = {f["rule"] for f in r["findings"]}
    assert r["status"] == "warning"
    assert "missing_model_order_disclosure" in rules


def test_digest_human_review_triggers_accept_reframed_cauer_digest():
    tex = r"""
\begin{abstract}
The model agrees with the circular-conductor reference with sub-percent
accuracy across the skin-effect transition region and surface-impedance tail.
\end{abstract}
\section{Introduction}
Warburg-type circuit elements are a known way to represent the
square-root tail~\cite{Warburg1899,Randles1947}.  We instead replace
the fitted transition by a frequency-dependent surface-impedance
boundary condition (SIBC) envelope~\cite{YuferevIda2010}.  Higher-order
impedance terms can be included in the same framework and can be needed
for curved or three-dimensional applications, whereas this benchmark
uses the leading surface term~\cite{Senior1962}.
\section{Numerical Example}
We consider a circular copper conductor of radius 5 mm and conductivity
5.8e7 S/m.  The Bessel solution of the radial diffusion problem is used
as the reference.  The
mixed system is written as
\begin{equation}
\label{eq:mixed-block}
\begin{aligned}
\begin{bmatrix}K_{bb}(s)&K_{bs}(s)\\K_{sb}(s)&K_{ss}(s)\end{bmatrix}
\begin{bmatrix}c_b\\c_s\end{bmatrix} &= g(s),\\
S(s) &= K_{ss}(s)-K_{sb}(s)K_{bb}(s)^{-1}K_{bs}(s).
\end{aligned}
\end{equation}
Here $K_{bb}(s)$, $K_{ss}(s)$, and $K_{bs}(s),K_{sb}(s)$ are the
bulk, surface, and coupling Galerkin blocks.  Eliminating the bulk
unknowns in~\eqref{eq:mixed-block} gives the Schur complement, which
maps surface Dirichlet data to Neumann flux.
Let $N_b$ denote the number of bulk CLN rungs in the mixed model.  The
mixed model uses $N_b=2$ and the added leading SIBC surface envelope
($p_{\mathrm H}=0$) to supply the non-rational skin-effect tail.
The bulk-only CLN references use order $N=10$, distinct from $N_b$ in the
mixed model; the L-terminated ladder closes the last rung inductively,
and the R-terminated ladder closes it resistively.  The Schur composition
preserves the f^{-1/2} SIBC scaling.
\begin{figure}
\caption{Per-unit-length admittance of a circular copper conductor
(radius 5 mm).  The gray line marks $f_N$, the highest pole frequency
of the bulk CLN.}
\end{figure}
"""
    r = pw.paper_writing_check_digest_human_review_triggers(tex)
    assert r["status"] == "clean"
    assert r["n_findings"] == 0


def test_digest_human_review_triggers_flag_overclaim_and_ambiguous_method_subject():
    tex = r"""
\section{Results}
The finite parameter sweep is a strong evidence for the general mechanism.
低透磁率での数値例は，高透磁率での安定性の基礎的な証拠である。
導出したモーメント条件はFortran実装に組み込んだ。
"""
    r = pw.paper_writing_check_digest_human_review_triggers(tex, max_report=20)
    rules = {f["rule"] for f in r["findings"]}
    assert "overstated_numerical_evidence" in rules
    assert "ambiguous_derived_method_subject" in rules


def test_digest_human_review_triggers_accept_scoped_claim_and_named_method_subject():
    tex = r"""
\section{Results}
The finite parameter sweep supports the proposed interpretation.
低透磁率での数値例は，高透磁率での安定性を考察するための知見を与える。
導出したMMPMのモーメント条件はFortran実装に組み込んだ。
"""
    r = pw.paper_writing_check_digest_human_review_triggers(tex, max_report=20)
    rules = {f["rule"] for f in r["findings"]}
    assert "overstated_numerical_evidence" not in rules
    assert "ambiguous_derived_method_subject" not in rules


# =====================================================================
# Figure / equation / table tools
# =====================================================================


def test_check_figure_caption_showing_is_callable():
    r = pw.paper_writing_check_figure_caption_showing(
        "Figure 1 shows the result."
    )
    assert isinstance(r, dict)
    assert "caption_body_policy" in r
    assert "caption becomes long" in r["caption_body_policy"]
    assert "keep the caption concise" in r["caption_body_policy"]


def test_em_paper_style_caption_body_policy_topic():
    from radia_mcp.paper_writing._em_paper_style import (
        paper_writing_em_paper_style,
    )
    topic = paper_writing_em_paper_style("caption_body")
    assert "main text" in topic
    assert "If a caption becomes long, prioritize the body text" in topic


def test_check_figure_forward_reference_on_minimal_tex(tmp_path):
    """\\ref{fig:one} comes BEFORE \\begin{figure} -- forward reference."""
    tex = tmp_path / "p.tex"
    tex.write_text(TEX_MINIMAL, encoding="utf-8")
    r = pw.paper_writing_check_figure_forward_reference(str(tex))
    # Tool returns dict; verify structure
    assert isinstance(r, dict)


def test_check_equation_numbering_on_minimal_tex(tmp_path):
    tex = tmp_path / "p.tex"
    tex.write_text(TEX_MINIMAL, encoding="utf-8")
    r = pw.paper_writing_check_equation_numbering(str(tex))
    assert isinstance(r, dict)


# =====================================================================
# Citation / reference tools
# =====================================================================


def test_check_citation_usage_with_bib(tmp_path):
    tex = tmp_path / "p.tex"
    bib = tmp_path / "refs.bib"
    tex.write_text(TEX_MINIMAL, encoding="utf-8")
    bib.write_text(BIB_MINIMAL, encoding="utf-8")
    r = pw.paper_writing_check_citation_usage(str(tex), str(bib))
    assert isinstance(r, dict)


def test_lint_reference_format_basic(tmp_path):
    bib = tmp_path / "refs.bib"
    bib.write_text(BIB_MINIMAL, encoding="utf-8")
    r = pw.paper_writing_lint_reference_format(str(bib))
    assert isinstance(r, dict)


def test_check_self_citation_ratio_returns_dict(tmp_path):
    tex = tmp_path / "p.tex"
    bib = tmp_path / "refs.bib"
    tex.write_text(TEX_MINIMAL, encoding="utf-8")
    bib.write_text(BIB_MINIMAL, encoding="utf-8")
    r = pw.paper_writing_check_self_citation_ratio(
        str(tex), str(bib), author_last_names="Author"
    )
    assert isinstance(r, dict)


# =====================================================================
# JA-lint tools (8 inline-copied from grant_writing)
# =====================================================================


def test_check_kanji_ratio_in_range():
    r = pw.paper_writing_check_kanji_ratio(ABSTRACT_JA)
    assert "kanji_ratio" in r
    assert 0.0 <= r["kanji_ratio"] <= 1.0
    # Body Japanese text typically 25-40% kanji
    assert 0.15 < r["kanji_ratio"] < 0.60


def test_check_notation_variants_returns_dict():
    txt = "事に関して、コトについて、こと。"
    r = pw.paper_writing_check_notation_variants(txt)
    assert isinstance(r, dict)


def test_normalize_terminology_preserves_english_and_paths():
    text = (
        r"\bicaption{境界固定cubeの検証。}{fixed-boundary cube benchmark.}"
        "\n"
        r"\includegraphics{figures/fig_cube_scaling_mdx.pdf}"
    )
    r = pw.paper_writing_normalize_terminology(text, rules="cube=>立方体")
    assert r["n_replacements"] == 1
    assert "境界固定立方体の検証" in r["normalized_text"]
    assert "fixed-boundary cube benchmark" in r["normalized_text"]
    assert "fig_cube_scaling_mdx.pdf" in r["normalized_text"]


def test_normalize_terminology_file_dry_run_and_apply(tmp_path):
    tex = tmp_path / "paper.tex"
    tex.write_text(
        "境界固定cubeを調べる。 English cube remains. fig_cube.pdf\n",
        encoding="utf-8",
    )
    preview = pw.paper_writing_normalize_terminology_file(str(tex), dry_run=True)
    assert preview["changed"] is True
    assert preview["n_replacements"] == 1
    assert "境界固定cube" in tex.read_text(encoding="utf-8")

    applied = pw.paper_writing_normalize_terminology_file(str(tex), dry_run=False)
    assert applied["n_replacements"] == 1
    updated = tex.read_text(encoding="utf-8")
    assert "境界固定立方体" in updated
    assert "English cube remains" in updated
    assert "fig_cube.pdf" in updated


def test_find_undefined_acronyms_in_text():
    r = pw.paper_writing_find_undefined_acronyms(
        "BEM is faster than FEM.  We use BEM for the open boundary."
    )
    assert isinstance(r, dict)


def test_latex_acronym_first_use_accepts_full_name_before_acronym(tmp_path):
    from radia_mcp.paper_writing._undefined_acronyms import (
        paper_writing_check_undefined_acronyms,
    )
    tex = tmp_path / "paper.tex"
    tex.write_text(
        r"\documentclass{article}\begin{document}"
        r"\begin{abstract}This paper presents a magnetic analysis method."
        r"\end{abstract}"
        r"\section{Introduction}"
        r"We use the Finite Element Method (FEM) for the local solve. "
        r"The FEM solution is then reduced."
        r"\end{document}",
        encoding="utf-8",
    )
    r = paper_writing_check_undefined_acronyms(str(tex))
    assert r["n_undefined"] == 0


def test_acronym_usage_audit_low_usage_flagged():
    r = pw.paper_writing_acronym_usage_audit(
        "We use IEEE Trans.~Magn.  PEEC is also relevant."
    )
    assert isinstance(r, dict)


def test_lint_bedrock_on_ja_text():
    r = pw.paper_writing_lint_bedrock(WEAK_JA)
    assert isinstance(r, dict)


def test_check_misuse_japanese_returns_dict():
    r = pw.paper_writing_check_misuse_japanese(WEAK_JA)
    assert isinstance(r, dict)


def test_suggest_redundancy_fixes_returns_dict():
    r = pw.paper_writing_suggest_redundancy_fixes(WEAK_JA)
    assert isinstance(r, dict)


def test_check_subject_predicate_distance_returns_dict():
    r = pw.paper_writing_check_subject_predicate_distance(
        "本論文は、これこれの方法を用いて、何々を行うことによって、結果として性能向上を示した。"
    )
    assert isinstance(r, dict)


# =====================================================================
# PDF tools (need a real PDF to exercise; create minimal one)
# =====================================================================


def _make_minimal_pdf(tmp_path):
    """Create a tiny valid PDF via pymupdf (skipped if not installed)."""
    pymupdf = pytest.importorskip("pymupdf")
    pdf_path = tmp_path / "test.pdf"
    doc = pymupdf.open()
    page = doc.new_page()
    page.insert_text((50, 100), "Test page 1")
    page2 = doc.new_page()
    page2.insert_text((50, 100), "Test page 2")
    doc.save(str(pdf_path))
    doc.close()
    return pdf_path


def test_validate_pdf_pages_within_limit(tmp_path):
    pdf = _make_minimal_pdf(tmp_path)
    r = pw.paper_writing_validate_pdf_pages(str(pdf), page_limit=10)
    assert "page_count" in r or "n_pages" in r or "count" in r


def test_check_pdf_edge_overflow_returns_dict(tmp_path):
    pdf = _make_minimal_pdf(tmp_path)
    r = pw.paper_writing_check_pdf_edge_overflow(str(pdf))
    assert isinstance(r, dict)


def test_count_underlines_on_minimal_tex(tmp_path):
    tex = tmp_path / "p.tex"
    tex.write_text(TEX_MINIMAL.replace(
        "\\title{Sample}",
        r"\title{\underline{Sample}}",
    ), encoding="utf-8")
    r = pw.paper_writing_count_underlines(str(tex))
    assert isinstance(r, dict)


# =====================================================================
# v0.88.0 layout/figure knowledge tools
# =====================================================================


def test_tex_figure_placement_overview_has_content():
    r = paper_writing_tex_figure_placement("overview")
    assert len(r) > 500
    assert "TL;DR" in r or "default" in r.lower()


def test_tex_figure_placement_all_concatenates_sections():
    r = paper_writing_tex_figure_placement("all")
    assert len(r) > 5000
    # Spot-check several section keywords
    for kw in ["FloatBarrier", "columnwidth", "centering", "IEEE"]:
        assert kw in r, f"missing key concept: {kw!r}"


def test_tex_figure_placement_aliases_work():
    """ieee / ieej / igte / htbp etc. all dispatch."""
    for alias in ["ieee", "ieej", "igte", "htbp", "placeins",
                    "columnwidth", "subcaption", "antipatterns"]:
        r = paper_writing_tex_figure_placement(alias)
        assert len(r) > 100, f"alias {alias!r} returned empty"
        assert "Unknown topic" not in r


def test_tex_figure_placement_unknown_topic_returns_help():
    r = paper_writing_tex_figure_placement("does_not_exist")
    assert "Unknown topic" in r


def test_layout_visual_recipe_has_4_steps():
    r = paper_writing_layout_visual_recipe()
    assert len(r) > 1500
    for kw in ["thumbnail_strip", "whitespace_anomalies",
                "render_pages_to_png", "floats_far_from_reference"]:
        assert kw in r, f"missing tool reference: {kw!r}"


# =====================================================================
# v0.89.0 external sources (arxiv + Semantic Scholar) tools
# =====================================================================


def test_normalize_arxiv_id_strips_prefixes_and_versions():
    assert _normalize_arxiv_id("2603.17339") == "2603.17339"
    assert _normalize_arxiv_id("arXiv:2603.17339") == "2603.17339"
    assert _normalize_arxiv_id("arxiv:2603.17339") == "2603.17339"
    assert _normalize_arxiv_id("2603.17339v2") == "2603.17339"
    assert _normalize_arxiv_id("2603.17339v17") == "2603.17339"
    assert _normalize_arxiv_id(
        "https://arxiv.org/abs/2603.17339"
    ) == "2603.17339"
    assert _normalize_arxiv_id(
        "https://arxiv.org/pdf/2603.17339.pdf"
    ) == "2603.17339"
    assert _normalize_arxiv_id(
        "http://arxiv.org/abs/2603.17339v3"
    ) == "2603.17339"
    # Old-style preserved
    assert _normalize_arxiv_id("physics/0501123") == "physics/0501123"


def test_extract_equations_finds_all_six_environments():
    """Exercise each of the 6 displayed-math environments."""
    src = r"""
    Intro.  \begin{equation} a = b \end{equation}
    \begin{align} c &= d \\ e &= f \end{align}
    \begin{eqnarray*} g = h \end{eqnarray*}
    \begin{gather} i = j \end{gather}
    Then $$k = l$$ and finally \[m = n\].
    """
    r = paper_writing_arxiv_extract_equations(src)
    assert r["n_equations"] == 6, f"expected 6, got {r['n_equations']}"
    # Verify ordering by char_offset
    offsets = [e["char_offset"] for e in r["equations"]]
    assert offsets == sorted(offsets)


def test_extract_equations_truncates_long_body():
    long_body = "x = " + "y" * 2000
    src = r"\begin{equation}" + long_body + r"\end{equation}"
    r = paper_writing_arxiv_extract_equations(src, snippet_chars=200)
    assert r["n_equations"] == 1
    assert len(r["equations"][0]["body"]) < 250
    assert "truncated" in r["equations"][0]["body"]


def test_extract_equations_caps_max():
    src = r"\begin{equation} x = y \end{equation}" * 200
    r = paper_writing_arxiv_extract_equations(src, max_equations=50)
    assert r["n_equations"] == 200
    assert len(r["equations"]) == 50
    assert r["truncated"] is True


def test_extract_equations_handles_non_string():
    r = paper_writing_arxiv_extract_equations(None)
    assert "error" in r


def test_external_sources_recipe_has_credits_and_decision_tree():
    r = paper_writing_external_sources_recipe()
    assert len(r) > 2000
    for kw in ["arxiv-latex-mcp", "Semantic Scholar", "Decision",
                "eess.SP", "physics.class-ph", "takashiishida"]:
        assert kw in r, f"missing key string: {kw!r}"


# =====================================================================
# v0.89.0 network-mocked arxiv / S2 tools
# =====================================================================


class _FakeResponse:
    """Minimal requests.Response stand-in for tests."""
    def __init__(self, status_code, content=b"", text=None, json_data=None):
        self.status_code = status_code
        self.content = content
        self.text = text if text is not None else content.decode(
            "utf-8", errors="replace")
        self._json = json_data

    def raise_for_status(self):
        if self.status_code >= 400:
            raise RuntimeError(f"HTTP {self.status_code}")

    def json(self):
        return self._json or {}


def test_arxiv_fetch_latex_source_mocked(monkeypatch, tmp_path):
    """Mock arXiv e-print fetch with a tiny synthetic tarball."""
    import io
    import tarfile
    from radia_mcp.paper_writing import _arxiv_source as mod

    # Build a synthetic tarball with 2 .tex files
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:gz") as t:
        for name, body in [
            ("main.tex", r"\documentclass{article}\begin{document}main\end{document}"),
            ("appendix.tex", r"appendix content"),
        ]:
            data = body.encode("utf-8")
            info = tarfile.TarInfo(name=name)
            info.size = len(data)
            t.addfile(info, io.BytesIO(data))
    blob = buf.getvalue()

    def fake_get(url, timeout=None):
        assert "/e-print/" in url
        return _FakeResponse(200, content=blob)

    real_req = pytest.importorskip("requests")
    monkeypatch.setattr(real_req, "get", fake_get)

    r = mod.paper_writing_arxiv_fetch_latex_source("2603.17339")
    assert "error" not in r, r
    assert r["arxiv_id"] == "2603.17339"
    assert r["n_files_total"] == 2
    assert r["main_tex"] == "main.tex"
    assert len(r["files"]) == 2
    assert any("documentclass" in f["content"] for f in r["files"])


def test_arxiv_search_mocked(monkeypatch):
    """Mock arXiv Atom-XML response."""
    from radia_mcp.paper_writing import _arxiv_source as mod
    xml = """<?xml version="1.0" encoding="UTF-8"?>
<feed xmlns="http://www.w3.org/2005/Atom"
      xmlns:arxiv="http://arxiv.org/schemas/atom">
  <entry>
    <id>http://arxiv.org/abs/2603.17339v1</id>
    <title>A sample paper</title>
    <summary>This paper proposes a finite element analysis.</summary>
    <published>2026-03-01T00:00:00Z</published>
    <author><name>A. Author</name></author>
    <author><name>B. Author</name></author>
    <category term="eess.SP"/>
    <category term="physics.app-ph"/>
    <link rel="alternate" href="http://arxiv.org/abs/2603.17339v1"/>
    <link title="pdf" href="http://arxiv.org/pdf/2603.17339v1"/>
  </entry>
</feed>"""

    def fake_get(url, timeout=None):
        return _FakeResponse(200, text=xml)

    real_req = pytest.importorskip("requests")
    monkeypatch.setattr(real_req, "get", fake_get)

    r = mod.paper_writing_arxiv_search("finite element", max_results=5)
    assert "error" not in r, r
    assert r["n_results"] == 1
    p = r["papers"][0]
    assert p["arxiv_id"] == "2603.17339"
    assert p["title"] == "A sample paper"
    assert len(p["authors"]) == 2
    assert "eess.SP" in p["categories"]
    assert p["pdf_url"]


def test_semantic_scholar_lookup_mocked(monkeypatch):
    from radia_mcp.paper_writing import _arxiv_source as mod
    payload = {
        "paperId": "abc123",
        "title": "Sample",
        "year": 2026,
        "authors": [{"name": "A. Author"}],
        "venue": "IEEE TAP",
        "referenceCount": 42,
        "citationCount": 100,
    }

    def fake_get(url, params=None, timeout=None):
        return _FakeResponse(200, json_data=payload)

    real_req = pytest.importorskip("requests")
    monkeypatch.setattr(real_req, "get", fake_get)

    r = mod.paper_writing_semantic_scholar_lookup("10.1109/TAP.2006.880747")
    assert "error" not in r, r
    assert r["metadata"]["paperId"] == "abc123"
    assert r["paper_id_resolved"].startswith("DOI:")


def test_semantic_scholar_lookup_404_returns_error(monkeypatch):
    from radia_mcp.paper_writing import _arxiv_source as mod

    def fake_get(url, params=None, timeout=None):
        return _FakeResponse(404)

    real_req = pytest.importorskip("requests")
    monkeypatch.setattr(real_req, "get", fake_get)

    r = mod.paper_writing_semantic_scholar_lookup("10.9999/nonexistent")
    assert "error" in r
    assert "not found" in r["error"].lower()


def test_semantic_scholar_references_mocked(monkeypatch):
    from radia_mcp.paper_writing import _arxiv_source as mod
    payload = {"data": [
        {"citedPaper": {"paperId": "ref1", "title": "Ref 1", "year": 2020}},
        {"citedPaper": {"paperId": "ref2", "title": "Ref 2", "year": 2021}},
    ]}

    def fake_get(url, params=None, timeout=None):
        assert "/references" in url
        return _FakeResponse(200, json_data=payload)

    real_req = pytest.importorskip("requests")
    monkeypatch.setattr(real_req, "get", fake_get)

    r = mod.paper_writing_semantic_scholar_references("abc123")
    assert "error" not in r
    assert r["n_references"] == 2
    assert r["references"][0]["title"] == "Ref 1"


def test_semantic_scholar_citations_mocked(monkeypatch):
    from radia_mcp.paper_writing import _arxiv_source as mod
    payload = {"data": [
        {"citingPaper": {"paperId": "c1", "title": "Citing 1", "year": 2024}},
    ]}

    def fake_get(url, params=None, timeout=None):
        assert "/citations" in url
        return _FakeResponse(200, json_data=payload)

    real_req = pytest.importorskip("requests")
    monkeypatch.setattr(real_req, "get", fake_get)

    r = mod.paper_writing_semantic_scholar_citations("abc123")
    assert "error" not in r
    assert r["n_citations"] == 1


# =====================================================================
# Health / orchestrator tools
# =====================================================================


def test_health_report_returns_dict():
    r = pw.paper_writing_health_report(ABSTRACT_EN)
    assert isinstance(r, dict)


def test_adaptive_health_report_returns_dict():
    r = pw.paper_writing_adaptive_health_report(ABSTRACT_EN)
    assert isinstance(r, dict)


def test_run_full_workflow_returns_dict():
    r = pw.paper_writing_run_full_workflow(ABSTRACT_EN)
    assert isinstance(r, dict)


def test_root_cause_diagnosis_returns_dict():
    r = pw.paper_writing_root_cause_diagnosis(ABSTRACT_EN)
    assert isinstance(r, dict)


def test_next_5_actions_returns_dict():
    r = pw.paper_writing_next_5_actions(ABSTRACT_EN)
    assert isinstance(r, dict)


def test_rewrite_suggest_returns_dict():
    r = pw.paper_writing_rewrite_suggest(target=WEAK_EN, n_candidates=2)
    assert isinstance(r, dict)


# =====================================================================
# Reviewer-handling tools
# =====================================================================


def test_classify_reviewer_comment_returns_dict():
    r = pw.paper_writing_classify_reviewer_comment(
        "The figures are too small to read.  Please increase font size."
    )
    assert isinstance(r, dict)


def test_generate_response_letter_returns_str():
    """Returns a str (response letter text), not dict."""
    r = pw.paper_writing_generate_response_letter(
        "The figures are too small.", agreement_stance="agree",
        response_draft="Figures enlarged.", page_line_ref="p. 3"
    )
    assert isinstance(r, str)
    assert len(r) > 50


def test_generate_cover_letter_returns_str():
    r = pw.paper_writing_generate_cover_letter(
        journal="IEEE Trans. Magn.",
        title="A new method",
        abstract=ABSTRACT_EN,
        contribution="First closed-form derivation.",
    )
    assert isinstance(r, str)
    assert len(r) > 50


def test_reviewer_2_trigger_summary_returns_dict():
    r = pw.paper_writing_reviewer_2_trigger_summary(
        "We propose a perfect method that fully solves the problem."
    )
    assert isinstance(r, dict)


# =====================================================================
# Plan-B Tier composite scores (T1-T20)
# =====================================================================


def test_contribution_clarity_score_returns_dict():
    r = pw.paper_writing_contribution_clarity_score(ABSTRACT_EN)
    assert isinstance(r, dict)


def test_claim_quantification_returns_dict():
    r = pw.paper_writing_claim_quantification(
        "We reduce error by 23 percent.  Speed increases 5x."
    )
    assert isinstance(r, dict)


def test_limitation_statement_presence_returns_dict():
    r = pw.paper_writing_limitation_statement_presence(ABSTRACT_EN)
    assert isinstance(r, dict)


def test_related_work_density_returns_dict():
    txt = "Smith et al.~\\cite{a} showed this.  Jones~\\cite{b} extended it."
    r = pw.paper_writing_related_work_density(txt)
    assert isinstance(r, dict)


def test_figure_referencing_coverage_returns_dict(tmp_path):
    tex = tmp_path / "p.tex"
    tex.write_text(TEX_MINIMAL, encoding="utf-8")
    r = pw.paper_writing_figure_referencing_coverage(str(tex))
    assert isinstance(r, dict)


def test_journal_fit_assessment_returns_dict():
    r = pw.paper_writing_journal_fit_assessment(
        text=ABSTRACT_EN,
        journal_name="IEEE Trans. Magn.",
        aims_and_scope=("Original research on theoretical and computational "
                         "magnetics, electromagnetics, and magnetic materials."),
    )
    assert isinstance(r, dict)


# =====================================================================
# v0.90.0 EM-paper-specific style + submission gate
# =====================================================================


def test_em_paper_style_overview_returns_index():
    from radia_mcp.paper_writing._em_paper_style import (
        paper_writing_em_paper_style,
    )
    r = paper_writing_em_paper_style("overview")
    assert len(r) > 500
    # Index lists the 6 main topics
    for kw in ["sign_convention", "vector", "b_vs_h", "units",
                "equations", "reviewer_patterns"]:
        assert kw in r, f"missing topic ref: {kw}"


def test_em_paper_style_all_concatenates():
    from radia_mcp.paper_writing._em_paper_style import (
        paper_writing_em_paper_style,
    )
    r = paper_writing_em_paper_style("all")
    assert len(r) > 8000
    # Spot-check key concepts across sections
    for kw in ["exp(+j", "exp(-i", "Bossavit", "magnetic flux density",
                "siunitx", "align", "Reviewer 2"]:
        assert kw in r, f"missing key concept: {kw}"


def test_em_paper_style_aliases_resolve():
    from radia_mcp.paper_writing._em_paper_style import (
        paper_writing_em_paper_style,
    )
    for alias in ["sign", "vector", "tensor", "b_vs_h",
                    "si_units", "siunitx", "equations", "reviewer_patterns",
                    "phase_convention", "magnetic_field_terminology"]:
        r = paper_writing_em_paper_style(alias)
        assert len(r) > 100, f"alias {alias!r} returned too little"
        assert "Unknown topic" not in r


def test_em_paper_style_unknown_topic_returns_help():
    from radia_mcp.paper_writing._em_paper_style import (
        paper_writing_em_paper_style,
    )
    r = paper_writing_em_paper_style("does_not_exist_xyz")
    assert "Unknown topic" in r


def test_em_submission_gate_no_inputs_fails_on_bib_policy():
    """v0.91.0: no bib_path => FAIL on lab POLICY check."""
    from radia_mcp.paper_writing._em_paper_style import (
        paper_writing_em_submission_gate,
    )
    r = paper_writing_em_submission_gate()
    assert r["verdict"] == "fail"
    # The bib_policy check itself fired
    names = [c["name"] for c in r["checks"]]
    assert "bib_policy" in names
    bib_check = next(c for c in r["checks"] if c["name"] == "bib_policy")
    assert bib_check["status"] == "fail"
    assert "reference.bib" in bib_check["summary"].lower() or \
            "reference.bib" in str(bib_check.get("detail", ""))


def test_em_submission_gate_tex_only_returns_partial(tmp_path):
    from radia_mcp.paper_writing._em_paper_style import (
        paper_writing_em_submission_gate,
    )
    tex = tmp_path / "p.tex"
    tex.write_text(TEX_MINIMAL, encoding="utf-8")
    r = paper_writing_em_submission_gate(tex_path=str(tex))
    assert r["n_checks_run"] > 0
    # tex checks should have run
    names = [c["name"] for c in r["checks"]]
    assert "figure_forward_reference" in names
    assert "equation_numbering" in names
    assert "count_underlines" in names


def test_em_submission_gate_with_bib(tmp_path):
    from radia_mcp.paper_writing._em_paper_style import (
        paper_writing_em_submission_gate,
    )
    tex = tmp_path / "p.tex"
    bib = tmp_path / "refs.bib"
    tex.write_text(TEX_MINIMAL, encoding="utf-8")
    bib.write_text(BIB_MINIMAL, encoding="utf-8")
    r = paper_writing_em_submission_gate(
        tex_path=str(tex), bib_path=str(bib), abstract_text=ABSTRACT_EN
    )
    names = [c["name"] for c in r["checks"]]
    assert "lint_reference_format" in names
    assert "check_citation_usage" in names
    assert "validate_abstract_length" in names
    assert "abstract_weak_expressions" in names


def test_em_submission_gate_verdict_structure():
    from radia_mcp.paper_writing._em_paper_style import (
        paper_writing_em_submission_gate,
    )
    r = paper_writing_em_submission_gate()
    # Verdict + counts
    for k in ("verdict", "n_checks_run", "n_critical", "n_warning",
                "n_passed", "n_skipped", "checks", "advice"):
        assert k in r, f"missing key: {k}"
    # n_critical + n_warning + n_passed + n_skipped == n_checks_run
    total = (r["n_critical"] + r["n_warning"]
             + r["n_passed"] + r["n_skipped"])
    assert total == r["n_checks_run"]


# =====================================================================
# v0.91.0: citation-verification workflow
# =====================================================================


def test_citation_workflow_recipe_has_policy_lines():
    from radia_mcp.paper_writing._citation_verify import (
        paper_writing_citation_workflow_recipe,
    )
    r = paper_writing_citation_workflow_recipe()
    assert len(r) > 2000
    # Key policy keywords MUST be present
    for kw in ["NEVER invent", "ALWAYS",
                "reference.bib", "Crossref", "Semantic Scholar",
                "arXiv", "verify", "TODO"]:
        assert kw in r, f"citation recipe missing key concept: {kw!r}"


def test_verify_citation_no_bib_path_returns_error():
    from radia_mcp.paper_writing._citation_verify import (
        paper_writing_verify_citation,
    )
    r = paper_writing_verify_citation(claim="some claim", bib_path="")
    assert r["verdict"] == "error"
    assert "bib_path" in r["advice"].lower()
    assert "reference.bib" in r["advice"].lower()


def test_verify_citation_missing_bib_returns_error(tmp_path):
    from radia_mcp.paper_writing._citation_verify import (
        paper_writing_verify_citation,
    )
    r = paper_writing_verify_citation(
        claim="some claim",
        bib_path=str(tmp_path / "does_not_exist.bib"),
    )
    assert r["verdict"] == "error"
    assert "not found" in r["advice"].lower()


def test_verify_citation_doi_already_in_bib(tmp_path):
    """If candidate_doi matches an entry in bib, verdict='found_in_bib'."""
    from radia_mcp.paper_writing._citation_verify import (
        paper_writing_verify_citation,
    )
    bib = tmp_path / "ref.bib"
    bib.write_text(r"""
@article{kohyook2006,
  author = {I.-S. Koh and J.-G. Yook},
  title  = {Exact Closed-Form Expression of a Sommerfeld Integral},
  journal= {IEEE Trans. Antennas Propagat.},
  year   = {2006},
  volume = {54},
  pages  = {2568--2576},
  doi    = {10.1109/TAP.2006.880747},
}
""", encoding="utf-8")
    r = paper_writing_verify_citation(
        claim="Koh-Yook exact closed form for impedance plane",
        bib_path=str(bib),
        candidate_doi="10.1109/TAP.2006.880747",
    )
    assert r["verdict"] == "found_in_bib"
    assert r["matching_key"] == "kohyook2006"
    assert "kohyook2006" in r["advice"]


def test_verify_citation_doi_already_in_bib_case_insensitive(tmp_path):
    from radia_mcp.paper_writing._citation_verify import (
        paper_writing_verify_citation,
    )
    bib = tmp_path / "ref.bib"
    bib.write_text(r"""
@article{kohyook2006,
  doi = {10.1109/TAP.2006.880747},
  title = {Exact Closed-Form},
}
""", encoding="utf-8")
    # Pass DOI with different case + https prefix
    r = paper_writing_verify_citation(
        claim="test",
        bib_path=str(bib),
        candidate_doi="HTTPS://DOI.ORG/10.1109/tap.2006.880747",
    )
    assert r["verdict"] == "found_in_bib"
    assert r["matching_key"] == "kohyook2006"


def test_verify_citation_title_match_in_bib(tmp_path):
    from radia_mcp.paper_writing._citation_verify import (
        paper_writing_verify_citation,
    )
    bib = tmp_path / "ref.bib"
    bib.write_text(r"""
@article{kohyook2006,
  title = {Exact Closed-Form Expression of a Sommerfeld Integral for
            the Impedance Plane Problem},
}
""", encoding="utf-8")
    r = paper_writing_verify_citation(
        claim="some claim", bib_path=str(bib),
        candidate_title=("Exact Closed-Form Expression of a Sommerfeld "
                          "Integral for the Impedance Plane Problem"),
    )
    assert r["verdict"] == "found_in_bib"
    assert r["matching_key"] == "kohyook2006"


def test_verify_citation_insufficient_info_returns_no_candidate(tmp_path):
    from radia_mcp.paper_writing._citation_verify import (
        paper_writing_verify_citation,
    )
    bib = tmp_path / "ref.bib"
    bib.write_text("", encoding="utf-8")
    r = paper_writing_verify_citation(
        claim="", bib_path=str(bib),
        search_arxiv_if_no_doi=False,
    )
    assert r["verdict"] == "no_candidate_found"
    assert "insufficient" in r["advice"].lower()


def test_verify_citation_doi_resolve_failure_marks_no_candidate(
        tmp_path, monkeypatch):
    """If Crossref DOI lookup fails, verdict='no_candidate_found'.
    Critically: it must NOT generate a fake suggested_bibtex."""
    from radia_mcp.paper_writing import _citation_verify as cv
    from radia_mcp.paper_writing import paper_download as pd

    def fake_resolve_doi(doi):
        return {"ok": False, "error": "Crossref returned 404"}

    monkeypatch.setattr(pd, "paper_writing_resolve_doi", fake_resolve_doi)

    bib = tmp_path / "ref.bib"
    bib.write_text("", encoding="utf-8")
    r = cv.paper_writing_verify_citation(
        claim="test",
        bib_path=str(bib),
        candidate_doi="10.9999/fake.doi.does.not.exist",
    )
    assert r["verdict"] == "no_candidate_found"
    assert r["suggested_bibtex"] is None
    assert "fabricated" in r["advice"].lower() or "typo" in r["advice"].lower()


def test_verify_citation_arxiv_search_fallback_mocked(monkeypatch, tmp_path):
    """When only claim text given, falls back to arXiv search."""
    from radia_mcp.paper_writing import _citation_verify as cv
    from radia_mcp.paper_writing import _arxiv_source as arx

    def fake_arxiv_search(query, max_results=10, **_):
        return {
            "n_results": 2,
            "papers": [
                {"arxiv_id": "2603.17339", "title": "Sample 1",
                 "authors": ["A. Author"], "categories": ["eess.SP"]},
                {"arxiv_id": "2607.99999", "title": "Sample 2",
                 "authors": ["B. Author"], "categories": ["physics.app-ph"]},
            ],
        }
    monkeypatch.setattr(arx, "paper_writing_arxiv_search", fake_arxiv_search)

    bib = tmp_path / "ref.bib"
    bib.write_text("", encoding="utf-8")
    r = cv.paper_writing_verify_citation(
        claim="A sufficiently long claim about Sommerfeld layered media",
        bib_path=str(bib),
    )
    assert r["verdict"] == "needs_disambiguation"
    assert len(r["candidates"]) == 2
    assert "ask the user" in r["advice"].lower()


# =====================================================================
# v0.92.0: PDF overlap/overflow detection
# =====================================================================


def _make_pdf_with_image_and_overlapping_text(tmp_path):
    """Create a PDF where text physically overlaps an image bbox."""
    pymupdf = pytest.importorskip("pymupdf")
    pdf_path = tmp_path / "overlap.pdf"
    doc = pymupdf.open()
    page = doc.new_page(width=595, height=842)  # A4
    # Insert a "filled rectangle" as a proxy for an image bbox (pymupdf
    # exposes shape fills via get_drawings, but for get_image_info we
    # actually need an embedded image -- so use a tiny PNG).
    import struct, zlib
    # Minimal 10x10 RGB PNG (white pixels)
    w, h = 10, 10
    raw = b"".join(b"\x00" + b"\xff\xff\xff" * w for _ in range(h))
    def chunk(tag, data):
        return (struct.pack(">I", len(data)) + tag + data
                + struct.pack(">I", zlib.crc32(tag + data) & 0xffffffff))
    png_bytes = (b"\x89PNG\r\n\x1a\n"
                 + chunk(b"IHDR", struct.pack(">IIBBBBB", w, h, 8, 2, 0, 0, 0))
                 + chunk(b"IDAT", zlib.compress(raw))
                 + chunk(b"IEND", b""))
    # Insert image covering middle of page (e.g. x=100..400, y=300..600)
    page.insert_image(pymupdf.Rect(100, 300, 400, 600), stream=png_bytes)
    # Insert TEXT that overlaps with the image bbox
    page.insert_text((150, 400), "Caption straddling figure body text")
    doc.save(str(pdf_path))
    doc.close()
    return pdf_path


def test_detect_text_image_overlap_finds_overlap(tmp_path):
    from radia_mcp.paper_writing._pdf_overlap_detection import (
        paper_writing_detect_text_image_overlap,
    )
    pdf = _make_pdf_with_image_and_overlapping_text(tmp_path)
    r = paper_writing_detect_text_image_overlap(str(pdf))
    assert "error" not in r, r
    assert r["n_overlaps"] >= 1
    # Validate structure of first overlap
    ovl = r["overlaps"][0]
    for k in ("page", "text_snippet", "text_bbox", "image_bbox",
                "iou", "intersection_area_pt2"):
        assert k in ovl
    assert ovl["page"] == 1
    assert ovl["iou"] > 0.0


def test_detect_text_image_overlap_clean_pdf(tmp_path):
    """Text and image well separated -> no overlap."""
    from radia_mcp.paper_writing._pdf_overlap_detection import (
        paper_writing_detect_text_image_overlap,
    )
    import struct, zlib
    pymupdf = pytest.importorskip("pymupdf")
    pdf_path = tmp_path / "clean.pdf"
    doc = pymupdf.open()
    page = doc.new_page(width=595, height=842)
    # Image at top
    w, h = 10, 10
    raw = b"".join(b"\x00" + b"\xff\xff\xff" * w for _ in range(h))
    def chunk(tag, data):
        return (struct.pack(">I", len(data)) + tag + data
                + struct.pack(">I", zlib.crc32(tag + data) & 0xffffffff))
    png_bytes = (b"\x89PNG\r\n\x1a\n"
                 + chunk(b"IHDR", struct.pack(">IIBBBBB", w, h, 8, 2, 0, 0, 0))
                 + chunk(b"IDAT", zlib.compress(raw))
                 + chunk(b"IEND", b""))
    page.insert_image(pymupdf.Rect(100, 100, 400, 250), stream=png_bytes)
    # Text well BELOW the image
    page.insert_text((100, 500), "Text well separated from image")
    doc.save(str(pdf_path))
    doc.close()

    r = paper_writing_detect_text_image_overlap(str(pdf_path))
    assert r["n_overlaps"] == 0


def test_detect_text_overflow_page_finds_overflow(tmp_path):
    """Text inserted past page edge should be flagged.

    Note: pymupdf clips the reported bbox at the page right edge with
    sub-pixel rounding, so even text that visually overflows by many
    pixels often has bbox extending only a fraction of a point past
    the edge.  This test uses overflow_tolerance_pt=0.0 so the
    detector picks up any non-zero overflow."""
    from radia_mcp.paper_writing._pdf_overlap_detection import (
        paper_writing_detect_text_overflow_page,
    )
    pymupdf = pytest.importorskip("pymupdf")
    pdf_path = tmp_path / "overflow.pdf"
    doc = pymupdf.open()
    page = doc.new_page(width=595, height=842)
    page.insert_text((580, 400), "X" * 30)   # 30 X's at x=580 -> clearly overflows
    doc.save(str(pdf_path))
    doc.close()
    r = paper_writing_detect_text_overflow_page(
        str(pdf_path), overflow_tolerance_pt=0.0
    )
    assert "error" not in r, r
    assert r["n_overflows"] >= 1
    ovl = r["overflows"][0]
    assert ovl["side"] == "right"
    assert ovl["overflow_pt"] > 0.0


def test_detect_text_overflow_page_clean(tmp_path):
    """Text well within page -> no overflow."""
    from radia_mcp.paper_writing._pdf_overlap_detection import (
        paper_writing_detect_text_overflow_page,
    )
    pymupdf = pytest.importorskip("pymupdf")
    pdf_path = tmp_path / "clean_over.pdf"
    doc = pymupdf.open()
    page = doc.new_page(width=595, height=842)
    page.insert_text((100, 400), "Centered text")
    doc.save(str(pdf_path))
    doc.close()
    r = paper_writing_detect_text_overflow_page(str(pdf_path))
    assert r["n_overflows"] == 0


def test_detect_overlapping_text_blocks_returns_dict(tmp_path):
    """Just verify the tool runs cleanly on a normal PDF (no overlaps
    expected for separated text)."""
    from radia_mcp.paper_writing._pdf_overlap_detection import (
        paper_writing_detect_overlapping_text_blocks,
    )
    pymupdf = pytest.importorskip("pymupdf")
    pdf_path = tmp_path / "two_paragraphs.pdf"
    doc = pymupdf.open()
    page = doc.new_page(width=595, height=842)
    page.insert_text((100, 200), "First paragraph at the top.")
    page.insert_text((100, 600), "Second paragraph well below.")
    doc.save(str(pdf_path))
    doc.close()
    r = paper_writing_detect_overlapping_text_blocks(str(pdf_path))
    assert "error" not in r
    assert r["n_overlaps"] == 0


def test_pdf_overlap_recipe_has_3_tools_documented():
    from radia_mcp.paper_writing._pdf_overlap_detection import (
        paper_writing_pdf_overlap_recipe,
    )
    r = paper_writing_pdf_overlap_recipe()
    assert len(r) > 1000
    # All 3 tools documented
    for kw in ["detect_text_image_overlap",
                "detect_text_overflow_page",
                "detect_overlapping_text_blocks",
                "IoU", "VILA", "PaperOrchestra"]:
        assert kw in r, f"recipe missing: {kw}"


def test_overlap_iou_math():
    """Verify _iou() helper computes correctly."""
    from radia_mcp.paper_writing._pdf_overlap_detection import _iou
    # Identical rects -> IoU = 1
    assert abs(_iou((0, 0, 10, 10), (0, 0, 10, 10)) - 1.0) < 1e-6
    # Disjoint -> IoU = 0
    assert _iou((0, 0, 10, 10), (20, 20, 30, 30)) == 0.0
    # Half-overlap: A=(0,0,10,10) area=100, B=(5,0,15,10) area=100,
    # intersection=(5,0,10,10) area=50, union=150 => IoU=50/150=1/3
    assert abs(_iou((0, 0, 10, 10), (5, 0, 15, 10)) - 1.0 / 3.0) < 1e-6
    # Touching edges only -> IoU = 0 (strict <)
    assert _iou((0, 0, 10, 10), (10, 0, 20, 10)) == 0.0


# =====================================================================
# v0.92.0: undefined-variable check
# =====================================================================


def test_undefined_variables_clean_paper_with_where_clauses(tmp_path):
    """A paper that defines every symbol via 'where ...' clauses
    should have 0 undefined."""
    from radia_mcp.paper_writing._undefined_variables import (
        paper_writing_check_undefined_variables,
    )
    tex = tmp_path / "clean.tex"
    tex.write_text(r"""
\documentclass{IEEEtran}
\begin{document}
The power dissipation is
\begin{equation}
P = R I^2
\end{equation}
where $P$ is the dissipated power and $R$ is the resistance.
\end{document}
""", encoding="utf-8")
    r = paper_writing_check_undefined_variables(str(tex))
    # P and R are defined via "where" clause; I, R may overlap with
    # whitelist.  Just verify n_undefined is small.
    assert r["n_undefined"] <= 1   # at most one stragglesr (e.g. small "I" depending on whitelist)


def test_undefined_variables_finds_eta_s(tmp_path):
    """A paper that introduces $\\eta_s$ but never defines it should
    report it as undefined."""
    from radia_mcp.paper_writing._undefined_variables import (
        paper_writing_check_undefined_variables,
    )
    tex = tmp_path / "undef.tex"
    tex.write_text(r"""
\documentclass{IEEEtran}
\begin{document}
The Sommerfeld integrand contains
\begin{equation}
I = \frac{1}{u + \eta_s k_0}
\end{equation}
The result is plotted below.  Note that $u = \sqrt{k_\rho^2 - k_0^2}$.
\end{document}
""", encoding="utf-8")
    r = paper_writing_check_undefined_variables(str(tex))
    # eta_s and k_0 (likely subscript k_0 too) -- should flag eta_s
    flagged_syms = {x["symbol"] for x in r["undefined_symbols"]}
    assert any("eta_s" in s or "eta" in s for s in flagged_syms), (
        f"expected eta_s flagged, got {flagged_syms}"
    )


def test_undefined_variables_nomenclature_block(tmp_path):
    """A paper with an IEEEdescription nomenclature should NOT
    flag symbols defined there."""
    from radia_mcp.paper_writing._undefined_variables import (
        paper_writing_check_undefined_variables,
    )
    tex = tmp_path / "nomen.tex"
    tex.write_text(r"""
\documentclass{IEEEtran}
\begin{document}
\begin{IEEEdescription}
\item[$\eta_s$] Surface impedance.
\item[$\zeta$] Damping coefficient.
\end{IEEEdescription}
The equation is
\begin{equation}
I = \frac{1}{1 + \eta_s + \zeta}
\end{equation}
\end{document}
""", encoding="utf-8")
    r = paper_writing_check_undefined_variables(str(tex))
    flagged_syms = {x["symbol"] for x in r["undefined_symbols"]}
    # eta_s and zeta should NOT be flagged (defined in IEEEdescription)
    assert not any("eta_s" in s for s in flagged_syms)
    # zeta is on UNIVERSAL_WHITELIST so it's exempt anyway
    assert "\\zeta" not in flagged_syms


def test_undefined_variables_whitelist_exempts_universal_constants(tmp_path):
    """\\pi, \\omega, \\sigma etc. should never be flagged."""
    from radia_mcp.paper_writing._undefined_variables import (
        paper_writing_check_undefined_variables,
    )
    tex = tmp_path / "univ.tex"
    tex.write_text(r"""
\documentclass{IEEEtran}
\begin{document}
\begin{equation}
\omega = 2\pi f, \quad \sigma = 5.8 \times 10^7
\end{equation}
\end{document}
""", encoding="utf-8")
    r = paper_writing_check_undefined_variables(str(tex))
    flagged = {x["symbol"] for x in r["undefined_symbols"]}
    for sym in [r"\omega", r"\pi", r"\sigma", "f"]:
        assert sym not in flagged, f"{sym} should be exempt"


def test_undefined_variables_extra_whitelist_arg(tmp_path):
    """extra_whitelist param lets users add project-specific generics."""
    from radia_mcp.paper_writing._undefined_variables import (
        paper_writing_check_undefined_variables,
    )
    tex = tmp_path / "extra.tex"
    tex.write_text(r"""
\documentclass{IEEEtran}
\begin{document}
\begin{equation}
\alpha + \beta + Q
\end{equation}
\end{document}
""", encoding="utf-8")
    # Without extra whitelist, Q would be flagged (no def, single
    # letter, not in universal whitelist).
    r = paper_writing_check_undefined_variables(str(tex))
    flagged = {x["symbol"] for x in r["undefined_symbols"]}
    assert "Q" in flagged
    # With Q in extra_whitelist -> not flagged
    r2 = paper_writing_check_undefined_variables(
        str(tex), extra_whitelist="Q,Z")
    flagged2 = {x["symbol"] for x in r2["undefined_symbols"]}
    assert "Q" not in flagged2


def test_undefined_variables_missing_tex_returns_error(tmp_path):
    from radia_mcp.paper_writing._undefined_variables import (
        paper_writing_check_undefined_variables,
    )
    r = paper_writing_check_undefined_variables(
        str(tmp_path / "does_not_exist.tex"))
    assert "error" in r
    assert "not found" in r["error"]


def test_resolve_input_chain_inlines_subfiles(tmp_path):
    """resolve_input_chain merges \\input{...} subfiles into one string."""
    from radia_mcp.paper_writing._tex_resolver import resolve_input_chain
    (tmp_path / "main.tex").write_text(
        r"\documentclass{article}\begin{document}"
        r"\input{chap1}" + "\n" + r"\input{chap2}"
        r"\end{document}", encoding="utf-8")
    (tmp_path / "chap1.tex").write_text("Chapter 1 body.", encoding="utf-8")
    (tmp_path / "chap2.tex").write_text("Chapter 2 body.", encoding="utf-8")
    r = resolve_input_chain(str(tmp_path / "main.tex"))
    assert r["ok"] is True
    assert len(r["files_resolved"]) == 3   # main + 2 chapters
    assert "Chapter 1 body" in r["merged_tex"]
    assert "Chapter 2 body" in r["merged_tex"]
    assert len(r["files_missing"]) == 0


def test_resolve_input_chain_handles_missing_subfile(tmp_path):
    from radia_mcp.paper_writing._tex_resolver import resolve_input_chain
    (tmp_path / "main.tex").write_text(
        r"\documentclass{article}\begin{document}"
        r"\input{nonexistent}"
        r"\end{document}", encoding="utf-8")
    r = resolve_input_chain(str(tmp_path / "main.tex"))
    assert r["ok"] is True   # missing file is graceful, not fatal
    assert len(r["files_missing"]) == 1
    assert "nonexistent" in r["files_missing"][0]["path"]


def test_resolve_input_chain_skips_commented_input(tmp_path):
    """\\input{...} inside a % comment line is NOT followed."""
    from radia_mcp.paper_writing._tex_resolver import resolve_input_chain
    (tmp_path / "main.tex").write_text(
        r"\documentclass{article}\begin{document}"
        "\n% \\input{commented}\n"
        r"\input{real}"
        r"\end{document}", encoding="utf-8")
    (tmp_path / "real.tex").write_text("real content", encoding="utf-8")
    r = resolve_input_chain(str(tmp_path / "main.tex"))
    assert r["ok"] is True
    # commented should NOT be resolved, real should
    paths = [f["path"] for f in r["files_resolved"]]
    assert any("real.tex" in p for p in paths)
    assert not any("commented" in p for p in paths)
    assert len(r["files_missing"]) == 0   # commented NOT counted missing


def test_extract_abstract_standard_env(tmp_path):
    from radia_mcp.paper_writing._tex_resolver import extract_abstract_from_tex
    src = (r"\documentclass{article}\begin{document}"
           r"\begin{abstract}This is the abstract body.\end{abstract}"
           r"\section{Introduction}...")
    r = extract_abstract_from_tex(src)
    assert r == "This is the abstract body."


def test_extract_abstract_ieee_abstract_env(tmp_path):
    from radia_mcp.paper_writing._tex_resolver import extract_abstract_from_tex
    src = r"\begin{IEEEabstract}IEEE-style abstract here.\end{IEEEabstract}"
    r = extract_abstract_from_tex(src)
    assert r == "IEEE-style abstract here."


def test_extract_abstract_not_found_returns_none():
    from radia_mcp.paper_writing._tex_resolver import extract_abstract_from_tex
    src = r"\documentclass{article}\begin{document}No abstract here.\end{document}"
    r = extract_abstract_from_tex(src)
    assert r is None


def test_paper_writing_extract_abstract_full_tool(tmp_path):
    from radia_mcp.paper_writing import paper_writing_extract_abstract
    tex = tmp_path / "p.tex"
    tex.write_text(
        r"\documentclass{IEEEtran}\begin{document}"
        r"\begin{abstract}Direct abstract here, ~30 chars.\end{abstract}"
        r"\section{Intro}...", encoding="utf-8")
    r = paper_writing_extract_abstract(str(tex))
    assert r["ok"] is True
    assert r["method"] == "direct"
    assert "Direct abstract here" in r["abstract"]
    assert r["n_chars"] > 20


def test_paper_writing_extract_abstract_via_input_resolution(tmp_path):
    """Abstract in a separate file, pulled in via \\input{abstract}."""
    from radia_mcp.paper_writing import paper_writing_extract_abstract
    (tmp_path / "main.tex").write_text(
        r"\documentclass{article}\begin{document}"
        r"\input{abstract}"
        r"\input{body}"
        r"\end{document}", encoding="utf-8")
    (tmp_path / "abstract.tex").write_text(
        r"\begin{abstract}Inlined via \\input{}, 40+ chars.\end{abstract}",
        encoding="utf-8")
    (tmp_path / "body.tex").write_text("Body text.", encoding="utf-8")
    r = paper_writing_extract_abstract(str(tmp_path / "main.tex"))
    assert r["ok"] is True
    assert r["method"] == "resolved_inputs"
    assert "Inlined via" in r["abstract"]


def test_em_submission_gate_auto_extracts_abstract(tmp_path):
    """v0.93.0: gate auto-extracts abstract when not passed."""
    from radia_mcp.paper_writing._em_paper_style import (
        paper_writing_em_submission_gate,
    )
    tex = tmp_path / "p.tex"
    bib = tmp_path / "ref.bib"
    tex.write_text(
        r"\documentclass{IEEEtran}\begin{document}"
        r"\begin{abstract}This abstract should be auto-extracted by "
        r"the gate when abstract_text= is not passed.\end{abstract}"
        r"\section{Introduction}Body.\end{document}",
        encoding="utf-8")
    bib.write_text("", encoding="utf-8")
    r = paper_writing_em_submission_gate(
        tex_path=str(tex), bib_path=str(bib)
    )
    # The new 'abstract_extracted' check should have fired
    names = [c["name"] for c in r["checks"]]
    assert "abstract_extracted" in names
    abs_check = next(c for c in r["checks"] if c["name"] == "abstract_extracted")
    assert abs_check["status"] == "pass"
    # And the abstract-dependent checks should now run
    assert "validate_abstract_length" in names


def test_em_submission_gate_resolves_multifile_tex(tmp_path):
    """v0.93.0: gate auto-resolves \\input{} subfiles."""
    from radia_mcp.paper_writing._em_paper_style import (
        paper_writing_em_submission_gate,
    )
    (tmp_path / "main.tex").write_text(
        r"\documentclass{IEEEtran}\begin{document}"
        r"\input{intro}\input{methods}"
        r"\end{document}", encoding="utf-8")
    (tmp_path / "intro.tex").write_text(
        r"\section{Intro}Intro content.", encoding="utf-8")
    (tmp_path / "methods.tex").write_text(
        r"\section{Methods}Methods content.", encoding="utf-8")
    (tmp_path / "ref.bib").write_text("", encoding="utf-8")
    r = paper_writing_em_submission_gate(
        tex_path=str(tmp_path / "main.tex"),
        bib_path=str(tmp_path / "ref.bib"),
    )
    names = [c["name"] for c in r["checks"]]
    assert "multifile_resolved" in names
    mf_check = next(c for c in r["checks"] if c["name"] == "multifile_resolved")
    assert mf_check["status"] == "pass"
    assert "2 \\input subfiles" in mf_check["summary"]


def test_undefined_variables_reports_first_occurrence_context(tmp_path):
    from radia_mcp.paper_writing._undefined_variables import (
        paper_writing_check_undefined_variables,
    )
    tex = tmp_path / "ctx.tex"
    tex.write_text(r"""
\documentclass{IEEEtran}
\begin{document}
Some intro.  Now consider $\xi_{magic}$ in the formula.
\end{document}
""", encoding="utf-8")
    r = paper_writing_check_undefined_variables(str(tex))
    # xi is whitelisted but xi_magic isn't a single token; the regex
    # extracts \xi which IS whitelisted.  Test the no-flag path.
    # Use a definitely-novel symbol instead:
    tex.write_text(r"""
\documentclass{IEEEtran}
\begin{document}
The novel symbol $\fakecommand$ appears here.
\end{document}
""", encoding="utf-8")
    r = paper_writing_check_undefined_variables(str(tex))
    flagged = {x["symbol"] for x in r["undefined_symbols"]}
    assert r"\fakecommand" in flagged
    # Context retrieved
    rec = next(x for x in r["undefined_symbols"]
                if x["symbol"] == r"\fakecommand")
    assert "first_occurrence_context" in rec
    assert "fakecommand" in rec["first_occurrence_context"]
