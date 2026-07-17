"""MCP Server: radia_mcp.paper_writing

Journal paper writing skill suite -- IMRaD lint, abstract / citation
/ figure / equation checks, image-based PDF layout verification
(v0.88.0), LaTeX figure placement knowledge (v0.88.0), and
PDF download (IEEE / ScienceDirect / Emerald with cookies).

Promoted 2026-05-26 from mcp-server-document.paper_writing
(LAB-private) to radia-mcp (public PyPI).

2026-07-17: the presentation server was MERGED into this one (Sugahara:
slide decks cannot yet be authored end-to-end by AI, so the slide lint /
PPTX toolset does not warrant a standalone server).  All presentation_*
tools are served here; the radia_mcp.presentation module remains the
implementation home.

Usage:
    mcp-server-paper-writing              # stdio
    mcp-server-paper-writing --selftest   # self-test
"""

import sys

from mcp.server.fastmcp import FastMCP

from ..common import register_status_tool
from ..presentation import register as _register_presentation

from . import tools as _tools
from ._pdf_layout_visual import (
    paper_writing_render_pages_to_png,
    paper_writing_detect_page_whitespace_anomalies,
    paper_writing_layout_thumbnail_strip,
    paper_writing_check_floats_far_from_reference,
    paper_writing_layout_visual_recipe,
)
from ._tex_figure_placement import paper_writing_tex_figure_placement
from ._arxiv_source import (
    paper_writing_arxiv_fetch_latex_source,
    paper_writing_arxiv_extract_equations,
    paper_writing_arxiv_search,
    paper_writing_semantic_scholar_lookup,
    paper_writing_semantic_scholar_references,
    paper_writing_semantic_scholar_citations,
    paper_writing_external_sources_recipe,
)
from ._em_paper_style import (
    paper_writing_em_paper_style,
    paper_writing_em_submission_gate,
)
from ._citation_verify import (
    paper_writing_verify_citation,
    paper_writing_citation_workflow_recipe,
    paper_writing_check_citation_keys_exist,
)
from ._undefined_acronyms import (
    paper_writing_check_undefined_acronyms,
)
from ._digest_lints import (
    paper_writing_check_ref_label_consistency,
    paper_writing_check_ieee_keywords,
    paper_writing_check_pdf_unresolved_markers,
)
from ._pdf_overlap_detection import (
    paper_writing_detect_text_image_overlap,
    paper_writing_detect_text_overflow_page,
    paper_writing_detect_overlapping_text_blocks,
    paper_writing_pdf_overlap_recipe,
)
from ._undefined_variables import (
    paper_writing_check_undefined_variables,
)
from ._tex_resolver import (
    paper_writing_resolve_input_chain,
    paper_writing_extract_abstract,
)


mcp = FastMCP("mcp-server-paper-writing")


# ============================================================
# Auto-register the 67 paper_writing_* tools from tools.py
# (matches the registration pattern previously used in
# mcp-server-document.paper_writing.__init__.register).
# ============================================================
_REGISTERED_TOOLS = []
for _name in sorted(dir(_tools)):
    if _name.startswith("paper_writing_"):
        _fn = getattr(_tools, _name)
        if callable(_fn):
            mcp.tool()(_fn)
            _REGISTERED_TOOLS.append(_name)


# ============================================================
# v0.88.0: image-based PDF layout verification (5 new tools)
# ============================================================
mcp.tool()(paper_writing_render_pages_to_png)
mcp.tool()(paper_writing_detect_page_whitespace_anomalies)
mcp.tool()(paper_writing_layout_thumbnail_strip)
mcp.tool()(paper_writing_check_floats_far_from_reference)
mcp.tool()(paper_writing_layout_visual_recipe)

# ============================================================
# v0.88.0: LaTeX figure placement knowledge
# ============================================================
mcp.tool()(paper_writing_tex_figure_placement)


# ============================================================
# v0.89.0: external paper sources (arXiv LaTeX + Semantic Scholar)
# GitHub survey 2026-05-26 absorption
# ============================================================
mcp.tool()(paper_writing_arxiv_fetch_latex_source)
mcp.tool()(paper_writing_arxiv_extract_equations)
mcp.tool()(paper_writing_arxiv_search)
mcp.tool()(paper_writing_semantic_scholar_lookup)
mcp.tool()(paper_writing_semantic_scholar_references)
mcp.tool()(paper_writing_semantic_scholar_citations)
mcp.tool()(paper_writing_external_sources_recipe)


# ============================================================
# v0.90.0: EM-paper-specific style knowledge + pre-submission gate
# ============================================================
mcp.tool()(paper_writing_em_paper_style)
mcp.tool()(paper_writing_em_submission_gate)


# ============================================================
# v0.91.0: citation verification workflow (reference.bib enforcement)
# ============================================================
mcp.tool()(paper_writing_verify_citation)
mcp.tool()(paper_writing_citation_workflow_recipe)


# ============================================================
# v0.92.0: PDF overlap/overflow + undefined-variable check
# ============================================================
mcp.tool()(paper_writing_detect_text_image_overlap)
mcp.tool()(paper_writing_detect_text_overflow_page)
mcp.tool()(paper_writing_detect_overlapping_text_blocks)
mcp.tool()(paper_writing_pdf_overlap_recipe)
mcp.tool()(paper_writing_check_undefined_variables)


# ============================================================
# v0.93.0: multi-file .tex resolver + abstract auto-extract
# ============================================================
mcp.tool()(paper_writing_resolve_input_chain)
mcp.tool()(paper_writing_extract_abstract)


# ============================================================
# 2026-07-17: presentation merged into paper-writing (Sugahara).
# AI cannot yet author slide decks end-to-end, so the slide lint /
# PPTX helper toolset rides this server instead of a standalone
# mcp-server-presentation (entry point retired).  Humans author the
# slides; these tools lint / extract / budget them.
# ============================================================
_N_PRESENTATION_TOOLS = _register_presentation(mcp)


# ============================================================
# Prompt: cite-with-verification reminder
# ============================================================

@mcp.prompt()
def cite_a_claim(claim: str, bib_path: str = "reference.bib") -> str:
    """Surface the lab POLICY for inserting any \\cite{} into a paper.

    Args:
        claim: the sentence / paragraph that needs a citation.
        bib_path: path to the user's reference.bib (default
            'reference.bib' in CWD).
    """
    return (
        f"Insert a citation for this claim:\n"
        f"  > {claim}\n\n"
        f"Lab POLICY (MANDATORY):\n"
        f"  1. Read the existing reference.bib at {bib_path!r}\n"
        f"     via paper_writing_lint_reference_format(bib_path).\n"
        f"  2. NEVER invent a DOI / author / year.  Verify via\n"
        f"     paper_writing_verify_citation(claim, bib_path, ...).\n"
        f"  3. If the verify tool returns verdict='found_in_bib',\n"
        f"     reuse the existing matching_key.\n"
        f"  4. If verdict='ready_to_insert', append the\n"
        f"     suggested_bibtex to reference.bib FIRST, then\n"
        f"     write \\cite{{citation_key}}.\n"
        f"  5. If verdict='needs_disambiguation', show the user the\n"
        f"     candidates and ask which one applies.\n"
        f"  6. If verdict='no_candidate_found', do NOT fabricate.\n"
        f"     Mark as \\cite{{TODO: verify -- ...}} and ask the user.\n\n"
        f"See paper_writing_citation_workflow_recipe() for the full\n"
        f"behavioral rule.  Failure to verify is the #1 reason\n"
        f"AI-assisted papers get rejected by reviewers."
    )


register_status_tool(
    mcp,
    server_name='mcp-server-paper-writing',
    description=(
        'Journal-paper writing skill suite: IMRaD/abstract/citation/'
        'figure lint (67 tools), image-based PDF layout verification '
        '(pymupdf), LaTeX figure placement knowledge (htbp/placeins/'
        'widths/anti-patterns), IEEE/ScienceDirect/Emerald PDF '
        'download with cookies. Also serves the merged presentation_* '
        'slide lint + PPTX toolset (2026-07-17: standalone '
        'presentation server retired).'
    ),
    subpackage='radia_mcp.paper_writing',
    related_servers=["literature-index", "graph", "chart2d", "figure",
                       "poster"],
    optional_deps=["pymupdf", "Pillow", "requests", "python-pptx"],
)


def main():
    """Entry point for mcp-server-paper-writing."""
    if "--selftest" in sys.argv:
        print("paper_writing MCP server self-test:")
        print(f"  registered paper_writing_* tools: {len(_REGISTERED_TOOLS)}")
        print(f"  registered presentation_* tools (merged): {_N_PRESENTATION_TOOLS}")
        assert _N_PRESENTATION_TOOLS > 60, (
            f"presentation merge lost tools: {_N_PRESENTATION_TOOLS}")
        # Smoke-test 3 representative tools
        r1 = _tools.paper_writing_check_kanji_ratio(
            'これは日本語の文章で、漢字の比率を計算します。'
        )
        assert 'kanji_ratio' in r1, f"kanji_ratio missing key: {r1}"
        print(f"    check_kanji_ratio: kanji_ratio={r1['kanji_ratio']:.3f}")

        r2 = _tools.paper_writing_count_weak_expressions(
            'It is possible that we may consider somewhat improving the result.'
        )
        assert 'total_weak_expressions' in r2
        print(f"    count_weak_expressions: total={r2['total_weak_expressions']}")

        r3 = _tools.paper_writing_analyze_sentences(
            'Short sentence. A medium one with several words. ' * 3
        )
        assert 'total_sentences' in r3
        print(f"    analyze_sentences: n_sentences={r3['total_sentences']}")

        # v0.88.0 layout knowledge tools
        full = paper_writing_tex_figure_placement('all')
        ov = paper_writing_tex_figure_placement('overview')
        unk = paper_writing_tex_figure_placement('does_not_exist')
        assert len(full) > 8000, f"tex_figure_placement too small: {len(full)}"
        assert len(ov) > 500, f"tex_figure_placement overview too small: {len(ov)}"
        assert "Unknown topic" in unk
        print(f"  tex_figure_placement: all={len(full)} chars, "
              f"overview={len(ov)} chars")

        rec = paper_writing_layout_visual_recipe()
        assert len(rec) > 1500, f"layout_visual_recipe too small: {len(rec)}"
        assert "thumbnail_strip" in rec
        print(f"  layout_visual_recipe: {len(rec)} chars")

        # v0.89.0 external sources tools
        rec2 = paper_writing_external_sources_recipe()
        assert len(rec2) > 1500, f"external_sources_recipe too small: {len(rec2)}"
        assert "arxiv-latex-mcp" in rec2 and "Semantic Scholar" in rec2
        print(f"  external_sources_recipe: {len(rec2)} chars")
        # Light smoke -- arxiv_id normalizer (offline, no network):
        from ._arxiv_source import _normalize_arxiv_id
        assert _normalize_arxiv_id("arXiv:2603.17339v3") == "2603.17339"
        assert _normalize_arxiv_id(
            "https://arxiv.org/abs/2603.17339"
        ) == "2603.17339"
        # Extract equations regex (offline)
        eqs = paper_writing_arxiv_extract_equations(
            r"Text. \begin{equation} E = mc^2 \end{equation} more. $$F = ma$$"
        )
        assert eqs["n_equations"] == 2, f"expected 2 eqs, got {eqs}"
        print(f"  arxiv_extract_equations: {eqs['n_equations']} eqs found")

        # v0.90.0 EM paper style + submission gate
        em_full = paper_writing_em_paper_style("all")
        em_ov = paper_writing_em_paper_style("overview")
        assert len(em_full) > 8000, f"em_paper_style all too small: {len(em_full)}"
        assert len(em_ov) > 500
        for kw in ["exp(+j", "magnetic flux density", "siunitx"]:
            assert kw in em_full, f"em_paper_style missing: {kw}"
        print(f"  em_paper_style: all={len(em_full)} chars, "
              f"overview={len(em_ov)} chars")
        gate = paper_writing_em_submission_gate()
        # v0.91.0: no bib_path => bib_policy gate fires => verdict=fail
        assert gate["verdict"] == "fail"
        names = [c["name"] for c in gate["checks"]]
        assert "bib_policy" in names
        print(f"  em_submission_gate (no inputs): verdict={gate['verdict']} "
              f"(bib_policy gate fired as designed), "
              f"{gate['n_checks_run']} checks")

        # v0.91.0: citation verification workflow
        rec3 = paper_writing_citation_workflow_recipe()
        assert len(rec3) > 2000, f"citation recipe too small: {len(rec3)}"
        for kw in ["NEVER invent", "reference.bib", "Crossref"]:
            assert kw in rec3, f"citation recipe missing: {kw}"
        print(f"  citation_workflow_recipe: {len(rec3)} chars")
        # Offline check: missing bib_path -> error verdict
        from ._citation_verify import paper_writing_verify_citation
        v = paper_writing_verify_citation(claim="x", bib_path="")
        assert v["verdict"] == "error"
        print(f"  verify_citation (no bib): verdict={v['verdict']} "
              f"(refuses to fabricate)")

        # v0.92.0: PDF overlap recipe + IoU math
        rec4 = paper_writing_pdf_overlap_recipe()
        assert len(rec4) > 1000
        for kw in ["IoU", "PaperOrchestra", "VILA"]:
            assert kw in rec4, f"overlap recipe missing: {kw}"
        from ._pdf_overlap_detection import _iou
        assert _iou((0, 0, 10, 10), (0, 0, 10, 10)) == 1.0
        assert _iou((0, 0, 10, 10), (20, 20, 30, 30)) == 0.0
        print(f"  pdf_overlap_recipe: {len(rec4)} chars, IoU math OK")

        # v0.92.0: undefined-variable check (synthetic .tex in memory
        # via a tmp NamedTemporaryFile; skip if can't write)
        import tempfile
        tf = tempfile.NamedTemporaryFile(
            mode="w", suffix=".tex", delete=False, encoding="utf-8")
        try:
            tf.write(r"""
\documentclass{IEEEtran}
\begin{document}
\begin{equation}
I = \frac{1}{1 + \fakecmd}
\end{equation}
\end{document}
""")
            tf.close()
            udv = paper_writing_check_undefined_variables(tf.name)
            assert udv["n_undefined"] >= 1
            assert any(r"\fakecmd" in s["symbol"]
                        for s in udv["undefined_symbols"])
            print(f"  check_undefined_variables: "
                  f"{udv['n_undefined']} undefined symbols "
                  f"(fakecmd flagged as designed)")
        finally:
            import os as _os
            _os.unlink(tf.name)

        # v0.93.0: multi-file resolver + abstract auto-extract
        from ._tex_resolver import (
            resolve_input_chain,
            extract_abstract_from_tex,
            paper_writing_resolve_input_chain,
            paper_writing_extract_abstract,
        )
        # abstract extractor offline check
        src = (r"\documentclass{IEEEtran}"
               r"\begin{document}"
               r"\begin{abstract}"
               r"This paper studies the impedance of an axisymmetric inductor."
               r"\end{abstract}"
               r"\end{document}")
        ab = extract_abstract_from_tex(src)
        assert ab is not None and "impedance" in ab, \
            f"extract_abstract failed: {ab}"
        print(f"  extract_abstract_from_tex: {len(ab)} chars extracted")

        # resolve_input_chain on a 2-file synthetic project
        import tempfile as _tf, os as _os2
        tmpdir = _tf.mkdtemp(prefix="radia_mcp_selftest_")
        try:
            main_tex = _os2.path.join(tmpdir, "main.tex")
            ch1_tex = _os2.path.join(tmpdir, "ch1.tex")
            with open(main_tex, "w", encoding="utf-8") as fh:
                fh.write(r"\documentclass{IEEEtran}"
                         "\n\\input{ch1}\n")
            with open(ch1_tex, "w", encoding="utf-8") as fh:
                fh.write(r"\section{Intro} subfile content.")
            r = paper_writing_resolve_input_chain(
                main_tex, return_merged_text=True)
            assert r["ok"] and len(r["files_resolved"]) == 2
            assert "subfile content" in r["merged_tex"]
            print(f"  resolve_input_chain: "
                  f"{len(r['files_resolved'])} files merged")

            # paper_writing_extract_abstract end-to-end via input
            main2 = _os2.path.join(tmpdir, "paper.tex")
            abs2 = _os2.path.join(tmpdir, "abstract.tex")
            with open(main2, "w", encoding="utf-8") as fh:
                fh.write(r"\documentclass{IEEEtran}"
                         "\n\\input{abstract}\n")
            with open(abs2, "w", encoding="utf-8") as fh:
                fh.write(r"\begin{abstract}AC eddy current."
                         r"\end{abstract}")
            r2 = paper_writing_extract_abstract(main2)
            assert r2["ok"] and r2["method"] == "resolved_inputs"
            assert "AC eddy current" in r2["abstract"]
            print(f"  extract_abstract (via \\input): "
                  f"{r2['n_chars']} chars, method={r2['method']}")
        finally:
            import shutil as _sh
            _sh.rmtree(tmpdir, ignore_errors=True)

        print("OK")
        return
    mcp.run()


if __name__ == "__main__":
    main()
