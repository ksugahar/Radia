from __future__ import annotations

import shutil

import pytest

from radia_mcp.bibliography.plans.T14_canonical import (
    _keys_in_order,
    bibliography_canonical_path,
    bibliography_make_bbl,
)


def test_canonical_path_describes_single_source_and_bbl_delivery():
    result = bibliography_canonical_path()
    assert "references.bib" in result
    assert "single source of truth" in result
    assert "generated .bbl" in result


def test_citation_keys_preserve_first_appearance_order():
    tex = r"\citep[see][p.~2]{beta, alpha} \citet{alpha} \nocite{gamma}"
    assert _keys_in_order(tex) == ["beta", "alpha", "gamma"]


def test_make_bbl_rejects_unknown_key_without_partial_output(tmp_path, monkeypatch):
    monkeypatch.setattr(shutil, "which", lambda _name: None)
    tex = tmp_path / "paper.tex"
    tex.write_text(r"\cite{not_in_the_canonical_bibliography}", encoding="utf-8")
    out = tmp_path / "paper.bbl"
    out.write_text("previous verified bibliography", encoding="utf-8")

    result = bibliography_make_bbl(str(tex), style="plain", out_path=str(out))

    assert result.startswith("Error:")
    assert "absent from canonical references.bib" in result
    assert out.read_text(encoding="utf-8") == "previous verified bibliography"


@pytest.mark.skipif(shutil.which("bibtex") is None, reason="BibTeX is unavailable")
def test_make_bbl_exports_only_the_generated_bbl(tmp_path):
    tex = tmp_path / "paper.tex"
    tex.write_text(
        r"\documentclass{article}\begin{document}\cite{abe2017passive}"
        r"\bibliographystyle{plain}\bibliography{references}\end{document}",
        encoding="utf-8",
    )

    result = bibliography_make_bbl(str(tex))

    bbl = tmp_path / "paper.bbl"
    assert result.startswith("bibliography_make_bbl:")
    assert bbl.is_file()
    assert r"\bibitem{abe2017passive}" in bbl.read_text(encoding="utf-8")
    assert not (tmp_path / "references.bib").exists()


@pytest.mark.skipif(shutil.which("bibtex") is None, reason="BibTeX is unavailable")
def test_make_bbl_collects_citations_from_input_files(tmp_path):
    (tmp_path / "body.tex").write_text(r"\cite{abe2017passive}", encoding="utf-8")
    tex = tmp_path / "paper.tex"
    tex.write_text(
        r"\documentclass{article}\begin{document}\input{body}"
        r"\bibliographystyle{plain}\bibliography{references}\end{document}",
        encoding="utf-8",
    )

    result = bibliography_make_bbl(str(tex))

    assert result.startswith("bibliography_make_bbl:")
    assert r"\bibitem{abe2017passive}" in (tmp_path / "paper.bbl").read_text(
        encoding="utf-8"
    )
