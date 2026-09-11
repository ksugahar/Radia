from __future__ import annotations

import json
import shutil

import pytest
from radia_mcp.bibliography.plans.T14_canonical import (
    _keys_in_order,
    bibliography_canonical_path,
    bibliography_get_entries,
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


def test_get_entries_returns_canonical_records_in_requested_order():
    result = json.loads(bibliography_get_entries("YuferevIda2009, Kameari2018"))

    assert result["ok"] is True
    assert result["canonical_path"].endswith("references.bib")
    assert len(result["canonical_sha256"]) == 64
    assert [entry["key"] for entry in result["entries"]] == [
        "YuferevIda2009",
        "Kameari2018",
    ]
    assert result["entries"][0]["fields"]["year"] == "2009"


def test_get_entries_fails_closed_for_missing_or_duplicate_keys():
    missing = json.loads(bibliography_get_entries("not-a-key"))
    duplicate = json.loads(bibliography_get_entries("Kameari2018, Kameari2018"))

    assert missing == {
        "ok": False,
        "error": "citation keys are absent from canonical references.bib",
        "missing": ["not-a-key"],
    }
    assert duplicate == {
        "ok": False,
        "error": "duplicate citation keys",
        "keys": ["Kameari2018"],
    }


def test_igte_cauer_sibc_reference_set_is_canonical():
    keys = (
        "Kameari2018, Kuriyama2019, kuriyama2021multiport, Matsuo2026jmmm, "
        "Senior1962, Mitzner1967, YuferevIda2009, QuarteroniValli1999, "
        "gautschi2004orthogonal, HirumaXFEM2023, vandyke1975perturbation, "
        "gallivan1996rational, oh1995efficient, deeley1990surface, "
        "JingguoLavers1993, warne1994eddy, yuferev2001surface, "
        "proekt2002overlapping, dauge2014corner, HirumaXFEMexact2024"
    )
    result = json.loads(bibliography_get_entries(keys))

    assert result["ok"] is True
    assert len(result["entries"]) == 20
    senior = next(entry for entry in result["entries"] if entry["key"] == "Senior1962")
    assert senior["fields"]["year"] == "1960"
    assert senior["fields"]["doi"] == "10.1007/BF02920074"


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


@pytest.mark.parametrize("text,expected", [
    ("% \\cite{ignored}\n\\cite{real}", ["real"]),
    (r"\citet*[see][p.~2]{beta,alpha} \parencite{gamma}", ["beta", "alpha", "gamma"]),
    (r"\verb|\cite{literal}| \cite{real}", ["real"]),
    (r"\begin{verbatim}\cite{literal}\end{verbatim}\cite{real}", ["real"]),
    (r"\begin{lstlisting}\cite{literal}\end{lstlisting}\cite{real}", ["real"]),
    (r"50\% \cite{real}", ["real"]),
    ("\\\\% \\cite{ignored}\n\\cite{real}", ["real"]),
    (r"\\cite{literal} \cite{real}", ["real"]),
])
def test_citation_scanner_ignores_comments_and_literal_code(text, expected):
    assert _keys_in_order(text) == expected


@pytest.mark.parametrize("text", [r"\parencites{one}{two}", r"\InputIfFileExists{body}{}{}",
                                  r"\includeonly{body}", r"\cite{*}"])
def test_unsupported_static_citation_syntax_fails_explicitly(text):
    with pytest.raises(ValueError):
        _keys_in_order(text)


@pytest.mark.skipif(shutil.which("bibtex") is None, reason="BibTeX is unavailable")
def test_make_bbl_nocite_all_uses_canonical_fixture_and_ignores_comments(tmp_path, monkeypatch):
    from radia_mcp.bibliography.plans import T14_canonical as canonical
    bib = tmp_path / "fixture.bib"
    bib.write_text("@misc{one,title={First},author={Doe, Jane},year=2024}\n"
                   "@misc{two,title={Second},author={Roe, John},year=2025}", encoding="utf-8")
    monkeypatch.setattr(canonical, "CANONICAL", bib)
    tex = tmp_path / "paper.tex"
    tex.write_text("% \\cite{not_real}\n\\nocite{*}\\bibliographystyle{plain}", encoding="utf-8")
    result = canonical.bibliography_make_bbl(str(tex))
    assert result.startswith("bibliography_make_bbl:"), result
    bbl = tex.with_suffix(".bbl").read_text(encoding="utf-8")
    assert r"\bibitem{one}" in bbl and r"\bibitem{two}" in bbl
