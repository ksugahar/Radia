"""Cross-check and bbl paths must use the same complete citation source."""
import pytest
from radia_mcp.bibliography.plans.T8_cite_validation import bibliography_cite_validation


def test_nested_input_uses_compile_root_and_explicit_pgf_suffix(tmp_path):
    (tmp_path / "sections").mkdir()
    (tmp_path / "sections" / "body.tex").write_text(r"\input{plot.pgf}", encoding="utf-8")
    (tmp_path / "plot.pgf").write_text(r"\citep*[see][2]{one}", encoding="utf-8")
    tex = tmp_path / "paper.tex"
    tex.write_text("% \\cite{not_real}\n\\input{sections/body}", encoding="utf-8")
    bib = tmp_path / "fixture.bib"
    bib.write_text("@misc{one,title={Test}}", encoding="utf-8")
    assert "PASS" in bibliography_cite_validation(str(tex), str(bib))


@pytest.mark.parametrize("mode", ["missing", "encoding", "directory"])
def test_failed_input_is_not_ignored(tmp_path, mode):
    tex = tmp_path / "paper.tex"
    tex.write_text(r"\cite{one}\input{body}", encoding="utf-8")
    child = tmp_path / "body.tex"
    if mode == "encoding":
        child.write_bytes(b"\x81")
    elif mode == "directory":
        child.mkdir()
    bib = tmp_path / "fixture.bib"
    bib.write_text("@misc{one,title={Test}}", encoding="utf-8")
    result = bibliography_cite_validation(str(tex), str(bib))
    assert result.startswith("Error:") and "PASS" not in result


def test_nocite_all_does_not_report_star_as_missing(tmp_path):
    tex = tmp_path / "paper.tex"
    tex.write_text(r"\nocite{*}", encoding="utf-8")
    bib = tmp_path / "fixture.bib"
    bib.write_text("@misc{one,title={Test}}\n@misc{two,title={Second}}", encoding="utf-8")
    assert "PASS" in bibliography_cite_validation(str(tex), str(bib))


def test_duplicate_bib_keys_are_ambiguous(tmp_path):
    tex = tmp_path / "paper.tex"
    tex.write_text(r"\cite{one}", encoding="utf-8")
    bib = tmp_path / "fixture.bib"
    bib.write_text("@misc{one,title={Test}}\n@misc{one,title={Second}}", encoding="utf-8")
    assert bibliography_cite_validation(str(tex), str(bib)).startswith("Error:")
