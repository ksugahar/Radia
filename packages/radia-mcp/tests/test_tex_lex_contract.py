"""Shared static TeX lexical boundaries and input-chain failure contracts."""
import pytest

from radia_mcp.paper_writing._tex_lex import mask_tex_noncode
from radia_mcp.paper_writing._tex_resolver import resolve_input_chain
from radia_mcp.bibliography.plans.T14_canonical import _keys_in_order


@pytest.mark.parametrize("literal", [
    r"\verb|\input{missing}\cite{fake}|",
    r"\verb*+\input{missing}\cite{fake}+",
    r"\begin{verbatim}\input{missing}\cite{fake}\end{verbatim}",
    r"\begin{minted}{tex}\input{missing}\cite{fake}\end{minted}",
    r"\begin{lstlisting}\input{missing}\cite{fake}\end{lstlisting}",
    "% \\input{missing} \\cite{fake}\n",
    "\\\\% \\input{missing} \\cite{fake}\n",
    r"\\input{missing} \\cite{fake}",
])
def test_literal_inputs_and_citations_do_not_reach_resolver(tmp_path, literal):
    tex = tmp_path / "paper.tex"
    source = literal + "\n" + r"\input{body}"
    tex.write_text(source, encoding="utf-8")
    (tmp_path / "body.tex").write_text(r"\cite{real}", encoding="utf-8")
    masked = mask_tex_noncode(source)
    assert len(masked) == len(source)
    assert [i for i, c in enumerate(masked) if c == "\n"] == [i for i, c in enumerate(source) if c == "\n"]
    result = resolve_input_chain(str(tex))
    assert result["ok"], result
    assert len(result["files_resolved"]) == 2
    assert _keys_in_order(result["merged_tex"]) == ["real"]


def test_commented_literal_start_does_not_hide_active_citations():
    assert _keys_in_order("% \\begin{verbatim}\n\\cite{real}") == ["real"]


def test_escaped_percent_does_not_hide_real_input(tmp_path):
    tex = tmp_path / "paper.tex"
    tex.write_text(r"50\% \input{body}", encoding="utf-8")
    (tmp_path / "body.tex").write_text("body", encoding="utf-8")
    result = resolve_input_chain(str(tex))
    assert result["ok"] and "body" in result["merged_tex"]


@pytest.mark.parametrize("body", [r"\input{paper}", r"\input{sub/../paper.tex}"])
def test_cyclic_inputs_fail_instead_of_succeeding_as_duplicate(tmp_path, body):
    (tmp_path / "sub").mkdir()
    tex = tmp_path / "paper.tex"
    tex.write_text(body, encoding="utf-8")
    result = resolve_input_chain(str(tex))
    assert result["ok"] is False
    assert "cyclic" in result["error"]


def test_repeated_noncyclic_input_keeps_deduplicated_static_contract(tmp_path):
    tex = tmp_path / "paper.tex"
    tex.write_text(r"\input{body}\input{body}", encoding="utf-8")
    (tmp_path / "body.tex").write_text("unique body", encoding="utf-8")
    result = resolve_input_chain(str(tex))
    assert result["ok"]
    assert len(result["files_duplicate"]) == 1
    assert result["merged_tex"].count("unique body") == 1


@pytest.mark.parametrize("source", [r"\begin{verbatim}\input{missing}", r"\verb|not closed"])
def test_unterminated_literal_is_explicit_error(tmp_path, source):
    tex = tmp_path / "paper.tex"
    tex.write_text(source, encoding="utf-8")
    result = resolve_input_chain(str(tex))
    assert not result["ok"]
    assert "unterminated" in result["error"]
