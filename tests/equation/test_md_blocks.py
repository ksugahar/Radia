"""Block structure: what a viewer lays out and a writer styles.

`MarkdownDoc` says which spans are maths; this says what the file is made of.
It lives in C++ because the viewer cannot afford Python -- importing the package
costs about 1.4 s against 4 ms for the module itself -- and because keeping one
implementation is the point: the Word writer calls this same scanner rather than
carrying a copy that would drift from it.

The guarantee that makes it safe to build an editor on is the same one
`MarkdownDoc` gives: concatenating every block's source rebuilds the file
exactly.
"""

from __future__ import annotations

import pytest

equation = pytest.importorskip("radia.equation")
md_blocks = equation.md_blocks
Kind = equation.MdBlock


def kinds(markdown):
    return [b.kind for b in md_blocks(markdown)]


def test_a_heading_carries_its_level_and_loses_its_hashes():
    (b,) = md_blocks("### Title")
    assert b.kind == Kind.Heading
    assert b.level == 3
    assert b.text == "Title"


def test_hashes_without_a_space_are_not_a_heading():
    (b,) = md_blocks("#hashtag")
    assert b.kind == Kind.Paragraph


@pytest.mark.parametrize("marker", ["-", "*", "+"])
def test_bullet_markers(marker):
    (b,) = md_blocks(f"{marker} item")
    assert b.kind == Kind.Bullet
    assert b.text == "item"


@pytest.mark.parametrize("src", ["1. item", "2) item", "  10. item"])
def test_numbered_items(src):
    (b,) = md_blocks(src)
    assert b.kind == Kind.Numbered
    assert b.text == "item"


def test_a_paragraph_runs_until_something_else_starts():
    blocks = md_blocks("one\ntwo\n# heading\n")
    assert [b.kind for b in blocks] == [Kind.Paragraph, Kind.Heading]
    assert blocks[0].text == "one\ntwo"


def test_a_code_fence_loses_its_markers_and_keeps_its_language():
    (b,) = md_blocks("```python\nx = 1\n```\n")
    assert b.kind == Kind.Code
    assert b.info == "python"
    assert b.text == "x = 1"


def test_markdown_inside_a_fence_is_not_markdown():
    """The whole reason fences are handled first."""
    (b,) = md_blocks("```\n# not a heading\n- not a bullet\n```\n")
    assert b.kind == Kind.Code
    assert "# not a heading" in b.text


def test_an_unclosed_fence_still_ends_somewhere():
    blocks = md_blocks("```\nx = 1\n")
    assert blocks[0].kind == Kind.Code
    assert "x = 1" in blocks[0].text


def test_tildes_fence_too_and_do_not_close_a_backtick_fence():
    (b,) = md_blocks("```\n~~~\nstill code\n```\n")
    assert b.kind == Kind.Code
    assert "~~~" in b.text


ROUND_TRIP = [
    "",
    "one line",
    "no trailing newline",
    "trailing newline\n",
    "\n\n\n",
    "# h\n\ntext\n\n- a\n- b\n\n```py\nx = 1\n```\n\nend\n",
    "para one\npara one cont\n\npara two\n",
    "```\nunclosed\n",
    "\r\nwindows\r\nline endings\r\n",
    "# 見出し\n\n数式 $x^{2}$ を含む段落。\n",
]


@pytest.mark.parametrize("markdown", ROUND_TRIP)
def test_the_sources_rebuild_the_file_exactly(markdown):
    assert "".join(b.source for b in md_blocks(markdown)) == markdown


def test_blank_runs_are_kept_so_the_round_trip_can_be_exact():
    blocks = md_blocks("a\n\n\n\nb\n")
    assert [b.kind for b in blocks] == [Kind.Paragraph, Kind.Blank, Kind.Paragraph]
    assert blocks[1].source == "\n\n\n"


def test_the_word_writer_and_the_viewer_see_the_same_blocks():
    """One scanner, so a document cannot mean two different things."""
    src = "# T\n\ntext with $x$\n\n- item\n\n```sh\necho hi\n```\n"
    blocks = md_blocks(src)
    assert [b.kind for b in blocks if b.kind != Kind.Blank] == [
        Kind.Heading, Kind.Paragraph, Kind.Bullet, Kind.Code]
    # and the maths inside a block is found by the span scanner
    para = next(b for b in blocks if b.kind == Kind.Paragraph)
    assert [p.latex for p in equation.split_math(para.text) if p.is_math] == ["x"]
