"""The document, laid out: what a viewer draws and what a click means.

Two things here are worth more than the rest.

A click must be able to tell an equation from prose, because that is the whole
interaction: clicking prose opens the text, clicking an equation opens the
equation editor with Equation Editor 3.0's chords.  If `math_at` is wrong the
editor opens the wrong thing.

And Japanese must wrap.  Japanese has no spaces, so breaking only at spaces
would put an entire paragraph on one line -- with kinsoku, so no line begins
with a full stop or a closing bracket.
"""

from __future__ import annotations

import pytest

equation = pytest.importorskip("radia.equation")
layout_markdown = equation.layout_markdown
DocStyle = equation.DocStyle
MdBlock = equation.MdBlock

WIDTH = 340.0

SAMPLE = """\
# Heading

Prose with $x^{2}$ inside it.

$$E = mc^{2}$$

- bullet with $y$
1. one
7. seven

```python
x = 1
```
"""

JAPANESE = (
    "日本語の段落です。"
    "行の折り返しが、空白のない"
    "日本語でも正しく起きるかを"
    "確かめます。数式 $E = mc^{2}$ "
    "も混ぜます。句読点が行頭に"
    "来てはいけません。\n"
)

# Characters that may not begin a line.
FORBIDDEN_AT_LINE_START = (
    "。、）」！？ーっゃゅょ"
    ".,)]}!?"
)


def lines_of(doc):
    """Group runs into lines by their baseline, left to right."""
    out = {}
    for r in doc.runs:
        out.setdefault(r.baseline, []).append(r)
    return [sorted(v, key=lambda r: r.x) for _, v in sorted(out.items())]


# ---- what a click means ----------------------------------------------------

def test_a_click_inside_an_equation_finds_that_equation():
    doc = layout_markdown(SAMPLE, WIDTH)
    assert doc.maths
    for m in doc.maths:
        assert doc.math_at(m.x + m.width / 2, m.baseline - m.ascent / 2) == m.index


def test_a_click_on_prose_finds_no_equation():
    doc = layout_markdown(SAMPLE, WIDTH)
    assert doc.math_at(1.0, doc.height - 1.0) == -1


def test_a_click_just_past_an_equation_misses_it():
    doc = layout_markdown("Text $x$ text\n", WIDTH)
    (m,) = doc.maths
    assert doc.math_at(m.x + m.width + 4.0, m.baseline - 1.0) == -1


def test_a_click_finds_the_block_it_landed_in():
    doc = layout_markdown(SAMPLE, WIDTH)
    for b in doc.blocks:
        mid = (b.top + b.bottom) / 2
        assert doc.block_at(WIDTH / 2, mid) == b.block


def test_blocks_are_in_order_and_do_not_overlap():
    doc = layout_markdown(SAMPLE, WIDTH)
    for a, b in zip(doc.blocks, doc.blocks[1:]):
        assert a.bottom <= b.top
        assert a.block < b.block


# ---- Japanese --------------------------------------------------------------

def test_japanese_wraps_instead_of_running_off_the_page():
    doc = layout_markdown(JAPANESE, WIDTH)
    assert len(lines_of(doc)) > 1


def test_no_line_begins_with_a_forbidden_character():
    doc = layout_markdown(JAPANESE, WIDTH)
    starts = [line[0].text[:1] for line in lines_of(doc)]
    assert starts, "the assertion below would be vacuous"
    assert not [c for c in starts if c in FORBIDDEN_AT_LINE_START]


def test_nothing_is_laid_out_past_the_page():
    doc = layout_markdown(JAPANESE + SAMPLE, WIDTH)
    for r in doc.runs:
        assert r.x <= doc.width
    for m in doc.maths:
        assert m.x + m.width <= doc.width


# ---- the blocks ------------------------------------------------------------

def test_a_numbered_list_shows_the_number_the_file_has():
    """The file says 7, so the viewer says 7 -- it does not renumber."""
    doc = layout_markdown("1. one\n7. seven\n", WIDTH)
    markers = [r.text for r in doc.runs if r.text in ("1.", "7.")]
    assert markers == ["1.", "7."]


def test_code_keeps_one_line_per_source_line():
    doc = layout_markdown("```\na\nb\nc\n```\n", WIDTH)
    baselines = {r.baseline for r in doc.runs if r.mono}
    assert len(baselines) == 3


def test_code_is_not_wrapped_even_when_it_is_long():
    """A wrapped line of code says something the file does not."""
    long_line = "x = " + " + ".join(str(i) for i in range(80))
    doc = layout_markdown("```\n" + long_line + "\n```\n", WIDTH)
    assert len([r for r in doc.runs if r.mono]) == 1


def test_an_inline_code_span_is_set_in_the_monospace_font():
    doc = layout_markdown("Set `x = 1` here.\n", WIDTH)
    mono = [r for r in doc.runs if r.mono]
    assert mono and "".join(r.text for r in mono).replace(" ", "") == "x=1"


def test_a_heading_is_larger_than_body_text():
    doc = layout_markdown("# Title\n\nbody\n", WIDTH)
    st = DocStyle()
    sizes = {r.size for r in doc.runs}
    assert st.heading[0] in sizes
    assert st.body in sizes


def test_display_math_gets_its_own_centred_line():
    doc = layout_markdown("$$x = 1$$\n", WIDTH)
    (m,) = doc.maths
    assert m.display
    assert abs((m.x + m.width / 2) - doc.width / 2) < 1.0


def test_blank_blocks_add_no_geometry_of_their_own():
    a = layout_markdown("one\n\ntwo\n", WIDTH)
    b = layout_markdown("one\n\n\n\n\ntwo\n", WIDTH)
    assert abs(a.height - b.height) < 1e-9


# ---- degenerate input ------------------------------------------------------

@pytest.mark.parametrize("markdown", ["", "\n", "   ", "\n\n\n"])
def test_an_empty_document_lays_out_without_complaint(markdown):
    doc = layout_markdown(markdown, WIDTH)
    assert doc.height > 0
    assert doc.math_at(0, 0) == -1
    assert doc.block_at(0, 0) == -1


def test_a_narrow_page_still_terminates():
    """One item wider than the line must not loop forever."""
    doc = layout_markdown("aaaaaaaaaaaaaaaaaaaaaaaa bb\n", 12.0)
    assert doc.runs


# ---- the equations are the ones that get pasted ----------------------------

def test_the_equations_carry_the_latex_they_were_written_with():
    """The viewer hands this same LaTeX to the Office writers, so what is on
    screen and what is pasted are the same equation."""
    doc = layout_markdown(SAMPLE, WIDTH)
    assert [m.latex for m in doc.maths] == ["x^{2}", "E = mc^{2}", "y"]


def test_every_equation_names_the_block_it_came_from():
    doc = layout_markdown(SAMPLE, WIDTH)
    blocks = {b.block for b in doc.blocks}
    for m in doc.maths:
        assert m.block in blocks


def test_equation_order_follows_the_document():
    doc = layout_markdown(SAMPLE, WIDTH)
    assert [m.index for m in doc.maths] == list(range(len(doc.maths)))


def test_the_style_is_a_parameter_not_a_constant():
    small = layout_markdown(SAMPLE, WIDTH)
    st = DocStyle()
    st.body = 22.0
    big = layout_markdown(SAMPLE, WIDTH, st)
    assert big.height > small.height
