"""Markdown -> PowerPoint, checked at the boundary that matters.

The same contract as the Word side: the equations arrive as EQUATIONS, one
``<m:oMath>`` per equation in the source, inside the ``<a14:m>`` wrapper a
slide uses.  A picture would follow neither the theme font nor the text size,
and could not be edited in PowerPoint's own editor.

The two layout rules are here because both were wrong in the first deck built
from a real note: every display equation carried a bullet, and a dense slide
ran off the bottom of the page.
"""

from __future__ import annotations

import re
import zipfile

import pytest

equation = pytest.importorskip("radia.equation")
pytest.importorskip("pptx")
pytest.importorskip("lxml")

B = chr(92)

SOURCE = """# First slide

Text with $a_{j}$ inline.

$$\\int_{0}^{1} f(x)\\, dx$$

- bullet with $\\frac{a}{b}$

# Second slide

Just words.

$$E = mc^{2}$$
"""


def slides(path):
    with zipfile.ZipFile(path) as z:
        return [z.read(n).decode("utf-8") for n in sorted(z.namelist())
                if re.match(r"ppt/slides/slide\d+\.xml$", n)]


def test_every_equation_arrives_as_an_equation(tmp_path):
    out = str(tmp_path / "deck.pptx")
    equation.markdown_to_pptx(SOURCE, out)
    from radia.equation.office import count_equations
    assert count_equations(out) == 4


def test_a_heading_starts_a_slide(tmp_path):
    out = str(tmp_path / "deck.pptx")
    equation.markdown_to_pptx(SOURCE, out)
    xml = slides(out)
    assert len(xml) == 2, len(xml)
    assert "First slide" in xml[0]
    assert "Second slide" in xml[1]


def test_a_title_makes_a_cover(tmp_path):
    out = str(tmp_path / "deck.pptx")
    equation.markdown_to_pptx(SOURCE, out, title="A note")
    xml = slides(out)
    assert len(xml) == 3, len(xml)
    assert "A note" in xml[0]


def test_the_equations_are_inside_the_powerpoint_wrapper(tmp_path):
    r"""<a14:m> is what makes it a PowerPoint equation rather than loose
    OMML that PowerPoint ignores."""
    out = str(tmp_path / "deck.pptx")
    equation.markdown_to_pptx(SOURCE, out)
    xml = "".join(slides(out))
    assert "a14:m" in xml, xml[:400]


def test_a_display_equation_is_not_a_bullet_point(tmp_path):
    r"""PowerPoint bullets every paragraph until told otherwise, so the first
    deck put a dot in front of each $$...$$."""
    out = str(tmp_path / "deck.pptx")
    equation.markdown_to_pptx("# S\n\n$$E = mc^{2}$$\n", out)
    xml = slides(out)[0]
    assert "buNone" in xml, xml[:400]
    assert 'algn="ctr"' in xml, xml[:400]


def test_the_body_asks_to_shrink_rather_than_overflow(tmp_path):
    r"""A slide whose last line is cut in half by the edge of the page has not
    been made."""
    out = str(tmp_path / "deck.pptx")
    equation.markdown_to_pptx(SOURCE, out)
    assert "normAutofit" in slides(out)[0]


def test_literal_dollars_are_not_equations(tmp_path):
    out = str(tmp_path / "deck.pptx")
    equation.markdown_to_pptx("# S\n\nPrices are $5 and $6.\n", out)
    from radia.equation.office import count_equations
    assert count_equations(out) == 0


def test_a_note_with_no_heading_still_makes_a_slide(tmp_path):
    out = str(tmp_path / "deck.pptx")
    equation.markdown_to_pptx("Just a paragraph with $x$.\n", out)
    from radia.equation.office import count_equations
    assert count_equations(out) == 1
    assert len(slides(out)) == 1
