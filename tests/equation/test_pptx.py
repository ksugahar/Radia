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


SCRIPT_SOURCE = """# Kelvin transformation

> Here I explain how the exterior is mapped onto a finite region.

- exterior $\\Omega^{E}$ maps to $\\Omega'$

# What it costs

> One extra region, and a matched interface.

Nothing else.
"""


def slides(path):
    with zipfile.ZipFile(path) as z:
        return [z.read(n).decode("utf-8") for n in sorted(z.namelist())
                if re.match(r"ppt/slides/slide\d+\.xml$", n)]


def notes(path):
    """The speaker notes, one string per notes slide that exists."""
    with zipfile.ZipFile(path) as z:
        return {n: z.read(n).decode("utf-8") for n in sorted(z.namelist())
                if re.match(r"ppt/notesSlides/notesSlide\d+\.xml$", n)}


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


def test_the_script_lands_in_the_notes_not_on_the_slide(tmp_path):
    r"""A blockquote is what you SAY; the slide holds what is SHOWN.

    Writing the script onto the slide is the failure this rule exists to
    prevent -- the audience then reads the talk instead of hearing it.
    """
    out = str(tmp_path / "deck.pptx")
    equation.markdown_to_pptx(SCRIPT_SOURCE, out)
    spoken = "Here I explain how the exterior is mapped"
    assert spoken not in "".join(slides(out))
    assert any(spoken in xml for xml in notes(out).values())


def test_each_script_belongs_to_the_slide_it_follows(tmp_path):
    out = str(tmp_path / "deck.pptx")
    equation.markdown_to_pptx(SCRIPT_SOURCE, out)
    got = notes(out)
    assert len(got) == 2, sorted(got)
    first, second = [got[k] for k in sorted(got)]
    assert "exterior is mapped" in first
    assert "One extra region" in second
    assert "One extra region" not in first


def test_the_quote_marker_is_not_spoken(tmp_path):
    r"""``>`` marks the script; it is not part of it."""
    out = str(tmp_path / "deck.pptx")
    equation.markdown_to_pptx("# S\n\n> Say this.\n", out)
    text = "".join(notes(out).values())
    assert "Say this." in text
    assert "&gt; Say this." not in text


def test_a_slide_with_no_script_has_no_notes(tmp_path):
    r"""An empty notes page is not the same as a slide nobody scripted --
    the coverage check has to be able to see the difference."""
    out = str(tmp_path / "deck.pptx")
    equation.markdown_to_pptx("# S\n\nJust shown text.\n", out)
    assert notes(out) == {}


def test_the_script_may_carry_equations(tmp_path):
    r"""You say the formula out loud too; it should read as one."""
    out = str(tmp_path / "deck.pptx")
    equation.markdown_to_pptx("# S\n\n> recall that $E = mc^{2}$ holds.\n", out)
    assert "m:oMath" in "".join(notes(out).values())


def test_a_note_with_no_heading_still_makes_a_slide(tmp_path):
    out = str(tmp_path / "deck.pptx")
    equation.markdown_to_pptx("Just a paragraph with $x$.\n", out)
    from radia.equation.office import count_equations
    assert count_equations(out) == 1
    assert len(slides(out)) == 1
