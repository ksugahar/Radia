"""Empty slots: the dotted boxes that make a half-typed equation legible.

A template nobody has typed into has no extent, so a renderer that draws only
what is there shows a fresh fraction as a bar floating in space, with no sign
that there are two places to type or where Tab is about to go.  Equation Editor
drew these boxes, and that is most of what made its structure readable.

The part worth guarding is that the two callers want opposite things.  The
editor MUST show them.  A picture on its way to a slide must NOT -- dotted
rectangles in a finished equation would be a defect a reader sees and the
author does not.  So the layout reports the boxes and the drawing is chosen by
whoever is drawing.

This also explains three palette buttons that looked broken: a cell for a
script, a matrix or an accent is ALL empty slots, so with nothing drawn for
them the button was blank.
"""

from __future__ import annotations

import pytest

equation = pytest.importorskip("radia.equation")
slots = equation.tex_empty_slots


def _style(em):
    st = equation.SvgStyle()
    st.empty_slot_em = em
    return st


EDITING = _style(0.55)      # on screen: an untyped slot must be visible
PICTURE = equation.SvgStyle()   # on a slide: it must not take up room


# ---- what has an empty slot ------------------------------------------------

@pytest.mark.parametrize("latex,count", [
    (r"\dfrac{}{}",   2),
    (r"\dfrac{a}{}",  1),
    (r"\dfrac{a}{b}", 0),
    (r"\sqrt{}",      1),
    (r"x^{2}",        0),
    ("",              1),      # the whole equation is one empty slot
    ("abc",           0),
])
def test_the_layout_reports_the_slots_that_are_empty(latex, count):
    assert len(slots(latex)) == count


def test_an_empty_slot_has_room_to_type_in():
    """Zero width would put the caret on top of the fraction bar.

    Asked for explicitly, because the two callers want opposite things here as
    well.  A PICTURE reserves nothing for an empty slot -- TeX has no such
    notion, and reserving it made a summation with no operand measure 6.6 pt
    wider than the summation -- so that is the default.  The editor asks for
    the room it needs."""
    for x, y, w, h in slots(r"\dfrac{}{}", EDITING):
        assert w > 0
        assert h > 0


def test_a_picture_reserves_nothing_for_an_empty_slot():
    """The box is still reported, so the editor can find it; it just has no
    width.  An equation on a slide is the equation, not the equation plus room
    for something nobody typed."""
    for x, y, w, h in slots(r"\dfrac{}{}"):
        assert w == 0
    empty = equation.tex_metrics(r"\dfrac{a}{}", PICTURE)[0]
    filled = equation.tex_metrics(r"\dfrac{a}{a}", PICTURE)[0]
    assert empty == pytest.approx(filled), (
        "an empty denominator widened the fraction by %.3f pt" % (empty - filled))


def test_the_two_slots_of_a_fraction_are_at_different_heights():
    boxes = slots(r"\dfrac{}{}")
    assert len(boxes) == 2
    ys = sorted(b[1] for b in boxes)
    assert ys[1] > ys[0], "numerator and denominator drew in the same place"


def test_a_filled_template_reports_nothing():
    assert slots(r"\dfrac{\vec{B}}{\mu_{0}}") == []


def test_slots_are_reported_from_nested_templates_too():
    assert len(slots(r"\dfrac{\sqrt{}}{}")) == 2


# ---- the picture stays clean -----------------------------------------------

def test_the_svg_of_a_half_typed_equation_has_no_boxes():
    """A picture must not carry the editor's scaffolding into a slide."""
    svg = equation.tex_to_svg(r"\dfrac{}{}")
    assert "dash" not in svg.lower()
    # The only rect a fraction needs is its bar.
    assert svg.count("<rect") <= 1


def test_a_finished_equation_and_a_half_typed_one_differ_in_the_layout_only():
    """The boxes exist in the layout, not in the picture: the SVG of an empty
    fraction has the bar and nothing else, while the layout has two boxes."""
    assert len(slots(r"\dfrac{}{}")) == 2
    svg = equation.tex_to_svg(r"\dfrac{}{}")
    assert "<text" not in svg          # nothing to typeset


# ---- what fixed the blank buttons ------------------------------------------

TEMPLATE_CELLS = ["sub", "sup", "subsup", "matrix2x2", "matrix3x3",
                  "hat", "vec", "bar", "frac", "sqrt"]


@pytest.mark.parametrize("name", TEMPLATE_CELLS)
def test_a_freshly_inserted_template_has_something_to_show(name):
    """A palette cell renders the result of actually inserting the template.
    With no boxes these were empty, and the button looked broken."""
    e = equation.Equation()
    assert e.insert_template(name)
    assert len(slots(e.latex())) > 0, f"{name} would draw as nothing"
