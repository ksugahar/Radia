"""Caret geometry: the bridge between the editing model and the layout.

The editing model addresses the caret as a path of (child, slot) steps; the
layout knows where things are.  Joining them is what lets an editor draw the
caret and turn a click into a position, and it is the one part of an equation
editor that cannot be checked by looking at the output -- a caret drawn in the
wrong place still produces correct LaTeX.

The invariant worth testing is the round trip: ask where the caret is, click
there, and land back in the same place.  A slot-numbering mismatch between the
layout and node_slots() breaks exactly that, and breaks nothing else.
"""

from __future__ import annotations

import pytest

equation = pytest.importorskip("radia.equation")
Equation = equation.Equation


def build(latex):
    eq = Equation()
    eq.load_latex(latex)
    return eq


CASES = [
    "abc",
    r"\frac{a}{b}",
    r"x_{i}",
    r"x_{i}^{2}",
    r"\sqrt{x+1}",
    r"\sqrt[3]{x}",
    r"\left(a+b\right)",
    r"\sum_{i}^{n} a",
    r"\int_{0}^{1} f",
    r"\frac{\sqrt{a}}{b+c}",
    r"\sigma_{f}(x) = \sum_{j} B_{fj}(x) a_{j}",
]


@pytest.mark.parametrize("latex", CASES)
def test_the_caret_has_a_position_everywhere_it_can_go(latex):
    """Walk the caret through the whole equation; every stop must be drawable."""
    eq = build(latex)
    eq.move_home()
    seen = 0
    for _ in range(200):
        found, x, top, bottom = eq.caret_geometry()
        assert found, f"no position for caret {eq.caret()} in {latex}"
        assert bottom > top, "the caret has no height"
        seen += 1
        if not eq.move_right():
            break
    assert seen > 1


@pytest.mark.parametrize("latex", CASES)
def test_clicking_where_the_caret_is_puts_it_back(latex):
    """The round trip that a slot-numbering mismatch breaks and nothing else
    does.

    Not "the same path": several positions genuinely share a point -- the end of
    an integrand and the end of the integral are the same place -- and a tall
    caret's middle can lie inside a subscript, where a click *should* go.  What
    must hold is that the click lands at the same horizontal position, in a slot
    that covers the height clicked at.  Those two together are what a mismatched
    slot number breaks.
    """
    eq = build(latex)
    eq.move_home()
    for _ in range(200):
        before = eq.caret()
        found, x, top, bottom = eq.caret_geometry()
        assert found
        y = (top + bottom) / 2.0
        assert eq.move_to_point(x, y)
        found2, x2, top2, bottom2 = eq.caret_geometry()
        assert found2
        assert abs(x2 - x) < 0.01, (
            f"{latex}: clicking at {before} ({x:.2f}) landed at "
            f"{eq.caret()} ({x2:.2f})")
        assert top2 - 0.01 <= y <= bottom2 + 0.01, (
            f"{latex}: clicking at y={y:.2f} landed in a slot spanning "
            f"{top2:.2f}..{bottom2:.2f} ({eq.caret()})")
        if not eq.move_right():
            break


def test_the_caret_moves_right_as_it_moves_through_a_line():
    eq = build("abcd")
    eq.move_home()
    xs = []
    while True:
        found, x, _t, _b = eq.caret_geometry()
        assert found
        xs.append(x)
        if not eq.move_right():
            break
    assert xs == sorted(xs), f"positions are not monotone: {xs}"
    assert xs[-1] > xs[0]


def test_a_numerator_sits_above_a_denominator():
    """Which is the whole reason the caret needs a vertical extent at all."""
    eq = build(r"\frac{a}{b}")
    eq.move_home()
    eq.move_right()                      # into the fraction
    num = eq.caret_geometry()
    assert num[0]
    eq.next_slot()                       # into the denominator
    den = eq.caret_geometry()
    assert den[0]
    assert num[2] < den[2], "the numerator's caret is not above the denominator's"


def test_an_empty_slot_still_has_a_visible_caret():
    """A fresh template is exactly where the caret matters most, and an empty
    slot has no extent of its own to borrow."""
    eq = Equation()
    eq.insert_template("frac")
    found, _x, top, bottom = eq.caret_geometry()
    assert found
    assert bottom - top > 1.0, "the caret in an empty numerator has no height"


def test_a_click_far_below_lands_in_the_denominator():
    eq = build(r"\frac{a}{b}")
    w, asc, desc = eq.extents()
    assert eq.move_to_point(w / 2.0, desc * 0.9)
    caret = eq.caret()
    eq2 = build(r"\frac{a}{b}")
    eq2.move_home()
    eq2.move_right()
    eq2.next_slot()                      # the denominator, by construction
    assert caret.split(":")[0] == eq2.caret().split(":")[0]


def test_extents_match_the_rendered_picture():
    import re

    latex = r"\sigma_{f}(x) = \sum_{j} B_{fj}(x) a_{j}"
    eq = build(latex)
    w, asc, desc = eq.extents()
    svg = equation.tex_to_svg(latex)
    m = re.search(r'viewBox="0 0 ([\d.]+) ([\d.]+)"', svg)
    pad = equation.SvgStyle().padding
    assert abs(float(m.group(1)) - (w + 2 * pad)) < 0.01
    assert abs(float(m.group(2)) - (asc + desc + 2 * pad)) < 0.01
