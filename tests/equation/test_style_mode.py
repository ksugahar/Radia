"""Style is a mode, the way Equation Editor 3.0 has it.

Its status bar reads "Style: Math" and keeps reading it, so picking Greek and
carrying on typing gives Greek.  Here it used to apply only to a selection,
which meant one Greek letter cost three operations -- type it, select it back,
restyle it -- against Equation Editor's one.  A lab that writes alpha, mu and
sigma all day pays that three times a line.

The mode also has to be visible.  A mode you cannot see is a trap: after
Ctrl+Shift+G everything typed comes out Greek, and with nowhere saying so the
only way to find out is to type something and be surprised.  Equation Editor
answers that with the status bar; so does this.
"""

from __future__ import annotations

import pytest

equation = pytest.importorskip("radia.equation")
Equation = equation.Equation

if not hasattr(Equation(), "style"):
    pytest.skip("built before Style became a mode", allow_module_level=True)


def fresh():
    return Equation()


# ---- it is a mode -----------------------------------------------------------

def test_it_starts_in_math():
    assert fresh().style() == "math"


def test_setting_it_without_a_selection_is_allowed():
    """This is the whole change: it used to refuse, so the style could only
    ever be applied to something already typed."""
    e = fresh()
    assert e.set_style("greek")
    assert e.style() == "greek"


def test_typing_then_follows_it():
    e = fresh()
    e.set_style("greek")
    e.insert_text("abg")
    out = e.latex()
    assert r"\alpha" in out and r"\beta" in out and r"\gamma" in out, out


def test_one_greek_letter_costs_one_operation():
    """Against three before: type, select, restyle."""
    e = fresh()
    e.set_style("greek")
    e.insert_text("a")
    assert e.latex().strip() == r"\alpha"


def test_it_stays_until_changed():
    e = fresh()
    e.set_style("greek")
    e.insert_text("a")
    e.insert_text("b")
    out = e.latex()
    assert r"\alpha" in out and r"\beta" in out


def test_going_back_to_math_stops_it():
    e = fresh()
    e.set_style("greek")
    e.insert_text("a")
    assert e.set_style("math")
    e.insert_text("b")
    out = e.latex()
    assert r"\alpha" in out
    assert "b" in out and r"\beta" not in out


def test_an_unknown_style_changes_nothing():
    e = fresh()
    e.set_style("greek")
    assert not e.set_style("bogus")
    assert e.style() == "greek"


# ---- and it still does what it did to a selection ---------------------------

def test_a_selection_is_restyled_as_well():
    e = fresh()
    e.insert_text("a")
    e.extend_left()
    assert e.set_style("greek")
    assert e.latex().strip() == r"\alpha"


def test_restyling_a_selection_also_leaves_the_mode_set():
    """Both halves of what the menu does: change what is highlighted, and
    change what comes next."""
    e = fresh()
    e.insert_text("a")
    e.extend_left()
    e.set_style("vector")
    assert e.style() == "vector"


@pytest.mark.parametrize("style,typed,want", [
    ("vector", "B", r"\bm{B}"),
    ("text", "B", r"\text{B}"),
    ("greek", "q", r"\theta"),
])
def test_every_style_reaches_what_is_typed(style, typed, want):
    e = fresh()
    e.set_style(style)
    e.insert_text(typed)
    assert e.latex().strip() == want
