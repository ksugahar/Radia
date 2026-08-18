"""Two things Equation Editor does that this did not.

Both come from its own key table, read out of the binary rather than guessed
at, and both are about the difference between an operation that is right for
correcting a letter and one that is right for getting on with the equation.
"""

from __future__ import annotations

import pytest

equation = pytest.importorskip("radia.equation")

if not hasattr(equation.Equation, "shortcuts"):
    pytest.skip("built without the shortcut table", allow_module_level=True)


def chords():
    return {b[0]: b[1] for b in equation.Equation.shortcuts()}


# ---- one Greek letter, without changing the style --------------------------
#
# Equation Editor has BOTH Ctrl+Shift+G, which switches the style so that
# everything typed after it is Greek, and Ctrl+G, which is a prefix for a
# single letter.  They are not interchangeable: a line of algebra with one mu
# in it wants the second, and only the first existed here.

def test_the_style_and_the_single_letter_are_different_chords():
    c = chords()
    assert c["Ctrl+Shift+G"] == "style.greek"
    assert c["Ctrl+G, M"] == "greek.m"


def test_a_single_greek_letter_leaves_the_style_alone():
    e = equation.Equation()
    e.insert_text("B")
    assert e.style() == "math"
    assert e.command("greek.m")
    assert e.latex().endswith(r"\mu")
    assert e.style() == "math", "one letter should not switch the keyboard"
    e.insert_text("c")
    assert e.latex().endswith("c"), "typing after it is still Latin"


@pytest.mark.parametrize("latin,name", [
    ("a", "alpha"), ("b", "beta"), ("m", "mu"), ("w", "omega"),
    ("f", "phi"), ("q", "theta"),
])
def test_the_letters_map_the_way_the_style_maps_them(latin, name):
    """The same table the Greek STYLE uses, so the two cannot disagree."""
    e = equation.Equation()
    assert e.command("greek." + latin)
    assert name in e.latex()


def test_every_greek_chord_in_the_table_actually_inserts_something():
    e = equation.Equation()
    n = 0
    for chord, cmd in chords().items():
        if not chord.startswith("Ctrl+G,"):
            continue
        n += 1
        one = equation.Equation()
        assert one.command(cmd), "%s -> %s inserted nothing" % (chord, cmd)
    assert n >= 40, "expected the Greek alphabet, got %d chords" % n


# ---- moving past a template rather than into it ----------------------------

def test_plain_right_walks_into_a_fraction():
    e = equation.Equation()
    e.load_latex(r"\frac{a}{b}c")
    e.move_home()
    assert e.command("caret.right")
    assert e.caret() != ":1", "plain Right should step inside"


def test_ctrl_right_steps_over_the_whole_fraction():
    """Four presses to get past \frac{a}{b} becomes one."""
    e = equation.Equation()
    e.load_latex(r"\frac{a}{b}c")
    e.move_home()
    assert e.command("caret.right_item")
    assert e.caret() == ":1"


def test_ctrl_shift_right_selects_the_whole_item():
    e = equation.Equation()
    e.load_latex(r"\frac{a}{b}c")
    e.move_home()
    assert e.command("select.right_item")
    assert e.has_selection()
    assert e.selected_latex() == r"\frac{a}{b}"


def test_item_motion_stops_at_the_edge_rather_than_leaving_the_slot():
    """A motion that silently changes which slot you are in is worse than one
    that stops -- the same reasoning the selection already follows."""
    e = equation.Equation()
    e.load_latex(r"\frac{a}{b}")
    e.move_home()
    e.command("caret.right")          # into the numerator
    inside = e.caret()
    e.command("caret.end")
    assert not e.command("caret.right_item")
    assert e.caret().split(":")[0] == inside.split(":")[0]


def test_the_four_item_chords_are_bound():
    c = chords()
    assert c["Ctrl+Left"] == "caret.left_item"
    assert c["Ctrl+Right"] == "caret.right_item"
    assert c["Ctrl+Shift+Left"] == "select.left_item"
    assert c["Ctrl+Shift+Right"] == "select.right_item"
