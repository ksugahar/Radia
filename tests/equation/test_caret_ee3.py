"""The caret, as Equation Editor 3.0's own key table describes it.

That table is not in its code, it is a resource -- twelve-byte records the
dispatcher walks, matching a context, a key and a set of modifiers.  Read out,
its caret section says:

    plain      Left Right Up Down     one command group
    Shift+     the same four          same indices, a different group:
                                      the same motion, extending a selection
    Ctrl+      the same four          a DIFFERENT group and different indices,
                                      so these are not the same operation
    Tab, Insert  next slot            the same command; Insert really is a
    Shift+Tab    previous slot        second Tab there
    Home End PageUp PageDown Enter    their own groups

Up and down were simply missing here, which is the gap this closes: no arrow
reached from a denominator to its numerator.
"""

from __future__ import annotations

import pytest

equation = pytest.importorskip("radia.equation")
Equation = equation.Equation


def at(latex, *moves):
    e = Equation()
    e.load_latex(latex)
    e.move_end()
    for m in moves:
        e.command(m)
    return e


def slot(e):
    """The caret's slot, without the index -- which part of the structure."""
    c = e.caret()
    return c.rsplit(":", 1)[0]


# ---- up and down move between the parts of a structure ----------------------

def test_up_from_a_denominator_reaches_the_numerator():
    """The one that was missing.  A numerator is not above a denominator in
    the tree -- only on the page -- so this is decided from the layout."""
    e = at(r"\frac{ab}{cd}", "caret.left", "caret.left")
    below = slot(e)
    assert e.command("caret.up")
    assert slot(e) != below


def test_and_down_comes_back():
    e = at(r"\frac{ab}{cd}", "caret.left", "caret.left")
    start = slot(e)
    e.command("caret.up")
    e.command("caret.down")
    assert slot(e) == start


def test_up_from_a_subscript_reaches_the_superscript():
    e = at(r"x_{i}^{2}", "caret.left")
    before = slot(e)
    if e.command("caret.up"):
        assert slot(e) != before


def test_they_report_failure_rather_than_moving_nowhere():
    """With nothing above, Up has to say so: an editor that silently does
    nothing and one that silently moves somewhere odd look the same."""
    e = at("abc")
    assert not e.command("caret.up")
    assert not e.command("caret.down")


def test_up_does_not_leave_a_selection_behind():
    e = at(r"\frac{ab}{cd}", "caret.left")
    e.extend_left()
    assert e.has_selection()
    e.command("caret.up")
    assert not e.has_selection()


# ---- the chords are the ones Equation Editor publishes ----------------------

def chords():
    return {c: cmd for c, cmd, _label in Equation.shortcuts()}


@pytest.mark.parametrize("chord,command", [
    ("Up", "caret.up"),
    ("Down", "caret.down"),
    ("Left", "caret.left"),
    ("Right", "caret.right"),
    ("Tab", "caret.next_slot"),
    ("Insert", "caret.next_slot"),      # the same command in its table
    ("Shift+Tab", "caret.prev_slot"),
    ("Home", "caret.home"),
    ("End", "caret.end"),
])
def test_the_caret_chord_is_bound(chord, command):
    assert chords().get(chord) == command, chord


def test_insert_is_a_second_tab():
    """Equation Editor's table gives Tab and Insert the same command.  Guessing
    would not have produced that."""
    c = chords()
    assert c["Insert"] == c["Tab"]


@pytest.mark.parametrize("chord,command", [
    ("Ctrl+Shift+=", "style.math"),
    ("Ctrl+Shift+E", "style.text"),
    ("Ctrl+Shift+F", "style.function"),
    ("Ctrl+Shift+I", "style.variable"),
    ("Ctrl+Shift+B", "style.vector"),
    ("Ctrl+Shift+G", "style.greek"),
])
def test_the_style_chords_are_equation_editors(chord, command):
    """It names the style by its EFFECT -- I for italic gives Variable, B for
    bold gives Matrix-Vector -- and "=" returns to plain Math.  Four of these
    were guessed here and four of the guesses were wrong."""
    assert chords().get(chord) == command, chord


@pytest.mark.parametrize("gone", [
    "Ctrl+Shift+M", "Ctrl+Shift+T", "Ctrl+Shift+V", "Ctrl+Shift+X",
])
def test_the_guessed_style_chords_are_gone(gone):
    assert gone not in chords()


def test_the_skewed_fraction_has_its_chord():
    assert chords().get("Ctrl+/") == "template.frac"
