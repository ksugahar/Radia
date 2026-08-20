"""A selection belongs to the slot it was made in.

This is the model-level half of a crash the window self-test found first: a
random editing walk died with an access violation, and AddressSanitizer named
``take_selection`` reading past the end of a NodeList.  The anchor was a bare
index with no record of WHICH slot it indexed, so Tab, Ctrl+Up, a click, undo
-- anything that moved to another slot or rebuilt the tree -- left it pointing
into a list it no longer belonged to.  A shorter list then made the next
insertion read past the end and erase a range that was not there, wrecking the
heap; the editor died seconds later somewhere unrelated.

Fuzzing found it, and fuzzing must not be what keeps it fixed: these run in
milliseconds and say what the rule IS.  The rule is that a stale anchor means
NO selection -- not a selection quietly shrunk to fit, which would delete a
different range than the one the user highlighted.
"""

from __future__ import annotations

import pytest

equation = pytest.importorskip("radia.equation")
Equation = equation.Equation


def _frac_with_text(text="abcdef"):
    """A fraction whose numerator holds `text` and whose denominator is empty."""
    e = Equation()
    e.insert_template("frac")
    e.insert_text(text)
    return e


# ---- the crash itself -------------------------------------------------------

def test_selecting_then_tabbing_then_inserting_does_not_read_past_the_slot():
    """The exact shape of the crash: a long numerator selected, Tab into the
    empty denominator, insert.  The anchor said 6; the denominator holds 0."""
    e = _frac_with_text()
    e.select_all()
    e.next_slot()
    e.insert_template("sqrt")           # take_selection() ran here
    assert r"\sqrt" in e.latex()


def test_the_selection_does_not_survive_into_another_slot():
    e = _frac_with_text()
    e.select_all()
    assert e.has_selection()
    e.next_slot()
    assert not e.has_selection(), (
        "the anchor counted the numerator, so it means nothing here")


def test_undo_ends_a_selection_rather_than_leaving_it_stale():
    """Undo replaces the whole tree under a live anchor."""
    e = Equation()
    e.insert_text("abcdef")
    e.select_all()
    assert e.has_selection()
    e.undo()
    assert not e.has_selection()
    e.insert_template("frac")           # would have indexed the old length
    assert r"\dfrac" in e.latex()


def test_load_latex_ends_a_selection():
    e = Equation()
    e.insert_text("abcdef")
    e.select_all()
    e.load_latex("x")
    assert not e.has_selection()
    e.insert_template("sqrt")
    assert r"\sqrt" in e.latex()


# ---- the rule still lets ordinary selection work ---------------------------

def test_a_selection_in_its_own_slot_is_still_wrapped():
    """The fix must not cost the feature: a template WRAPS what is selected."""
    e = Equation()
    e.insert_text("B")
    e.extend_left()
    assert e.has_selection()
    e.insert_template("sqrt")
    assert e.latex().strip() in (r"\sqrt{B}", r"\sqrt{B }"), e.latex()


def test_select_all_then_delete_clears_the_slot():
    e = Equation()
    e.insert_text("abc")
    e.select_all()
    assert e.delete_selection()
    assert e.latex().strip() == ""


def test_returning_to_the_slot_does_not_resurrect_the_selection():
    """Coming back must not make a stale anchor live again -- the tree may
    have changed while away."""
    e = _frac_with_text()
    e.select_all()
    e.next_slot()
    e.prev_slot()
    assert not e.has_selection()


def test_the_style_menu_ignores_a_stale_selection():
    """set_style() indexed the slot with the same lo/hi."""
    e = _frac_with_text()
    e.select_all()
    e.next_slot()
    assert e.set_style("vector")
    assert e.style() == "vector"


def test_selected_latex_is_empty_when_the_anchor_is_stale():
    e = _frac_with_text()
    e.select_all()
    assert e.selected_latex() != ""
    e.next_slot()
    assert e.selected_latex() == ""
