"""Selecting a range, and what happens to it.

Equation Editor's Edit menu is Select All, Cut, Copy, Paste and Clear, and all
of them rest on this.  Without it a range cannot be deleted or replaced, which
is most of editing -- you can only ever back out the way you came in.

The model is a range within ONE slot.  A template is a single item in its
parent's slot, so selecting a whole fraction is the same operation as selecting
a run of characters, and nothing needs a model spanning levels of the tree.
Shift at a slot edge stops there rather than jumping out: a selection that
silently changes what it is anchored to is worse than one that will not grow.
"""

from __future__ import annotations

import pytest

equation = pytest.importorskip("radia.equation")
Equation = equation.Equation


def at_end(latex):
    e = Equation()
    e.load_latex(latex)
    e.move_end()
    return e


# ---- extending --------------------------------------------------------------

def test_nothing_is_selected_to_begin_with():
    assert not at_end("abc").has_selection()


def test_shift_left_selects():
    e = at_end("abc")
    assert e.extend_left()
    assert e.has_selection()


def test_shift_left_twice_selects_two():
    e = at_end("abc")
    e.extend_left()
    e.extend_left()
    assert e.selected_latex() == "bc"


def test_shift_right_from_the_start_selects_forwards():
    e = Equation()
    e.load_latex("abc")
    e.move_home()
    e.extend_right()
    assert e.selected_latex() == "a"


def test_shifting_back_to_where_it_started_selects_nothing():
    e = at_end("abc")
    e.extend_left()
    e.extend_right()
    assert not e.has_selection()


def test_shift_stops_at_the_slot_edge():
    e = Equation()
    e.load_latex("ab")
    e.move_home()
    assert not e.extend_left()
    assert not e.has_selection()


def test_select_all_takes_the_whole_slot():
    e = at_end("abc")
    e.select_all()
    assert e.selected_latex() == "abc"


def test_shift_home_and_end():
    e = at_end("abc")
    e.extend_home()
    assert e.selected_latex() == "abc"
    e2 = Equation()
    e2.load_latex("abc")
    e2.move_home()
    e2.extend_end()
    assert e2.selected_latex() == "abc"


# ---- a plain move drops it --------------------------------------------------

@pytest.mark.parametrize("move", ["move_left", "move_right", "move_home",
                                  "move_end", "next_slot", "prev_slot"])
def test_moving_without_shift_clears_the_selection(move):
    """A caret that wandered off while a range stayed highlighted would delete
    something other than what is shown."""
    e = at_end(r"\dfrac{ab}{c}")
    e.extend_left()
    getattr(e, move)()
    assert not e.has_selection()


# ---- what replaces it -------------------------------------------------------

def test_typing_replaces_the_selection():
    e = at_end("abc")
    e.extend_left()
    e.extend_left()
    e.insert_text("X")
    assert e.latex() == "aX"


def test_a_symbol_replaces_the_selection():
    e = at_end("abc")
    e.select_all()
    e.insert_symbol(r"\alpha")
    assert e.latex() == r"\alpha"


def test_delete_and_backspace_remove_the_selection():
    for method in ("backspace", "erase"):
        e = at_end("abcd")
        e.extend_left()
        e.extend_left()
        assert getattr(e, method)()
        assert e.latex() == "ab"


def test_replacing_a_selection_is_one_undo_not_two():
    e = at_end("abc")
    e.select_all()
    e.insert_text("X")
    assert e.latex() == "X"
    assert e.undo()
    assert e.latex() == "abc"


# ---- a template wraps it ----------------------------------------------------

def test_a_template_wraps_what_was_selected():
    """Select B, press the vector chord, get a vector.  This is how a vector
    actually gets written, which is why selection came before the Style work."""
    e = at_end("B")
    e.select_all()
    assert e.insert_template("vec")
    assert e.latex() == r"\vec{B}"


@pytest.mark.parametrize("kind,opening", [
    ("vec",   r"\vec{"),
    ("hat",   r"\hat{"),
    ("sqrt",  r"\sqrt{"),
    ("paren", r"\left("),
])
def test_wrapping_covers_the_usual_templates(kind, opening):
    e = at_end("AB")
    e.select_all()
    assert e.insert_template(kind)
    out = e.latex()
    assert out.startswith(opening), out
    assert "AB" in out, out


def test_a_fraction_takes_the_selection_as_its_numerator():
    e = at_end("ab")
    e.select_all()
    e.insert_template("frac")
    out = e.latex()
    assert out.startswith(r"\frac{ab}") or out.startswith(r"\dfrac{ab}"), out


def test_a_script_takes_the_selection_as_its_base():
    e = at_end("xy")
    e.select_all()
    e.insert_template("sup")
    assert e.latex().startswith("xy^"), e.latex()


def test_without_a_selection_a_script_still_takes_the_previous_character():
    """The older behaviour has to survive: Ctrl+L after x subscripts the x."""
    e = at_end("xy")
    e.insert_template("sub")
    assert e.latex().startswith("xy_"), e.latex()


# ---- copying part of an equation --------------------------------------------

def test_the_selection_can_be_read_as_latex():
    e = at_end(r"a + \dfrac{b}{c}")
    e.extend_left()
    assert e.selected_latex() == r"\dfrac{b}{c}"


def test_reading_the_selection_does_not_disturb_the_equation():
    """It borrows the nodes out to write them; they have to come back."""
    e = at_end(r"a + \dfrac{b}{c}")
    before = e.latex()
    e.select_all()
    e.selected_latex()
    assert e.latex() == before
    assert e.selected_latex() == before


def test_nothing_selected_reads_as_nothing():
    assert at_end("abc").selected_latex() == ""


# ---- what the editor highlights ---------------------------------------------

def test_the_selection_has_a_box_to_draw():
    e = at_end("abc")
    e.extend_left()
    box = e.selection_geometry()
    assert box.found
    assert box.x1 > box.x0
    assert box.bottom > box.top


def test_no_selection_means_no_box():
    assert not at_end("abc").selection_geometry().found


def test_a_wider_selection_has_a_wider_box():
    one, two = at_end("abc"), at_end("abc")
    one.extend_left()
    two.extend_left()
    two.extend_left()
    a, b = one.selection_geometry(), two.selection_geometry()
    assert (b.x1 - b.x0) > (a.x1 - a.x0)


# ---- the chords are published -----------------------------------------------

# ---- paste ------------------------------------------------------------------
#
# Without this the editor was a one-way door: an equation could be written and
# sent out, and nothing could be brought back in.  It goes through the model
# rather than replacing the document, so it lands where the caret is.


def test_paste_lands_at_the_caret():
    e = at_end("a + ")
    assert e.insert_latex(r"\dfrac{b}{c}")
    assert e.latex() == r"a+\dfrac{b}{c}"


def test_paste_strips_delimiters():
    """An equation copied out of a Markdown file arrives as an equation."""
    e = at_end("x")
    assert e.insert_latex("$y$")
    assert e.latex() == "xy"


def test_paste_replaces_a_selection_in_one_undo():
    e = at_end("abc")
    e.select_all()
    assert e.insert_latex(r"\alpha")
    assert e.latex() == r"\alpha"
    assert e.undo()
    assert e.latex() == "abc"


def test_paste_goes_inside_the_slot_the_caret_is_in():
    """A whole-document load could not do this."""
    e = Equation()
    e.load_latex(r"\dfrac{}{}")
    e.move_home()
    e.move_right()                      # into the numerator
    assert e.insert_latex("q")
    assert "q" in e.latex()


def test_pasting_nothing_changes_nothing():
    e = at_end("abc")
    assert not e.insert_latex("")
    assert e.latex() == "abc"


def test_what_was_cut_can_be_pasted_back():
    """Ctrl+X then Ctrl+V has to be a no-op, which is the whole contract."""
    e = at_end(r"a + \dfrac{b}{c}")
    e.extend_left()
    taken = e.selected_latex()
    e.delete_selection()
    assert e.insert_latex(taken)
    assert e.latex() == r"a+\dfrac{b}{c}"


def test_the_selection_chords_are_in_the_table():
    """The window builds its keys from this table, so a chord missing here is
    a chord that does not exist."""
    commands = {command for _chord, command, _label in Equation.shortcuts()}
    for name in ("select.left", "select.right", "select.home",
                 "select.end", "select.all"):
        assert name in commands


@pytest.mark.parametrize("command", ["select.left", "select.right",
                                     "select.home", "select.end", "select.all"])
def test_every_selection_command_dispatches(command):
    """From the middle, where every direction has somewhere to go -- at a slot
    edge extending outwards correctly reports that it did nothing."""
    e = Equation()
    e.load_latex("abc")
    e.move_home()
    e.move_right()
    assert e.command(command)
    assert e.has_selection()


@pytest.mark.parametrize("command,at_edge", [
    ("select.left", "home"), ("select.right", "end"),
])
def test_extending_past_a_slot_edge_reports_that_it_did_nothing(command, at_edge):
    e = Equation()
    e.load_latex("abc")
    getattr(e, "move_" + at_edge)()
    assert not e.command(command)
    assert not e.has_selection()
