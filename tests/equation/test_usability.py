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


# ---- one symbol, from Equation Editor's own table ---------------------------
#
# Ctrl+K and a letter.  The pairing was read out of the binary, not recalled:
# the key resource stores these as records whose command number IS the code
# point -- 8706 is U+2202, 8734 is U+221E -- which is the only part of that
# resource still legible, since its label column is mis-joined.

@pytest.mark.parametrize("key,cp,latex", [
    ("T", 0x00D7, r"\times"),
    ("A", 0x2192, r"\rightarrow"),
    ("D", 0x2202, r"\partial"),
    ("E", 0x2208, r"\in"),
    ("I", 0x221E, r"\infty"),
    ("C", 0x2282, r"\subset"),
])
def test_ctrl_k_inserts_the_symbol_the_table_names(key, cp, latex):
    c = chords()
    assert c["Ctrl+K, " + key] == "symbol.%d" % cp
    e = equation.Equation()
    assert e.command("symbol.%d" % cp)
    assert latex in e.latex()


def test_the_shifted_pairs_are_the_negated_ones():
    """Shift+E is "not an element of" and Shift+C "not a subset of" -- the
    table pairs them with the plain letters deliberately."""
    c = chords()
    assert c["Ctrl+K, Shift+E"] == "symbol.%d" % 0x2209
    assert c["Ctrl+K, Shift+C"] == "symbol.%d" % 0x2284


def test_ctrl_k_and_ctrl_g_are_separate_prefixes():
    c = chords()
    assert c["Ctrl+K, D"].startswith("symbol.")
    assert c["Ctrl+G, D"].startswith("greek.")


# ---- Enter breaks the line --------------------------------------------------
#
# Only possible now that a stack has a layout.  Before, a second line would
# have been stored, written out, read back -- and drawn as nothing.

def test_enter_makes_a_second_line():
    e = equation.Equation()
    e.insert_text("a+b")
    assert e.command("edit.newline")
    e.insert_text("c=d")
    latex = e.latex()
    assert "gathered" in latex
    assert "a+b" in latex and "c = d" in latex


def test_enter_takes_the_rest_of_the_line_with_it():
    """It splits at the caret, which is what Enter does everywhere else."""
    e = equation.Equation()
    e.insert_text("xy")
    e.move_home()
    e.move_right()
    e.command("edit.newline")
    latex = e.latex()
    assert "x" in latex.split(chr(92) * 2)[0]
    assert "y" in latex.split(chr(92) * 2)[1]


def test_a_third_line_goes_after_the_second_not_around_it():
    e = equation.Equation()
    e.insert_text("a")
    e.command("edit.newline"); e.insert_text("b")
    e.command("edit.newline"); e.insert_text("c")
    assert e.latex().count(chr(92) * 2) == 2, e.latex()


def test_the_stack_is_taller_than_one_line():
    """The layout draws it -- which is the thing that was missing."""
    st = equation.SvgStyle()
    st.padding = 0.0
    one = equation.Equation(); one.insert_text("a")
    two = equation.Equation(); two.insert_text("a")
    two.command("edit.newline"); two.insert_text("b")
    h1 = sum(equation.tex_metrics(one.latex(), st)[1:])
    h2 = sum(equation.tex_metrics(two.latex(), st)[1:])
    assert h2 > h1 * 1.5, "two lines drew no taller than one (%.2f vs %.2f)" % (h2, h1)


def test_enter_is_bound():
    assert chords()["Enter"] == "edit.newline"
