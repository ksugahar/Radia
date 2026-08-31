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
    assert e.selected_latex() == r"\dfrac{a}{b}"


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


def test_the_item_chords_are_bound():
    """All four, and Ctrl+Shift with them selects.

    These were briefly given away to the Format menu on a misreading: the key
    table's records group by command pair, which merges records of different
    KINDS, and read that way Ctrl+Shift+Left looked like Align Center.  Record
    by record it is kind 6 -- navigation -- in group 3, beside group 1 (bare
    arrows, move by character), group 2 (Shift+arrows, select by character)
    and group 8 (Ctrl+arrows, move by item).  Align lives on Ctrl+Shift+L, C
    and R, which are kind 1."""
    c = chords()
    assert c["Ctrl+Left"] == "caret.left_item"
    assert c["Ctrl+Right"] == "caret.right_item"
    assert c["Ctrl+Shift+Left"] == "select.left_item"
    assert c["Ctrl+Shift+Right"] == "select.right_item"


def test_selecting_by_item_still_works_from_its_command():
    e = equation.Equation()
    e.load_latex(r"\frac{a}{b}c")
    e.move_home()
    assert e.command("select.right_item")
    assert e.selected_latex() == r"\dfrac{a}{b}"


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


# ---- Format > Align ---------------------------------------------------------
#
# The chords are Equation Editor's own, read off the same key table: command
# 3,0 is Ctrl+Shift+L, 3,1 is Ctrl+Shift+C, 3,2 is Ctrl+Shift+R, and group 3
# is the Format menu.  The group numbering was confirmed independently -- the
# Style menu matches all six of its chords and the Edit menu all four of its.


def two_lines():
    e = equation.Equation()
    e.insert_text("a+b")
    e.command("edit.newline")
    e.insert_text("c")
    return e


def first_x_of_each_line(latex):
    import re
    st = equation.SvgStyle()
    st.padding = 0.0
    svg = equation.tex_to_svg(latex, st)
    return [float(m) for m in re.findall(r'<text[^>]*x="([-0-9.]+)"', svg)]


def test_the_three_align_chords_are_bound():
    c = chords()
    assert c["Ctrl+Shift+L"] == "format.left"
    assert c["Ctrl+Shift+C"] == "format.center"
    assert c["Ctrl+Shift+R"] == "format.right"


@pytest.mark.parametrize("how,env", [
    ("center", "gathered"), ("left", "aligned"), ("right", "aligned"),
])
def test_aligning_a_stack_reaches_the_latex(how, env):
    e = two_lines()
    assert e.command("format." + how)
    assert env in e.latex()


def test_left_and_right_survive_being_saved_and_read_back():
    """Flush left and flush right have no environment of their own, so they
    are written as an `aligned` with the content in one column.  A `gathered`
    would come back centred, which is the whole reason for doing it that way."""
    for how in ("left", "center", "right"):
        e = two_lines()
        e.command("format." + how)
        text = e.latex()
        back = equation.Equation()
        back.load_latex(text)
        assert back.latex() == text, how


def test_the_short_line_actually_moves():
    """Not just the markup -- the drawing."""
    pos = {}
    for how in ("left", "center", "right"):
        e = two_lines()
        e.command("format." + how)
        pos[how] = first_x_of_each_line(e.latex())[-1]
    assert pos["left"] < pos["center"] < pos["right"], pos
    assert pos["left"] == 0.0


def test_aligning_outside_a_stack_says_so():
    """There is nothing to align, and returning false is better than silently
    doing nothing to the whole equation."""
    e = equation.Equation()
    e.insert_text("a+b")
    assert not e.command("format.left")


# ---- Align at = -------------------------------------------------------------
#
# Command 3,3 in the key table, on Ctrl+Shift+Right.  It splits each line at
# its first relation into two cells; the layout then sets the first flush right
# and the second flush left, which puts every equals sign on one vertical line.
# That is what \begin{align} does with its &, and the LaTeX written out is
# exactly that.


def two_equations():
    e = equation.Equation()
    e.insert_text("a+b=c")
    e.command("edit.newline")
    e.insert_text("x=y+z")
    return e


def test_align_at_equals_has_a_chord_of_our_own():
    """Equation Editor binds none: the Format menu's fourth item is command
    3,3 and no record in the key table names it.  With no menu bar here, a
    command with no chord is a command nobody can reach."""
    assert chords()["Ctrl+Shift+A"] == "format.at_eq"


def test_align_at_equals_puts_the_signs_on_one_line():
    import re
    e = two_equations()
    assert e.command("format.at_eq")
    st = equation.SvgStyle()
    st.padding = 0.0
    svg = equation.tex_to_svg(e.latex(), st)
    xs = [float(x) for x, t in
          re.findall(r'<text[^>]*x="([-0-9.]+)"[^>]*>([^<]*)</text>', svg)
          if t == "="]
    assert len(xs) == 2, xs
    assert abs(xs[0] - xs[1]) < 0.01, "the two = did not line up: %s" % xs


def test_align_at_equals_writes_a_real_aligned():
    e = two_equations()
    e.command("format.at_eq")
    latex = e.latex()
    assert "aligned" in latex
    assert "&" in latex


def test_align_at_equals_survives_a_save():
    e = two_equations()
    e.command("format.at_eq")
    text = e.latex()
    back = equation.Equation()
    back.load_latex(text)
    assert back.latex() == text


def test_going_back_to_centre_rejoins_the_cells():
    e = two_equations()
    before = e.latex()
    e.command("format.at_eq")
    e.command("format.center")
    assert e.latex() == before


# ---- Size -------------------------------------------------------------------
#
# Equation Editor's Size menu, minus the two entries that no longer mean
# anything here: its Symbol and Sub-Symbol sizes exist to say how big a
# summation sign should be, and that is now read from the font.  A size marker
# is a SWITCH, which is exactly what TeX's \scriptstyle is, so the two map
# onto each other and the setting survives a save.


@pytest.mark.parametrize("cmd,chord", [
    ("size.full", "Ctrl+Shift+1"),
    ("size.sub", "Ctrl+Shift+2"),
    ("size.sub2", "Ctrl+Shift+3"),
])
def test_the_size_chords_are_bound(cmd, chord):
    assert chords()[chord] == cmd


def test_a_size_is_a_mode_for_what_follows():
    st = equation.SvgStyle()
    st.padding = 0.0
    e = equation.Equation()
    e.insert_text("A")
    assert e.command("size.sub")
    e.insert_text("bc")
    small = equation.tex_metrics(e.latex(), st)[0]
    plain = equation.tex_metrics("Abc", st)[0]
    assert small < plain, "the marker changed nothing (%.2f vs %.2f)" % (small, plain)


def test_a_size_on_a_selection_stops_at_its_end():
    e = equation.Equation()
    e.insert_text("abc")
    e.select_all()
    assert e.command("size.sub")
    latex = e.latex()
    assert latex.startswith(chr(92) + "scriptstyle")
    assert latex.endswith(chr(92) + "displaystyle")


def test_a_size_survives_a_save():
    for cmd in ("size.sub", "size.sub2"):
        e = equation.Equation()
        e.insert_text("A")
        e.command(cmd)
        e.insert_text("b")
        text = e.latex()
        back = equation.Equation()
        back.load_latex(text)
        assert back.latex() == text, cmd


# ---- one slot up or down ----------------------------------------------------

def test_ctrl_down_goes_from_numerator_to_denominator():
    """Not the same as Down, which reads the drawing: this walks the
    structure and lands at the start of the next slot."""
    e = equation.Equation()
    e.insert_template("frac")
    e.insert_text("p")
    assert e.command("caret.slot_down")
    e.insert_text("q")
    assert e.latex() == chr(92) + "dfrac{p}{q}"


def test_ctrl_up_comes_back():
    e = equation.Equation()
    e.insert_template("frac")
    e.command("caret.slot_down")
    assert e.command("caret.slot_up")
    e.insert_text("n")
    assert "frac{n}" in e.latex()


def test_it_stops_at_the_last_slot_rather_than_leaving():
    e = equation.Equation()
    e.insert_template("frac")
    assert e.command("caret.slot_down")
    assert not e.command("caret.slot_down")


def test_the_slot_chords_are_bound():
    c = chords()
    assert c["Ctrl+Up"] == "caret.slot_up"
    assert c["Ctrl+Down"] == "caret.slot_down"


# ---- explicit space ---------------------------------------------------------
#
# These were being DROPPED: a\\,b came out the same width as ab.  For a lab that
# writes "5\\,mm" all day that is not a small thing.
#
# They are carried as the Unicode spaces that mean the same thing, so they are
# ordinary characters everywhere -- the tree holds one, the writers emit one,
# Word receives one -- rather than a node type every writer would need a case
# for to say what a character already says.

import pytest as _pytest


def width(latex):
    st = equation.SvgStyle()
    st.padding = 0.0
    return equation.tex_metrics(latex, st)[0]


@_pytest.mark.parametrize("latex,em", [
    (r"a\,b", 3.0 / 18.0),
    (r"a\:b", 4.0 / 18.0),
    (r"a\;b", 5.0 / 18.0),
    (r"a\quad b", 1.0),
])
def test_a_space_is_as_wide_as_tex_makes_it(latex, em):
    """At 12 point: 2, 2.667, 3.333 and 12 pt."""
    assert width(latex) - width("ab") == _pytest.approx(em * 12.0, abs=0.01)


@_pytest.mark.parametrize("latex", [r"a\,b", r"a\:b", r"a\;b", r"a\quad b"])
def test_a_space_reads_back_as_what_was_typed(latex):
    e = equation.Equation()
    e.load_latex(latex)
    back = equation.Equation()
    back.load_latex(e.latex())
    assert back.latex() == e.latex()
    assert chr(92) in e.latex(), "the space came back as a raw character"


def test_the_four_space_chords_are_bound():
    """Equation Editor keeps them in one command group, on these four keys."""
    c = chords()
    assert c["Shift+Space"] == "symbol.8198"        # U+2006, thin
    assert c["Ctrl+Space"] == "symbol.8287"         # U+205F, medium
    assert c["Ctrl+Shift+Space"] == "symbol.8197"   # U+2005, thick
    assert c["Ctrl+Alt+Space"] == "symbol.8195"     # U+2003, quad


def test_typing_a_unit_the_way_a_physicist_does():
    e = equation.Equation()
    e.insert_text("5")
    assert e.command("symbol.8198")
    e.insert_text("mm")
    assert chr(92) + "," in e.latex()
    assert width(e.latex()) > width("5mm")


# ---- the templates that have two ways in -----------------------------------

@_pytest.mark.parametrize("second,cmd", [
    ("F", "template.frac"), ("R", "template.sqrt"), ("H", "template.sup"),
    ("L", "template.sub"), ("J", "template.subsup"), ("I", "template.int"),
    ("N", "template.nthroot"), ("S", "template.sum"), ("P", "template.prod"),
    ("M", "template.matrix2x2"),
])
def test_ctrl_t_reaches_the_template_too(second, cmd):
    """Equation Editor carries each of these twice -- once as a direct chord
    and once as a second key after the prefix, with the same command pair
    either way -- so F reaches the fraction whether it follows Ctrl+T or Ctrl.
    Only one half of that was here."""
    assert chords()["Ctrl+T, " + second] == cmd


def test_both_bracket_keys_reach_parentheses():
    c = chords()
    assert c["Ctrl+9"] == "template.paren"
    assert c["Ctrl+0"] == "template.paren"


# ---- what an unknown command does ------------------------------------------
#
# It used to drop the backslash and set the rest as text, so inom{n}{k} drew
# as the word "binom" followed by an n and a k and looked like a finished
# equation.  That is the quiet-wrong-answer failure this repository has a
# policy against.  It keeps the backslash now, which says what has happened.


def test_an_unknown_command_fails_loudly():
    e = equation.Equation()
    with _pytest.raises(ValueError, match="unsupported TeX control sequence"):
        e.load_latex(chr(92) + "xhookrightarrow{f}")


def test_unknown_command_arguments_are_not_misrepresented_as_valid_output():
    e = equation.Equation()
    with _pytest.raises(ValueError, match="unsupported TeX control sequence"):
        e.load_latex(chr(92) + "boxed{abc}")


# ---- the constructs added here ---------------------------------------------

def test_binom_is_a_ruleless_fraction_in_parentheses():
    """It is built as one -- a FracNode with the rule turned off, inside a
    paren fence -- and written back out as \binom, because the parts alone
    would go out as {n \atop k} and \atop is infix: the parser has no way to
    take it, so the equation came back with the word "atop" in it."""
    e = equation.Equation()
    e.load_latex(chr(92) + "binom{n}{k}")
    assert e.latex() == chr(92) + "binom{n}{k}"
    assert width(chr(92) + "binom{n}{k}") > width("nk")


def test_a_binomial_has_no_rule_and_a_fraction_does():
    import re
    st = equation.SvgStyle()
    st.padding = 0.0
    rules = lambda t: equation.tex_to_svg(t, st).count("<rect")
    assert rules(chr(92) + "frac{n}{k}") > rules(chr(92) + "binom{n}{k}")


@_pytest.mark.parametrize("latex,taller", [
    (chr(92) + "overbrace{a+b}", True),
    (chr(92) + "overbrace{a+b}^{n}", True),
])
def test_an_overbrace_is_drawn(latex, taller):
    """The node existed and its layout returned the CONTENT and nothing else,
    so the brace was simply absent with no sign that anything was missing."""
    st = equation.SvgStyle()
    st.padding = 0.0
    plain = equation.tex_metrics("a+b", st)[1]
    braced = equation.tex_metrics(latex, st)[1]
    assert braced > plain + 2.0, "no brace was drawn (%.2f vs %.2f)" % (braced, plain)


def test_an_underbrace_goes_below():
    st = equation.SvgStyle()
    st.padding = 0.0
    plain = equation.tex_metrics("x+y", st)[2]
    braced = equation.tex_metrics(chr(92) + "underbrace{x+y}_{m}", st)[2]
    assert braced > plain + 2.0


@_pytest.mark.parametrize("latex", [
    chr(92) + "overbrace{a+b}^{n}",
    chr(92) + "underbrace{x+y}_{m}",
    chr(92) + "binom{n}{k}",
])
def test_the_new_constructs_survive_a_save(latex):
    e = equation.Equation()
    e.load_latex(latex)
    once = e.latex()
    back = equation.Equation()
    back.load_latex(once)
    assert back.latex() == once


@_pytest.mark.parametrize("kind", ["bigcup", "bigcap", "coprod"])
def test_the_remaining_large_operators_are_offered(kind):
    """The layout could always draw them; there was no way in."""
    assert kind in equation.Equation.templates()
    e = equation.Equation()
    assert e.insert_template(kind)


def test_ctrl_t_u_reaches_the_union():
    """Command 39,1 in the key table, which had no template to reach."""
    assert chords()["Ctrl+T, U"] == "template.bigcup"
