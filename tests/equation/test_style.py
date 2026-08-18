"""Equation Editor's Style menu: Math, Text, Function, Variable, Matrix-Vector, Greek.

The tree carried these typefaces all along and nothing could set them, so a
bold vector -- most of what this lab writes -- was unreachable from the editor.

Greek is not a typeface but a keyboard.  Equation Editor's Greek style maps the
Latin letters onto the Adobe Symbol layout, so `a` gives alpha and `q` gives
theta; the character itself changes.  Restyling a Latin `a` instead would leave
a Latin `a` in the OMML and Word would set it in a Latin font.
"""

from __future__ import annotations

import pytest

equation = pytest.importorskip("radia.equation")
Equation = equation.Equation


def styled(latex, style):
    e = Equation()
    e.load_latex(latex)
    e.move_end()
    e.select_all()
    assert e.set_style(style)
    return e


# ---- the styles exist and apply ---------------------------------------------

def test_the_style_list_is_the_old_editors_menu():
    assert set(Equation.styles()) == {
        "math", "text", "function", "variable", "vector", "greek"}


@pytest.mark.parametrize("style", ["math", "text", "function",
                                   "variable", "vector", "greek"])
def test_every_style_applies_to_a_selection(style):
    assert styled("B", style).latex()


def test_an_unknown_style_is_refused():
    e = Equation()
    e.load_latex("B")
    e.select_all()
    assert not e.set_style("bogus")


def test_without_a_selection_it_sets_the_mode_instead():
    """It used to refuse, so a style could only ever be applied to something
    already typed -- one Greek letter for three operations.  Equation Editor's
    Style menu is a mode, so with nothing highlighted this changes what comes
    NEXT rather than doing nothing.  What it must never do is silently restyle
    the whole equation."""
    e = Equation()
    e.load_latex("B")
    e.move_end()
    assert e.set_style("vector")
    assert e.style() == "vector"
    assert e.latex() == "B"          # what was already there is untouched


# ---- what the lab writes ----------------------------------------------------

def test_matrix_vector_is_bold():
    assert styled("B", "vector").latex() == r"\mathbf{B}"


def test_bold_reaches_office_as_bold():
    e = styled("B", "vector")
    assert 'm:val="bi"' in equation.tex_to_omml(e.latex())


def test_text_style_marks_it_as_text():
    assert styled("B", "text").latex() == r"\text{B}"


def test_math_style_puts_it_back():
    e = styled("B", "vector")
    e.select_all()
    assert e.set_style("math")
    assert e.latex() == "B"


# ---- greek is a keyboard, not a font ----------------------------------------

@pytest.mark.parametrize("latin,expected", [
    ("a", r"\alpha"), ("b", r"\beta"), ("g", r"\gamma"),
    ("q", r"\theta"), ("p", r"\pi"), ("w", r"\omega"),
    ("D", r"\Delta"), ("S", r"\Sigma"), ("W", r"\Omega"),
])
def test_the_symbol_keyboard(latin, expected):
    assert styled(latin, "greek").latex().strip() == expected


UPPERCASE = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"


@pytest.mark.parametrize("letter", list(UPPERCASE))
def test_no_uppercase_greek_letter_is_lost(letter):
    """Beta, Eta, Rho and the rest look like Latin letters, so LaTeX has no
    \\Beta to write them with -- and they used to come out as '?'.  The
    fallback has to be the character itself, the same fix the accents and the
    kanji needed."""
    out = styled(letter, "greek").latex().strip()
    assert out
    assert "?" not in out


@pytest.mark.parametrize("letter", list(UPPERCASE))
def test_uppercase_greek_survives_a_round_trip(letter):
    once = styled(letter, "greek").latex()
    again = Equation()
    again.load_latex(once)
    assert again.latex() == once


def test_greek_leaves_non_letters_alone():
    e = styled("a+b", "greek")
    out = e.latex()
    assert r"\alpha" in out and r"\beta" in out and "+" in out


# ---- it reaches inside structures -------------------------------------------

def test_style_reaches_into_a_fraction():
    """A selection may hold a whole fraction; restyling only its top level
    would leave the numerator in the old face."""
    e = Equation()
    e.load_latex(r"\dfrac{ab}{cd}")
    e.move_end()
    e.select_all()
    assert e.set_style("greek")
    out = e.latex()
    for cmd in (r"\alpha", r"\beta", r"\chi", r"\delta"):
        assert cmd in out, out


def test_style_is_one_undo():
    e = styled("abc", "greek")
    assert e.undo()
    assert e.latex() == "abc"


# ---- the chords are published -----------------------------------------------

def test_the_style_chords_are_in_the_table():
    commands = {c for _chord, c, _label in Equation.shortcuts()}
    for s in Equation.styles():
        assert f"style.{s}" in commands


@pytest.mark.parametrize("style", ["math", "text", "function",
                                   "variable", "vector", "greek"])
def test_every_style_dispatches_as_a_command(style):
    e = Equation()
    e.load_latex("B")
    e.move_end()
    e.select_all()
    assert e.command(f"style.{style}")
