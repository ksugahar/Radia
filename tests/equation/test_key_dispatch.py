"""What a key press does -- checked by pressing keys, without a keyboard.

This is the layer that was once entirely dead: every Ctrl+letter chord in the
editor did nothing, and six hundred tests stayed green, because the decision
lived inside a WIN32 window and nothing below the Python binding could reach
it.  The answer at the time was a probe that took over the real keyboard and
the foreground for a minute at a time -- which is what you do when there is no
seam, and an imposition on whoever is at the machine.

There is a seam now.  The window's whole part is reading the modifier state,
which only a window can; everything that decides anything takes the modifiers
as arguments.  So the chords, the two-step prefixes and the rule about not
typing a chord's own key are all ordinary function calls.

Virtual-key codes, so a test says what a person's fingers do.
"""

from __future__ import annotations

import pytest

equation = pytest.importorskip("radia.equation")
Equation = equation.Equation

if not hasattr(Equation(), "press"):
    pytest.skip("built before the key dispatch had a seam",
                allow_module_level=True)

VK = {c: ord(c) for c in "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"}
VK.update({"Tab": 0x09, "Left": 0x25, "Up": 0x26, "Right": 0x27, "Down": 0x28,
           "Home": 0x24, "End": 0x23, "Insert": 0x2D, "Delete": 0x2E,
           "Backspace": 0x08})


def press(e, key, ctrl=False, shift=False, alt=False):
    return e.press(VK[key], ctrl, shift, alt)


def fresh(latex=""):
    e = Equation()
    if latex:
        e.load_latex(latex)
        e.move_end()
    return e


# ---- the templates ----------------------------------------------------------

@pytest.mark.parametrize("key,fragment", [
    ("F", r"\dfrac"), ("R", r"\sqrt"), ("H", "^"), ("L", "_"),
    ("I", r"\int"),
])
def test_ctrl_letter_inserts_its_template(key, fragment):
    """All of these once did nothing at all."""
    e = fresh()
    assert press(e, key, ctrl=True) == "consumed", key
    assert fragment in e.latex(), (key, e.latex())


def test_a_plain_letter_is_not_a_command():
    """It has to fall through to typing, or nothing could be written."""
    assert press(fresh(), "F") == "ignored"


def test_ctrl_and_ctrl_shift_are_different_keys():
    """Ctrl+F is a fraction and Ctrl+Shift+F is the Function style.  Reading
    the letter as a character rather than a key made every plain Ctrl+letter
    chord unreachable, because typing a capital F does take shift."""
    a, b = fresh(), fresh()
    press(a, "F", ctrl=True)
    press(b, "F", ctrl=True, shift=True)
    assert r"\dfrac" in a.latex()
    assert r"\frac" not in b.latex()


# ---- two-step chords --------------------------------------------------------

def test_the_prefix_waits_for_its_second_key():
    e = fresh()
    assert press(e, "T", ctrl=True) == "pending"
    assert e.latex() == ""


def test_and_then_builds_the_template():
    e = fresh()
    press(e, "T", ctrl=True)
    assert press(e, "S") == "consumed"
    assert r"\sum" in e.latex()


def test_the_second_key_is_not_also_typed():
    """The chord consumed the key and WM_CHAR delivered it again, so summation
    came out as \\sum_{sn}.  Whatever answers the press has to say it swallowed
    the key, which is what "consumed" is for."""
    e = fresh()
    press(e, "T", ctrl=True)
    press(e, "S")
    assert "s" not in e.latex().replace(r"\sum", "")


def test_a_prefix_followed_by_nonsense_is_spent_not_typed():
    e = fresh()
    press(e, "T", ctrl=True)
    assert press(e, "Q") == "consumed"      # eaten, not left to be typed
    assert "q" not in e.latex().lower()


# ---- the caret --------------------------------------------------------------

def test_up_reaches_the_numerator():
    e = fresh(r"\frac{ab}{cd}")
    press(e, "Left")
    press(e, "Left")
    below = e.caret().rsplit(":", 1)[0]
    assert press(e, "Up") == "consumed"
    assert e.caret().rsplit(":", 1)[0] != below


def test_insert_moves_like_tab():
    a, b = fresh(), fresh()
    press(a, "F", ctrl=True)
    press(b, "F", ctrl=True)
    press(a, "Tab")
    press(b, "Insert")
    assert a.caret() == b.caret()


def test_shift_left_extends_rather_than_moves():
    e = fresh("ab")
    assert press(e, "Left", shift=True) == "consumed"
    assert e.has_selection()


def test_left_alone_does_not_select():
    e = fresh("ab")
    press(e, "Left")
    assert not e.has_selection()


# ---- the styles -------------------------------------------------------------

@pytest.mark.parametrize("key,style", [
    ("E", "text"), ("F", "function"), ("I", "variable"),
    ("B", "vector"), ("G", "greek"),
])
def test_ctrl_shift_letter_sets_equation_editors_style(key, style):
    e = fresh()
    assert press(e, key, ctrl=True, shift=True) == "consumed", key
    assert e.style() == style


def test_greek_then_typing_gives_greek():
    """The whole point of the mode, exercised the way a person does it."""
    e = fresh()
    press(e, "G", ctrl=True, shift=True)
    e.insert_text("a")
    assert e.latex().strip() == r"\alpha"
