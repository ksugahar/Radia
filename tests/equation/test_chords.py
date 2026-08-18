"""What the window actually watches for when you press a key.

Equation::shortcuts() publishes the chords as text so one table serves the
window, the documentation and the user.  Everything in this file is about the
step in between: turning "Ctrl+F" into a virtual-key code and a set of
modifiers.  It had no test, and a mistake there cost the editor every keyboard
shortcut it has -- Ctrl+F, Ctrl+A, Ctrl+I, all of them dead -- while the whole
suite stayed green, because the suite drives the model directly and never asks
what a key press does.
"""

from __future__ import annotations

import pytest

equation = pytest.importorskip("radia.equation")
Equation = equation.Equation

VK = {name: ord(name) for name in "ABCDEFGHIJKLMNOPQRSTUVWXYZ"}


def steps(chord):
    return [tuple(s) for s in Equation.chord_steps(chord)]


# ---- a letter names a key; it does not ask for shift -------------------------

def test_ctrl_f_is_the_f_key_without_shift():
    """The bug, stated once.  Asking the keyboard layout what character "F" is
    answers "shift plus the F key", because that is how you type a capital F --
    so Ctrl+F silently became Ctrl+Shift+F and nothing matched a plain Ctrl+F."""
    assert steps("Ctrl+F") == [(VK["F"], True, False, False)]


def test_ctrl_a_is_the_a_key_without_shift():
    assert steps("Ctrl+A") == [(VK["A"], True, False, False)]


@pytest.mark.parametrize("chord,command,label", Equation.shortcuts())
def test_no_press_asks_for_a_modifier_it_did_not_name(chord, command, label):
    """The invariant the fix restores: what the window waits for is what the
    table says, and nothing more.

    Per press, not per chord -- the second half of "Ctrl+T, S" is a bare S.
    Shift is checked only where the key is a letter or a name; punctuation is
    written as the character, so there the layout decides and "{" carries a
    shift the table never spells out."""
    parts = chord.split(", ")
    got = steps(chord)
    assert len(got) == len(parts), chord

    for part, (vk, ctrl, shift, alt) in zip(parts, got):
        assert ctrl == ("Ctrl+" in part), part
        assert alt == ("Alt+" in part), part
        key = part.rsplit("+", 1)[-1] if part.endswith("+") is False else part
        if key.isalpha():
            assert shift == ("Shift+" in part), part


@pytest.mark.parametrize("chord,command,label", Equation.shortcuts())
def test_every_published_chord_resolves(chord, command, label):
    """A chord that fails to parse is dropped, so the shortcut just goes
    missing -- silently, which is how this class of bug survives."""
    assert steps(chord), chord


# ---- shift, when it is asked for --------------------------------------------

def test_shift_is_honoured_when_the_table_says_so():
    assert steps("Ctrl+Shift+F") == [(VK["F"], True, True, False)]


def test_ctrl_f_and_ctrl_shift_f_are_different_keys():
    """They are two different commands -- fraction and function style -- so if
    they resolved the same way one of them would be unreachable."""
    assert steps("Ctrl+F") != steps("Ctrl+Shift+F")


def test_shift_tab_is_tab_with_shift():
    assert steps("Shift+Tab") == [(0x09, False, True, False)]


# ---- punctuation is a character, and the layout decides ----------------------

def test_brackets_and_braces_are_one_key_told_apart_by_shift():
    """"[" and "{" are the same physical key on this layout; the shift flag is
    the whole difference, and here it genuinely comes from the character."""
    (vk_open, _, shift_open, _), = steps("Ctrl+[")
    (vk_brace, _, shift_brace, _), = steps("Ctrl+{")
    assert vk_open == vk_brace
    assert shift_open is False
    assert shift_brace is True


def test_a_digit_needs_no_shift():
    (_vk, ctrl, shift, alt), = steps("Ctrl+9")
    assert (ctrl, shift, alt) == (True, False, False)


# ---- two-step chords ---------------------------------------------------------

def test_a_two_step_chord_is_two_presses():
    assert steps("Ctrl+T, S") == [(VK["T"], True, False, False),
                                  (VK["S"], False, False, False)]


def test_the_two_step_chords_hang_off_a_small_number_of_prefixes():
    """Templates hang off Ctrl+T, the Greek alphabet off Ctrl+G and the symbols
    off Ctrl+K -- Equation Editor's own three -- so the window has to hold any
    of those first presses and wait rather than act on it.

    Two prefixes, not one: this test used to pin Ctrl+T as the only one, which
    was true when the only two-step chords were summation, product and matrix.
    Keeping the set small is the thing worth guarding -- every prefix is a key
    that does nothing on its own, and a person who presses it by accident is
    left waiting."""
    prefixes = {steps(c)[0] for c, _cmd, _lbl in Equation.shortcuts()
                if ", " in c}
    assert prefixes == {(VK["T"], True, False, False),
                        (VK["G"], True, False, False),
                        (VK["K"], True, False, False)}


def test_the_second_key_of_a_greek_chord_carries_the_shift():
    """Ctrl+G then A is alpha; Ctrl+G then Shift+A is capital alpha.  The
    shift belongs to the SECOND press, which is the part a chord table can get
    wrong without anything else noticing."""
    assert steps("Ctrl+G, A") == [(VK["G"], True, False, False),
                                  (VK["A"], False, False, False)]
    assert steps("Ctrl+G, Shift+A") == [(VK["G"], True, False, False),
                                        (VK["A"], False, True, False)]


# ---- named keys --------------------------------------------------------------

@pytest.mark.parametrize("name,vk", [
    ("Tab", 0x09), ("Left", 0x25), ("Right", 0x27), ("Home", 0x24),
    ("End", 0x23), ("Backspace", 0x08), ("Delete", 0x2E),
])
def test_named_keys(name, vk):
    assert steps(name) == [(vk, False, False, False)]


def test_an_unknown_key_name_is_refused():
    assert steps("Ctrl+Nonsense") == []
