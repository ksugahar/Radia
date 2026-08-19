r"""The display fraction Equation Editor writes in pieces.

An INLINE fraction carries both parts in the template.  A DISPLAY one does
not: the template holds the first line of the numerator, and the rest of it --
and the whole denominator -- follow as siblings, separated by "back to full
size" records.  Read straight that gives \dfrac{a}{} with "+ b" and "c" beside
it, which looks like a finished equation and is not one.

The streams here are assembled byte by byte so the case can be stated without
shipping a lecture file, and so the shape under test is visible in the test
rather than in a fixture.
"""

from __future__ import annotations

import pytest

equation = pytest.importorskip("radia.equation")

B = chr(92)

HEADER = bytes([3, 1, 1, 3, 10])

REC_END, REC_LINE, REC_CHAR, REC_TMPL = 0, 1, 2, 3
REC_FULL = 10
TM_FRACT = 14

TFW_VAR = 128 + 3          # variable
TFW_SYM = 128 + 6          # symbol
TFW_NUM = 128 + 8          # number


def char(ch, tf=TFW_VAR):
    cp = ord(ch)
    return bytes([REC_CHAR, tf, cp & 0xFF, (cp >> 8) & 0xFF])


def line(*records):
    return bytes([REC_LINE]) + b"".join(records) + bytes([REC_END])


def end():
    return bytes([REC_END])


def full():
    return bytes([REC_FULL])


def frac_first_line(numer_first):
    """A fraction template with an empty slot 0 and one line in slot 1 --
    which is how the display form starts."""
    return (bytes([REC_TMPL, TM_FRACT, 0])
            + end()                       # slot 0: empty
            + line(numer_first) + end())  # slot 1: one line, then close


def test_the_pieces_are_put_back_together():
    data = (HEADER + full()
            + frac_first_line(char("a"))
            + char("+", TFW_SYM) + char("b")
            + full()                       # separator
            + char("c")
            + full()                       # terminator
            + end())
    out = equation.mtef_to_latex(data)
    assert B + "dfrac{" in out, out
    # a + b over c, whatever the spacing
    body = out.replace(" ", "")
    assert "a+b" in body, out
    assert "{c}" in body, out


def test_an_inline_fraction_is_left_alone():
    r"""Both parts in the template: nothing to recover, nothing to touch."""
    data = (HEADER
            + bytes([REC_TMPL, TM_FRACT, 0])
            + end()
            + line(char("a")) + line(char("b")) + end()
            + char("+", TFW_SYM) + char("z")
            + end())
    out = equation.mtef_to_latex(data).replace(" ", "")
    assert "{a}{b}" in out, out
    assert "+z" in out, out          # the z stayed outside, as it should


def test_no_terminator_means_no_repair():
    r"""Refuse rather than guess.  A fraction that is merely not repaired is
    visibly wrong; one repaired WRONGLY reads as finished."""
    data = (HEADER + full()
            + frac_first_line(char("a"))
            + char("+", TFW_SYM) + char("b")
            + full()                       # separator, but nothing closes it
            + char("c")
            + end())
    out = equation.mtef_to_latex(data)     # must not hang, must not throw
    assert "c" in out


def test_an_empty_denominator_is_not_invented():
    data = (HEADER + full()
            + frac_first_line(char("a"))
            + char("+", TFW_SYM) + char("b")
            + full()
            + full()                       # nothing between the separators
            + end())
    out = equation.mtef_to_latex(data)
    assert "a" in out and "b" in out


def test_the_stream_still_reads_whole():
    r"""The repair runs after the whole stream is read, not instead of it."""
    data = (HEADER + full()
            + frac_first_line(char("a"))
            + char("+", TFW_SYM) + char("b")
            + full() + char("c") + full()
            + end()
            + char("z") + end())           # a second top-level list
    out = equation.mtef_to_latex(data)
    assert "z" in out, out
