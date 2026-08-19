r"""A big operator's limits, written away from the operator.

Equation Editor puts a large operator's limits after it rather than in it:
the template, then its body, then the limit lines, then a switch to symbol
size and the operator's own sign.  Flat, that block sits beside the operator
it belongs to and Pass 2 finds it by scanning back.

Nested, it does not.  A double product writes the OUTER operator's block past
the end of the enclosing list, where nothing precedes it but a line -- so the
scan gave up and the limits came out as loose text at the end of the equation:

    \prod \limits_{j=1}^{n} \prod (1+x_{ij})\scriptstyle i = 1m∏

Reading into that line finds the operator.  Which one it is takes care: the
deepest one still without limits, EXCEPT where a block is already coming for
it in its own list -- an inline double sum writes the inner operator's block
one level down, and giving the outer's limits to the inner made them vanish
from the output altogether.  Loose text says it is wrong; a silent
disappearance does not.

The streams are built here rather than harvested, so the shape under test is
the thing in the file.
"""

from __future__ import annotations

import pytest

equation = pytest.importorskip("radia.equation")

B = chr(92)

HEADER = bytes([3, 1, 1, 3, 10])

REC_END, REC_LINE, REC_CHAR, REC_TMPL = 0, 1, 2, 3
REC_FULL, REC_SUB, REC_SYM = 10, 11, 13
TM_SUM = 29                 # summation, the operator with limits

TFW_VAR = 128 + 3
TFW_SYM = 128 + 6
TFW_NUM = 128 + 8


def char(ch, tf=TFW_VAR):
    cp = ord(ch)
    return bytes([REC_CHAR, tf, cp & 0xFF, (cp >> 8) & 0xFF])


def line(*records):
    return bytes([REC_LINE]) + b"".join(records) + bytes([REC_END])


def null_line():
    """A LINE record whose options say it is empty."""
    return bytes([REC_LINE | (1 << 4)])


def end():
    return bytes([REC_END])


def bigop():
    """A summation template with an empty body slot."""
    return bytes([REC_TMPL, TM_SUM, 1]) + end()


def sigma():
    return char("∑", TFW_SYM)


def block(lower, upper=None):
    """The display data: SUB, the limits, SYM, the sign."""
    out = bytes([REC_SUB]) + line(lower)
    out += line(upper) if upper is not None else null_line()
    return out + bytes([REC_SYM]) + sigma()


def test_a_flat_operator_keeps_its_limits():
    data = (HEADER + bytes([REC_FULL])
            + line(bigop(), line(char("a")), line(char("i")), null_line(),
                   bytes([REC_SYM]), sigma())
            + end())
    out = equation.mtef_to_latex(data)
    assert B + "sum" in out, out
    assert "i" in out, out


def test_a_nested_operator_gets_the_block_past_the_end():
    r"""The outer sum's block is outside the line that holds it, which is
    where the limits used to be dropped as loose text."""
    inner = line(bigop(), line(char("a")), line(char("j")), null_line(),
                 bytes([REC_SYM]), sigma())
    data = (HEADER + bytes([REC_FULL])
            + line(bigop(), inner)
            + block(char("i"))
            + end())
    out = equation.mtef_to_latex(data)
    assert out.count(B + "sum") == 2, out
    assert out.count(B + "limits") == 2, out
    # both indices survive, and neither is loose text at the end
    assert "i" in out and "j" in out, out
    assert B + "scriptstyle" not in out, out


def test_the_inner_operator_keeps_its_own_block():
    r"""The outer sum's limits must not be handed to the inner one, which has
    its own coming: that is how they went missing."""
    inner = line(bigop(), line(char("a")), line(char("j")), null_line(),
                 bytes([REC_SYM]), sigma())
    data = (HEADER + bytes([REC_FULL])
            + line(bigop(), inner)
            + block(char("i"))
            + end())
    out = equation.mtef_to_latex(data)
    # i belongs to the first sum, j to the second
    assert out.index("i") < out.index("j"), out


def test_a_block_with_nothing_to_own_it_is_left_alone():
    r"""No operator anywhere: refuse rather than attach it to something."""
    data = (HEADER + bytes([REC_FULL])
            + line(char("x"))
            + block(char("i"))
            + end())
    out = equation.mtef_to_latex(data)     # must not throw
    assert "x" in out
