r"""An integral's body and limits, which arrive in one slot together.

Equation Editor writes a native integral as a template whose slot holds the
integrand AND the display block -- limits, a switch to symbol size, the
integral signs.  The slot count comes from the variation, so that whole slot
was taken as a limit and the integrand became a superscript:

    p = \iint ^{r\sigma ds\scriptstyle \int \int }

which is not an integral at all.  59 of the corpus's 779 documents were some
version of this.

The other half of the same story is that no pass could see into a slot: the
pipeline ran on a LINE's children, and a slot is a plain list.  An integral
inside a bracket kept its display block for exactly that reason.
"""

from __future__ import annotations

import pytest

equation = pytest.importorskip("radia.equation")

B = chr(92)

HEADER = bytes([3, 1, 1, 3, 10])

REC_END, REC_LINE, REC_CHAR, REC_TMPL = 0, 1, 2, 3
REC_FULL, REC_SUB, REC_SYM = 10, 11, 13

TM_PAREN = 1                # a fence, to put an integral inside
TM_SINT = 21                # single integral

VAR_LOWER, VAR_UPPER = 1, 2

TFW_VAR = 128 + 3
TFW_SYM = 128 + 6
TFW_NUM = 128 + 8


def char(ch, tf=TFW_VAR):
    cp = ord(ch)
    return bytes([REC_CHAR, tf, cp & 0xFF, (cp >> 8) & 0xFF])


def line(*records):
    return bytes([REC_LINE]) + b"".join(records) + bytes([REC_END])


def null_line():
    return bytes([REC_LINE | (1 << 4)])


def end():
    return bytes([REC_END])


def integral(body, lower=None, upper=None, sel=TM_SINT):
    """The template, then one slot holding the body and the display block."""
    # One bit, always: the real files say sel=21 var=2 whether or not there
    # are limits, because the limits are in the block and not in slots.  Two
    # bits would make it a three-slot template and the parser would read an
    # object list that is not there.
    var = VAR_UPPER
    slot = line(body)
    slot += bytes([REC_SUB])
    slot += line(lower) if lower is not None else null_line()
    slot += line(upper) if upper is not None else null_line()
    slot += bytes([REC_SYM]) + char("∫", TFW_SYM)
    return bytes([REC_TMPL, sel, var]) + end() + slot + end()


def test_the_integrand_is_the_body_and_not_a_superscript():
    data = HEADER + bytes([REC_FULL]) + line(integral(char("f"))) + end()
    out = equation.mtef_to_latex(data)
    assert B + "int" in out, out
    assert "^{" not in out, out          # the body used to land up here
    assert "f" in out, out


def test_the_limits_come_back():
    data = (HEADER + bytes([REC_FULL])
            + line(integral(char("f"), lower=char("a"), upper=char("b")))
            + end())
    out = equation.mtef_to_latex(data)
    assert "_{a}" in out.replace(" ", ""), out
    assert "^{b}" in out.replace(" ", ""), out


def test_an_integral_with_no_limits_gets_none_invented():
    data = HEADER + bytes([REC_FULL]) + line(integral(char("f"))) + end()
    out = equation.mtef_to_latex(data)
    assert "_{}" not in out.replace(" ", ""), out
    assert "^{}" not in out.replace(" ", ""), out


def test_one_nested_inside_another_is_split_too():
    r"""The inner integral is in the outer one's body, which is a slot -- so
    this only works if the passes walk into slots."""
    inner = integral(char("g"), lower=char("c"), upper=char("d"))
    data = (HEADER + bytes([REC_FULL])
            + line(integral(inner, lower=char("a"), upper=char("b")))
            + end())
    out = equation.mtef_to_latex(data).replace(" ", "")
    assert out.count(B + "int") == 2, out
    for limit in ("_{a}", "^{b}", "_{c}", "^{d}"):
        assert limit in out, (limit, out)


def test_an_integral_inside_a_bracket_is_reached():
    r"""A fence's content is a slot as well; the display block sat in it,
    untouched by every pass, until the pipeline walked in."""
    body = bytes([REC_TMPL, TM_PAREN, 0]) + end() \
        + line(integral(char("f"), lower=char("a"), upper=char("b")))
    data = HEADER + bytes([REC_FULL]) + line(body) + end()
    out = equation.mtef_to_latex(data)
    assert B + "scriptstyle" not in out, out
    assert "a" in out and "b" in out, out


def test_a_slot_that_is_not_that_shape_is_left_alone():
    r"""No size marker in the slot: nothing to split, and nothing is done."""
    data = (HEADER + bytes([REC_FULL])
            + line(bytes([REC_TMPL, TM_SINT, 0]) + end() + line(char("f")) + end())
            + end())
    out = equation.mtef_to_latex(data)     # must not throw
    assert "f" in out, out
