r"""Reading a whole MTEF stream, not just its first object list.

Equation Editor does not wrap a document in one list: it closes a list with
END and carries on.  The editor's reader took the first list and stopped, so a
real 1805-byte lecture file arrived as 86 bytes of it -- seven items, drawn
without complaint, and looking like a short equation rather than a truncated
one.

The streams here are built byte by byte, because that is the only way to state
the case without shipping someone's lecture: a header, a character, the END
that used to finish the job, and a character after it.
"""

from __future__ import annotations

import pytest

equation = pytest.importorskip("radia.equation")

B = chr(92)

# MTEF v3 header: version, platform, product, product major, product minor.
HEADER = bytes([3, 1, 1, 3, 10])

REC_END = 0
REC_CHAR = 2
TFW_VAR = 128 + 3          # TF_VARIABLE with a 16-bit code


def char(ch):
    """One CHAR record: tag, typeface, then the code little-endian."""
    cp = ord(ch)
    return bytes([REC_CHAR, TFW_VAR, cp & 0xFF, (cp >> 8) & 0xFF])


def end():
    return bytes([REC_END])


def test_one_list_reads():
    data = HEADER + char("a") + end()
    assert "a" in equation.mtef_to_latex(data)


def test_what_follows_the_end_is_not_thrown_away():
    """This is the whole bug: b lived past an END and was never read."""
    data = HEADER + char("a") + end() + char("b") + end()
    out = equation.mtef_to_latex(data)
    assert "a" in out
    assert "b" in out, "everything after the first END was dropped: " + out


def test_several_lists_in_a_row():
    data = HEADER + b"".join(char(c) + end() for c in "abcde")
    out = equation.mtef_to_latex(data)
    for c in "abcde":
        assert c in out, (c, out)


def test_a_truncated_stream_stops_rather_than_spinning():
    """A record cut off mid-way must end the read, not loop on no progress."""
    data = HEADER + char("a") + end() + bytes([REC_CHAR, TFW_VAR])
    out = equation.mtef_to_latex(data)          # must return, not hang
    assert "a" in out


def test_the_editor_route_and_the_legacy_route_are_both_available():
    r"""mtef_to_tex answers with the legacy standalone converter; mtef_to_latex
    answers with the editor's own parser and emitter.  Which reading you get
    is now something the caller states rather than something they inherit."""
    data = equation.tex_to_mtef(B + "frac{a}{b}")
    assert "a" in equation.mtef_to_tex(data)
    assert "a" in equation.mtef_to_latex(data)
