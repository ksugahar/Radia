r"""A matrix's cells: how many bytes the header takes, and where the rows go.

Two faults, one shape.

The row and column partition arrays are TWO BITS per line, packed.  Read as
one byte per line, a 2x1 matrix consumed five bytes where the file has two,
and the three it took past them were its own first cell -- so the cells came
out empty, their contents stayed outside, and characters went missing: the
x of x_1 was eaten and the vector read

    [ _{1} \\ x_{2} ]

Then, with the cells read correctly, some documents turn out to put EVERY row
in the first cell and leave the others empty.  Handing the lines back out is
the second half.

Both are pinned here from bytes, because both are byte-level and the number
five-against-two is the whole argument.
"""

from __future__ import annotations

import pytest

equation = pytest.importorskip("radia.equation")

B = chr(92)

HEADER = bytes([3, 1, 1, 3, 10])

REC_END, REC_LINE, REC_CHAR = 0, 1, 2
REC_MATRIX, REC_FULL = 5, 10

TFW_VAR = 128 + 3


def char(ch):
    cp = ord(ch)
    return bytes([REC_CHAR, TFW_VAR, cp & 0xFF, (cp >> 8) & 0xFF])


def line(*records):
    return bytes([REC_LINE]) + b"".join(records) + bytes([REC_END])


def end():
    return bytes([REC_END])


def parts(n):
    """(n + 1) partition lines at two bits each, packed into bytes."""
    return bytes([0]) * (((n + 1) * 2 + 7) // 8)


def matrix(rows, cols, cells):
    """valign, hjust, vjust, rows, cols, partitions, then one list per cell."""
    out = bytes([REC_MATRIX, 1, 1, 1, rows, cols]) + parts(rows) + parts(cols)
    for c in cells:
        out += c + end()
    return out


def test_a_two_by_one_matrix_keeps_both_cells():
    data = (HEADER + bytes([REC_FULL])
            + line(matrix(2, 1, [line(char("a")), line(char("b"))]))
            + end())
    out = equation.mtef_to_latex(data).replace(" ", "").replace("\n", "")
    assert "a" + B + B + "b" in out, out


def test_the_first_character_of_the_first_cell_survives():
    r"""The old header arithmetic ate it: x_1 arrived as _1."""
    data = (HEADER + bytes([REC_FULL])
            + line(matrix(2, 1, [line(char("x"), char("1")), line(char("y"))]))
            + end())
    out = equation.mtef_to_latex(data)
    assert "x" in out, out
    assert "1" in out, out


def test_every_row_in_the_first_cell_is_handed_back_out():
    r"""Some documents write all the rows into cell one and leave the rest
    empty; the column then came out as a single run."""
    crammed = line(char("a")) + line(char("b"))
    data = (HEADER + bytes([REC_FULL])
            + line(matrix(2, 1, [crammed, b""]))
            + end())
    out = equation.mtef_to_latex(data).replace(" ", "").replace("\n", "")
    assert "a" + B + B + "b" in out, out


def body_of(out):
    """What is between \\begin{matrix} and \\end{matrix}."""
    out = out.replace(" ", "").replace("\n", "")
    i = out.index("{matrix}") + len("{matrix}")
    j = out.index(B + "end{matrix}")
    return out[i:j]


def test_a_full_matrix_is_not_redistributed():
    data = (HEADER + bytes([REC_FULL])
            + line(matrix(2, 1, [line(char("a")), line(char("b"))]))
            + end())
    # counting letters would count the ones in \\begin{matrix}
    assert body_of(equation.mtef_to_latex(data)) == "a" + B + B + "b"


def test_a_cell_without_line_boundaries_is_left_alone():
    r"""Nothing to distribute on: two characters in one cell could be one row
    or two, and guessing would be worse than the run it already is."""
    crammed = line(char("a"), char("b"))
    data = (HEADER + bytes([REC_FULL])
            + line(matrix(2, 1, [crammed, b""]))
            + end())
    out = equation.mtef_to_latex(data)      # must not throw, must not guess
    assert "a" in out and "b" in out, out
