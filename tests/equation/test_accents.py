"""Vectors and accents: \\vec{B}, \\hat{n}, \\dot{x}, \\overline{A}.

These were silently lost on the way back out to LaTeX.  The picture was right
and the Office paste was right, so nothing looked wrong -- but saving an
equation with a vector in it and reopening the file gave an equation without
one.  In electromagnetics that is most of the equations.

The cause is worth remembering: the emitter dropped a standalone embellishment
node on the reasoning that an embellishment is "usually attached to CHAR".
True of a tree read from MTEF, where Equation Editor hangs the accent off the
character; false of a tree from the LaTeX parser, which builds a node.  An
assumption that held for one producer and not the other, in a path the MTEF
corpus tests never reached.
"""

from __future__ import annotations

import pytest

equation = pytest.importorskip("radia.equation")
Equation = equation.Equation

# (written, read back) -- a few normalise to an equivalent spelling.
ROUND_TRIP = [
    (r"\vec{B}",        r"\vec{B}"),
    (r"\hat{n}",        r"\hat{n}"),
    (r"\tilde{f}",      r"\tilde{f}"),
    (r"\dot{x}",        r"\dot{x}"),
    (r"\ddot{x}",       r"\ddot{x}"),
    (r"\overline{A}",   r"\overline{A}"),
    (r"\bar{A}",        r"\overline{A}"),      # one bar in the tree
    (r"\mathbf{B}",     r"\mathbf{B}"),
    # A vector is bold ITALIC here, so \bm keeps its own name on the way
    # out; \mathbf is upright bold and stays a separate thing.
    (r"\bm{H}",         r"\bm{H}"),
    (r"\boldsymbol{J}", r"\bm{J}"),
]


@pytest.mark.parametrize("written,expected", ROUND_TRIP)
def test_an_accent_survives_being_written_back(written, expected):
    e = Equation()
    e.load_latex(written)
    assert e.latex() == expected


@pytest.mark.parametrize("written,expected", ROUND_TRIP)
def test_and_survives_a_second_pass(written, expected):
    """Reopening a saved file must not erode it further."""
    e = Equation()
    e.load_latex(written)
    once = e.latex()
    e2 = Equation()
    e2.load_latex(once)
    assert e2.latex() == expected


def test_an_accent_over_an_expression_keeps_the_expression():
    e = Equation()
    e.load_latex(r"\vec{AB}")
    assert e.latex() == r"\vec{AB}"


def test_a_vector_inside_a_fraction_survives():
    e = Equation()
    e.load_latex(r"\dfrac{\vec{B}}{\mu}")
    out = e.latex()
    assert r"\vec{B}" in out
    assert r"\mu" in out


def test_an_accent_reaches_office_as_an_accent():
    """It always did; this is here so a fix to the LaTeX side cannot break it."""
    xml = equation.tex_to_omml(r"\vec{B}")
    assert "<m:acc>" in xml
    assert "<m:t>B</m:t>" in xml


def test_an_accent_draws_something():
    """An accent whose base is present must not render as blank paper."""
    st = equation.SvgStyle()
    with_accent = equation.tex_to_png(r"\vec{B}", st, 1.0)
    plain = equation.tex_to_png("B", st, 1.0)
    assert len(with_accent) > len(plain)


# ---- what the lab actually writes -------------------------------------------

FIELD_EQUATIONS = [
    r"\nabla \times \vec{H} = \vec{J}",
    r"\nabla \cdot \vec{B} = 0",
    r"\vec{B} = \mu \vec{H}",
    r"\vec{F} = \vec{J} \times \vec{B}",
    r"\dfrac{\partial \vec{B}}{\partial t}",
]


@pytest.mark.parametrize("latex", FIELD_EQUATIONS)
def test_maxwell_survives_a_save_and_reopen(latex):
    e = Equation()
    e.load_latex(latex)
    once = e.latex()
    for vec in ("\\vec{H}", "\\vec{J}", "\\vec{B}", "\\vec{F}"):
        if vec in latex:
            assert vec in once, once
