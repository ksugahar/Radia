"""MathML: the one clipboard format the whole of Office reads as maths.

Measured, by offering Word's own formats to PowerPoint one at a time: with only
`MathML` present PowerPoint produces a native equation; with only
`Rich Text Format` it produces a text box.  RTF is therefore Word's route and
MathML is PowerPoint's, and a single Copy carries both.

The spellings asserted here are Word's own output for the same equations.
"""

from __future__ import annotations

import pytest

equation = pytest.importorskip("radia.equation")

# LaTeX in, the fragments Word writes for it.
SPELLING = [
    (r"\frac{a}{b}",       ["<mfrac>", "<mi>a</mi>", "<mi>b</mi>"]),
    (r"a^{2}",             ["<msup>", "<mn>2</mn>"]),
    (r"a_{i}",             ["<msub>", "<mi>i</mi>"]),
    (r"a_{i}^{2}",         ["<msubsup>"]),
    (r"\sqrt{x}",          ["<msqrt>"]),
    (r"\left(x\right)",    ["<mo>(</mo>", "<mo>)</mo>"]),
    (r"\left[x\right]",    ["<mo>[</mo>", "<mo>]</mo>"]),
    (r"\sum_{i}^{n} a",    ["<munderover>", '<mo stretchy="false">∑</mo>']),
    (r"\int_{0}^{1} a",    ["<msubsup>", '<mo stretchy="false">∫</mo>']),
    (r"\oint_{C} a",       ['<mo stretchy="false">∮</mo>']),
    (r"\begin{pmatrix} a & b \\ c & d \end{pmatrix}",
                           ["<mtable>", "<mtr>", "<mtd>"]),
    (r"\hat{x}",           ['<mover accent="true">']),
    (r"\overline{x}",      ["<mover"]),
    (r"\sigma + \mu",      ["<mi>σ</mi>", "<mo>+</mo>", "<mi>μ</mi>"]),
]


@pytest.mark.parametrize("latex,expected", SPELLING, ids=[c[0] for c in SPELLING])
def test_spelling(latex, expected):
    mml = equation.tex_to_mathml(latex)
    for frag in expected:
        assert frag in mml, f"{frag!r} missing from the MathML for {latex}"


def test_root_takes_its_arguments_the_other_way_round():
    """<mroot> is body first, index second -- the reverse of the tree's order,
    and of every other output here."""
    mml = equation.tex_to_mathml(r"\sqrt[3]{x}")
    assert "<mroot>" in mml
    assert mml.index("<mi>x</mi>") < mml.index("<mn>3</mn>")


def test_runs_are_split_three_ways():
    """The split carries the spacing: <mo> gets operator spacing, <mi> is
    italicised as a variable, <mn> is a number."""
    mml = equation.tex_to_mathml("x+1")
    assert "<mi>x</mi>" in mml
    assert "<mo>+</mo>" in mml
    assert "<mn>1</mn>" in mml


def test_is_a_complete_math_element():
    mml = equation.tex_to_mathml("a+b")
    assert mml.startswith("<math ")
    assert mml.endswith("</math>")
    assert 'xmlns="http://www.w3.org/1998/Math/MathML"' in mml


def test_display_and_inline_are_distinguished():
    opts = equation.MathMLOptions()
    opts.display = False
    assert 'display="inline"' in equation.tex_to_mathml("a", opts)
    assert 'display="block"' in equation.tex_to_mathml("a")


@pytest.mark.parametrize("latex", [c[0] for c in SPELLING])
def test_every_output_describes_the_same_equation(latex):
    """OMML, RTF and MathML come from one walk.  A construct that appears in
    one and not the others means the walk was duplicated somewhere."""
    import re

    omml = equation.tex_to_omml(latex)
    rtf = equation.tex_to_rtf(latex)
    mml = equation.tex_to_mathml(latex)
    for omml_el, rtf_cw, mml_els in [
        ("m:f", "mf", ["<mfrac"]),
        ("m:rad", "mrad", ["<msqrt", "<mroot"]),
        ("m:nary", "mnary", ["<munderover", "<msubsup", "<munder", "<msub",
                             "<mover", "<msup", "stretchy"]),
        ("m:m", "mm", ["<mtable"]),
        ("m:acc", "macc", ["<mover"]),
    ]:
        in_omml = f"<{omml_el}>" in omml
        in_rtf = re.search(r"\{\\" + rtf_cw + r"(?![A-Za-z])", rtf) is not None
        assert in_omml == in_rtf, f"{latex}: {omml_el} differs between OMML and RTF"
        if in_omml:
            assert any(e in mml for e in mml_els), \
                f"{latex}: {omml_el} has no MathML counterpart"
