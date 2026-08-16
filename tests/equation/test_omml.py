"""LaTeX -> OMML: the path an equation takes into Word and PowerPoint.

Two properties are checked, and the second is the one a LaTeX round trip
structurally cannot see.  LaTeX binds `x` and `_{c}` implicitly, so a script
with an empty base round-trips clean; OMML has no implicit binding, and an
empty ``<m:e>`` renders as a dotted placeholder box next to the operator.  The
same applies to ``<m:nary>`` limits and bodies.
"""

from __future__ import annotations

import re

import pytest

equation = pytest.importorskip("radia.equation")

_EMPTY_SLOT = re.compile(r"<m:(e|num|den)></m:\1>")

# Everything a technical note actually contains.
EQUATIONS = [
    r"q_{0} = 0",
    r"H(x_{c})",
    r"\nabla H(x_{c})",
    r"\sigma_{f}(x) = \sum_{j} B_{fj}(x) a_{j}",
    r"M_{ij} = \sum_{f} \int_{S_{f}} \phi_{i}(d) B_{fj}(x) dS",
    r"\frac{a+b}{c}",
    r"\sqrt{x^{2}+y^{2}}",
    r"\sqrt[3]{x}",
    r"\left(\frac{1}{2}\right)",
    r"\left[\frac{a}{b}\right]",
    r"\int_{0}^{\infty} e^{-x} dx",
    r"\oint_{C} F \cdot dr",
    r"\sum_{i=1}^{n} a_{i}",
    r"\prod_{k} p_{k}",
    r"\overline{AB}",
    r"\hat{x}",
    r"\vec{B}",
    r"\mathbf{A}",
    r"\sin(\theta)+\cos(\theta)",
    r"\begin{pmatrix} a & b \\ c & d \end{pmatrix}",
    r"\begin{bmatrix} 1 & 0 \\ 0 & 1 \end{bmatrix}",
    r"\mu_{0} \nabla \times H = J",
    r"\frac{\mu_{0}}{4\pi} \oint_{C} \frac{I dl \times r}{|r|^{3}}",
]


@pytest.mark.parametrize("latex", EQUATIONS)
def test_produces_omml(latex):
    omml = equation.tex_to_omml(latex)
    assert "<m:oMath" in omml


@pytest.mark.parametrize("latex", EQUATIONS)
def test_no_empty_slot(latex):
    """An empty slot is not a crash -- it is a placeholder box on the slide."""
    omml = equation.tex_to_omml(latex)
    assert not _EMPTY_SLOT.search(omml), (
        f"{latex} leaves an empty OMML slot, which Office draws as a box")


def test_display_wraps_in_omath_para():
    opts = equation.OmmlOptions()
    opts.display = True
    assert "<m:oMathPara" in equation.tex_to_omml(r"a+b", opts)
    assert "<m:oMathPara" not in equation.tex_to_omml(r"a+b")


def test_variables_are_italic_and_operators_upright():
    omml = equation.tex_to_omml("x+1")
    # OMML marks style per run; digits and operators are upright, letters are
    # left alone so Office italicises them as variables.
    assert omml.count('<m:sty m:val="p"/>') == 2


def test_nary_carries_its_operator_and_limits():
    omml = equation.tex_to_omml(r"\oint_{C} F")
    assert 'm:chr m:val="\u222e"' in omml
    assert "<m:sub><m:r><m:t>C</m:t></m:r></m:sub>" in omml
    assert "<m:supHide" in omml           # no upper limit was given


def test_unparsable_command_does_not_lose_the_equation():
    """A typo costs one glyph, not the whole equation."""
    omml = equation.tex_to_omml(r"a + \notacommand + b")
    assert "<m:oMath" in omml
    assert "<m:t>a</m:t>" in omml
    assert "<m:t>b</m:t>" in omml
