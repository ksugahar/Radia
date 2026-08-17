"""The layout engine, checked by invariants rather than by eye.

Layout defects are the kind a person notices instantly and a test never does,
so the checks here are the numeric shadows of what actually went wrong:

  * a subscript that sank a third of a line below its base, because the glyph
    box came from the font's global descent instead of the glyph's own ink --
    Cambria Math reserves room for extensible brackets, Times does not, so a
    Greek base broke while a Latin one looked fine;
  * a parenthesis measured through the legacy Symbol code page, which put a gap
    after every "(" and let relations overlap a fraction bar;
  * a big operator whose operand was joined without the inter-atom space, so a
    sigma touched the symbol after it.

The layout is also what the editor draws with and positions its caret from, so
these numbers matter beyond the SVG they are read out of here.
"""

from __future__ import annotations

import re

import pytest

equation = pytest.importorskip("radia.equation")

_VIEWBOX = re.compile(r'viewBox="0 0 ([\d.]+) ([\d.]+)"')
_TEXT = re.compile(r'<text x="([-\d.]+)" y="([-\d.]+)"[^>]*'
                   r'font-size="([\d.]+)"[^>]*>([^<]*)</text>')
_RECT = re.compile(r'<rect x="([-\d.]+)" y="([-\d.]+)" '
                   r'width="([\d.]+)" height="([\d.]+)"')


class Rendered:
    def __init__(self, latex, style=None):
        self.style = style or equation.SvgStyle()
        self.svg = equation.tex_to_svg(latex, self.style)
        m = _VIEWBOX.search(self.svg)
        self.width, self.height = float(m.group(1)), float(m.group(2))
        self.glyphs = [(float(x), float(y), float(s), t)
                       for x, y, s, t in _TEXT.findall(self.svg)]
        self.rects = [tuple(float(v) for v in r) for r in _RECT.findall(self.svg)]

    def glyph(self, ch):
        for g in self.glyphs:
            if g[3] == ch:
                return g
        raise AssertionError(f"{ch!r} was not drawn")

    @property
    def ink_width(self):
        return self.width - 2 * self.style.padding


@pytest.mark.parametrize("latex", [
    r"\sigma_{f}", r"\sum_{i=1}^{n} a_{i}", r"\frac{a}{b}", r"\sqrt{x^{2}}",
    r"\left(\frac{1}{2}\right)", r"\int_{0}^{\infty} e^{-x} dx",
    r"\nabla \times H = J",
])
def test_nothing_is_drawn_outside_the_declared_box(latex):
    r = Rendered(latex)
    for x, y, _size, text in r.glyphs:
        assert -0.01 <= x <= r.width + 0.01, f"{text!r} at x={x}"
        assert -0.01 <= y <= r.height + 0.01, f"{text!r} at y={y}"
    for x, _y, w, _h in r.rects:
        assert -0.01 <= x and x + w <= r.width + 0.01


@pytest.mark.parametrize("base,sub", [
    (r"\sigma", "f"), ("B", "j"), (r"\mu", "0"), ("x", "i"),
])
def test_a_subscript_sits_just_below_the_baseline(base, sub):
    """Not a line below it: the regression used the font's global descent."""
    r = Rendered(f"{base}_{{{sub}}}")
    drop = r.glyphs[-1][1] - r.glyphs[0][1]      # base is drawn first
    assert 0.05 * r.style.full <= drop <= 0.45 * r.style.full


def test_a_parenthesis_measures_like_a_parenthesis():
    """The regression measured it through SYMBOL_CHARSET."""
    style = equation.SvgStyle()
    per_paren = (Rendered("(x)", style).width - Rendered("x", style).width) / 2.0
    assert 0.20 * style.full <= per_paren <= 0.55 * style.full


def test_a_binary_operator_gets_two_medium_spaces():
    """Four eighteenths on each side.  Equation Editor 3.1 measures 0.156 em
    there, which is three; appearance follows TeX."""
    style = equation.SvgStyle()
    mu = style.full / 18.0
    plus = Rendered("+", style).ink_width
    gained = Rendered("a+b", style).width - Rendered("ab", style).width - plus
    assert abs(gained - 2 * 4 * mu) <= 1.5 * mu


def test_a_leading_minus_is_unary():
    """TeX's rule: a binary operator with nothing to bind on its left is not
    binary, so `-x` opens without a gap."""
    style = equation.SvgStyle()
    mu = style.full / 18.0
    minus = Rendered("{-}", style).ink_width
    gap = Rendered("-x", style).width - Rendered("x", style).width - minus
    assert abs(gap) <= 1.5 * mu


def test_a_relation_is_spaced_more_widely_than_a_binary_operator():
    style = equation.SvgStyle()
    assert Rendered("a=b", style).width > Rendered("a+b", style).width


def test_a_big_operator_does_not_touch_its_operand():
    """The operand is joined inside the operator's own layout, so the space has
    to be applied there rather than by the atom loop."""
    r = Rendered(r"\sum_{j} B")
    sigma, b = r.glyph("∑"), r.glyph("B")
    assert b[0] - sigma[0] >= 0.6 * sigma[2]


def test_a_fraction_rule_spans_both_parts():
    r = Rendered(r"\frac{a+b}{c}")
    assert r.rects, "a fraction drew no rule"
    bar = r.rects[0]
    for x, _y, _size, text in r.glyphs:
        assert bar[0] - 0.01 <= x <= bar[0] + bar[2] + 0.01, f"{text!r} outside"


def test_type_sizes_come_from_the_style():
    small, big = equation.SvgStyle(), equation.SvgStyle()
    big.full, big.sub, big.sym = 24.0, 14.0, 36.0
    assert Rendered("x", big).width > Rendered("x", small).width
