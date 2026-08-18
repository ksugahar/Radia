r"""What the editor sets, measured against what TeX sets.

The rule for this editor is that appearance follows TeX and usability follows
Equation Editor 3.  This file is the appearance half of it, and it is a
comparison against numbers rather than against a picture of numbers: TeX will
state the width, height and depth of any box it makes, so the reference is
`\the\wd0` and not a screenshot.

The reference was produced by XeLaTeX with unicode-math and Latin Modern Math
-- the same OpenType file this renders with, read through the same MATH table.
That matters more than it sounds.  `\usepackage{lmodern}` looks like the same
typeface and is not: it sends TeX to the Type1 lmex10 extension font, whose
radicals and large operators are built from different pieces at different
sizes.  Comparing against it once reported an 18.5 % error in the radical that
belonged to the reference, and a "correction" made from that number had to be
undone.  The generator is validation_test/equation/tex_reference.tex.

The tolerances below are not a target to creep upward.  Where a case is exact
it is pinned exact, and the three that are not carry the reason.
"""

from __future__ import annotations

import pytest

equation = pytest.importorskip("radia.equation")

if not hasattr(equation, "tex_metrics"):
    pytest.skip("built without tex_metrics", allow_module_level=True)


def measure(latex):
    style = equation.SvgStyle()
    style.padding = 0.0
    return equation.tex_metrics(latex, style)


# (latex, TeX's width/height/depth in points at 12 pt, tolerance in points)
#
# 0.01 pt is "the same number": TeX rounds to scaled points and prints five
# figures, so the last digit is its own rounding, not a disagreement.
EXACT = 0.01

CASES = [
    # Plain text, which is where the fonts have to agree at all.
    ("x",                    (6.8640,  5.3040, 0.1320), EXACT),
    ("ab",                   (11.6640, 8.3280, 0.1320), EXACT),
    ("abc",                  (17.1600, 8.3280, 0.1320), EXACT),

    # Spacing between atom classes.
    ("a+b",                  (26.3339, 8.3280, 0.9960), EXACT),
    ("a=b",                  (27.6673, 8.3280, 0.1320), EXACT),
    ("a,b",                  (17.0002, 8.3280, 2.3160), EXACT),
    ("a(b)",                 (21.0000, 8.9760, 2.9760), 0.03),

    # Fractions: the rule is as wide as the parts, the box adds a null
    # delimiter on each side, and a fraction inside one steps down a size.
    (r"\frac{a}{b}",         (8.7480, 13.4289, 8.3649), EXACT),
    (r"\frac{abc}{d}",       (19.5600, 16.4529, 8.3649), EXACT),

    # Radicals: the variant is chosen by ink, and the slack is split the way
    # make_radical splits it.
    (r"\sqrt{\frac{a}{b}}",  (20.7480, 18.5400, 10.7400), EXACT),
    (r"\sqrt{2}",            (15.9960, 11.6040, 0.8760), 0.30),

    # Scripts sit at the shifts the font states; the widths still want the
    # font's `ssty` alternates, which are not wired up -- see below.
    ("x^{2}",                (12.1436, 9.9340, 0.1320), 0.60),
    ("x_{i}",                (10.7576, 5.3040, 3.0483), 0.60),
    ("x_{i}^{2}",            (12.1436, 9.9340, 3.2167), 0.60),
]


@pytest.mark.parametrize("latex,tex,tol", CASES)
def test_the_box_is_the_size_tex_makes_it(latex, tex, tol):
    got = measure(latex)
    for name, g, t in zip(("width", "height", "depth"), got, tex):
        assert abs(g - t) <= tol, (
            "%s of %r: TeX %.4f, ours %.4f (off by %+.4f, allowed %.2f)"
            % (name, latex, t, g, g - t, tol))


# ---------------------------------------------------------------------------
# Known open, with the cause named.  These are xfail rather than a loose
# tolerance so that fixing one is reported rather than passing silently.
# ---------------------------------------------------------------------------

OPEN = [
    # A font ships `ssty` alternates -- letters redrawn for script and
    # scriptscript sizes, slightly wider so they hold up small.  TeX applies
    # them; this does not, so every script is a few per cent narrow, and it
    # compounds when scripts nest.  Latin Modern Math has no ssty for "=",
    # which is why that one atom scales at exactly 0.7 and the letters do not.
    ("x^{y^{z}}",            (16.9372, 10.0580, 0.1320), 0.05, "ssty alternates"),

    # A large operator with limits at its side: TeX takes the italic
    # correction OUT of the operator's own width and uses it to shift the
    # superscript right (Appendix G, make_op).  Here it is still in the width,
    # so the operator is about 7 pt too wide before the limits are placed.
    (r"\int_{0}^{T}",        (19.3592, 19.9105, 12.1799), 0.05, "operator italic correction"),
    (r"\oint_{C}",           (12.7124, 16.3323, 12.1799), 0.05, "operator italic correction"),

    # A summation is not getting the display-size variant the font offers, and
    # an empty operand slot still reserves the editor's caret placeholder.
    (r"\sum_{n=1}^{N}",      (17.3280, 20.9771, 14.5049), 0.05, "display variant + empty-slot placeholder"),
]


@pytest.mark.parametrize("latex,tex,tol,why", OPEN)
def test_known_disagreements_with_tex(latex, tex, tol, why):
    pytest.xfail(why)
    got = measure(latex)
    for g, t in zip(got, tex):
        assert abs(g - t) <= tol
