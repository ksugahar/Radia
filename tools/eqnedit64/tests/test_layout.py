"""The layout engine, checked by invariants rather than by eye.

Layout defects are the kind that a person notices instantly and a test never
does, so the checks here are the numeric shadows of what went wrong in
practice:

  * a subscript that sank a third of a line below its base, because the glyph
    box came from the font's global descent instead of the glyph's own ink --
    Cambria Math reserves room for extensible brackets, Times does not, so a
    Greek base broke while a Latin one looked fine;
  * a parenthesis measured through the legacy Symbol code page, which put a
    gap after every "(" and let relations overlap a fraction bar;
  * a big operator whose operand was joined without the inter-atom space, so a
    sigma touched the symbol after it.

Each of those has a number attached to it below.

Run:  python tests\\test_layout.py
"""
from __future__ import annotations

import os
import re
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
sys.path.insert(0, os.path.join(_ROOT, "build"))

from eqnedit_core import SvgStyle, tex_to_svg, math_font_loaded  # noqa: E402

_VIEWBOX = re.compile(r'viewBox="0 0 ([\d.]+) ([\d.]+)"')
_TEXT = re.compile(r'<text x="([-\d.]+)" y="([-\d.]+)"[^>]*font-size="([\d.]+)"[^>]*>([^<]*)</text>')
_RECT = re.compile(r'<rect x="([-\d.]+)" y="([-\d.]+)" width="([\d.]+)" height="([\d.]+)"')
_PATHGLYPH = re.compile(
    r'<path data-char="([^"]*)" data-size="([\d.]+)" '
    r'transform="translate\(([-\d.]+),([-\d.]+)\)')


def fold_math_italic(text: str) -> str:
    """Math-italic code points back to the ASCII letters they stand for.

    A variable is drawn as its Unicode math-italic character, the way TeX
    sets it -- Latin Modern Math has no italic face, so asking GDI to slant
    the upright glyph would give the wrong shape.  Tests still name letters
    as "x" and "B".
    """
    out = []
    for ch in text:
        cp = ord(ch)
        if 0x1D44E <= cp <= 0x1D467:
            out.append(chr(ord("a") + cp - 0x1D44E))
        elif 0x1D434 <= cp <= 0x1D44D:
            out.append(chr(ord("A") + cp - 0x1D434))
        elif cp == 0x210E:                       # no MATHEMATICAL ITALIC H
            out.append("h")
        elif 0x1D6FC <= cp <= 0x1D714:
            out.append(chr(0x03B1 + cp - 0x1D6FC))
        else:
            out.append(ch)
    return "".join(out)


class Rendered:
    def __init__(self, latex: str, style=None):
        self.svg = tex_to_svg(latex, style or SvgStyle())
        m = _VIEWBOX.search(self.svg)
        self.width, self.height = float(m.group(1)), float(m.group(2))
        self.glyphs = [(float(x), float(y), float(s), t)
                       for x, y, s, t in _TEXT.findall(self.svg)]
        # Size variants are drawn as outlines and are glyphs all the same.
        self.glyphs += [(float(m.group(3)), float(m.group(4)),
                         float(m.group(2)), m.group(1))
                        for m in _PATHGLYPH.finditer(self.svg)]
        self.glyphs.sort(key=lambda g: g[0])
        self.rects = [tuple(float(v) for v in r) for r in _RECT.findall(self.svg)]

    def glyph(self, ch: str):
        for g in self.glyphs:
            if fold_math_italic(g[3]) == ch:
                return g
        raise AssertionError(f"{ch!r} was not drawn in: {self.svg[:200]}")


def approx(a: float, b: float, tol: float) -> bool:
    return abs(a - b) <= tol


def main() -> int:
    # Every number below is a font measurement.  When the embedded math font
    # does not load, GDI has nothing to measure with and every glyph comes
    # back zero wide -- which reads as three unrelated layout failures rather
    # than as one missing font.  Say which it is.
    if not math_font_loaded():
        print("FAIL  the math font did not load, so every measurement here "
              "would be of a substituted face")
        return 1

    failures = []
    style = SvgStyle()          # full 12, sub 7, sym 18

    # --- nothing may be drawn outside the box the SVG declares -------------
    for latex in [r"\sigma_{f}", r"\sum_{i=1}^{n} a_{i}", r"\frac{a}{b}",
                  r"\sqrt{x^{2}}", r"\left(\frac{1}{2}\right)",
                  r"\int_{0}^{\infty} e^{-x} dx", r"\nabla \times H = J"]:
        r = Rendered(latex, style)
        for x, y, size, text in r.glyphs:
            if x < -0.01 or x > r.width + 0.01 or y < -0.01 or y > r.height + 0.01:
                failures.append(f"{latex}: glyph {text!r} at ({x},{y}) is outside "
                                f"the {r.width}x{r.height} box")
        for x, y, w, h in r.rects:
            if x < -0.01 or x + w > r.width + 0.01:
                failures.append(f"{latex}: rule at x={x} w={w} overflows the box")

    # Every palette structure must remain visible in the native renderer.
    for latex, required in [
        (r"\begin{bmatrix}a&b\\c&d\end{bmatrix}", "abcd"),
        (r"\begin{cases}x&x>0\\-x&x\leq0\end{cases}", "xx0"),
        (r"\hat{x}+\vec{B}+\overline{y}", "xBy"),
        (r"\operatorname{curl} H", "curlH"),
    ]:
        r = Rendered(latex, style)
        shown = fold_math_italic("".join(g[3] for g in r.glyphs))
        for ch in required:
            if ch not in shown:
                failures.append(f"{latex}: {ch!r} vanished from the renderer ({shown!r})")

    # --- a subscript sits just below the baseline, not a line below --------
    # The regression: the font's global descent was used as the base's depth.
    # The base is drawn first, the subscript last, so their baselines can be
    # compared without knowing which glyph a LaTeX name produced.
    for base, sub in [(r"\sigma", "f"), ("B", "j"), (r"\mu", "0"), ("x", "i")]:
        r = Rendered(f"{base}_{{{sub}}}", style)
        drop = r.glyphs[-1][1] - r.glyphs[0][1]
        if not 0.05 * style.full <= drop <= 0.45 * style.full:
            failures.append(f"{base}_{{{sub}}}: subscript dropped {drop:.2f}pt, "
                            f"expected between {0.05*style.full:.2f} and "
                            f"{0.45*style.full:.2f}")

    # --- parentheses cost about what a parenthesis costs -------------------
    # The regression: measured through SYMBOL_CHARSET, "(" came back with a
    # wildly wrong advance.
    w_x = Rendered("x", style).width
    w_px = Rendered("(x)", style).width
    per_paren = (w_px - w_x) / 2.0
    if not 0.20 * style.full <= per_paren <= 0.55 * style.full:
        failures.append(f"a parenthesis measures {per_paren:.2f}pt at "
                        f"{style.full}pt, which is not a parenthesis")

    # --- TeX inter-atom spacing is actually applied ------------------------
    mu = style.full / 18.0
    w_ab = Rendered("ab", style).width
    w_apb = Rendered("a+b", style).width
    plus = Rendered("+", style).width - 2 * style.padding
    gained = w_apb - w_ab - plus
    if not approx(gained, 2 * 4 * mu, 1.5 * mu):
        failures.append(f"a+b gained {gained:.2f}pt of space, expected "
                        f"{2*4*mu:.2f}pt (two medium spaces)")

    # A leading minus is unary: no space at all around it.
    w_mx = Rendered("-x", style).width
    w_minus = Rendered("{-}", style).width - 2 * style.padding
    unary_gap = w_mx - w_x - w_minus
    if abs(unary_gap) > 1.5 * mu:
        failures.append(f"a leading minus was spaced as a binary operator "
                        f"({unary_gap:.2f}pt of extra space)")

    # A relation gets more room than a binary operator.
    w_aeqb = Rendered("a=b", style).width
    if not w_aeqb - w_apb > 0:
        failures.append("a relation is not spaced more widely than a binary "
                        "operator")

    # --- a big operator does not touch its operand -------------------------
    # The operand is joined inside the operator's own layout, so the space has
    # to be applied there rather than by the atom loop.
    r = Rendered(r"\sum_{j} B", style)
    sigma = r.glyph("∑")
    b = r.glyph("B")
    if b[0] - sigma[0] < 0.6 * sigma[2]:
        failures.append(f"the operand starts {b[0]-sigma[0]:.2f}pt after a "
                        f"{sigma[2]:.0f}pt operator -- they overlap")

    # --- integral limits follow TeX's display-integral script positions ---
    r = Rendered(r"\int_{1}^{2}x", style)
    main = r.glyph("x")
    upper = r.glyph("2")
    lower = r.glyph("1")
    upper_rise = main[1] - upper[1]
    lower_drop = lower[1] - main[1]
    upper_correction = upper[0] - lower[0]
    if not 0.95 * style.full <= upper_rise <= 1.20 * style.full:
        failures.append(
            f"integral upper limit rose {upper_rise:.2f}pt; TeX reference is "
            f"about {1.10*style.full:.2f}pt")
    if not 0.75 * style.full <= lower_drop <= 1.05 * style.full:
        failures.append(
            f"integral lower limit dropped {lower_drop:.2f}pt; TeX reference is "
            f"about {0.90*style.full:.2f}pt")
    if not 0.30 * style.full <= upper_correction <= 0.55 * style.full:
        failures.append(
            f"integral upper-limit italic correction is {upper_correction:.2f}pt; "
            f"TeX reference is about {0.44*style.full:.2f}pt")

    # --- TeX's style chain sets nested parts smaller -----------------------
    # A fraction sets its parts one style down, so the size drops once the
    # chain reaches script: pdfLaTeX gives 12, 8 and 6 pt as fractions nest.
    # Every part stayed at 12 pt here, which made a nested fraction a third
    # taller than TeX's and dragged the fences around it with it.  The
    # expected sizes are pdfLaTeX's own, measured with tools/tex_geometry.py.
    for latex, expected in ((r"\frac{a}{b}", [12]),
                            (r"\frac{\frac{a}{b}}{c}", [8, 12]),
                            (r"\frac{\frac{\frac{a}{b}}{c}}{d}", [6, 8, 12]),
                            (r"x^{\frac{a}{b}}", [6, 12]),
                            (r"\frac{a}{b}^{2}", [8, 12]),
                            (r"\sqrt{\frac{a}{b}}", [12])):
        sizes = sorted({round(g[2]) for g in Rendered(latex, style).glyphs})
        if sizes != expected:
            failures.append(
                f"{latex}: glyph sizes {sizes}, pdfLaTeX uses {expected} -- "
                f"the style chain is not being followed")

    # --- a fraction's rule spans both parts --------------------------------
    r = Rendered(r"\frac{a+b}{c}", style)
    if not r.rects:
        failures.append("a fraction drew no rule")
    else:
        bar = r.rects[0]
        for x, y, size, text in r.glyphs:
            if x < bar[0] - 0.01 or x > bar[0] + bar[2] + 0.01:
                failures.append(f"{text!r} sits outside the fraction rule")

    # A bra-ket assembled from plain characters gives itself away with a
    # one-line bar and one-line angle brackets around a fraction.  It has to
    # grow with its content exactly as \left...\right does, and stay flat
    # when the content is flat.
    # Measured as the height the delimiter ends up, not as a scale factor:
    # fences take the font's designed larger glyph now and only stretch when
    # even the largest falls short, so looking for a transform asked about
    # the mechanism rather than the result and failed on a fence that had in
    # fact grown correctly.
    def fence_height(latex):
        r = Rendered(latex, style)
        return r.height

    tall = fence_height(r"\left\langle \frac{p}{q} \middle| \psi \right\rangle")
    flat = fence_height(r"\left\langle a \middle| b \right\rangle")
    fence = fence_height(r"\left( \frac{p}{q} \right)")
    plain = fence_height(r"\left( a \right)")
    if not tall > flat * 1.4:
        failures.append(
            f"a bra-ket around a fraction is {tall:.1f}pt tall against "
            f"{flat:.1f}pt around plain atoms -- it did not grow")
    if not fence > plain * 1.4:
        failures.append(
            f"a fence around a fraction is {fence:.1f}pt tall against "
            f"{plain:.1f}pt around a letter -- it did not grow")
    if abs(tall - fence) > 0.51:
        failures.append(
            f"bra-ket grew to {tall:.2f}pt, an ordinary fence to {fence:.2f}pt "
            f"around the same fraction")

    # cases left-aligns its columns; matrix centres them.  Sharing one layout
    # made the canvas disagree with the PDF it was about to produce.
    # The two cells of the first column have very different widths, so a
    # left-aligned column starts them at the same x and a centred one does
    # not.  Compare the glyphs themselves: the brace a cases environment
    # draws is not part of either row.
    def glyph_x(latex, wanted):
        xs = [x for x, y, size, text in Rendered(latex).glyphs
              if fold_math_italic(text) == wanted]
        return min(xs) if xs else None

    body = r"a & bbbbbb \\ cccccc & d"
    for env, aligned in ((r"cases", True), (r"matrix", False)):
        latex = r"\begin{%s}%s\end{%s}" % (env, body, env)
        top, bottom = glyph_x(latex, "a"), glyph_x(latex, "c")
        if top is None or bottom is None:
            failures.append(f"{env}: could not find both column-one cells")
        elif aligned and abs(top - bottom) > 0.01:
            failures.append(
                f"cases rows do not share a left edge: {top} vs {bottom}")
        elif not aligned and abs(top - bottom) <= 0.01:
            failures.append(
                f"matrix rows should be centred, not aligned: {top}")

    if failures:
        print(f"FAIL  {len(failures)}")
        for f in failures:
            print("  " + f)
        return 1
    print("ok    layout: boxes, script depth, parens, TeX spacing, operators, "
          "bra-ket stretch")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
