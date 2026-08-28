"""Look at the pixels the canvas actually draws.

Every other check in this suite reads numbers.  The defects that got reported
here were seen, not measured: a vinculum that did not meet its radical, limits
sitting away from their operator.  Geometry tests missed them because the
numbers were self-consistent, and the SVG export could not stand in for the
canvas -- it goes through different code, and a viewer that lacks Cambria Math
substitutes a differently shaped radical, which is how one "fixed" render was
read as still broken.

So this renders through the real GDI path, via `Eqnedit64.exe --render-png`,
and asserts properties of the ink that do not depend on the font: a rule that
continues a glyph must not have a gap in it, a fraction bar must reach across
both of its parts, and nothing may be clipped at the edge of the image.

Run:  python tests\\test_render_image.py
"""
from __future__ import annotations

import os
import subprocess
import sys
import tempfile

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)

try:
    import fitz  # PyMuPDF, used only to decode the PNG
except ImportError:  # pragma: no cover
    print("skip  PyMuPDF is not installed; cannot decode the rendered PNG")
    raise SystemExit(0)

INK = 128          # a pixel darker than this counts as ink
# build/ first, deliberately.  dist/ is a copy made after the link step and
# it can be stale -- a running instance locks it, and the copy then fails
# while the build still reports success.  Testing the stale copy is how a
# regression once looked like a pass.  `dist == build` is checked separately
# by the release gate, which is where that belongs.
EXE_CANDIDATES = [os.path.join(_ROOT, "build", "Eqnedit64.exe"),
                  os.path.join(_ROOT, "dist", "Eqnedit64.exe")]


def executable():
    for path in EXE_CANDIDATES:
        if os.path.exists(path):
            return path
    return None


class Render:
    """The rendered equation as a grid of ink flags."""

    def __init__(self, exe, latex, workdir):
        self.latex = latex
        png = os.path.join(workdir, "render.png")
        if os.path.exists(png):
            os.remove(png)
        result = subprocess.run([exe, "--render-png", latex, png],
                                capture_output=True, text=True)
        if result.returncode != 0 or not os.path.exists(png):
            raise RuntimeError("render failed for %r (exit %d)"
                               % (latex, result.returncode))
        pix = fitz.Pixmap(png)
        self.w, self.h, n = pix.width, pix.height, pix.n
        data = pix.samples
        self.rows = []
        for y in range(self.h):
            base = y * self.w * n
            self.rows.append([data[base + x * n] < INK for x in range(self.w)])

    def runs(self, y):
        """Contiguous ink spans in one row, as (first, last)."""
        out, start = [], None
        for x, on in enumerate(self.rows[y]):
            if on and start is None:
                start = x
            elif not on and start is not None:
                out.append((start, x - 1))
                start = None
        if start is not None:
            out.append((start, self.w - 1))
        return out

    def first_ink_row(self):
        for y in range(self.h):
            if any(self.rows[y]):
                return y
        return None

    def ink_columns(self, y0, y1):
        cols = set()
        for y in range(y0, y1 + 1):
            for x, on in enumerate(self.rows[y]):
                if on:
                    cols.add(x)
        return cols


def main() -> int:
    exe = executable()
    if exe is None:
        print("skip  Eqnedit64.exe has not been built")
        return 0

    failures = []
    with tempfile.TemporaryDirectory() as work:
        def render(latex):
            return Render(exe, latex, work)

        # A radical's vinculum continues the flag of the sign.  Drawn at the
        # advance width instead of the flag's end, or at the nominal rule
        # thickness instead of the flag's, it detaches -- which is what was
        # reported.  The top row of ink must therefore be one unbroken run.
        for latex in (r"\sqrt{x}", r"\sqrt{x^{2}}", r"\sqrt{\frac{a}{b}}",
                      r"\sqrt{x_{i}^{2}}"):
            try:
                r = render(latex)
            except RuntimeError as exc:
                failures.append(str(exc))
                continue
            top = r.first_ink_row()
            if top is None:
                failures.append("%s: drew nothing" % latex)
                continue
            runs = r.runs(top)
            if len(runs) != 1:
                failures.append(
                    "%s: the top of the radical is in %d pieces %s -- the "
                    "vinculum does not meet the sign"
                    % (latex, len(runs), runs))

        # The vinculum is one thickness for its whole length.  Cambria's surd
        # ends in a flat flag thicker than TeX's rule, and drawing the glyph
        # whole left that flag protruding under the left end of the bar: a
        # short deeper stub that reads as a chipped bar.  Measure the depth of
        # the topmost ink run column by column across the bar; it has to be
        # constant.  The geometry checks cannot see this -- they read where
        # the rule is, not what the glyph draws next to it.
        # `x+` in front of two of these is not decoration: it moves the
        # radical off the origin, which is the only way the clip-does-not-
        # travel-with-the-glyph defect shows up.  Every radical case here
        # used to start at x=0, so the sign could vanish everywhere else
        # while these all passed.
        for latex in (r"\sqrt{x}", r"\sqrt{\left\langle a \right\rangle }",
                      r"\sqrt{\frac{a}{b}}", r"x+\sqrt{a}",
                      r"x+\sqrt{a^{2}}"):
            try:
                r = render(latex)
            except RuntimeError as exc:
                failures.append(str(exc))
                continue
            top = r.first_ink_row()
            if top is None:
                continue
            first, last = r.runs(top)[0]
            depths = []
            for x in range(first, last + 1):
                depth = 0
                for y in range(top, r.h):
                    if not r.rows[y][x]:
                        break
                    depth += 1
                depths.append((x, depth))
            # The profile runs: the sign's diagonal (much deeper than the
            # rule), then the rule.  A flag protruding under the bar shows up
            # as a plateau in between -- deeper than the rule but nothing like
            # the diagonal.  The defect drew nine such columns; a clean join
            # has none, give or take an antialiased one.
            bar = min(d for _, d in depths)
            peak = max(d for _, d in depths)

            # The sign itself has to be drawn.  Clipping the glyph to hide its
            # flag once removed the whole radical: the clip boundary did not
            # move with the glyph, so anywhere but the origin it swallowed
            # everything and left a bare bar floating over the content.  The
            # top run then has one uniform depth -- which the flag check below
            # is perfectly happy with.  The sign's peak joins the rule at the
            # top, so the deepest column of that run must be far deeper than
            # the rule; with no sign, peak and rule are the same number.
            # A multiple of the rule would be wrong here: how much deeper the
            # peak is depends on how heavy the font's radical is, and Latin
            # Modern's is light enough that 9 px meets a 5 px rule.  What
            # says "no sign" is the run being uniform -- peak and rule the
            # same number -- so a small absolute margin is the right test.
            if peak < bar + 3:
                failures.append(
                    "%s: the top run is %d px deep everywhere and the rule is "
                    "%d px -- the radical sign is not being drawn"
                    % (latex, peak, bar))

            between = [x for x, d in depths if bar < d < (bar + peak) / 2.0]
            if len(between) > 2:
                failures.append(
                    "%s: %d columns %s sit between the %d px rule and the "
                    "%d px diagonal -- the sign's flag is protruding under "
                    "the left end of the bar"
                    % (latex, len(between),
                       "%d..%d" % (between[0], between[-1]), bar, peak))

        # A display operator straddles the baseline: TeX centres it on the
        # math axis.  Setting it on the baseline instead stood the integral
        # about 0.24 em too high and its top curl ran into the upper limit --
        # while every geometry probe still passed, because they measured
        # where the limits were, not where the sign was.
        #
        # This cannot be probed from glyph origins: a big operator's origin is
        # a font's own choice, so cmex and Latin Modern disagree by 1.6 em
        # while their ink agrees.  Measured as ink, from a 12 pt pdfLaTeX
        # \displaystyle run: 61.5% of the integral and 68.3% of the sum sit
        # above the baseline.  `x` shares the baseline and gives it away.
        for latex, share, name in ((r"\int x", 61.5, "integral"),
                                   (r"\sum x", 68.3, "sum")):
            try:
                r = render(latex)
            except RuntimeError as exc:
                failures.append(str(exc))
                continue
            # split into ink groups by blank columns: first is the sign, last
            # is the letter, whose ink bottom is the baseline
            spans, run = [], []
            for x in range(r.w):
                ink = [y for y in range(r.h) if r.rows[y][x]]
                if ink:
                    run.append((min(ink), max(ink)))
                elif run:
                    spans.append(run)
                    run = []
            if run:
                spans.append(run)
            if len(spans) < 2:
                failures.append("%s: %d ink group(s); cannot find the baseline"
                                % (latex, len(spans)))
                continue
            sign, letter = spans[0], spans[-1]
            sign_top = min(t for t, _ in sign)
            sign_bottom = max(b for _, b in sign)
            baseline = max(b for _, b in letter)
            height = float(sign_bottom - sign_top)
            above = 100.0 * (baseline - sign_top) / height
            if abs(above - share) > 2.0:
                failures.append(
                    "%s: %.1f%% of the %s sign is above the baseline, TeX puts "
                    "%.1f%% -- it is not centred on the math axis"
                    % (latex, above, name, share))

        # A display operator is the font's designed larger glyph, not the text
        # one scaled up.  Scaling widens a sign as much as it heightens it:
        # the base integral is 0.498 wide per unit of height where TeX's
        # display one is 0.400, and that extra width is what a reader sees as
        # the limits overlapping the sign.  Measured as ink, since it is a
        # property of the shape and not of any coordinate.
        for latex, ratio, name in ((r"\int", 0.400, "integral"),
                                   (r"\sum", 0.950, "sum"),
                                   (r"\prod", 0.832, "product")):
            try:
                r = render(latex)
            except RuntimeError as exc:
                failures.append(str(exc))
                continue
            cols = [x for x in range(r.w) if any(r.rows[y][x] for y in range(r.h))]
            rows = [y for y in range(r.h) if any(r.rows[y])]
            if not cols or not rows:
                failures.append("%s: drew nothing" % latex)
                continue
            got = float(cols[-1] - cols[0] + 1) / float(rows[-1] - rows[0] + 1)
            if abs(got / ratio - 1.0) > 0.06:
                failures.append(
                    "%s: the %s sign is %.3f wide per unit of height, TeX's is "
                    "%.3f (%+.0f%%) -- it is the base glyph scaled up rather "
                    "than the font's designed size"
                    % (latex, name, got, ratio, 100.0 * (got / ratio - 1.0)))
        # A fraction rule reaches across both parts.  Its row is the widest
        # single run in the image.
        for latex in (r"\frac{a}{b}", r"\frac{a+b}{c}", r"\frac{1}{x^{2}}"):
            try:
                r = render(latex)
            except RuntimeError as exc:
                failures.append(str(exc))
                continue
            widest, at = None, None
            for y in range(r.h):
                for a, b in r.runs(y):
                    if widest is None or b - a > widest[1] - widest[0]:
                        widest, at = (a, b), y
            if widest is None:
                failures.append("%s: drew nothing" % latex)
                continue
            above = r.ink_columns(0, max(0, at - 2))
            below = r.ink_columns(min(r.h - 1, at + 2), r.h - 1)
            # pdfLaTeX's rule contains its parts' ink with room to spare
            # (measured: rule 60..154, denominator 63..146).  The canvas is
            # built from advance widths, and a glyph whose ink starts left of
            # its origin -- italic x is one -- reaches a little past the rule.
            # Closing this needs the layout to track ink on the left as well
            # as on the right; until then the allowance is the size of the
            # effect measured at 300 dpi, not a number picked to pass.
            SIDE_BEARING = 4       # px at 300 dpi, about 0.04 em
            for name, cols in (("numerator", above), ("denominator", below)):
                if not cols:
                    failures.append("%s: no %s ink" % (latex, name))
                elif (min(cols) < widest[0] - SIDE_BEARING
                      or max(cols) > widest[1] + SIDE_BEARING):
                    failures.append(
                        "%s: the rule spans %s but the %s spans %d..%d"
                        % (latex, widest, name, min(cols), max(cols)))

        # Nothing may be clipped: the reported size has to contain the ink.
        for latex in (r"\sqrt{x^{2}}", r"\frac{a}{b}", r"\int_{a}^{b}f",
                      r"\sum_{k=1}^{n}x^{k}", r"E = mc^{2}"):
            try:
                r = render(latex)
            except RuntimeError as exc:
                failures.append(str(exc))
                continue
            edges = []
            if any(r.rows[0]):
                edges.append("top")
            if any(r.rows[r.h - 1]):
                edges.append("bottom")
            if any(row[0] for row in r.rows):
                edges.append("left")
            if any(row[r.w - 1] for row in r.rows):
                edges.append("right")
            if edges:
                failures.append("%s: ink touches the %s edge of the image"
                                % (latex, "/".join(edges)))

    if failures:
        print("FAIL  %d" % len(failures))
        for f in failures:
            print("  " + f)
        return 1
    print("ok    rendered images: radical join, fraction rule span, no "
          "clipping")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
