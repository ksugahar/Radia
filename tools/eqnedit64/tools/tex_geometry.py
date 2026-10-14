"""Measure how a construct is laid out, by pdfLaTeX and by Eqnedit64.

The renderer's layout constants were chosen by eye and every one of them was
wrong: the fraction rule was 0.39 em too wide, the numerator 0.24 em too
close to it, exponents 0.08 em too small.  Nothing caught that, because the
only checks compared the renderer with itself.

This module renders the same equation twice -- once with pdflatex, once with
the editor -- and reports both geometries in em.  `refresh` writes the
pdfLaTeX side to `tests/tex_reference.json` so the test can run on a machine
without a TeX installation; re-run it when a construct is added.

    python tools\\tex_geometry.py refresh      # needs pdflatex + PyMuPDF
    python tools\\tex_geometry.py report       # show both sides now
"""
from __future__ import annotations

import json
import os
import re
import subprocess
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
sys.path.insert(0, os.path.join(_ROOT, "build"))

REFERENCE = os.path.join(_ROOT, "tests", "tex_reference.json")
EM = 12.0

# name -> (latex, [probes]).  A probe is (label, kind, *glyphs).
CASES = [
    ("superscript", r"x^{2}", [
        ("superscript raise", "raise", "x", "2"),
        ("script size ratio", "size", "2"),
    ]),
    ("subscript", r"x_{2}", [
        ("subscript drop", "drop", "x", "2"),
    ]),
    # A superscript starts at the base's width plus the base's italic
    # correction.  Without it the 2 of f^2 sat on the f's hook -- 0.097 em
    # short for f and 0.202 for V, while x and a were fine, so it read as a
    # problem with one letter rather than a missing rule.  Three bases with
    # very different leans, so leaving the rule out cannot pass.
    ("superscript on f", r"f^{2}", [
        ("superscript offset", "dx", "f", "2"),
    ]),
    ("superscript on V", r"V^{2}", [
        ("superscript offset", "dx", "V", "2"),
    ]),
    ("superscript on x", r"x^{2}", [
        ("superscript offset", "dx", "x", "2"),
    ]),
    # A function name is an operator and takes a thin space before its
    # argument.  \sin arrives as a run of Function-styled letters whose first
    # atom is an ordinary "s", so the whole word was classed as ordinary and
    # the space went missing: sin omega t set as sinωt.
    ("function then argument", r"\sin x", [
        ("space after the name", "dx", "n", "x"),
    ]),
    ("fraction", r"\frac{a}{b}", [
        ("numerator above bar", "above_rule", "a"),
        ("denominator below bar", "below_rule", "b"),
        ("rule width", "rule_width"),
        ("rule thickness", "rule_thickness"),
    ]),
    ("radical", r"\sqrt{x}", [
        ("bar above baseline", "above_baseline", "x"),
        ("bar thickness", "rule_thickness"),
    ]),
    # The vinculum is default_rule_thickness whatever it covers.  Measuring
    # only \sqrt{x} missed that the bar used to be scaled by the radical
    # glyph's stretch, so it grew as the content grew.
    ("tall radical", r"\sqrt{\langle a\rangle}", [
        ("bar thickness", "rule_thickness"),
    ]),
    ("integral limits", r"\int_a^b", [
        ("upper right of lower", "dx", "a", "b"),
        ("limit vertical span", "dy", "b", "a"),
        # How far each limit sits from the integral sign itself.  Without
        # these two the limits could drift arbitrarily far right while the
        # checks above still passed, which is exactly what happened.
        ("lower limit from sign", "dx_first", "a"),
        ("upper limit from sign", "dx_first", "b"),
    ]),
    ("sum limits", r"\sum_a^b", [
        ("limit vertical span", "dy", "b", "a"),
    ]),
    # Where a display operator sits relative to the baseline is deliberately
    # NOT probed here.  A big operator's glyph origin is a font's own choice --
    # cmex puts it near the bottom of the sign, Latin Modern at the letter
    # baseline -- so comparing origins compares two different conventions and
    # reports a 1.6 em disagreement where the ink actually agrees.  That is the
    # same apples-to-oranges reading that once made the integral look too big
    # when it was too small.  It is checked as rendered ink instead, in
    # tests/test_render_image.py.
]

_DOC = (r"\documentclass[preview,border=6pt,12pt]{standalone}" "\n"
        r"\usepackage{amsmath}" "\n" r"\begin{document}" "\n"
        r"$\displaystyle %s$" "\n" r"\end{document}" "\n")

_TEXT = re.compile(r'<text x="([-\d.]+)" y="([-\d.]+)"[^>]*'
                   r'font-size="([\d.]+)"[^>]*>([^<]*)</text>')
_RECT = re.compile(r'<rect x="([-\d.]+)" y="([-\d.]+)" '
                   r'width="([\d.]+)" height="([\d.]+)"')


_PATHGLYPH = re.compile(
    r'<path data-char="([^"]*)" data-size="([\d.]+)" '
    r'transform="translate\(([-\d.]+),([-\d.]+)\)')

def _fold(ch):
    """Math-italic letters back to ASCII.

    The canvas draws a variable as its Unicode math-italic code point, the way
    TeX sets it, while pdfLaTeX reports a plain "x".  Comparing the raw
    characters made every probe that names a letter come back empty -- and an
    empty probe is silently skipped, so the comparison table quietly shrank
    from fifteen rows to four instead of failing.
    """
    cp = ord(ch)
    if 0x1D44E <= cp <= 0x1D467:
        return chr(ord("a") + cp - 0x1D44E)
    if 0x1D434 <= cp <= 0x1D44D:
        return chr(ord("A") + cp - 0x1D434)
    if cp == 0x210E:
        return "h"
    if 0x1D6FC <= cp <= 0x1D714:
        return chr(0x03B1 + cp - 0x1D6FC)
    return ch


def _probe(glyphs, rules, probe):
    """glyphs: [(char, x, y, size)], rules: [(x0, y0, x1, y1)]."""
    def find(c):
        for g in glyphs:
            if _fold(g[0]) == c:
                return g
        return None

    kind = probe[1]
    if kind == "size":
        g = find(probe[2])
        return g[3] / EM if g else None
    if kind == "rule_thickness":
        if not rules:
            return None
        r = max(rules, key=lambda r: r[2] - r[0])
        return (r[3] - r[1]) / EM
    if kind in ("dx_first", "dy_first"):
        # Offset from the leftmost glyph -- the big operator -- without
        # naming it, because cmex encodes the integral as "Z" and Cambria
        # as U+222B.  This is the distance the user sees between an
        # integral and its limits.
        g = find(probe[2])
        if not g or not glyphs:
            return None
        first = min(glyphs, key=lambda q: q[1])
        return ((g[1] - first[1]) / EM if kind == "dx_first"
                else (g[2] - first[2]) / EM)
    if kind == "rule_width":
        if not rules:
            return None
        r = max(rules, key=lambda r: r[2] - r[0])
        return (r[2] - r[0]) / EM
    if kind in ("above_rule", "below_rule", "above_baseline"):
        g = find(probe[2])
        if not g or not rules:
            return None
        if kind == "above_baseline":
            r = min(rules, key=lambda r: r[1])
            return (g[2] - r[1]) / EM
        r = max(rules, key=lambda r: r[2] - r[0])
        return ((r[1] - g[2]) / EM if kind == "above_rule"
                else (g[2] - r[3]) / EM)
    a, b = find(probe[2]), find(probe[3])
    if not a or not b:
        return None
    return ((b[1] - a[1]) / EM if kind == "dx" else (b[2] - a[2]) / EM)


def tex_geometry(latex, name, workdir):
    import fitz
    os.makedirs(workdir, exist_ok=True)
    path = os.path.join(workdir, name + ".tex")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(_DOC % latex)
    # pdfTeX refuses to start when the working directory has non-ASCII
    # characters, and this project lives under a Japanese path, so run it
    # from the scratch directory with a bare file name.
    subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error",
                    name + ".tex"], capture_output=True, cwd=workdir)
    pdf = os.path.join(workdir, name + ".pdf")
    if not os.path.exists(pdf):
        raise RuntimeError("pdflatex could not typeset %r" % latex)
    doc = fitz.open(pdf)
    page = doc[0]
    glyphs = []
    for block in page.get_text("rawdict")["blocks"]:
        for line in block.get("lines", []):
            for span in line.get("spans", []):
                for ch in span["chars"]:
                    if ch["c"].strip():
                        glyphs.append((ch["c"], ch["origin"][0],
                                       ch["origin"][1], span["size"]))
    rules = []
    for d in page.get_drawings():
        r = d["rect"]
        # A TeX rule is a stroked line: its rect is degenerate and the
        # thickness lives in the stroke width.  Reading only the rect
        # reported every pdfLaTeX bar as zero thick, which is why nothing
        # here ever noticed that every bar in the canvas was 25% too thick.
        height = r.y1 - r.y0
        stroke = d.get("width") or 0.0
        if height < stroke:
            mid = (r.y0 + r.y1) / 2.0
            rules.append((r.x0, mid - stroke / 2.0, r.x1, mid + stroke / 2.0))
        else:
            rules.append((r.x0, r.y0, r.x1, r.y1))
    doc.close()
    return glyphs, rules


def eqnedit_geometry(latex):
    from eqnedit_core import SvgStyle, tex_to_svg
    svg = tex_to_svg(latex, SvgStyle())        # shipped defaults, untouched
    glyphs = [(m.group(4), float(m.group(1)), float(m.group(2)),
               float(m.group(3))) for m in _TEXT.finditer(svg)]
    # A designed size variant is drawn as an outline, because no viewer can
    # name it.  It is still a glyph and still has to be measured; skipping it
    # made the leftmost glyph the first limit rather than the operator, and
    # the limit offsets came back measured from the wrong thing.
    glyphs += [(m.group(1), float(m.group(3)), float(m.group(4)),
                float(m.group(2))) for m in _PATHGLYPH.finditer(svg)]
    glyphs.sort(key=lambda g: g[1])
    rules = [(float(m.group(1)), float(m.group(2)),
              float(m.group(1)) + float(m.group(3)),
              float(m.group(2)) + float(m.group(4)))
             for m in _RECT.finditer(svg)]
    return glyphs, rules


def measure_eqnedit():
    out = {}
    for name, latex, probes in CASES:
        glyphs, rules = eqnedit_geometry(latex)
        for probe in probes:
            value = _probe(glyphs, rules, probe)
            if value is not None:
                out["%s / %s" % (name, probe[0])] = round(value, 4)
    return out


def measure_tex(workdir):
    out = {}
    for name, latex, probes in CASES:
        glyphs, rules = tex_geometry(latex, re.sub(r"\W", "_", name), workdir)
        for probe in probes:
            value = _probe(glyphs, rules, probe)
            if value is not None:
                out["%s / %s" % (name, probe[0])] = round(value, 4)
    return out


def load_reference():
    with open(REFERENCE, encoding="utf-8") as fh:
        return json.load(fh)


def main(argv):
    action = argv[1] if len(argv) > 1 else "report"
    if action == "refresh":
        work = os.path.join(os.environ.get("TEMP", "."), "eqnedit_texref")
        ref = measure_tex(work)
        with open(REFERENCE, "w", encoding="utf-8") as fh:
            json.dump(ref, fh, indent=2, sort_keys=True)
            fh.write("\n")
        print("wrote %d reference measurements to %s" % (len(ref), REFERENCE))
        return 0
    ref = load_reference()
    mine = measure_eqnedit()
    print("%-46s %9s %10s %8s" % ("measurement", "pdfLaTeX", "Eqnedit64",
                                  "diff"))
    print("-" * 78)
    for key in sorted(set(ref) | set(mine)):
        t, v = ref.get(key), mine.get(key)
        if t is None or v is None:
            print("%-46s %9s %10s" % (key, t if t is not None else "-",
                                      v if v is not None else "-"))
            continue
        print("%-46s %9.3f %10.3f %+8.3f%s"
              % (key, t, v, v - t, "  <-- off" if abs(v - t) > 0.03 else ""))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
