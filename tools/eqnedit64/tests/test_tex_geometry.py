"""The canvas is checked against LaTeX, not against itself.

Every layout constant in the renderer was once picked by eye, and every one
of them was wrong: the fraction rule was 0.39 em wider than TeX's, the
numerator sat 0.24 em too close to it, exponents were set 0.08 em too small,
and a superscript rode 0.04 em high.  The layout tests in place at the time
all passed, because they compared the renderer with its own output.

`tests/tex_reference.json` holds geometry measured from real pdfLaTeX output
by `tools/tex_geometry.py refresh`.  This test compares the renderer against
those numbers, so a constant that drifts away from TeX fails here rather than
being noticed in a screenshot.

Run:  python tests\\test_tex_geometry.py
"""
from __future__ import annotations

import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
sys.path.insert(0, os.path.join(_ROOT, "build"))
sys.path.insert(0, _ROOT)

from tools.tex_geometry import load_reference, measure_eqnedit  # noqa: E402

# How far from TeX a measurement may sit, in em.  0.03 em is under half a
# point at 12 pt -- below what the eye picks up, above the noise in reading
# glyph origins out of a PDF.
TOLERANCE = 0.03

# Measurements known to disagree, with the reason.  An entry here is a debt,
# not a licence: it says the difference is understood and still open.
# Empty: every probed construct is now within tolerance of pdfLaTeX.
KNOWN = {}


def main() -> int:
    reference = load_reference()
    measured = measure_eqnedit()

    failures, allowed = [], []
    for key in sorted(reference):
        if key not in measured:
            failures.append("%s: the renderer produced no measurement" % key)
            continue
        want, got = reference[key], measured[key]
        delta = got - want
        if abs(delta) <= TOLERANCE:
            continue
        if key in KNOWN:
            allowed.append("%s: %+.3f em -- %s" % (key, delta, KNOWN[key]))
        else:
            failures.append("%s: TeX %.3f, canvas %.3f, off by %+.3f em"
                            % (key, want, got, delta))

    for key in KNOWN:
        if key in measured and key in reference:
            if abs(measured[key] - reference[key]) <= TOLERANCE:
                failures.append(
                    "%s now matches TeX -- remove it from KNOWN" % key)

    if failures:
        print("FAIL  %d" % len(failures))
        for f in failures:
            print("  " + f)
        return 1
    print("ok    %d measurements within %.3f em of pdfLaTeX (%d known "
          "differences)" % (len(reference) - len(allowed), TOLERANCE,
                            len(allowed)))
    for a in allowed:
        print("      known: " + a)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
