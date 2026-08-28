"""Compare many constructs against pdfLaTeX at once, and list what differs.

`tex_geometry.py` measures a handful of named distances very precisely.  This
is the other half: it renders a broad corpus both ways and compares the ink,
so a construct nobody thought to probe still gets looked at.  That gap is not
hypothetical -- the accents floated a full base-height above their letters,
and every geometry probe passed while they did, because none of them measured
an accent.

Two scale-free numbers per construct:

  aspect    the ink box's width over its height.  Catches anything set too
            wide, too tall, too far apart or too close.
  weight    the fraction of that box which is ink.  Catches a glyph drawn at
            the wrong size or a rule at the wrong thickness, which can leave
            the aspect untouched.

Both sides are rasterised and measured identically -- comparing a PDF's own
metrics against rendered ink is the apples-to-oranges reading that once made
the integral look too big when it was too small.

    python tools\\tex_sweep.py            # report, worst first
    python tools\\tex_sweep.py --all      # every construct, not just the bad
"""
from __future__ import annotations

import os
import subprocess
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)

try:
    import fitz
    from PIL import Image
except ImportError:  # pragma: no cover
    print("skip  PyMuPDF and Pillow are needed for the sweep")
    raise SystemExit(0)

DPI = 900
DOC = ("\\documentclass[preview,border=2pt,12pt]{standalone}\n"
       "\\usepackage{amsmath}\n"
       "\\begin{document}$\\displaystyle %s$\\end{document}\n")

CORPUS = [
    # accents -- the family that was reported
    (r"\tilde{a}", "tilde"), (r"\hat{a}", "hat"), (r"\bar{a}", "bar"),
    (r"\dot{a}", "dot"), (r"\ddot{a}", "ddot"), (r"\vec{a}", "vec"),
    (r"\tilde{M}", "tilde on a capital"),
    # scripts
    (r"x^{2}", "superscript"), (r"x_{i}", "subscript"),
    (r"x_{i}^{2}", "both scripts"), (r"x^{2^{2}}", "nested superscript"),
    (r"\alpha^{\beta}", "greek script"),
    # fractions and radicals
    (r"\frac{a}{b}", "fraction"), (r"\frac{a+b}{c+d}", "wider fraction"),
    (r"\frac{\frac{a}{b}}{c}", "nested fraction"),
    (r"\sqrt{x}", "radical"), (r"\sqrt{x^{2}+y^{2}}", "wide radical"),
    (r"\sqrt{\frac{a}{b}}", "radical over a fraction"),
    # operators
    (r"\int_a^b", "integral"), (r"\sum_{i=1}^{n}", "sum"),
    (r"\prod_{i=1}^{n}", "product"), (r"\oint_C", "contour integral"),
    (r"\lim_{x \to 0}", "limit"),
    (r"\int_a^b f(x)dx", "integral with a body"),
    # fences
    (r"\left( \frac{a}{b} \right)", "grown parentheses"),
    (r"\left[ x \right]", "brackets"),
    (r"\left\{ x \right\}", "braces"),
    (r"\left| x \right|", "bars"),
    (r"\left\langle a \right\rangle", "angle brackets"),
    (r"\left( \frac{\frac{a}{b}}{c} \right)", "tall parentheses"),
    # decorations
    (r"\overline{abc}", "overline"), (r"\underline{abc}", "underline"),
    (r"\overrightarrow{AB}", "over arrow"),
    (r"\overbrace{a+b}", "over brace"),
    # spacing and atoms
    (r"a+b", "binary operator"), (r"a=b", "relation"),
    (r"-x", "unary minus"), (r"f(x)", "function call"),
    (r"\sin x", "named function"), (r"a \quad b", "quad space"),
    (r"x \in A", "set relation"), (r"a \times b", "times"),
    # matrices
    (r"\begin{matrix} a & b \\ c & d \end{matrix}", "matrix"),
    (r"\begin{cases} x & x>0 \\ -x & x<0 \end{cases}", "cases"),
]


NORMAL_HEIGHT = 400   # both sides are scaled to this before ink is counted


def ink_box(path):
    """Aspect of the ink box, and how much of it is ink.

    Both sides are scaled to one height first.  Without that the ink fraction
    is really a measure of rasterisation resolution -- pdfLaTeX rendered at
    900 dpi against a canvas PNG a third that size showed every glyph as
    "thinner", which is a property of the comparison and not of the layout.
    """
    image = Image.open(path).convert("L")
    box = image.point(lambda v: 255 if v < 200 else 0).getbbox()
    if not box:
        return None
    crop = image.crop(box)
    width = max(1, int(round(crop.width * NORMAL_HEIGHT / crop.height)))
    crop = crop.resize((width, NORMAL_HEIGHT), Image.LANCZOS)
    mask = crop.point(lambda v: 255 if v < 200 else 0)
    ink = sum(1 for v in mask.get_flattened_data() if v)
    return width / float(NORMAL_HEIGHT), ink / float(width * NORMAL_HEIGHT)


def tex_side(latex, name, workdir):
    stem = "sw%d" % (abs(hash(name)) % 100000)
    with open(os.path.join(workdir, stem + ".tex"), "w", encoding="utf-8") as fh:
        fh.write(DOC % latex)
    subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error",
                    stem + ".tex"], capture_output=True, cwd=workdir)
    pdf = os.path.join(workdir, stem + ".pdf")
    if not os.path.exists(pdf):
        return None
    doc = fitz.open(pdf)
    png = os.path.join(workdir, stem + ".png")
    doc[0].get_pixmap(dpi=DPI).save(png)
    doc.close()
    return ink_box(png)


def canvas_side(latex, name, workdir):
    exe = os.path.join(_ROOT, "build", "Eqnedit64.exe")
    if not os.path.exists(exe):
        exe = os.path.join(_ROOT, "dist", "Eqnedit64.exe")
    png = os.path.join(workdir, "canvas%d.png" % (abs(hash(name)) % 100000))
    subprocess.run([exe, "--render-png", latex, png], capture_output=True)
    if not os.path.exists(png):
        return None
    return ink_box(png)


def main(show_all):
    work = os.path.join(os.environ.get("TEMP", "."), "eqnsweep")
    os.makedirs(work, exist_ok=True)
    rows = []
    for latex, name in CORPUS:
        tex = tex_side(latex, name, work)
        ours = canvas_side(latex, name, work)
        if not tex or not ours:
            rows.append((99.0, name, latex, None, None))
            continue
        aspect = abs(ours[0] / tex[0] - 1.0)
        weight = abs(ours[1] / tex[1] - 1.0)
        rows.append((max(aspect, weight), name, latex, tex, ours))
    rows.sort(reverse=True)

    print("%-26s %9s %9s   %9s %9s" %
          ("construct", "aspect", "(TeX)", "weight", "(TeX)"))
    print("-" * 72)
    shown = 0
    for worst, name, latex, tex, ours in rows:
        if not show_all and worst < 0.12:
            continue
        shown += 1
        if tex is None:
            print("%-26s  could not be rendered on one side" % name)
            continue
        print("%-26s %9.3f %9.3f   %9.3f %9.3f   %s"
              % (name, ours[0], tex[0], ours[1], tex[1],
                 "<-- %.0f%%" % (100 * worst)))
    if not shown:
        print("nothing over 12%% across %d constructs" % len(rows))
    else:
        print("\n%d of %d constructs differ by more than 12%%"
              % (shown, len(rows)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main("--all" in sys.argv))
