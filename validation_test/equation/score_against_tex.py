"""How close is this editor's setting of an equation to TeX's?

Appearance follows TeX here, so "how close to TeX" is the quality measure, and
it can be computed rather than judged: set the equation twice, once with
pdflatex and once with this editor, and compare the two pictures.

Dimensions alone miss things.  A bar that stops short, a limit lying across an
integral, a script at the wrong height -- all can leave the overall box the
right size.  Ink does not.

Run it by hand after touching the layout:

    python validation_test/equation/score_against_tex.py [pairs-dir]

Give it a directory and it writes each pair as one picture, TeX above and this
below, which is what actually tells you WHAT is wrong.

READ THE SCORE RELATIVELY.  A single letter -- the same letter, from the same
font, laid out by nobody -- scores about 73, because two thin strokes a pixel
apart lose points that have nothing to do with layout.  73 is therefore the
ceiling, not 100, and what matters is which equations fall far below it.

NOTHING IS RESIZED.  Both are set at 12 point and rendered at the same
resolution, so they are already the same physical size, and scaling either to
match the other would hide exactly what this is looking for.  A first version
scaled both to a common WIDTH, which stretched a tall narrow fraction to six
thousand pixels and reported 27 out of 100 for two pictures that agree.  They
are aligned on their ink instead, and the score is read as an indicator, not a
verdict: thin strokes a pixel apart cost real points even when nothing is
wrong.
"""
import os
import subprocess
import sys
import tempfile

sys.path.insert(0, r"S:\Radia\01_GitHub\.claude\worktrees\radia-equation\src")
import radia.equation as eq          # noqa: E402
from PIL import Image, ImageChops    # noqa: E402
import fitz                          # noqa: E402  (PyMuPDF)

PDFLATEX = r"C:\texlive\2026\bin\windows\pdflatex.exe"
WORK = "C:/temp/ee3/_score"          # NOT the system temp: its short form
DPI = 300                            # holds a "~", which TeX reads as active


def tex_png(latex, out):
    doc = (r"\documentclass[12pt]{article}"
           r"\usepackage[paperwidth=40cm,paperheight=14cm,margin=1cm]{geometry}"
           r"\usepackage{amsmath}\pagestyle{empty}"
           r"\begin{document}\noindent$\displaystyle " + latex +
           r"$\end{document}")
    os.makedirs(WORK, exist_ok=True)
    src = os.path.join(WORK, "e.tex")
    with open(src, "w", encoding="utf-8") as f:
        f.write(doc)
    r = subprocess.run([PDFLATEX, "-interaction=nonstopmode", "-halt-on-error",
                        "e.tex"], cwd=WORK, capture_output=True, text=True)
    pdf = os.path.join(WORK, "e.pdf")
    if not os.path.exists(pdf):
        raise RuntimeError("pdflatex failed for %r\n%s" % (latex, r.stdout[-600:]))
    fitz.open(pdf).load_page(0).get_pixmap(dpi=DPI).save(out)
    return out


def our_png(latex, out):
    st = eq.SvgStyle()
    st.padding = 4.0
    with open(out, "wb") as f:
        f.write(eq.tex_to_png(latex, st, DPI / 72.0))
    return out


def ink(path):
    im = Image.open(path).convert("L")
    bw = im.point(lambda v: 255 if v < 160 else 0, mode="L")
    box = bw.getbbox()
    if not box:
        raise RuntimeError("nothing drawn in " + path)
    return bw.crop(box)


def best_overlap(a, b):
    """Slide one over the other by a few pixels and keep the best fit.

    A one-pixel offset between two thin strokes costs more than most real
    layout errors, and it says nothing about the layout -- so let each picture
    find its place before the difference is counted."""
    W = max(a.width, b.width) + 8
    H = max(a.height, b.height) + 8
    pa = Image.new("L", (W, H), 0)
    pa.paste(a, (4, 4))
    A = pa.point(lambda v: 255 if v > 96 else 0)
    na = sum(1 for p in A.getdata() if p)

    best = None
    for dy in range(-4, 5):
        for dx in range(-4, 5):
            pb = Image.new("L", (W, H), 0)
            pb.paste(b, (4 + dx, 4 + dy))
            B = pb.point(lambda v: 255 if v > 96 else 0)
            inter = sum(1 for p in ImageChops.darker(A, B).getdata() if p)
            nb = sum(1 for p in B.getdata() if p)
            union = na + nb - inter
            iou = inter / union if union else 0.0
            if best is None or iou > best[0]:
                best = (iou, inter / na if na else 0.0,
                        (nb - inter) / nb if nb else 0.0, A, B)
    return best


CASES = [
    r"x", r"a+b", r"a=b", r"x^{2}", r"x_{i}", r"x_{i}^{2}",
    r"\frac{a}{b}", r"\frac{abc}{d}", r"\frac{\frac{p}{q}}{c}",
    r"\sqrt{2}", r"\sqrt{\frac{a}{b}}", r"\sqrt[3]{x}",
    r"\int_{0}^{T}f", r"\sum_{n=1}^{N}a_{n}", r"\oint_{C}g",
    r"\left(\frac{1}{2}\right)", r"\alpha+\beta=\gamma",
    r"\nabla\times H=J",
]


def main():
    keepdir = sys.argv[1] if len(sys.argv) > 1 else None
    if keepdir:
        os.makedirs(keepdir, exist_ok=True)
    print("%-28s %6s %8s %9s %s"
          % ("", "score", "covered", "spurious", "size vs TeX"))
    total, n = 0.0, 0
    with tempfile.TemporaryDirectory() as d:
        for i, latex in enumerate(CASES):
            try:
                a = ink(tex_png(latex, os.path.join(d, "t.png")))
                b = ink(our_png(latex, os.path.join(d, "o.png")))
                iou, cov, spur, A, B = best_overlap(a, b)
            except Exception as ex:
                print("%-28s  %s" % (latex, str(ex)[:70]))
                continue
            wr = b.width / a.width
            hr = b.height / a.height
            print("%-28s %6.1f %7.1f%% %8.1f%%   %+5.1f%% wide %+5.1f%% tall"
                  % (latex, iou * 100, cov * 100, spur * 100,
                     (wr - 1) * 100, (hr - 1) * 100))
            if keepdir:
                W, H = A.size
                side = Image.new("L", (W, H * 2 + 6), 0)
                side.paste(A, (0, 0))
                side.paste(B, (0, H + 6))
                side.point(lambda v: 255 - v).save(
                    os.path.join(keepdir, "%02d.png" % i))
            total += iou
            n += 1
    if n:
        print("\nmean %.1f over %d equations "
              "(an indicator, not a verdict -- see the header)" % (total / n * 100, n))


if __name__ == "__main__":
    main()
