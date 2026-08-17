"""End-to-end: does each Office application take our clipboard as maths?

The unit tests assert the markup; this asserts the only thing that finally
matters -- put it on the clipboard, paste, save, and read back what the
application made of it.  An equation that arrives as text runs instead of
<m:oMath> is a picture of an equation at best, and no amount of correct-looking
markup would reveal that.

One Copy carries both formats and each application takes what it understands:
Word reads the RTF as maths and ignores the MathML, PowerPoint does the
opposite.  Measured by offering them one at a time.

Excel is not covered.  Its equations live in shapes, and pasting into a shape's
text is a UI operation its object model does not expose -- all three COM routes
(Shape.Select + Paste, TextRange.Select + Paste, PasteSpecial) are refused.  The
clipboard payload is the same one PowerPoint accepts, so it is likely to work,
but that is not measured and is not claimed here.

Needs Office installed and drives it through COM, so it lives here rather than
in tests/: it is slow, it is a real application, and it cannot run on CI.

Run:  python -m pytest validation_test/equation -q
"""

from __future__ import annotations

import re
import time
import zipfile

import pytest

equation = pytest.importorskip("radia.equation")
cb = pytest.importorskip("win32clipboard")
com = pytest.importorskip("win32com.client")

CASES = [
    r"\frac{a}{b}",
    r"a^{2}",
    r"a_{i}^{2}",
    r"\sqrt{x}",
    r"\sqrt[3]{x}",
    r"\left[\frac{a}{b}\right]",
    r"\sum_{i}^{n} a",
    r"\oint_{C} F \cdot dr",
    r"\begin{pmatrix} a & b \\ c & d \end{pmatrix}",
    r"\hat{x}",
    r"\overline{x}",
    r"\sigma_{f}(x) = \sum_{j} B_{fj}(x) a_{j}",
    r"\frac{\mu_{0}}{4\pi} \oint_{C} \frac{I dl \times r}{|r|^{3}}",
]


def _put(latex: str) -> None:
    """Everything an Office application might want, in one go."""
    rtf_fmt = cb.RegisterClipboardFormat("Rich Text Format")
    mml_fmt = cb.RegisterClipboardFormat("MathML")
    for _ in range(12):
        try:
            cb.OpenClipboard()
            break
        except Exception:                   # another process holds it briefly
            time.sleep(0.25)
    else:
        pytest.skip("the clipboard stayed busy")
    try:
        cb.EmptyClipboard()
        cb.SetClipboardData(rtf_fmt,
                            equation.tex_to_rtf(latex).encode("latin-1", "replace"))
        cb.SetClipboardData(mml_fmt, equation.tex_to_mathml(latex).encode("utf-8"))
    finally:
        cb.CloseClipboard()


@pytest.fixture(scope="module")
def word():
    app = com.Dispatch("Word.Application")
    app.Visible = False
    yield app
    app.Quit()


@pytest.fixture(scope="module")
def powerpoint():
    app = com.Dispatch("PowerPoint.Application")
    app.Visible = True          # PowerPoint refuses its paste path when hidden
    yield app
    app.Quit()


@pytest.mark.validation
@pytest.mark.parametrize("latex", CASES)
def test_word_pastes_it_as_an_equation(word, latex, tmp_path):
    _put(latex)
    doc = word.Documents.Add()
    try:
        doc.Content.Paste()
        out = str(tmp_path / "pasted.docx")
        doc.SaveAs2(out, 16)                # wdFormatDocumentDefault
    finally:
        doc.Close(0)

    with zipfile.ZipFile(out) as z:
        xml = z.read("word/document.xml").decode("utf-8")

    n = len(re.findall(r"<m:oMath[ >]", xml))
    if not n:
        text = "".join(re.findall(r"<w:t[^>]*>([^<]*)</w:t>", xml))
        pytest.fail(f"{latex} arrived as text, not maths: {text[:80]!r}")
    assert n == 1


@pytest.mark.validation
@pytest.mark.parametrize("latex", CASES)
def test_powerpoint_pastes_it_as_an_equation(powerpoint, latex, tmp_path):
    _put(latex)
    pres = powerpoint.Presentations.Add(True)
    slide = pres.Slides.Add(1, 12)              # ppLayoutBlank
    time.sleep(0.35)
    try:
        slide.Shapes.Paste()
        out = str(tmp_path / "pasted.pptx")
        pres.SaveAs(out)
    finally:
        pres.Close()

    with zipfile.ZipFile(out) as z:
        xml = z.read("ppt/slides/slide1.xml").decode("utf-8")

    n = len(re.findall(r"<m:oMath[ >]", xml))
    if not n:
        pytest.fail(f"{latex} arrived in PowerPoint as something other than maths")
    assert n == 1
