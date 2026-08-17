"""End-to-end: does Word accept our clipboard RTF as a real equation?

The unit tests assert the control words; this asserts the only thing that
finally matters -- put the RTF on the clipboard, paste it into Word, save as
.docx, and read back what Word made of it.  An equation that arrives as text
runs instead of <m:oMath> is a picture of an equation at best, and no amount of
correct-looking markup would reveal that.

Needs Word installed and drives it through COM, so it lives here rather than in
tests/: it is slow, it is a real application, and it cannot run on CI.

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


def _put_rtf(rtf: str) -> None:
    fmt = cb.RegisterClipboardFormat("Rich Text Format")
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
        cb.SetClipboardData(fmt, rtf.encode("latin-1", "replace"))
    finally:
        cb.CloseClipboard()


@pytest.fixture(scope="module")
def word():
    app = com.Dispatch("Word.Application")
    app.Visible = False
    yield app
    app.Quit()


@pytest.mark.validation
@pytest.mark.parametrize("latex", CASES)
def test_word_pastes_it_as_an_equation(word, latex, tmp_path):
    _put_rtf(equation.tex_to_rtf(latex))
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
