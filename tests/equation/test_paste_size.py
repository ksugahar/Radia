"""A pasted equation is 24 pt, because we say so and not the destination box.

Sugahara, 2026-08-20: *powerpointに貼り付けたときは、24ptにしてほしいよ、18ptでは
小さい*.

The clipboard could not say it.  PowerPoint takes the MathML and turns it into
a native equation, but it ignores MathML sizing -- `mathsize` on `<math>` and
on an `<mstyle>` were both pasted into a real slide and both came out at the
placeholder's own size.  So an equation dropped into 18 pt body text was 18 pt,
and there was nothing on the clipboard that disagreed.

PowerPoint's own copy does say it: `Art::GVML ClipFormat`, an OPC package whose
runs carry an explicit `sz`.  `tex_to_gvml` writes that, the same way the RTF
here was written by transcribing what Word puts on the clipboard.

Most of this checks the package without Office, so it runs anywhere.  The last
test pastes into a real slide through PowerPoint's object model and reads the
size back, and skips where PowerPoint is not installed -- the claim is about
what PowerPoint does, so something has to actually ask it.
"""

from __future__ import annotations

import io
import zipfile

import pytest

_eq = pytest.importorskip("radia._equation")

if not hasattr(_eq, "tex_to_gvml"):
    pytest.skip("built before the GVML clipboard format",
                allow_module_level=True)

tex_to_gvml = _eq.tex_to_gvml
PASTE_SIZE_PT = _eq.PASTE_SIZE_PT

LATEX = r"E = mc^{2}"


def parts(latex=LATEX, **kw):
    with zipfile.ZipFile(io.BytesIO(tex_to_gvml(latex, **kw))) as z:
        assert z.testzip() is None, "a bad CRC in the package"
        return {n: z.read(n).decode("utf-8") for n in z.namelist()}


# ---- the package -----------------------------------------------------------

def test_the_default_paste_size_is_24_pt():
    assert PASTE_SIZE_PT == 24.0


def test_it_is_an_opc_package():
    """PowerPoint reads this format as a ZIP; anything else is rejected in
    silence, which looks exactly like the format not being offered."""
    pkg = tex_to_gvml(LATEX)
    assert pkg[:4] == b"PK\x03\x04", pkg[:8]
    assert set(parts()) == {
        "[Content_Types].xml",
        "_rels/.rels",
        "clipboard/drawings/drawing1.xml",
    }


def test_the_size_is_in_the_drawing():
    drawing = parts()["clipboard/drawings/drawing1.xml"]
    assert 'sz="2400"' in drawing, drawing[:400]


def test_the_size_follows_the_argument():
    drawing = parts(size_pt=18.0)["clipboard/drawings/drawing1.xml"]
    assert 'sz="1800"' in drawing
    assert 'sz="2400"' not in drawing


def test_it_carries_an_equation_not_a_picture():
    """<a14:m> is what makes PowerPoint treat the OMML as maths; without it
    the same XML is ignored, and a picture would follow neither the theme nor
    the equation editor."""
    drawing = parts()["clipboard/drawings/drawing1.xml"]
    assert "<a14:m" in drawing
    assert "<m:oMath" in drawing
    assert "xmlns:m=" in drawing, "the a14:m wrapper does not declare m:"


def test_the_content_types_declare_the_drawing():
    ct = parts()["[Content_Types].xml"]
    assert "/clipboard/drawings/drawing1.xml" in ct
    assert "drawing+xml" in ct


def test_the_relationship_points_at_the_drawing():
    rels = parts()["_rels/.rels"]
    assert "clipboard/drawings/drawing1.xml" in rels


@pytest.mark.parametrize("latex", [
    r"\frac{\partial B}{\partial t}",
    r"\sqrt{\alpha^{2}+1}",
    r"\int_{0}^{1} f(x)\,dx",
    r"日本語 x",
])
def test_every_equation_makes_a_readable_package(latex):
    drawing = parts(latex)["clipboard/drawings/drawing1.xml"]
    assert "<m:oMath" in drawing
    assert 'sz="2400"' in drawing


# ---- and what PowerPoint actually does with it ------------------------------

@pytest.mark.slow
def test_powerpoint_pastes_it_at_24_pt():
    """The claim is about PowerPoint, so ask PowerPoint."""
    win32com = pytest.importorskip("win32com.client")
    import ctypes
    import ctypes.wintypes as wt

    try:
        ppt = win32com.Dispatch("PowerPoint.Application")
    except Exception as exc:                       # noqa: BLE001
        pytest.skip(f"PowerPoint is not available: {exc}")

    u32 = ctypes.WinDLL("user32", use_last_error=True)
    k32 = ctypes.WinDLL("kernel32", use_last_error=True)
    u32.OpenClipboard.argtypes = [ctypes.c_void_p]
    u32.RegisterClipboardFormatW.argtypes = [ctypes.c_wchar_p]
    u32.RegisterClipboardFormatW.restype = wt.UINT
    u32.SetClipboardData.argtypes = [wt.UINT, ctypes.c_void_p]
    u32.SetClipboardData.restype = ctypes.c_void_p
    k32.GlobalAlloc.argtypes = [wt.UINT, ctypes.c_size_t]
    k32.GlobalAlloc.restype = ctypes.c_void_p
    k32.GlobalLock.argtypes = [ctypes.c_void_p]
    k32.GlobalLock.restype = ctypes.c_void_p
    k32.GlobalUnlock.argtypes = [ctypes.c_void_p]

    data = tex_to_gvml(LATEX)
    if not u32.OpenClipboard(None):
        pytest.skip("the clipboard is held by another process")
    try:
        u32.EmptyClipboard()
        h = k32.GlobalAlloc(0x0002, len(data))
        ctypes.memmove(k32.GlobalLock(h), data, len(data))
        k32.GlobalUnlock(h)
        u32.SetClipboardData(u32.RegisterClipboardFormatW(
            "Art::GVML ClipFormat"), h)
    finally:
        u32.CloseClipboard()

    pres = ppt.Presentations.Add(WithWindow=True)
    try:
        slide = pres.Slides.Add(1, 12)             # ppLayoutBlank
        slide.Shapes.Paste()
        shape = slide.Shapes(slide.Shapes.Count)
        text = shape.TextFrame.TextRange
        # Runs(i) off the range itself; the collection object it returns with
        # no argument is not indexable through late binding.
        sizes = {round(text.Runs(i).Font.Size, 1)
                 for i in range(1, text.Runs().Count + 1)}
        assert sizes == {24.0}, sizes
    finally:
        pres.Close()
        ppt.Quit()
