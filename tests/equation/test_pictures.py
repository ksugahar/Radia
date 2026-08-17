"""EMF and PNG: the equation as a picture, for targets with no equation object.

Which targets those are was measured.  Word and PowerPoint take a native
equation and need no picture at all.  Google Slides has no equation object,
rejects SVG on upload, and accepts raster -- so PNG is the one that always
works, and EMF is how a vector gets there through Google Drawings.  Excel keeps
the metafile beside the bitmap when it takes a picture, so EMF is not lossy
there either.

Both come from the same GDI drawing routine over the same layout, which is also
what an editor draws on screen with; that is why there is no separate encoder
here and no second set of metrics to drift.
"""

from __future__ import annotations

import struct

import pytest

equation = pytest.importorskip("radia.equation")

CASES = [
    r"\frac{a}{b}",
    r"a_{i}^{2}",
    r"\sqrt[3]{x}",
    r"\left[\frac{a}{b}\right]",
    r"\sum_{i}^{n} a",
    r"\sigma_{f}(x) = \sum_{j} B_{fj}(x) a_{j}",
]


def _png_size(data: bytes) -> tuple[int, int]:
    """Width and height from the IHDR, which is always the first chunk."""
    assert data[:8] == b"\x89PNG\r\n\x1a\n", "not a PNG"
    assert data[12:16] == b"IHDR"
    return struct.unpack(">II", data[16:24])


@pytest.mark.parametrize("latex", CASES)
def test_png_is_a_png(latex):
    data = equation.tex_to_png(latex)
    w, h = _png_size(data)
    assert w > 0 and h > 0


@pytest.mark.parametrize("latex", CASES)
def test_emf_is_an_enhanced_metafile(latex):
    data = equation.tex_to_emf(latex)
    # EMR_HEADER: record type 1, then its size, then the bounds.
    rec_type, rec_size = struct.unpack("<II", data[:8])
    assert rec_type == 1, "first record is not EMR_HEADER"
    assert rec_size >= 88
    assert data[40:44] == b" EMF", "missing the EMF signature"


@pytest.mark.parametrize("latex", CASES)
def test_a_picture_carries_ink(latex):
    """A metafile with no drawing records is a blank rectangle -- the failure
    mode that looks fine until it is pasted."""
    emf = equation.tex_to_emf(latex)
    assert b"\x00" in emf and len(emf) > 200
    png = equation.tex_to_png(latex)
    assert len(png) > 200


def test_the_picture_matches_the_layout():
    """Both the PNG and the SVG come from one layout, so their aspect ratios
    have to agree; if they drift, two sets of metrics have appeared."""
    import re

    latex = r"\sigma_{f}(x) = \sum_{j} B_{fj}(x) a_{j}"
    svg = equation.tex_to_svg(latex)
    m = re.search(r'viewBox="0 0 ([\d.]+) ([\d.]+)"', svg)
    svg_ratio = float(m.group(1)) / float(m.group(2))

    w, h = _png_size(equation.tex_to_png(latex))
    assert abs(w / h - svg_ratio) / svg_ratio < 0.05


def test_scale_makes_a_bigger_bitmap_not_a_different_shape():
    latex = r"\frac{a}{b}"
    w1, h1 = _png_size(equation.tex_to_png(latex, equation.SvgStyle(), 1.0))
    w4, h4 = _png_size(equation.tex_to_png(latex, equation.SvgStyle(), 4.0))
    assert w4 > w1 and h4 > h1
    assert abs((w4 / h4) - (w1 / h1)) / (w1 / h1) < 0.1


def test_type_sizes_reach_the_picture():
    small, big = equation.SvgStyle(), equation.SvgStyle()
    big.full, big.sub, big.sym = 24.0, 14.0, 36.0
    w_small, _ = _png_size(equation.tex_to_png("x", small, 1.0))
    w_big, _ = _png_size(equation.tex_to_png("x", big, 1.0))
    assert w_big > w_small


# ---- how big the picture says it is ----------------------------------------
#
# Sampling finely and being enormous are different things, and only the second
# is a bug.  Without a declared resolution an application has nothing but pixels
# to go on and assumes screen dots, so a 600 dpi equation pastes at 600/96 =
# 6.25 times its proper size.  The metafile has always carried its physical size
# -- that is what a vector format is -- so only the raster needed telling.


def _png_phys(data):
    """(x, y) pixels per metre from the pHYs chunk, or None."""
    i = 8
    while i + 8 <= len(data):
        length = struct.unpack(">I", data[i:i + 4])[0]
        if data[i + 4:i + 8] == b"pHYs":
            x, y, _unit = struct.unpack(">IIB", data[i + 8:i + 17])
            return x, y
        i += 12 + length
    return None


def _dib_header(data):
    keys = ("size", "w", "h", "planes", "bits", "compression",
            "image_size", "xppm", "yppm", "used", "important")
    return dict(zip(keys, struct.unpack("<IiiHHIIiiII", data[:40])))


@pytest.mark.parametrize("dpi", [96.0, 288.0, 600.0])
def test_the_png_declares_the_resolution_it_was_rendered_at(dpi):
    png = equation.tex_to_png("x", equation.SvgStyle(), dpi / 72.0)
    phys = _png_phys(png)
    assert phys is not None, "no pHYs chunk: the picture pastes oversized"
    assert abs(phys[0] * 0.0254 - dpi) < 1.0
    assert phys[0] == phys[1]


@pytest.mark.parametrize("dpi", [96.0, 288.0, 600.0])
def test_the_dib_declares_the_resolution_it_was_rendered_at(dpi):
    h = _dib_header(equation.tex_to_dib("x", equation.SvgStyle(), dpi / 72.0))
    assert abs(h["xppm"] * 0.0254 - dpi) < 1.0
    assert h["xppm"] == h["yppm"]


def test_raising_the_resolution_does_not_change_the_declared_size():
    """The whole point: more pixels, same equation."""
    st = equation.SvgStyle()
    sizes = []
    for dpi in (96.0, 288.0, 600.0):
        h = _dib_header(equation.tex_to_dib(r"E = mc^{2}", st, dpi / 72.0))
        sizes.append(h["w"] / h["xppm"])          # metres across
    for s in sizes[1:]:
        assert abs(s - sizes[0]) / sizes[0] < 0.02


def test_the_dib_is_bottom_up_as_a_pasted_bitmap_must_be():
    h = _dib_header(equation.tex_to_dib("x"))
    assert h["h"] > 0
    assert h["bits"] == 32
    assert h["compression"] == 0                  # BI_RGB


def test_the_dib_carries_all_of_its_pixels():
    data = equation.tex_to_dib("x")
    h = _dib_header(data)
    assert len(data) == 40 + h["w"] * h["h"] * 4


def test_the_two_rasters_are_pictures_of_the_same_equation():
    """PNG and DIB come from one rasterisation, so they cannot disagree."""
    st = equation.SvgStyle()
    png_w, png_h = _png_size(equation.tex_to_png(r"\frac{a}{b}", st, 4.0))
    h = _dib_header(equation.tex_to_dib(r"\frac{a}{b}", st, 4.0))
    assert (h["w"], h["h"]) == (png_w, png_h)
