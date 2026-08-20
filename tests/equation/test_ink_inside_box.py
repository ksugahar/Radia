r"""What is drawn has to be inside the box the layout reports.

A layout returns a width, an ascent and a descent, and every consumer trusts
them: the picture path sizes its bitmap from them, the window centres the
equation with them, and a palette button shrinks to fit them.  A construct
placed OUTSIDE that box therefore does not merely look wrong -- it is cropped
away by the picture, and it paints over whatever is next to it on screen.

That is not hypothetical.  `\overbrace` shifted its brace by the brace's
ASCENT where it should have used its descent, putting it a whole brace-height
too high: the pasted picture lost the brace completely (a blank strip of paper
where the annotation should be), and the editor's palette drew it into the row
of buttons above.  Nothing failed; the drawing was simply somewhere else.

Rendering and counting ink is the only check that sees this, because both the
geometry and the drawing agree with each other -- they are just both outside
the frame.
"""

from __future__ import annotations

import io
import struct
import zlib

import pytest

equation = pytest.importorskip("radia.equation")


def ink_rows(png: bytes):
    """Rows of the PNG that contain a non-white pixel, and the height."""
    # Minimal PNG reader: enough for the greyscale/RGB(A) images tex_to_png
    # writes, and it keeps this test free of an image dependency.
    assert png[:8] == b"\x89PNG\r\n\x1a\n", "not a PNG"
    pos, idat, w, h, depth, colour = 8, b"", 0, 0, 0, 0
    while pos < len(png):
        (length,) = struct.unpack(">I", png[pos:pos + 4])
        kind = png[pos + 4:pos + 8]
        body = png[pos + 8:pos + 8 + length]
        if kind == b"IHDR":
            w, h, depth, colour = struct.unpack(">IIBB", body[:10])
        elif kind == b"IDAT":
            idat += body
        pos += 12 + length
    assert depth == 8, f"unexpected bit depth {depth}"
    channels = {0: 1, 2: 3, 4: 2, 6: 4}[colour]
    raw = zlib.decompress(idat)
    stride = w * channels
    prev = bytearray(stride)
    rows, at = [], 0
    for y in range(h):
        f = raw[at]
        line = bytearray(raw[at + 1:at + 1 + stride])
        at += 1 + stride
        # undo the PNG row filter
        for i in range(stride):
            a = line[i - channels] if i >= channels else 0
            b = prev[i]
            c = prev[i - channels] if i >= channels else 0
            x = line[i]
            if f == 1:
                x += a
            elif f == 2:
                x += b
            elif f == 3:
                x += (a + b) // 2
            elif f == 4:
                p = a + b - c
                pa, pb, pc = abs(p - a), abs(p - b), abs(p - c)
                x += a if (pa <= pb and pa <= pc) else (b if pb <= pc else c)
            line[i] = x & 0xFF
        prev = line
        if any(line[i] < 200 for i in range(0, stride, channels)):
            rows.append(y)
    return rows, h


def png_of(tex, scale=6.0):
    return equation.tex_to_png(tex, equation.SvgStyle(), scale)


@pytest.mark.parametrize("tex", [
    r"\overbrace{ab}",
    r"\overbrace{ab}^{c}",
    r"\underbrace{ab}",
    r"\underbrace{ab}_{c}",
    r"\overline{ab}",
    r"\vec{a}",
    r"\hat{a}",
    r"\sqrt{a}",
    r"\sum_{i}^{n}",
    r"\frac{a}{b}",
])
def test_the_picture_is_not_blank(tex):
    """A construct that draws nothing at all has left its own box."""
    rows, h = ink_rows(png_of(tex))
    assert rows, f"{tex} rendered a blank picture ({h} rows of nothing)"


@pytest.mark.parametrize("tex,where", [
    (r"\overbrace{ab}", "top"),
    (r"\overbrace{ab}^{c}", "top"),
    (r"\overline{ab}", "top"),
    (r"\vec{a}", "top"),
    (r"\underbrace{ab}", "bottom"),
    (r"\underbrace{ab}_{c}", "bottom"),
    (r"\underline{ab}", "bottom"),
])
def test_the_decoration_is_where_it_belongs(tex, where):
    """An overbrace draws in the upper part of its own picture, and an
    underbrace in the lower one.  Both were cropped to nothing when the brace
    was placed by the wrong edge."""
    rows, h = ink_rows(png_of(tex))
    assert rows, f"{tex} rendered blank"
    if where == "top":
        assert min(rows) < h / 3, (
            f"{tex}: nothing in the top third of {h} rows (ink from "
            f"{min(rows)} to {max(rows)})")
    else:
        assert max(rows) > 2 * h / 3, (
            f"{tex}: nothing in the bottom third of {h} rows (ink from "
            f"{min(rows)} to {max(rows)})")
