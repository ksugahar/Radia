"""SVG-backed pictures in a PPTX -- the half python-pptx does not model.

PowerPoint stores a vector picture as an ``<asvg:svgBlip>`` extension on the
blip, with the raster fallback in the blip itself.  ``python-pptx`` only knows
the raster half and raises ``ValueError: no embedded image`` when a picture has
just the SVG -- which is exactly what PowerPoint's COM ``AddPicture`` leaves
behind.

The paste-scale audit used to swallow that exception, so a fully vectorised
deck reported "0 pictures, 0 flagged": eight figures, none of them measured,
and nothing saying so.  A checker that goes quiet when its input changes shape
is worse than one that fails, because the silence reads as a pass.

Measuring a vector picture is in fact easier than measuring a raster one: an
SVG states its own size, so the authored width needs no DPI metadata, and the
label sizes are written in the file rather than assumed from the lab's 24 pt
authoring convention.
"""
from __future__ import annotations

SVG_NS = "http://schemas.microsoft.com/office/drawing/2016/SVG/main"

# CSS units in points.  A unitless SVG length is CSS px (1/96 in).
_UNIT_PT = {"pt": 1.0, "px": 0.75, "in": 72.0, "cm": 72.0 / 2.54,
            "mm": 7.2 / 2.54, "pc": 12.0, "em": 12.0}

_NUMBER_CHARS = "0123456789."
_LEAD_CHARS = ": ='" + '"'

# No label on a slide is set at this size; a reading above it means the file's
# font-size values are not in the document's own units (see svg_picture_row).
_IMPLAUSIBLE_FONT_PT = 96.0
# A figure authored AT its paste width lands a hair under the floor through
# rounding alone (20.81 cm is 589.98 pt against an authored 590.0), and a check
# that fires on the third decimal teaches people to ignore it.
_FONT_PT_TOL = 0.25


def svg_blob(shape):
    """The SVG bytes behind a picture, or None if it is a plain raster."""
    from pptx.oxml.ns import qn
    try:
        blip = shape._element.blipFill.blip
    except Exception:
        return None
    ext_lst = blip.find(qn("a:extLst"))
    if ext_lst is None:
        return None
    for ext in ext_lst:
        for child in ext:
            if child.tag == "{%s}svgBlip" % SVG_NS:
                rid = child.get(qn("r:embed"))
                if rid:
                    try:
                        return shape.part.related_part(rid).blob
                    except Exception:
                        return None
    return None


def length_pt(text):
    """An SVG length in points, or None if it is not a length."""
    if text is None:
        return None
    s = str(text).strip()
    for unit, factor in _UNIT_PT.items():
        if s.endswith(unit):
            try:
                return float(s[:-len(unit)]) * factor
            except ValueError:
                return None
    try:
        return float(s) * _UNIT_PT["px"]
    except ValueError:
        return None


def _attr(tag_text, name):
    for quote in ('"', "'"):
        key = name + "=" + quote
        at = tag_text.find(key)
        if at < 0:
            continue
        at += len(key)
        end = tag_text.find(quote, at)
        if end > at:
            return tag_text[at:end]
    return None


def svg_size_pt(text):
    """(width_pt, height_pt) of an SVG document, from width/height or viewBox."""
    start = text.find("<svg")
    if start < 0:
        return None, None
    end = text.find(">", start)
    tag = text[start:end if end > start else start + 400]
    w = length_pt(_attr(tag, "width"))
    h = length_pt(_attr(tag, "height"))
    if w and h:
        return w, h
    box = _attr(tag, "viewBox")
    if box:
        parts = box.replace(",", " ").split()
        if len(parts) == 4:
            try:
                vw, vh = float(parts[2]), float(parts[3])
            except ValueError:
                return w, h
            # viewBox units are user units, i.e. CSS px
            return w or vw * _UNIT_PT["px"], h or vh * _UNIT_PT["px"]
    return w, h


def svg_viewbox(text):
    """The viewBox width/height in user units, or (None, None)."""
    start = text.find("<svg")
    if start < 0:
        return None, None
    end = text.find(">", start)
    tag = text[start:end if end > start else start + 400]
    box = _attr(tag, "viewBox")
    if not box:
        return None, None
    parts = box.replace(",", " ").split()
    if len(parts) != 4:
        return None, None
    try:
        return float(parts[2]), float(parts[3])
    except ValueError:
        return None, None


def _transform_scale(tag_text):
    """The uniform scale of a ``transform="matrix(...)"``, or 1.0.

    PyMuPDF puts the text's whole placement matrix on the ``<text>`` element,
    and for a figure that was itself placed at a reduction the matrix carries
    that reduction: ``font-size="173.9"`` inside ``matrix(.0733 ...)`` draws at
    12.7 units, not 173.9.  Reading the attribute alone reported a slide label
    at 237 pt.
    """
    at = tag_text.find("matrix(")
    if at < 0:
        for name in ("scale(",):
            s = tag_text.find(name)
            if s >= 0:
                try:
                    return abs(float(tag_text[s + len(name):
                                              tag_text.find(")", s)].split(",")[0]))
                except ValueError:
                    return 1.0
        return 1.0
    end = tag_text.find(")", at)
    parts = tag_text[at + 7:end].replace(",", " ").split()
    if len(parts) < 2:
        return 1.0
    try:
        a, b = float(parts[0]), float(parts[1])
    except ValueError:
        return 1.0
    scale = (a * a + b * b) ** 0.5
    return scale if scale > 0 else 1.0


def svg_font_sizes_units(text):
    """Every font-size in the document, in USER UNITS, smallest first.

    Deliberately not converted to points here.  Whether a document's user unit
    is a point or a CSS pixel is genuinely ambiguous -- matplotlib writes
    ``width="590.4pt" viewBox="0 0 590.4 295.2"`` (a point), PyMuPDF writes
    ``width="453.48"`` for a 453.48 pt PDF page (also a point, but unitless, so
    the CSS rule would call it a pixel and shrink every label by a quarter).
    Guessing that unit is how a 25 pt label came to be reported as 15 pt.

    The ambiguity disappears one step later: what a label measures ON THE SLIDE
    is ``units * pasted_pt / viewBox_width``, and the unit cancels.  So the
    sizes come out of the file in its own units and are converted against the
    pasted width, never against an assumed physical size.

    Only ``<text>`` elements that actually contain characters are measured, and
    each one's size is multiplied by the scale of its own transform.  Both
    matter: a ``font-size`` sitting on an empty element says nothing about what
    the audience sees, and one sitting inside a reduction matrix is not the
    size it is drawn at.

    Written as a scan rather than a regular expression on purpose: the value
    appears both as an attribute (``font-size="24px"``) and inside a style
    string (``style="font-size: 24px"``), and one scan handles both.
    """
    out = []
    pos = 0
    key = "font-size"
    while True:
        open_at = text.find("<text", pos)
        if open_at < 0:
            break
        close_at = text.find("</text>", open_at)
        if close_at < 0:
            break
        block = text[open_at:close_at]
        pos = close_at + 7
        body = "".join(part.split("<")[0] for part in block.split(">")[1:])
        if not any(c.isalnum() or ord(c) > 127 for c in body):
            continue                      # an element with nothing in it
        tag_end = block.find(">")
        scale = _transform_scale(block[:tag_end if tag_end > 0 else len(block)])
        at = block.find(key)
        if at < 0:
            continue
        idx = at + len(key)
        rest = block[idx:idx + 48]
        j = 0
        while j < len(rest) and rest[j] in _LEAD_CHARS:
            j += 1
        k = j
        while k < len(rest) and rest[k] in _NUMBER_CHARS:
            k += 1
        if k == j:
            continue
        unit = ""
        m = k
        while m < len(rest) and rest[m].isalpha():
            unit += rest[m]
            m += 1
        try:
            value = float(rest[j:k])
        except ValueError:
            value = None
        if value and unit not in ("", "px"):
            # An absolute unit is stated outright; it still has to be expressed
            # in user units, which only works when the document states its own
            # size.  Without that the value cannot be placed and is skipped.
            absolute = length_pt(rest[j:k] + unit)
            vw, _ = svg_viewbox(text)
            doc_w, _ = svg_size_pt(text)
            value = (absolute * vw / doc_w) if (absolute and vw and doc_w) else None
        if value:
            out.append(value * scale)
    return sorted(out)


def svg_body_font_pt(sizes):
    """The size most of the figure's labels are set in.

    Informational only.  It is NOT what the floor is checked against: which
    size is "the body" depends on the figure -- in an axis-heavy plot the mode
    is the axis labels, but in a diagram full of subscripted symbols the mode
    is the subscripts, and calling those the body would fail a figure whose
    actual labels are large.
    """
    if not sizes:
        return None
    counts = {}
    for s in sizes:
        key = round(s, 1)
        counts[key] = counts.get(key, 0) + 1
    best = max(counts.items(), key=lambda kv: (kv[1], kv[0]))
    return best[0]


def svg_picture_row(slide_no, name, shape, blob, slide_area_pt,
                    scale_tol, aspect_tol, min_area_fraction,
                    min_visible_font_pt):
    """One audit row for a vector picture.

    The authored size comes from the SVG itself, so there is no DPI to be
    missing and no resolution risk to report -- a vector picture cannot be
    pixelated.  What can still go wrong is what goes wrong for a raster: pasted
    at a size other than the authored one, which moves every label's on-page
    size, or stretched, which distorts them.
    """
    import hashlib
    from pptx.util import Emu

    text = blob.decode("utf-8", "replace")
    auth_w, auth_h = svg_size_pt(text)
    view_w, view_h = svg_viewbox(text)
    fonts_units = svg_font_sizes_units(text)
    body_units = svg_body_font_pt(fonts_units)

    disp_pt_w = float(Emu(shape.width).pt)
    disp_pt_h = float(Emu(shape.height).pt)
    disp_cm_w = float(Emu(shape.width).cm)
    disp_cm_h = float(Emu(shape.height).cm)
    area_fraction = (disp_pt_w * disp_pt_h / slide_area_pt) if slide_area_pt else 0.0
    minor = area_fraction < min_area_fraction

    risks = []
    # Scale is reported but NOT flagged on its own: resizing vector artwork is
    # lossless, so the only thing it can break is the on-page size of the text,
    # and that is measured directly below.  (A raster pasted at the wrong size
    # is a defect in itself -- it gets interpolated.)
    scale = (disp_pt_w / auth_w) if (auth_w and auth_h) else None
    aspect_err = 0.0
    nat_w, nat_h = (view_w, view_h) if (view_w and view_h) else (auth_w, auth_h)
    if nat_w and nat_h:
        native_aspect = nat_h / nat_w
        disp_aspect = (disp_pt_h / disp_pt_w) if disp_pt_w else 0.0
        aspect_err = (disp_aspect / native_aspect - 1.0) if native_aspect else 0.0
        if not minor and abs(aspect_err) > aspect_tol:
            risks.append(
                "ASPECT DISTORTED by %+.1f%% -- the picture was stretched; "
                "hold the aspect or re-author the figure." % (aspect_err * 100.0))
    else:
        risks.append(
            "SVG STATES NO SIZE -- authored width UNVERIFIABLE (no width/height "
            "and no viewBox); re-export the figure with an explicit size.")

    # The unit the figure's text is written in cancels here: a label of N user
    # units, in a document view_w units wide, pasted at disp_pt_w points, is
    # N * disp_pt_w / view_w points on the slide.
    per_unit = (disp_pt_w / view_w) if view_w else None
    body_on_page = None
    smallest_on_page = None
    text_check = "measured" if (fonts_units and per_unit) else "not-verifiable"
    if not fonts_units or not per_unit:
        # Reported, never silently passed -- and never flagged either, because
        # the file simply does not say, and a risk nobody can clear is noise.
        pass
    elif body_units * per_unit > _IMPLAUSIBLE_FONT_PT:
        # Nothing on a slide is set this large; the reading is not a size.
        text_check = "not-verifiable"
    else:
        smallest_on_page = fonts_units[0] * per_unit
        body_on_page = body_units * per_unit
        largest_on_page = fonts_units[-1] * per_unit
        # The floor is checked against the LARGEST text in the figure, which is
        # the one claim no one can argue with: if even the biggest label is
        # under the floor, the figure is too small for the room.  Anything
        # finer -- "the body", "the smallest" -- misfires on subscripts and on
        # exponent ticks, and a check that misfires gets ignored.
        if not minor and largest_on_page < min_visible_font_pt - _FONT_PT_TOL:
            risks.append(
                "FIGURE TEXT TOO SMALL -- its LARGEST label is %.1f pt on the "
                "slide (floor %.0f pt): everything in this figure is smaller "
                "than that.  Enlarge the picture (it is pasted %.2f cm wide) "
                "or re-author the figure with larger text."
                % (largest_on_page, min_visible_font_pt, disp_cm_w))

    return {
        "slide": slide_no,
        "shape": name,
        "kind": "svg",
        "pixels": None,
        "dpi": None,
        "authored_cm": round(auth_w / 72.0 * 2.54, 2) if auth_w else None,
        "authored_cm_height": round(auth_h / 72.0 * 2.54, 2) if auth_h else None,
        "displayed_cm": round(disp_cm_w, 2),
        "displayed_pt": round(disp_pt_w, 1),
        "displayed_cm_height": round(disp_cm_h, 2),
        "scale": round(scale, 4) if scale else None,
        "figure_font_units": sorted({round(f, 1) for f in fonts_units}),
        "body_font_units": body_units,
        "text_check": text_check,
        "body_font_on_slide_pt": (round(body_on_page, 1) if body_on_page else None),
        "smallest_font_on_slide_pt": (round(smallest_on_page, 1)
                                      if smallest_on_page else None),
        "largest_font_on_slide_pt": (round(fonts_units[-1] * per_unit, 1)
                                     if (fonts_units and per_unit) else None),
        "effective_dpi": None,
        "aspect_error": round(aspect_err, 4),
        "area_fraction": round(area_fraction, 4),
        "cropped": False,
        "sha1": hashlib.sha1(blob).hexdigest()[:10],
        "minor": bool(minor),
        "risks": risks,
    }
