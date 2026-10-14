"""Shared SVG geometry extraction for PowerPoint picture shapes.

PowerPoint stores a vector picture as an ``a:blip`` raster fallback plus an
``asvg:svgBlip`` extension pointing at the real SVG part.  Both the figure
paste-scale audit and the presentation figure-text audit need the vector
original's intrinsic size, and each grew its own copy of the extraction with
different unit handling -- one converting to points, the other returning raw
viewBox numbers.  Two copies of the same measurement can disagree about the
same file, so the parsing lives here once.
"""

from __future__ import annotations

import re
import xml.etree.ElementTree as ET


_RELATIONSHIP_EMBED = (
    "{http://schemas.openxmlformats.org/officeDocument/2006/"
    "relationships}embed"
)

_LENGTH = re.compile(
    r"^\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)"
    r"\s*(px|pt|pc|in|cm|mm)?\s*$",
    re.IGNORECASE,
)

# CSS/SVG user units are 96 per inch; PowerPoint measures in 72 per inch.
_UNIT_TO_PT = {
    "px": 72.0 / 96.0,
    "pt": 1.0,
    "pc": 12.0,
    "in": 72.0,
    "cm": 72.0 / 2.54,
    "mm": 72.0 / 25.4,
}


def svg_length_pt(value: str) -> float:
    """Convert an SVG length to points; a unitless value is CSS pixels."""
    match = _LENGTH.match(str(value or ""))
    if not match:
        return 0.0
    return float(match.group(1)) * _UNIT_TO_PT[(match.group(2) or "px").lower()]


def svg_geometry(blob: bytes) -> dict[str, float]:
    """Return the intrinsic geometry of an SVG document.

    Keys:
        ``width_pt`` / ``height_pt``: intrinsic physical size in points, taken
            from the ``width``/``height`` attributes and completed from the
            viewBox aspect when only one of them is given.  Both are 0.0 for a
            viewBox-only document, which declares an aspect but no physical
            size -- callers that need an authored width must reject that rather
            than assume one user unit is one point.
        ``view_width`` / ``view_height``: raw viewBox extents, or 0.0 when the
            document has no usable viewBox.
    """
    root = ET.fromstring(blob)
    view_box = root.attrib.get("viewBox") or root.attrib.get("viewbox")
    view_width = view_height = 0.0
    if view_box:
        values = [
            float(value)
            for value in re.split(r"[\s,]+", view_box.strip())
            if value
        ]
        if len(values) == 4 and values[2] > 0 and values[3] > 0:
            view_width, view_height = values[2], values[3]

    width_pt = svg_length_pt(root.attrib.get("width", ""))
    height_pt = svg_length_pt(root.attrib.get("height", ""))
    if width_pt <= 0 and height_pt > 0 and view_height > 0:
        width_pt = height_pt * view_width / view_height
    if height_pt <= 0 and width_pt > 0 and view_width > 0:
        height_pt = width_pt * view_height / view_width
    if width_pt <= 0 or height_pt <= 0:
        width_pt = height_pt = 0.0

    return {
        "width_pt": width_pt,
        "height_pt": height_pt,
        "view_width": view_width,
        "view_height": view_height,
    }


def picture_svg_blob(shape) -> bytes | None:
    """Return the SVG part behind a picture shape, or None if it has none."""
    for element in shape._element.iter():
        if not str(element.tag).endswith("}svgBlip"):
            continue
        relationship_id = element.get(_RELATIONSHIP_EMBED)
        if not relationship_id:
            continue
        return bytes(shape.part.related_part(relationship_id).blob)
    return None


__all__ = ["svg_length_pt", "svg_geometry", "picture_svg_blob"]
