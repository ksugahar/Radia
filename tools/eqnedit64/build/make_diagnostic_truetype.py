"""Make a private CFF/glyf A/B asset; not a distribution/release font.

Keep family names for the controlled GDI selection comparison. The output file
has a distinct name; preserve the GUST/LPPL notices and review derived naming
before any publication. This script never registers or renders a font.
"""
import argparse
import hashlib
import json
from pathlib import Path

import fontTools
from fontTools.fontBuilder import FontBuilder
from fontTools.pens.cu2quPen import Cu2QuPen
from fontTools.pens.ttGlyphPen import TTGlyphPen
from fontTools.ttLib import TTFont, newTable


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.output.exists() or args.source.resolve() == args.output.resolve():
        parser.error("output must be a new, separate diagnostic file")
    font = TTFont(args.source, recalcTimestamp=False)
    if "CFF " not in font or "glyf" in font or "CFF2" in font:
        parser.error("expected static CFF source")
    order = font.getGlyphOrder()
    original_tables = {tag: font.getTableData(tag) for tag in ("cmap", "MATH", "GSUB", "GPOS", "hmtx") if tag in font}
    source_metrics = dict(font["hmtx"].metrics)
    glyph_set = font.getGlyphSet()
    glyphs = {}
    for name in order:
        pen = TTGlyphPen(glyph_set)
        glyph_set[name].draw(Cu2QuPen(pen, max_err=0.5, reverse_direction=True))
        glyphs[name] = pen.glyph()
    post_values = {key: getattr(font["post"], key) for key in (
        "italicAngle", "underlinePosition", "underlineThickness", "isFixedPitch",
        "minMemType42", "maxMemType42", "minMemType1", "maxMemType1")}
    del font["CFF "]
    if "VORG" in font:
        del font["VORG"]
    font.sfntVersion = "\x00\x01\x00\x00"
    font["glyf"] = newTable("glyf")
    builder = FontBuilder(font=font)
    builder.setupGlyf(glyphs)
    builder.setupMaxp()
    builder.setupPost(keepGlyphNames=True, **post_values)
    font["head"].flags &= ~2  # Do not claim every rounded xMin equals hmtx LSB.
    args.output.parent.mkdir(parents=True, exist_ok=True)
    font.save(args.output)
    with TTFont(args.output, recalcTimestamp=False) as check:
        assert "glyf" in check and "CFF " not in check and "CFF2" not in check
        assert check.getGlyphOrder() == order
        assert check["hmtx"].metrics == source_metrics
        preserved = {tag: check.getTableData(tag) == value for tag, value in original_tables.items()}
        assert all(preserved.values()), preserved
        report = dict(diagnostic_only=True, source_sha256=digest(args.source),
                      output_sha256=digest(args.output), fonttools=fontTools.__version__,
                      glyphs=len(order), units_per_em=check["head"].unitsPerEm,
                      max_cubic_approximation_error_design_units=0.5,
                      note="Coordinate integer rounding is additional; no visual-equivalence claim.",
                      byte_identical_tables=preserved,
                      family_names_preserved_for_ab=True, license="GUST Font License / LPPL")
    args.output.with_suffix(".json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
