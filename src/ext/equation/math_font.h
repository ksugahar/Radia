/*
 * math_font.h -- the typesetting parameters a math font carries
 *
 * A math font is not just outlines.  Since OpenType 1.6 it also ships a MATH
 * table stating, in the designer's own numbers, how far above a radicand the
 * bar sits, how thick that bar is, where the fraction axis lies, and -- the
 * part no constant can substitute for -- a set of progressively taller radical
 * and brace glyphs, plus parts to assemble one taller still.
 *
 * The layout here used to guess all of that: fifty-nine hand-tuned multiples of
 * the point size, none of them from the font.  It shows.  A radical drawn as
 * the plain character with a rule tacked on cannot meet a tall radicand, and no
 * tuning of the gap reaches that, because the fix is to ask for a bigger
 * radical -- which the font has been offering all along.
 *
 * So this reads the table.  Everything is returned in em, so a caller multiplies
 * by the type size and is done; the font's design units never leak out.
 *
 * References: OpenType MATH table specification; MathML Core, which states the
 * layout rules in terms of exactly these constants; TeXbook Appendix G, whose
 * parameters they descend from.
 */
#ifndef MTEF_MATH_FONT_H
#define MTEF_MATH_FONT_H

#include <cstdint>
#include <string>
#include <vector>

namespace mtef {

/* The MathConstants, in the table's own order.  Named rather than numbered
 * because the layout reads them by name and a wrong index is silent. */
struct MathConstants {
    /* Percentages, as fractions: 0.71 means script is 71% of full size. */
    double scriptPercentScaleDown = 0.71;
    double scriptScriptPercentScaleDown = 0.5;
    /* Everything below is in em. */
    double delimitedSubFormulaMinHeight = 0;
    double displayOperatorMinHeight = 0;
    double mathLeading = 0;
    double axisHeight = 0;
    double accentBaseHeight = 0;
    double flattenedAccentBaseHeight = 0;
    double subscriptShiftDown = 0;
    double subscriptTopMax = 0;
    double subscriptBaselineDropMin = 0;
    double superscriptShiftUp = 0;
    double superscriptShiftUpCramped = 0;
    double superscriptBottomMin = 0;
    double superscriptBaselineDropMax = 0;
    double subSuperscriptGapMin = 0;
    double superscriptBottomMaxWithSubscript = 0;
    double spaceAfterScript = 0;
    double upperLimitGapMin = 0;
    double upperLimitBaselineRiseMin = 0;
    double lowerLimitGapMin = 0;
    double lowerLimitBaselineDropMin = 0;
    double stackTopShiftUp = 0;
    double stackTopDisplayStyleShiftUp = 0;
    double stackBottomShiftDown = 0;
    double stackBottomDisplayStyleShiftDown = 0;
    double stackGapMin = 0;
    double stackDisplayStyleGapMin = 0;
    double stretchStackTopShiftUp = 0;
    double stretchStackBottomShiftDown = 0;
    double stretchStackGapAboveMin = 0;
    double stretchStackGapBelowMin = 0;
    double fractionNumeratorShiftUp = 0;
    double fractionNumeratorDisplayStyleShiftUp = 0;
    double fractionDenominatorShiftDown = 0;
    double fractionDenominatorDisplayStyleShiftDown = 0;
    double fractionNumeratorGapMin = 0;
    double fractionNumDisplayStyleGapMin = 0;
    double fractionRuleThickness = 0;
    double fractionDenominatorGapMin = 0;
    double fractionDenomDisplayStyleGapMin = 0;
    double skewedFractionHorizontalGap = 0;
    double skewedFractionVerticalGap = 0;
    double overbarVerticalGap = 0;
    double overbarRuleThickness = 0;
    double overbarExtraAscender = 0;
    double underbarVerticalGap = 0;
    double underbarRuleThickness = 0;
    double underbarExtraDescender = 0;
    double radicalVerticalGap = 0;
    double radicalDisplayStyleVerticalGap = 0;
    double radicalRuleThickness = 0;
    double radicalExtraAscender = 0;
    double radicalKernBeforeDegree = 0;
    double radicalKernAfterDegree = 0;
    /* A fraction of the radical's height, not an em length. */
    double radicalDegreeBottomRaisePercent = 0.6;
};

/* One piece of an assembled glyph: the font's answer for heights past its
 * largest ready-made variant.  Extenders repeat; the others appear once. */
struct GlyphPart {
    uint16_t glyph = 0;
    double startConnector = 0;   /* em of overlap available at the near end */
    double endConnector = 0;
    double fullAdvance = 0;      /* em along the stretch axis */
    bool extender = false;
};

/* How one character grows.  A radical, brace or parenthesis is not scaled --
 * the font supplies a taller drawing, so the stem thickens as little as the
 * designer intended and the hook still meets the bar. */
struct Stretch {
    /* Ready-made sizes, smallest first; the advance is along the stretch
     * axis, so for a vertical stretch it is the height. */
    std::vector<std::pair<uint16_t, double>> variants;
    std::vector<GlyphPart> assembly;
    bool empty() const { return variants.empty() && assembly.empty(); }
};

class MathFont {
public:
    /* Reads a font file.  A .ttc holds several faces; `face` picks one. */
    static MathFont load(const std::string& path, int face = -1);

    /* The font this editor sets mathematics in, read once. */
    static const MathFont& math();

    bool ok() const { return ok_; }
    const std::string& why_not() const { return why_; }
    const std::string& family() const { return family_; }

    const MathConstants& constants() const { return c_; }

    /* 0 when the font has no glyph for it. */
    uint16_t glyph_for(uint32_t codepoint) const;

    /* Null when the glyph does not stretch. */
    const Stretch* vertical(uint16_t glyph) const;
    const Stretch* horizontal(uint16_t glyph) const;

    /* The smallest variant at least `target` em tall, or 0 with `grew` set
     * false when even the largest falls short -- then assemble instead. */
    uint16_t vertical_variant(uint16_t glyph, double target_em,
                              double* got_em) const;

    /* Extra advance a slanted glyph needs before a following script. */
    double italics_correction(uint16_t glyph) const;

    /* The glyph a font would rather draw at a small size.
     *
     * `ssty` is the OpenType feature every maths font ships for this: the
     * same letter redrawn a little wider and a little heavier so it holds up
     * as a subscript.  TeX applies it; it is why TeX's x in a superscript is
     * 5.44 pt where a simply-shrunk x would be 4.80.  `level` is 1 for script
     * and 2 for scriptscript; 0 comes back when the font offers nothing. */
    uint16_t script_variant(uint16_t glyph, int level) const;

    /* Where an accent hangs from, and where a letter wants one hung: the
     * point, measured from the glyph's origin, that the two are aligned by.
     * Returns `dflt` for a glyph the font says nothing about -- normally half
     * its width, which is what centring means when nothing is stated. */
    double top_accent_attachment(uint16_t glyph, double dflt) const;

    /* The minimum overlap the font asks for between assembled parts. */
    double min_connector_overlap() const { return minOverlap_; }

private:
    double m_upem() const { return upem_; }
    void parse_math(const std::vector<uint8_t>& b, uint32_t off);
    void parse_constants(const std::vector<uint8_t>& b, uint32_t off);
    void parse_variants(const std::vector<uint8_t>& b, uint32_t off);
    void parse_construction(const std::vector<uint8_t>& b, uint32_t off,
                            Stretch& s);
    void parse_italics(const std::vector<uint8_t>& b, uint32_t off);
    void parse_top_accent(const std::vector<uint8_t>& b, uint32_t off);
    void parse_ssty(const std::vector<uint8_t>& b, uint32_t off);
    void parse_cmap(const std::vector<uint8_t>& b, uint32_t off);
    void read_coverage(const std::vector<uint8_t>& b, uint32_t off,
                       std::vector<uint16_t>& out);

    bool ok_ = false;
    std::string why_, family_;
    MathConstants c_;
    double upem_ = 1000.0;
    double minOverlap_ = 0.0;
    std::vector<uint32_t> cmapCodes_;      /* sorted, parallel to cmapGlyphs_ */
    std::vector<uint16_t> cmapGlyphs_;
    std::vector<std::pair<uint16_t, Stretch>> vert_, horiz_;   /* by glyph */
    std::vector<std::pair<uint16_t, double>> italics_;
    std::vector<std::pair<uint16_t, double>> topAccent_;
    /* glyph -> { script, scriptscript } */
    std::vector<std::pair<uint16_t, std::pair<uint16_t, uint16_t>>> ssty_;
};

/* Where that font lives, empty when it could not be found. */
const std::string& math_font_path();

}  // namespace mtef

#endif  /* MTEF_MATH_FONT_H */
