#include "math_font.h"

#include <algorithm>
#include <cstdio>
#include <cstring>

namespace mtef {
namespace {

/* sfnt is big-endian throughout. */
uint16_t rd16(const uint8_t* p) { return uint16_t(p[0] << 8 | p[1]); }
int16_t  rs16(const uint8_t* p) { return int16_t(rd16(p)); }
uint32_t rd32(const uint8_t* p) {
    return uint32_t(p[0]) << 24 | uint32_t(p[1]) << 16 |
           uint32_t(p[2]) << 8  | uint32_t(p[3]);
}

bool at(const std::vector<uint8_t>& b, size_t off, size_t need) {
    return off + need <= b.size();
}

uint32_t tag(const char* s) {
    return uint32_t(uint8_t(s[0])) << 24 | uint32_t(uint8_t(s[1])) << 16 |
           uint32_t(uint8_t(s[2])) << 8  | uint32_t(uint8_t(s[3]));
}

std::vector<uint8_t> slurp(const std::string& path) {
    std::vector<uint8_t> out;
    FILE* f = nullptr;
    if (fopen_s(&f, path.c_str(), "rb") != 0 || !f) return out;
    std::fseek(f, 0, SEEK_END);
    const long n = std::ftell(f);
    std::fseek(f, 0, SEEK_SET);
    if (n > 0) {
        out.resize(size_t(n));
        if (std::fread(out.data(), 1, out.size(), f) != out.size()) out.clear();
    }
    std::fclose(f);
    return out;
}

/* A MathValueRecord is a design-unit value plus an optional device table; the
 * device table only matters for hinting at a specific ppem, which is not what
 * this layout works in. */
struct Reader {
    const std::vector<uint8_t>& b;
    size_t p = 0;
    bool bad = false;

    int16_t i16() {
        if (!at(b, p, 2)) { bad = true; return 0; }
        const int16_t v = rs16(&b[p]);
        p += 2;
        return v;
    }
    uint16_t u16() {
        if (!at(b, p, 2)) { bad = true; return 0; }
        const uint16_t v = rd16(&b[p]);
        p += 2;
        return v;
    }
    int16_t value() { const int16_t v = i16(); u16(); return v; }  /* + device */
};

}  // namespace

/* ---- the file ----------------------------------------------------------- */

MathFont MathFont::load(const std::string& path, int face) {
    MathFont m;
    const std::vector<uint8_t> b = slurp(path);
    if (b.size() < 12) { m.why_ = "cannot read " + path; return m; }

    /* A collection holds several faces; walk them and take the one asked for,
     * or the first that carries a MATH table. */
    std::vector<uint32_t> offsets;
    if (rd32(&b[0]) == tag("ttcf")) {
        if (!at(b, 12, 4)) { m.why_ = "truncated collection header"; return m; }
        const uint32_t n = rd32(&b[8]);
        for (uint32_t i = 0; i < n && at(b, 12 + i * 4, 4); ++i)
            offsets.push_back(rd32(&b[12 + i * 4]));
    } else {
        offsets.push_back(0);
    }
    if (offsets.empty()) { m.why_ = "no faces"; return m; }

    for (size_t fi = 0; fi < offsets.size(); ++fi) {
        if (face >= 0 && int(fi) != face) continue;
        const uint32_t base = offsets[fi];
        if (!at(b, base + 6, 2)) continue;
        const uint16_t numTables = rd16(&b[base + 4]);

        uint32_t mathOff = 0, mathLen = 0, headOff = 0, cmapOff = 0;
        for (uint16_t i = 0; i < numTables; ++i) {
            const size_t rec = base + 12 + size_t(i) * 16;
            if (!at(b, rec, 16)) break;
            const uint32_t t = rd32(&b[rec]);
            const uint32_t off = rd32(&b[rec + 8]), len = rd32(&b[rec + 12]);
            if (t == tag("MATH")) { mathOff = off; mathLen = len; }
            else if (t == tag("head")) headOff = off;
            else if (t == tag("cmap")) cmapOff = off;
        }
        if (!mathOff && face < 0) continue;      /* keep looking for a math face */

        if (headOff && at(b, headOff + 20, 2)) {
            const uint16_t upem = rd16(&b[headOff + 18]);
            if (upem) m.upem_ = double(upem);
        }
        if (!mathOff || !at(b, mathOff, mathLen ? 10 : 10)) {
            m.why_ = "no MATH table in " + path;
            return m;
        }
        m.parse_math(b, mathOff);
        if (cmapOff) m.parse_cmap(b, cmapOff);
        m.ok_ = true;
        return m;
    }
    m.why_ = "no face with a MATH table in " + path;
    return m;
}

/* ---- MATH --------------------------------------------------------------- */

void MathFont::parse_math(const std::vector<uint8_t>& b, uint32_t off) {
    if (!at(b, off, 10)) return;
    const uint16_t constantsOff = rd16(&b[off + 4]);
    const uint16_t glyphInfoOff = rd16(&b[off + 6]);
    const uint16_t variantsOff  = rd16(&b[off + 8]);

    if (constantsOff) parse_constants(b, off + constantsOff);
    if (variantsOff)  parse_variants(b, off + variantsOff);
    if (glyphInfoOff) parse_italics(b, off + glyphInfoOff);
}

void MathFont::parse_constants(const std::vector<uint8_t>& b, uint32_t off) {
    Reader r{b, off};
    const double em = 1.0 / m_upem();

    c_.scriptPercentScaleDown       = r.i16() / 100.0;
    c_.scriptScriptPercentScaleDown = r.i16() / 100.0;
    c_.delimitedSubFormulaMinHeight = r.u16() * em;
    c_.displayOperatorMinHeight     = r.u16() * em;

    /* The fifty-one MathValueRecords, in the table's order.  Writing them out
     * rather than indexing an array keeps the name at the point of use, which
     * is where a mistake would otherwise be invisible. */
    double* const v[] = {
        &c_.mathLeading, &c_.axisHeight, &c_.accentBaseHeight,
        &c_.flattenedAccentBaseHeight, &c_.subscriptShiftDown,
        &c_.subscriptTopMax, &c_.subscriptBaselineDropMin,
        &c_.superscriptShiftUp, &c_.superscriptShiftUpCramped,
        &c_.superscriptBottomMin, &c_.superscriptBaselineDropMax,
        &c_.subSuperscriptGapMin, &c_.superscriptBottomMaxWithSubscript,
        &c_.spaceAfterScript, &c_.upperLimitGapMin,
        &c_.upperLimitBaselineRiseMin, &c_.lowerLimitGapMin,
        &c_.lowerLimitBaselineDropMin, &c_.stackTopShiftUp,
        &c_.stackTopDisplayStyleShiftUp, &c_.stackBottomShiftDown,
        &c_.stackBottomDisplayStyleShiftDown, &c_.stackGapMin,
        &c_.stackDisplayStyleGapMin, &c_.stretchStackTopShiftUp,
        &c_.stretchStackBottomShiftDown, &c_.stretchStackGapAboveMin,
        &c_.stretchStackGapBelowMin, &c_.fractionNumeratorShiftUp,
        &c_.fractionNumeratorDisplayStyleShiftUp,
        &c_.fractionDenominatorShiftDown,
        &c_.fractionDenominatorDisplayStyleShiftDown,
        &c_.fractionNumeratorGapMin, &c_.fractionNumDisplayStyleGapMin,
        &c_.fractionRuleThickness, &c_.fractionDenominatorGapMin,
        &c_.fractionDenomDisplayStyleGapMin, &c_.skewedFractionHorizontalGap,
        &c_.skewedFractionVerticalGap, &c_.overbarVerticalGap,
        &c_.overbarRuleThickness, &c_.overbarExtraAscender,
        &c_.underbarVerticalGap, &c_.underbarRuleThickness,
        &c_.underbarExtraDescender, &c_.radicalVerticalGap,
        &c_.radicalDisplayStyleVerticalGap, &c_.radicalRuleThickness,
        &c_.radicalExtraAscender, &c_.radicalKernBeforeDegree,
        &c_.radicalKernAfterDegree,
    };
    for (double* d : v) *d = r.value() * em;
    c_.radicalDegreeBottomRaisePercent = r.i16() / 100.0;
}

/* Coverage tables list the glyphs a subtable speaks about, in glyph order. */
void MathFont::read_coverage(const std::vector<uint8_t>& b, uint32_t off,
                             std::vector<uint16_t>& out) {
    if (!at(b, off, 4)) return;
    const uint16_t format = rd16(&b[off]);
    if (format == 1) {
        const uint16_t n = rd16(&b[off + 2]);
        for (uint16_t i = 0; i < n && at(b, off + 4 + size_t(i) * 2, 2); ++i)
            out.push_back(rd16(&b[off + 4 + size_t(i) * 2]));
    } else if (format == 2) {
        const uint16_t n = rd16(&b[off + 2]);
        for (uint16_t i = 0; i < n; ++i) {
            const size_t rec = off + 4 + size_t(i) * 6;
            if (!at(b, rec, 6)) break;
            const uint16_t first = rd16(&b[rec]), last = rd16(&b[rec + 2]);
            for (uint32_t g = first; g <= last; ++g) out.push_back(uint16_t(g));
        }
    }
}

void MathFont::parse_variants(const std::vector<uint8_t>& b, uint32_t off) {
    if (!at(b, off, 10)) return;
    minOverlap_ = rd16(&b[off]) / m_upem();
    const uint16_t vertCov = rd16(&b[off + 2]);
    const uint16_t horizCov = rd16(&b[off + 4]);
    const uint16_t vertCount = rd16(&b[off + 6]);
    const uint16_t horizCount = rd16(&b[off + 8]);

    struct Axis { uint16_t cov, count; size_t list; std::vector<std::pair<uint16_t, Stretch>>* into; };
    const Axis axes[2] = {
        {vertCov,  vertCount,  off + 10,                          &vert_},
        {horizCov, horizCount, off + 10 + size_t(vertCount) * 2,  &horiz_},
    };

    for (const Axis& a : axes) {
        if (!a.cov || !a.count) continue;
        std::vector<uint16_t> glyphs;
        read_coverage(b, off + a.cov, glyphs);
        for (uint16_t i = 0; i < a.count && i < glyphs.size(); ++i) {
            if (!at(b, a.list + size_t(i) * 2, 2)) break;
            const uint16_t rel = rd16(&b[a.list + size_t(i) * 2]);
            if (!rel) continue;
            Stretch s;
            parse_construction(b, off + rel, s);
            if (!s.empty()) a.into->emplace_back(glyphs[i], std::move(s));
        }
        std::sort(a.into->begin(), a.into->end(),
                  [](const std::pair<uint16_t, Stretch>& x,
                     const std::pair<uint16_t, Stretch>& y) {
                      return x.first < y.first;
                  });
    }
}

void MathFont::parse_construction(const std::vector<uint8_t>& b, uint32_t off,
                                  Stretch& s) {
    if (!at(b, off, 4)) return;
    const uint16_t assemblyOff = rd16(&b[off]);
    const uint16_t nVariants = rd16(&b[off + 2]);
    for (uint16_t i = 0; i < nVariants; ++i) {
        const size_t rec = off + 4 + size_t(i) * 4;
        if (!at(b, rec, 4)) break;
        s.variants.emplace_back(rd16(&b[rec]), rd16(&b[rec + 2]) / m_upem());
    }
    if (!assemblyOff) return;

    const uint32_t a = off + assemblyOff;
    if (!at(b, a, 6)) return;
    const uint16_t nParts = rd16(&b[a + 4]);
    for (uint16_t i = 0; i < nParts; ++i) {
        const size_t rec = a + 6 + size_t(i) * 10;
        if (!at(b, rec, 10)) break;
        GlyphPart p;
        p.glyph          = rd16(&b[rec]);
        p.startConnector = rd16(&b[rec + 2]) / m_upem();
        p.endConnector   = rd16(&b[rec + 4]) / m_upem();
        p.fullAdvance    = rd16(&b[rec + 6]) / m_upem();
        p.extender       = (rd16(&b[rec + 8]) & 1) != 0;
        s.assembly.push_back(p);
    }
}

void MathFont::parse_italics(const std::vector<uint8_t>& b, uint32_t off) {
    if (!at(b, off, 2)) return;
    const uint16_t italicOff = rd16(&b[off]);
    if (!italicOff) return;
    const uint32_t t = off + italicOff;
    if (!at(b, t, 4)) return;
    const uint16_t cov = rd16(&b[t]);
    const uint16_t n = rd16(&b[t + 2]);
    std::vector<uint16_t> glyphs;
    read_coverage(b, t + cov, glyphs);
    for (uint16_t i = 0; i < n && i < glyphs.size(); ++i) {
        const size_t rec = t + 4 + size_t(i) * 4;
        if (!at(b, rec, 4)) break;
        italics_.emplace_back(glyphs[i], rs16(&b[rec]) / m_upem());
    }
    std::sort(italics_.begin(), italics_.end());
}

/* ---- cmap --------------------------------------------------------------- */

void MathFont::parse_cmap(const std::vector<uint8_t>& b, uint32_t off) {
    if (!at(b, off, 4)) return;
    const uint16_t n = rd16(&b[off + 2]);
    uint32_t best = 0;
    int bestScore = -1;
    for (uint16_t i = 0; i < n; ++i) {
        const size_t rec = off + 4 + size_t(i) * 8;
        if (!at(b, rec, 8)) break;
        const uint32_t sub = off + rd32(&b[rec + 4]);
        if (!at(b, sub, 2)) continue;
        /* Score by the subtable's FORMAT, not by its platform: Cambria Math
         * advertises a Unicode subtable that is format 14 -- variation
         * sequences, no plain mappings at all -- and choosing by platform
         * picks that empty table over the real one. */
        const uint16_t fmt = rd16(&b[sub]);
        const int score = (fmt == 12) ? 3 : (fmt == 4) ? 2 : -1;
        if (score > bestScore) { bestScore = score; best = sub; }
    }
    if (!best || !at(b, best, 4)) return;

    const uint16_t format = rd16(&b[best]);
    if (format == 4) {
        const uint16_t segX2 = rd16(&b[best + 6]);
        const uint16_t seg = segX2 / 2;
        const size_t endP = best + 14, startP = endP + segX2 + 2;
        const size_t deltaP = startP + segX2, rangeP = deltaP + segX2;
        for (uint16_t s = 0; s < seg; ++s) {
            if (!at(b, rangeP + size_t(s) * 2, 2)) break;
            const uint16_t end = rd16(&b[endP + size_t(s) * 2]);
            const uint16_t start = rd16(&b[startP + size_t(s) * 2]);
            const int16_t delta = rs16(&b[deltaP + size_t(s) * 2]);
            const uint16_t range = rd16(&b[rangeP + size_t(s) * 2]);
            if (start > end || end == 0xFFFF) continue;
            for (uint32_t c = start; c <= end; ++c) {
                uint16_t g = 0;
                if (range == 0) {
                    g = uint16_t(c + delta);
                } else {
                    const size_t gp = rangeP + size_t(s) * 2 + range
                                    + (c - start) * 2;
                    if (!at(b, gp, 2)) continue;
                    g = rd16(&b[gp]);
                    if (g) g = uint16_t(g + delta);
                }
                if (g) { cmapCodes_.push_back(c); cmapGlyphs_.push_back(g); }
            }
        }
    } else if (format == 12) {
        if (!at(b, best + 16, 4)) return;
        const uint32_t groups = rd32(&b[best + 12]);
        for (uint32_t i = 0; i < groups; ++i) {
            const size_t rec = best + 16 + size_t(i) * 12;
            if (!at(b, rec, 12)) break;
            const uint32_t first = rd32(&b[rec]), last = rd32(&b[rec + 4]);
            const uint32_t gid = rd32(&b[rec + 8]);
            if (last < first || last - first > 0x10000) continue;
            for (uint32_t c = first; c <= last; ++c) {
                cmapCodes_.push_back(c);
                cmapGlyphs_.push_back(uint16_t(gid + (c - first)));
            }
        }
    }
    /* Sorted so a lookup is a binary search. */
    std::vector<size_t> idx(cmapCodes_.size());
    for (size_t i = 0; i < idx.size(); ++i) idx[i] = i;
    std::sort(idx.begin(), idx.end(), [&](size_t x, size_t y) {
        return cmapCodes_[x] < cmapCodes_[y];
    });
    std::vector<uint32_t> cs;
    std::vector<uint16_t> gs;
    cs.reserve(idx.size());
    gs.reserve(idx.size());
    for (size_t i : idx) { cs.push_back(cmapCodes_[i]); gs.push_back(cmapGlyphs_[i]); }
    cmapCodes_.swap(cs);
    cmapGlyphs_.swap(gs);
}

/* ---- lookups ------------------------------------------------------------ */

uint16_t MathFont::glyph_for(uint32_t cp) const {
    const auto it = std::lower_bound(cmapCodes_.begin(), cmapCodes_.end(), cp);
    if (it == cmapCodes_.end() || *it != cp) return 0;
    return cmapGlyphs_[size_t(it - cmapCodes_.begin())];
}

static const Stretch* find(const std::vector<std::pair<uint16_t, Stretch>>& v,
                           uint16_t g) {
    const auto it = std::lower_bound(
        v.begin(), v.end(), g,
        [](const std::pair<uint16_t, Stretch>& p, uint16_t x) {
            return p.first < x;
        });
    return (it != v.end() && it->first == g) ? &it->second : nullptr;
}

const Stretch* MathFont::vertical(uint16_t g) const { return find(vert_, g); }
const Stretch* MathFont::horizontal(uint16_t g) const { return find(horiz_, g); }

uint16_t MathFont::vertical_variant(uint16_t g, double target, double* got) const {
    const Stretch* s = vertical(g);
    if (got) *got = 0;
    if (!s) return 0;
    for (const auto& v : s->variants) {
        if (v.second >= target) {
            if (got) *got = v.second;
            return v.first;
        }
    }
    /* Past the largest ready-made size the caller must assemble; report the
     * biggest so a caller that cannot assemble still gets the best fit. */
    if (!s->variants.empty()) {
        if (got) *got = s->variants.back().second;
        return s->variants.back().first;
    }
    return 0;
}

double MathFont::italics_correction(uint16_t g) const {
    const auto it = std::lower_bound(
        italics_.begin(), italics_.end(), g,
        [](const std::pair<uint16_t, double>& p, uint16_t x) {
            return p.first < x;
        });
    return (it != italics_.end() && it->first == g) ? it->second : 0.0;
}

const MathFont& MathFont::cambria() {
    static const MathFont f = [] {
        const char* const roots[] = {
            "C:\\Windows\\Fonts\\cambria.ttc",
            "C:\\Windows\\Fonts\\CAMBRIA.TTC",
        };
        for (const char* p : roots) {
            MathFont m = MathFont::load(p);
            if (m.ok()) return m;
        }
        MathFont bad;
        bad.why_ = "Cambria Math not found";
        return bad;
    }();
    return f;
}

}  // namespace mtef
