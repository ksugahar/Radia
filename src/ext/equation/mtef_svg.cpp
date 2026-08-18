/*
 * mtef_svg.cpp -- MTEF node tree -> SVG renderer (milestone 1)
 *
 * Layout model: every node produces a Layout -- a display list of glyphs and
 * rules positioned relative to its own origin, with the origin ON THE BASELINE
 * at the left edge, plus the box extents (width / ascent / descent).  Parents
 * translate their children.  Nothing is emitted until the whole tree is laid
 * out, so the final viewBox is exact and no second measuring pass is needed.
 *
 * Covered in this milestone: LINE, CHAR, SIZE, SCRIPT (sub/sup/subsup), FENCE,
 * FRACT, ROOT, and the integral / big-operator families (operator glyph plus
 * limits, inline or stacked).  Remaining template classes recurse into their
 * content so nothing silently disappears.
 */
#include "mtef_svg.h"

#include "math_font.h"
#include "mtef_gdi.h"
#include "math_layout.h"
#include "math_writer.h"      /* accent_drawing_char: one table, two spellings */
#include "mtef_parser.h"
#include "tex_parser.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <sstream>
#include <vector>

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX          /* else windows.h's max/min macros eat std::max */
#include <windows.h>
#endif

namespace mtef {

/* The first code point of a UTF-8 string, 0 for an empty one.  Wanted where a
 * glyph is known only by the text it draws. */
uint32_t first_code_of(const std::string& s) {
    if (s.empty()) return 0;
    const unsigned char c0 = (unsigned char)s[0];
    auto tail = [&](size_t i) -> uint32_t {
        return i < s.size() ? ((unsigned char)s[i] & 0x3F) : 0;
    };
    if (c0 < 0x80) return c0;
    if ((c0 & 0xE0) == 0xC0) return ((c0 & 0x1Fu) << 6) | tail(1);
    if ((c0 & 0xF0) == 0xE0) return ((c0 & 0x0Fu) << 12) | (tail(1) << 6) | tail(2);
    return ((c0 & 0x07u) << 18) | (tail(1) << 12) | (tail(2) << 6) | tail(3);
}

/* Pure arithmetic, so it sits outside the Windows-only metric layer: every
 * backend needs the same answer about which face a character belongs to. */
bool is_cjk(uint32_t c) {
    return (c >= 0x3000 && c <= 0x30FF) ||    /* punctuation, kana        */
           (c >= 0x3400 && c <= 0x4DBF) ||    /* ideographs ext A         */
           (c >= 0x4E00 && c <= 0x9FFF) ||    /* ideographs               */
           (c >= 0xF900 && c <= 0xFAFF) ||    /* compatibility            */
           (c >= 0xFF00 && c <= 0xFF60) ||    /* fullwidth forms          */
           (c >= 0xFFE0 && c <= 0xFFE6);
}

namespace {

/* ------------------------------------------------------------------ */
/* UTF-8 / UTF-16 helpers                                              */
/* ------------------------------------------------------------------ */
std::string utf8_of(uint32_t cp) { return mtef_utf8_of(cp); }

std::string xml_escape(const std::string& s) {
    std::string out;
    for (char c : s) {
        switch (c) {
            case '&': out += "&amp;"; break;
            case '<': out += "&lt;"; break;
            case '>': out += "&gt;"; break;
            default: out += c;
        }
    }
    return out;
}

/* ------------------------------------------------------------------ */
/* Font metrics                                                        */
/* ------------------------------------------------------------------ */
#ifdef _WIN32
/* Measure with GDI at a fixed em so the result scales linearly.  Widths come
 * back in 1/1000 em, which is what the SVG consumer will reproduce as long as
 * it resolves the same family. */
constexpr int kEm = 1000;

HFONT make_font(bool italic, bool symbol, bool cjk, int optical) {
    LOGFONTW lf = {};
    lf.lfHeight = -kEm;
    lf.lfItalic = italic ? TRUE : FALSE;
    /* DEFAULT_CHARSET, never SYMBOL_CHARSET.  Latin Modern Math is a Unicode font;
     * asking for the symbol charset makes GDI apply the legacy Symbol code
     * page, so U+0028 is measured as whatever sits at 0x28 in that page and
     * U+2264 is not found at all.  That single flag was behind the spurious
     * gap after "(" and the relation overlapping the fraction bar. */
    lf.lfCharSet = DEFAULT_CHARSET;
    /* Yu Mincho is the serif that sits beside Times without looking borrowed;
     * neither Latin Modern face has a single kana. */
    if (cjk) {
        wcscpy_s(lf.lfFaceName, L"Yu Mincho");
    } else if (symbol) {
        wcscpy_s(lf.lfFaceName, L"Latin Modern Math");
    } else {
        /* The optical cut the size asks for -- "LM Roman 12" at 12 point, "LM
         * Roman 8" for a script.  They are separate drawings and TeX picks
         * between them the same way; setting everything from the 10 pt cut
         * left an operator name 0.3 pt wide. */
        const std::string face = mtef_roman_face(double(optical));
        std::wstring w(face.begin(), face.end());
        wcscpy_s(lf.lfFaceName, w.c_str());
    }
    return CreateFontIndirectW(&lf);
}

struct MetricCache {
    HDC hdc = nullptr;
    std::map<int, HFONT> fonts;              /* key: italic*2 + symbol */
    std::map<std::pair<int, uint32_t>, double> widths;
    std::map<int, std::pair<double, double>> vmetrics;  /* asc, desc per em */

    /* The faces are loaded into this process rather than installed, so they
     * have to be there before anything is measured with them. */
    MetricCache() {
        load_private_fonts();
        hdc = CreateCompatibleDC(nullptr);
    }
    ~MetricCache() {
        for (auto& kv : fonts) DeleteObject(kv.second);
        if (hdc) DeleteDC(hdc);
    }
    /* The key carries the optical size as well as the face, because the text
     * face is EIGHT drawings and which one is wanted depends on the size. */
    HFONT font(int key) {
        auto it = fonts.find(key);
        if (it != fonts.end()) return it->second;
        HFONT f = make_font((key & 1) != 0, (key & 2) != 0, (key & 4) != 0,
                            key >> 3);
        fonts[key] = f;
        return f;
    }
    /* The face follows the character, so a caller cannot measure with one font
     * and have the renderer draw with another. */
    static int face_key(uint32_t cp, bool italic, bool symbol, double pt) {
        return (italic ? 1 : 0) | (symbol ? 2 : 0) | (is_cjk(cp) ? 4 : 0) |
               (mtef_optical_size(pt) << 3);
    }
    double width_em(uint32_t cp, bool italic, bool symbol, double pt) {
        int key = face_key(cp, italic, symbol, pt);
        auto k = std::make_pair(key, cp);
        auto it = widths.find(k);
        if (it != widths.end()) return it->second;
        HGDIOBJ old = SelectObject(hdc, font(key));
        wchar_t buf[3];
        int n = 0;
        if (cp < 0x10000) {
            buf[n++] = wchar_t(cp);
        } else {
            uint32_t v = cp - 0x10000;
            buf[n++] = wchar_t(0xD800 + (v >> 10));
            buf[n++] = wchar_t(0xDC00 + (v & 0x3FF));
        }
        SIZE sz = {};
        GetTextExtentPoint32W(hdc, buf, n, &sz);
        SelectObject(hdc, old);
        double w = double(sz.cx) / kEm;
        widths[k] = w;
        return w;
    }
    std::pair<double, double> vmetric(bool italic, bool symbol, double pt) {
        int key = (italic ? 1 : 0) | (symbol ? 2 : 0) | (mtef_optical_size(pt) << 3);
        auto it = vmetrics.find(key);
        if (it != vmetrics.end()) return it->second;
        HGDIOBJ old = SelectObject(hdc, font(key));
        TEXTMETRICW tm = {};
        GetTextMetricsW(hdc, &tm);
        SelectObject(hdc, old);
        auto v = std::make_pair(double(tm.tmAscent) / kEm,
                                double(tm.tmDescent) / kEm);
        vmetrics[key] = v;
        return v;
    }

    /* Per-glyph ink box, which is what a script has to clear.  The font's own
     * ascent and descent are the wrong measure: Latin Modern Math reserves room for
     * extensible brackets and integral signs, so using its global descent puts
     * the subscript of a sigma a third of a line too low, while the subscript
     * of a Times "B" sits correctly.  TeX has always used per-glyph height and
     * depth for exactly this reason. */
    struct Box { double asc = 0, desc = 0, ink_w = 0; };

    Box glyph_box(uint32_t cp, bool italic, bool symbol, double pt) {
        int key = face_key(cp, italic, symbol, pt);
        auto k = std::make_pair(key, cp);
        auto it = boxes.find(k);
        if (it != boxes.end()) return it->second;

        auto v = vmetric(italic, symbol, pt);
        Box b;
        b.asc = v.first;
        b.desc = v.second;
        b.ink_w = width_em(cp, italic, symbol, pt);
        /* Past U+FFFF, GetGlyphOutlineW takes a glyph index rather than a
         * character -- and the maths font is exactly where the characters
         * past U+FFFF are, since that is the block an italic alphabet lives
         * in.  Skipping the per-glyph box for them left every letter carrying
         * the font's own ascent and descent, which Latin Modern Math sizes for
         * an extensible integral: a plain x came out three and a half ems
         * tall. */
        uint32_t which = cp;
        UINT ggo = GGO_METRICS;
        if (cp >= 0x10000) {
            const mtef::MathFont& mf = mtef::MathFont::math();
            const uint16_t gid = mf.ok() ? mf.glyph_for(cp) : 0;
            if (gid) { which = gid; ggo = GGO_METRICS | GGO_GLYPH_INDEX; }
        }
        if (which < 0x10000) {
            HGDIOBJ old = SelectObject(hdc, font(key));
            GLYPHMETRICS gm = {};
            MAT2 id = {{0, 1}, {0, 0}, {0, 0}, {0, 1}};
            DWORD r = GetGlyphOutlineW(hdc, which, ggo, &gm, 0, nullptr, &id);
            SelectObject(hdc, old);
            if (r != GDI_ERROR && gm.gmBlackBoxY > 0) {
                b.asc = std::max(double(gm.gmptGlyphOrigin.y) / kEm, 0.0);
                b.desc = std::max(
                    double(int(gm.gmBlackBoxY) - gm.gmptGlyphOrigin.y) / kEm, 0.0);
                b.ink_w = std::max(
                    b.ink_w, double(gm.gmptGlyphOrigin.x + int(gm.gmBlackBoxX)) / kEm);
            }
        }
        boxes[k] = b;
        return b;
    }

    /* The same measurement for a glyph the font names rather than a character
     * it maps: a radical variant has no codepoint of its own. */
    Box glyph_box_index(uint16_t gid) {
        auto it = index_boxes.find(gid);
        if (it != index_boxes.end()) return it->second;

        const int key = 2;                       /* the maths face */
        HGDIOBJ old = SelectObject(hdc, font(key));
        Box b;
        int adv = 0;
        /* double(adv), not adv: both were int, so this was integer division.
         * A glyph wider than one em came back as exactly 1.0 -- which is why a
         * summation's display variant measured 12.0 pt against the 17.328 the
         * font gives it -- and anything narrower came back as 0, fell through
         * to the black box below, and got its INK where its ADVANCE was
         * wanted.  Every glyph measured by index went through here: the
         * radical variants, the fences, the large operators. */
        if (GetCharWidthI(hdc, gid, 1, nullptr, &adv))
            b.ink_w = double(adv) / kEm;
        GLYPHMETRICS gm = {};
        MAT2 id = {{0, 1}, {0, 0}, {0, 0}, {0, 1}};
        DWORD r = GetGlyphOutlineW(hdc, gid, GGO_METRICS | GGO_GLYPH_INDEX,
                                   &gm, 0, nullptr, &id);
        SelectObject(hdc, old);
        if (r != GDI_ERROR && gm.gmBlackBoxY > 0) {
            b.asc = double(gm.gmptGlyphOrigin.y) / kEm;
            b.desc = double(int(gm.gmBlackBoxY) - gm.gmptGlyphOrigin.y) / kEm;
            if (b.ink_w <= 0)
                b.ink_w = double(gm.gmptGlyphOrigin.x + int(gm.gmBlackBoxX)) / kEm;
        }
        index_boxes[gid] = b;
        return b;
    }

    std::map<std::pair<int, uint32_t>, Box> boxes;
    std::map<uint16_t, Box> index_boxes;
};

MetricCache& metrics() {
    static MetricCache c;
    return c;
}

double char_width(uint32_t cp, double sizePt, bool italic, bool symbol) {
    return metrics().width_em(cp, italic, symbol, sizePt) * sizePt;
}
/* The font's own extent, for things sized against the face rather than
 * against one glyph (fence stretching, the fallback line height). */
void char_vmetrics(double sizePt, bool italic, bool symbol,
                   double& asc, double& desc) {
    auto v = metrics().vmetric(italic, symbol, sizePt);
    asc = v.first * sizePt;
    desc = v.second * sizePt;
}
/* The ink box of one glyph, which is what neighbours and scripts must clear. */
void glyph_vmetrics(uint32_t cp, double sizePt, bool italic, bool symbol,
                    double& asc, double& desc) {
    auto b = metrics().glyph_box(cp, italic, symbol, sizePt);
    asc = b.asc * sizePt;
    desc = b.desc * sizePt;
}
/* How far the drawing actually reaches.  A large operator is drawn wider than
 * it advances, so laying the next atom out at the advance alone lets a sigma
 * touch the symbol after it. */
double glyph_ink_width(uint32_t cp, double sizePt, bool italic, bool symbol) {
    return metrics().glyph_box(cp, italic, symbol, sizePt).ink_w * sizePt;
}
#else
/* Portable fallback: enough to keep the tree walking and the tests honest
 * about being unmeasured, not enough for production output. */
double char_width(uint32_t, double sizePt, bool, bool) { return 0.5 * sizePt; }
void char_vmetrics(double sizePt, bool, bool, double& asc, double& desc) {
    asc = 0.75 * sizePt;
    desc = 0.25 * sizePt;
}
void glyph_vmetrics(uint32_t, double sizePt, bool, bool,
                    double& asc, double& desc) {
    asc = 0.70 * sizePt;
    desc = 0.05 * sizePt;
}
double glyph_ink_width(uint32_t, double sizePt, bool, bool) { return 0.5 * sizePt; }
#endif

/* ------------------------------------------------------------------ */
/* Typeface mapping                                                    */
/* ------------------------------------------------------------------ */
bool typeface_is_italic(int tf) { return tf == 3 || tf == 4; }   /* VARIABLE, LCGREEK */

/* Which family draws a code point.  This is a property of the character, not
 * of the MTEF typeface: a TF_SYMBOL "(" is still an ordinary parenthesis and
 * belongs in the text face, while anything past Latin-1 needs the math face. */
bool needs_math_face(uint32_t cp) { return cp > 0xFF; }

/* Which face a character belongs to, given what it IS as well as what it is.
 *
 * Everything in an equation comes from the maths font -- digits and upright
 * function names included -- because that is where TeX takes them from, and
 * because it is the only face that carries the `ssty` alternates a script
 * size wants.  Only TF_TEXT, which is prose the author asked for with \text,
 * belongs in the text face.
 *
 * A digit was the last thing still coming from the text face, and it showed:
 * the i in x_{i} picked up its script alternate and the 2 in x^{2} did not,
 * leaving that one a half point narrow where its neighbour was exact. */
bool math_face_for(uint32_t cp, int typeface) {
    /* Prose the author asked for with \text. */
    if (typeface == TF_TEXT) return false;
    /* An operator NAME is set in the text roman, not in the maths font's own
     * upright.  LaTeX spells this \operator@font and it reads as \mathrm: TeX
     * sets \sin at 14.424 pt, which is \mathrm{sin} exactly, where the maths
     * font's own upright -- \symup{sin} -- comes to 14.820.  Sending names to
     * the maths face along with the digits made every function name that much
     * wide, and it looked like kerning inside the name until the two were
     * measured side by side. */
    if (typeface == TF_FUNCTION) return false;
    if (is_cjk(cp)) return false;
    return true;
}

/* The Mathematical Alphanumeric Symbols block: where a maths font keeps its
 * italic letters.  Returns 0 for anything that is not one, which leaves
 * digits, punctuation and text upright alone.
 *
 * U+1D455 -- italic small h -- was never assigned, because Planck's constant
 * already had U+210E.  A font has the glyph there and nowhere else, so a
 * naive a + (c - 'a') draws nothing for h. */
uint32_t math_italic_of(uint32_t cp) {
    if (cp >= 'A' && cp <= 'Z') return 0x1D434 + (cp - 'A');
    if (cp == 'h')              return 0x210E;
    if (cp >= 'a' && cp <= 'z') return 0x1D44E + (cp - 'a');
    if (cp >= 0x3B1 && cp <= 0x3C9) return 0x1D6FC + (cp - 0x3B1);  /* alpha..omega */
    return 0;
}

/* ------------------------------------------------------------------ */
/* TeX atom classes and the spacing between them                       */
/* ------------------------------------------------------------------ */
enum AtomClass { kOrd, kOp, kBin, kRel, kOpen, kClose, kPunct, kInner };

AtomClass class_of_char(uint32_t cp) {
    switch (cp) {
        /* The ASCII forms as well as the typographic ones.  A person types
         * "a*b" and "a|b" on the keyboard; U+2217 and U+2225 arrive from the
         * palette.  Only the palette forms were listed, so an asterisk typed
         * at the keyboard was set as an ordinary letter with no space around
         * it at all -- measured against Equation Editor, 0.03 em where it
         * gives 0.16. */
        case '+': case '-': case '*': case 0x2212: case 0x00B1: case 0x2213:
        case 0x00D7: case 0x00F7: case 0x22C5: case 0x2217: case 0x2218:
        case 0x2229: case 0x222A: case 0x2227: case 0x2228: case 0x2295:
        case 0x2297: case 0x2299: case 0x2296: case 0x228E: case 0x2216:
            return kBin;
        case '=': case '<': case '>': case '|': case 0x2260: case 0x2264:
        case 0x2265:
        case 0x2248: case 0x2261: case 0x223C: case 0x2243: case 0x2245:
        case 0x221D: case 0x22A5: case 0x2225: case 0x2208: case 0x2209:
        case 0x2282: case 0x2283: case 0x2286: case 0x2287: case 0x2192:
        case 0x2190: case 0x2194: case 0x21D2: case 0x21D0: case 0x21D4:
        case 0x2262:
            return kRel;
        case '(': case '[': case '{': case 0x27E8: case 0x230A: case 0x2308:
            return kOpen;
        case ')': case ']': case '}': case 0x27E9: case 0x230B: case 0x2309:
            return kClose;
        case ',': case ';': case ':':
            return kPunct;
        default:
            return kOrd;
    }
}

/* TeX's table, in mu (18 mu = 1 em).  Entries marked "text only" in TeX are
 * kept unconditionally here: equations in a document are set in text style or
 * larger, and dropping them inside scripts costs more than it saves. */
int space_mu(AtomClass l, AtomClass r) {
    /* TeX's thin, medium and thick spaces.
     *
     * Equation Editor 3.1 sets tighter: measured off its screen it leaves
     * 0.219 em around a relation, 0.156 around a binary operator and 0.109
     * after a comma, which are 4, 3 and 2 eighteenths -- one less than TeX's
     * at every step, and a fit within three thousandths of an em.  This was
     * briefly set to those.
     *
     * It is back to TeX's because appearance follows TeX here and usability
     * follows Equation Editor.  Nearly everyone who reads an equation has
     * read TeX's spacing far more often than Equation Editor's, and the
     * openness is most of why TeX-set mathematics reads as it does. */
    static const int kThin = 3, kMed = 4, kThick = 5;
    switch (l) {
        case kOrd:
            if (r == kOp) return kThin;
            if (r == kBin) return kMed;
            if (r == kRel) return kThick;
            if (r == kInner) return kThin;
            return 0;
        case kOp:
            if (r == kOrd || r == kOp) return kThin;
            if (r == kBin) return kMed;
            if (r == kRel) return kThick;
            if (r == kInner) return kThin;
            return 0;
        case kBin:
            if (r == kClose || r == kPunct) return 0;
            return kMed;
        case kRel:
            if (r == kRel || r == kClose || r == kPunct) return 0;
            return kThick;
        case kOpen:
            return 0;
        case kClose:
            if (r == kOp) return kThin;
            if (r == kBin) return kMed;
            if (r == kRel) return kThick;
            if (r == kInner) return kThin;
            return 0;
        case kPunct:
            return kThin;
        case kInner:
            if (r == kBin) return kMed;
            if (r == kRel) return kThick;
            if (r == kClose) return 0;
            return kThin;
    }
    return 0;
}

/* Operator glyph for the integral / big-operator template families. */
uint32_t bigop_glyph(int selector) {
    switch (selector) {
        case tmSINT:  return 0x222B;   /* integral */
        case tmDINT:  return 0x222C;
        case tmTINT:  return 0x222D;
        case tmSSINT: return 0x222E;   /* contour */
        case tmDSINT: return 0x222F;
        case tmTSINT: return 0x2230;
        case tmSUM: case tmISUM:       return 0x2211;
        case tmPROD: case tmIPROD:     return 0x220F;
        case tmCOPROD: case tmICOPROD: return 0x2210;
        case tmUNION: case tmIUNION:   return 0x22C3;
        case tmINTER: case tmIINTER:   return 0x22C2;
        default: return 0x2211;
    }
}

/* Fence glyphs: {left, right} */
std::pair<uint32_t, uint32_t> fence_glyphs(int selector) {
    switch (selector) {
        case tmANGLE: return {0x27E8, 0x27E9};
        case tmPAREN: return {'(', ')'};
        case tmBRACE: return {'{', '}'};
        case tmBRACK: return {'[', ']'};
        case tmBAR:   return {'|', '|'};
        case tmDBAR:  return {0x2016, 0x2016};
        case tmFLOOR: return {0x230A, 0x230B};
        case tmCEIL:  return {0x2308, 0x2309};
        default:      return {'(', ')'};
    }
}

/* TeX's \scriptspace: the padding after a subscript or superscript.  A
 * length, not a proportion -- it does not scale with the type size, because
 * what it protects against is two pieces of ink touching. */
const double kScriptSpace = 0.5;

/* ------------------------------------------------------------------ */
/* Renderer                                                            */
/* ------------------------------------------------------------------ */
class Renderer {
public:
    explicit Renderer(const SvgStyle& s) : st_(s) {}

    Layout run(const LineNode& root) {
        return layout_list(root.children, st_.full, std::string());
    }

private:
    const SvgStyle& st_;

    /* A script is set at the percentage the FONT states, not at the size the
     * Equation Editor dialog happens to default to.  Latin Modern Math says 70
     * and 50; the dialog said 58 and 42, which is where a superscript came out
     * an eighth too small against TeX -- and TeX reads the same two numbers.
     * (Style follows the font; the dialog stays as it is for usability.) */
    double script_pt() const {
        const mtef::MathConstants& mc = mtef::MathFont::math().constants();
        return st_.full * mc.scriptPercentScaleDown;
    }
    double script_script_pt() const {
        const mtef::MathConstants& mc = mtef::MathFont::math().constants();
        return st_.full * mc.scriptScriptPercentScaleDown;
    }

    double size_of(int sizeType) const {
        switch (sizeType) {
            case SIZETYPE_SUB:    return script_pt();
            case SIZETYPE_SUB2:   return script_script_pt();
            case SIZETYPE_SYM:    return st_.sym;
            case SIZETYPE_SUBSYM: return st_.subsym;
            default:              return st_.full;
        }
    }
    /* One step down for scripts, floored at the sub-subscript size. */
    /* 0 for full size, 1 for script, 2 for scriptscript -- which of the
     * font's three drawings of a character is wanted. */
    int script_level(double pt) const {
        if (pt > script_pt() + 1e-9) return 0;
        if (pt > script_script_pt() + 1e-9) return 1;
        return 2;
    }

    double script_size(double cur) const {
        if (cur > script_pt() + 1e-9) return script_pt();
        return script_script_pt();
    }

    /* What a fraction sets its numerator and denominator in.
     *
     * TeX steps display and text style down to script, and script down to
     * script-script, so the first fraction stays full size and every fraction
     * inside one is set smaller.  Without that a fraction of a fraction just
     * grows taller and taller, which is what this editor did.
     *
     * Equation Editor 3.1 does NOT do this.  Its palette carries TWO fraction
     * templates side by side -- "Full-size vertical fraction" and
     * "Reduced-size vertical fraction" -- and the person writing picks which.
     * That works when the input is a template chosen by hand; it cannot work
     * here, where the input is LaTeX and rac is one command.  So the size
     * has to come from the nesting, and TeX's rule is the one that does. */
    double frac_child_size(double cur) const {
        return fracDepth_ == 0 ? cur : script_size(cur);
    }

    /* How many fractions enclose the one being laid out.  The step has to be
     * counted, not inferred from the type size: display and text style are the
     * SAME size in TeX, so an inner fraction at full size is indistinguishable
     * from the outer one by size alone -- which is how a first attempt at this
     * left every level full-sized. */
    int fracDepth_ = 0;

    Layout glyph_layout(uint32_t cp, double sizePt, bool italic, bool symbol) {
        /* A variable is set in the MATHS font's own italic alphabet, not in
         * the text italic that happens to look like it.  They are different
         * fonts with different widths, and the difference is not small: TeX
         * sets x at 6.864 pt where the text italic advances 5.570, because
         * TeX takes the letter from latinmodern-math.otf like everything else
         * and this took it from LM Roman 10 Italic.
         *
         * Digits and punctuation stay where they are -- upright, and the same
         * width in both faces. */
        if (italic) {
            const uint32_t m = math_italic_of(cp);
            if (m) { cp = m; italic = false; symbol = true; }
        }

        Layout L;
        Glyph g;
        g.x = 0; g.y = 0; g.size = sizePt;
        g.italic = italic; g.symbol = symbol; g.cjk = is_cjk(cp);
        g.text = utf8_of(cp);

        const mtef::MathFont& mf = mtef::MathFont::math();
        uint16_t gid = (symbol && mf.ok()) ? mf.glyph_for(cp) : 0;

        /* At a script size the font would rather draw a different letter.
         *
         * `ssty` is what every maths font ships for this: the same character
         * redrawn wider and heavier so it holds up small.  Simply shrinking
         * the full-size drawing, which is what this did, leaves every
         * subscript and superscript a few per cent narrow, and it compounds
         * when scripts nest -- TeX sets x^{y^{z}} a tenth wider than this
         * managed.  The tell was that "=" scaled at exactly 0.7 where every
         * letter scaled at 0.76 to 0.82: Latin Modern Math has no ssty
         * alternate for "=", so that one atom really is just shrunk. */
        if (gid) {
            const int level = script_level(sizePt);
            if (const uint16_t alt = mf.script_variant(gid, level)) gid = alt;
        }

        if (gid && gid != mf.glyph_for(cp)) {
            /* A named alternate has no character of its own, so it is drawn
             * and measured by index. */
            const MetricCache::Box b = metrics().glyph_box_index(gid);
            g.glyph_id = gid;
            L.w = b.ink_w * sizePt;
            L.asc = b.asc * sizePt;
            L.desc = b.desc * sizePt;
        } else {
            L.w = char_width(cp, sizePt, italic, symbol);
            glyph_vmetrics(cp, sizePt, italic, symbol, L.asc, L.desc);
        }

        /* The same table states the italic correction, and TeX appends it
         * after every maths character (Appendix G / TeX 755).  It is what
         * makes "ab" 11.664 pt rather than 11.496: 0.168 of it is the kern
         * after the b. */
        if (gid) L.w += mf.italics_correction(gid) * sizePt;

        L.glyphs.push_back(g);
        return L;
    }

    /* An accent over an expression: \vec{B}, \hat{n}, \dot{x}.
     *
     * This used to fall through to the default and produce nothing, so a
     * vector rendered as a bare letter -- the picture of \vec{B} was byte for
     * byte the picture of B.  The markup path was always right, which is why
     * it went unnoticed: Office drew the arrow, we did not. */
    Layout layout_embell(const EmbellNode& em, double sizePt,
                         const std::string& listPath, int child) {
        /* node_slots(kEmbell) = { content } */
        Layout body = layout_list(em.content, sizePt,
                                  slot_path(listPath, child, 0));
        const uint32_t acc = accent_drawing_char(em.embellType);

        Layout out;
        out.absorb(body, 0, 0);
        out.w = body.w;
        out.asc = body.asc;
        out.desc = body.desc;

        const mtef::MathConstants& mc = mtef::MathFont::math().constants();

        /* A bar is a rule the width of what it covers; a fixed-width macron
         * over a wide expression would look wrong. */
        if (!acc) {
            const double gap   = mc.overbarVerticalGap * sizePt;
            const double thick = mc.overbarRuleThickness * sizePt;
            const double extra = mc.overbarExtraAscender * sizePt;
            Rule r;
            r.x = 0;
            r.y = -(body.asc + gap + thick);
            r.w = body.w;
            r.h = thick;
            out.rules.push_back(r);
            out.asc = body.asc + gap + thick + extra;
            return out;
        }

        /* TeX's make_math_accent, which is not what was here.
         *
         *     delta = min(height(nucleus), accent_base_height)
         *     the accent is set at FULL size, its advance discarded, and
         *     lowered onto the nucleus by delta
         *
         * so the mark OVERLAPS what it sits on rather than standing clear
         * above it, and how far down it comes is capped at the font's
         * accentBaseHeight -- 0.45 em here.  That cap is the idea: over a
         * short letter the accent follows the letter, over a tall one it
         * stops, so a dot over an x and a dot over a B end up the same
         * distance apart as the letters are tall, up to a limit.
         *
         * Setting it at 0.62 of the size with a fixed 0.08 em of clear air, as
         * this did, made a dot over an x more than a third too tall.
         *
         * Checked against TeX: a dot over an x comes to 8.124 pt, which is the
         * dot's own 8.124 less delta 5.304 plus the x's 5.304. */
        Layout mark = glyph_layout(acc, sizePt, false, needs_math_face(acc));
        const double delta = std::min(body.asc, mc.accentBaseHeight * sizePt);

        /* Line the accent's attachment point up with the base's, which is
         * what the MATH table's two entries are for.  Half the width is only
         * the right answer when the font says nothing -- and it is never the
         * right answer for a combining mark, whose ink sits a quarter of an em
         * to the LEFT of an origin it advances nothing from.
         *
         * Using the base's own attachment rather than half its width is also
         * what leans the accent over an italic letter, which is the "skew" in
         * TeX's make_math_accent. */
        const mtef::MathFont& mf = mtef::MathFont::math();
        double baseAttach = body.w / 2.0;
        if (mf.ok() && body.glyphs.size() == 1 && body.glyphs.front().symbol) {
            const Glyph& only = body.glyphs.front();
            uint16_t g = only.glyph_id;
            if (!g) g = mf.glyph_for(first_code_of(only.text));
            if (g)
                baseAttach = mf.top_accent_attachment(g, body.w / 2.0 / sizePt)
                           * sizePt;
        }
        double markAttach = 0.0;
        if (mf.ok()) {
            if (const uint16_t g = mf.glyph_for(acc))
                markAttach = mf.top_accent_attachment(g, mark.w / 2.0 / sizePt)
                           * sizePt;
        }
        const double dx = baseAttach - markAttach;
        /* The mark's ink foot lands (body.asc - delta) above the baseline. */
        const double dy = -(body.asc - delta + mark.desc);
        out.absorb(mark, dx, dy);
        out.asc = std::max(out.asc, mark.asc + mark.desc - delta + body.asc);
        return out;
    }

    /* \overline / \underline: a rule above or below, spanning the content. */
    Layout layout_bar(const NodeList& content, bool above, double sizePt,
                      const std::string& listPath, int child) {
        Layout body = layout_list(content, sizePt,
                                  slot_path(listPath, child, 0));
        Layout out;
        out.absorb(body, 0, 0);
        out.w = body.w;
        out.asc = body.asc;
        out.desc = body.desc;

        /* The font states the gap, the rule and the band of clear space
         * outside it, separately for over and under.  These were 0.10 and 0.05
         * of the type size, picked by eye, and left the under case nearly a
         * fifth shallower than TeX sets it. */
        const mtef::MathConstants& mc = mtef::MathFont::math().constants();
        const double gap   = (above ? mc.overbarVerticalGap
                                    : mc.underbarVerticalGap) * sizePt;
        const double thick = (above ? mc.overbarRuleThickness
                                    : mc.underbarRuleThickness) * sizePt;
        const double extra = (above ? mc.overbarExtraAscender
                                    : mc.underbarExtraDescender) * sizePt;
        Rule r;
        r.x = 0;
        r.w = body.w;
        r.h = thick;
        if (above) {
            r.y = -(body.asc + gap + thick);
            out.asc = body.asc + gap + thick + extra;
        } else {
            r.y = body.desc + gap;
            out.desc = body.desc + gap + thick + extra;
        }
        out.rules.push_back(r);
        return out;
    }

    /* The atom class a node contributes to the spacing between its
     * neighbours.  A structure is Inner; a character carries its own class. */
    static AtomClass class_of(const Node& n) {
        switch (n.tag()) {
            case Node::kChar: {
                const auto& c = static_cast<const CharNode&>(n);
                uint32_t cp = c.charCode ? c.charCode : uint32_t(uint8_t(c.ch));
                return class_of_char(cp);
            }
            case Node::kIntegral:
            case Node::kBigOp:
                return kOp;
            case Node::kScript: {
                /* A scripted atom keeps the class of what is being scripted:
                 * \log_{2} is still an operator and still takes a thin space
                 * before its operand, and \sum_{i} beside its limits is still
                 * an operator.  Treating the whole thing as ordinary lost
                 * exactly one thin space from \log_{2} n. */
                const auto& sc = static_cast<const ScriptNode&>(n);
                for (const auto& b : sc.base)
                    if (b && b->tag() != Node::kSize) return class_of(*b);
                return kOrd;
            }
            case Node::kFence:
                return kInner;
            case Node::kFrac:
                return kInner;
            case Node::kLine: {
                /* A group takes the class of its first atom, so \sin(x)
                 * spaces like a function and not like a bare group. */
                const auto& l = static_cast<const LineNode&>(n);
                for (const auto& c : l.children) {
                    if (!c || c->tag() == Node::kSize) continue;
                    /* A function NAME is an operator, and the class belongs to
                     * the whole name rather than to its letters: TeX leaves a
                     * thin space after \sin and none between the s and the i.
                     * Classing the letters individually would space the name
                     * apart from the inside. */
                    if (c->tag() == Node::kChar &&
                        static_cast<const CharNode&>(*c).typeface == TF_FUNCTION)
                        return kOp;
                    return class_of(*c);
                }
                return kOrd;
            }
            default:
                return kOrd;
        }
    }

    /* One step down: into child `child` of the current slot, then into that
     * child's slot `slot` -- the same spelling Equation::caret() uses. */
    static std::string slot_path(const std::string& listPath, int child, int slot) {
        std::string step = std::to_string(child) + "." + std::to_string(slot);
        return listPath.empty() ? step : listPath + "/" + step;
    }

    Layout layout_list(const NodeList& list, double sizePt,
                       const std::string& path) {
        Layout out;
        double x = 0, cur = sizePt;
        bool have_prev = false;
        AtomClass prev = kOrd;
        int child = 0;
        size_t own_first = 0;               /* the stops this call records */

        out.stops.push_back({path, 0, 0.0, 0.0, 0.0});

        for (const auto& n : list) {
            if (!n) continue;
            if (n->tag() == Node::kSize) {
                cur = size_of(static_cast<const SizeNode*>(n.get())->sizeType);
                continue;
            }
            AtomClass cls = class_of(*n);
            /* TeX's rule: a binary operator with nothing to bind on its left
             * is not binary.  Without this, "-x" is set as if it were a
             * subtraction and opens with a gap. */
            if (cls == kBin && (!have_prev || prev == kBin || prev == kOp ||
                                prev == kRel || prev == kOpen || prev == kPunct))
                cls = kOrd;

            if (have_prev) {
                int mu = space_mu(prev, cls);
                /* TeX 766: the medium and thick spaces are inserted only in
                 * display and text styles -- in script and scriptscript they
                 * are dropped, and only the thin space survives.  Scaling
                 * them down with the type size instead, which is what this
                 * did, is not the same thing: it made the limits under a
                 * summation wider than the summation sign, where TeX sets
                 * them narrower, and put the whole construct a third too
                 * wide. */
                if (mu > 3 && cur < st_.full - 1e-9) mu = 0;
                x += mu * cur / 18.0;
            }

            Layout piece = layout_node(*n, cur, path, child);
            out.absorb(piece, x, 0);
            x += piece.w;
            out.asc = std::max(out.asc, piece.asc);
            out.desc = std::max(out.desc, piece.desc);
            prev = cls;
            have_prev = true;
            ++child;
            out.stops.push_back({path, child, x, 0.0, 0.0});
        }
        out.w = x;

        /* The caret spans the slot, which is only known now.  An empty slot has
         * no extent at all, so it gets the current type size -- otherwise the
         * caret in a fresh template would be invisible, which is exactly where
         * it is most needed. */
        double top = -out.asc, bottom = out.desc;
        if (out.asc <= 0 && out.desc <= 0) {
            top = -0.70 * sizePt;
            bottom = 0.20 * sizePt;

            /* An empty slot takes up room and says so.  Without a width a
             * fresh fraction is a bar with nothing above or below it, and the
             * editor gives no sign that there are two places to type. */
            out.w = st_.empty_slot_em * sizePt;
            out.asc = -top;
            out.desc = bottom;
            SlotBox b;
            b.x = 0;
            b.y = top;
            b.w = out.w;
            b.h = bottom - top;
            out.empty_slots.push_back(b);
        }
        for (size_t k = own_first; k < out.stops.size(); ++k)
            if (out.stops[k].path == path) {
                out.stops[k].top = top;
                out.stops[k].bottom = bottom;
            }
        return out;
    }

    Layout layout_node(const Node& n, double sizePt,
                       const std::string& listPath, int child) {
        switch (n.tag()) {
            case Node::kLine: {
                const auto& ln = static_cast<const LineNode&>(n);
                if (ln.isNull) return Layout();
                /* node_slots(kLine) = { children } */
                return layout_list(ln.children, sizePt,
                                   slot_path(listPath, child, 0));
            }
            case Node::kChar: {
                const auto& c = static_cast<const CharNode&>(n);
                uint32_t cp = c.charCode ? c.charCode : uint32_t(uint8_t(c.ch));
                if (!cp) return Layout();
                return glyph_layout(cp, sizePt, typeface_is_italic(c.typeface),
                                    math_face_for(cp, c.typeface));
            }
            case Node::kEmbell:
                return layout_embell(static_cast<const EmbellNode&>(n), sizePt,
                                     listPath, child);
            case Node::kMatrix: {
                const auto& m = static_cast<const MatrixNode&>(n);
                return layout_grid(m.elements, m.rows, m.cols, sizePt,
                                   listPath, child);
            }
            case Node::kPile: {
                /* A gathered or aligned stack is a grid one column wide (or
                 * ncols wide when it is aligned), so it lays out the same
                 * way. */
                const auto& pl = static_cast<const PileNode&>(n);
                const int cols = std::max(1, pl.ncols);
                const int rows = int((pl.lines.size() + cols - 1) / cols);
                return layout_grid(pl.lines, rows, cols, sizePt,
                                   listPath, child, pl.halign);
            }
            case Node::kDecoration: {
                const auto& d = static_cast<const DecorationNode&>(n);
                /* \overline / \underline: a rule the width of what it covers. */
                return layout_bar(d.content, d.selector == tmOBAR, sizePt,
                                  listPath, child);
            }
            case Node::kScript:
                return layout_script(static_cast<const ScriptNode&>(n), sizePt,
                                     listPath, child);
            case Node::kFence:
                return layout_fence(static_cast<const FenceNode&>(n), sizePt,
                                    listPath, child);
            case Node::kFrac:
                return layout_frac(static_cast<const FracNode&>(n), sizePt,
                                   listPath, child);
            case Node::kSqrt:
                return layout_sqrt(static_cast<const SqrtNode&>(n), sizePt,
                                   listPath, child);
            case Node::kIntegral: {
                const auto& i = static_cast<const IntegralNode&>(n);
                return layout_bigop(bigop_glyph(i.selector), i.body, i.lower, i.upper,
                                    i.hasLower, i.hasUpper, i.hasLimits, sizePt,
                                    listPath, child);
            }
            case Node::kBigOp: {
                const auto& b = static_cast<const BigOpNode&>(n);
                return layout_bigop(bigop_glyph(b.selector), b.body, b.lower, b.upper,
                                    b.hasLower, b.hasUpper, b.hasLimits, sizePt,
                                    listPath, child);
            }
            default:
                return layout_fallback(n, sizePt, listPath, child);
        }
    }

    /* node_slots(kScript) = { base } + (hasSub ? sub) + (hasSup ? sup) */
    /* Limits above and below, everything centred on the widest of the three.
     * The same arithmetic the large operators use -- Appendix G rule 13 --
     * so a limit under \lim sits where a limit under \sum does. */
    Layout stack_limits(const Layout& op, const NodeList& lower,
                        const NodeList& upper, bool hasLower, bool hasUpper,
                        double sizePt, const std::string& lp, int c,
                        int lowerSlot, int upperSlot) {
        const mtef::MathConstants& mc = mtef::MathFont::math().constants();
        const double ss = script_size(sizePt);
        Layout up, lo;
        if (hasUpper) up = layout_list(upper, ss, slot_path(lp, c, upperSlot));
        if (hasLower) lo = layout_list(lower, ss, slot_path(lp, c, lowerSlot));

        double w = op.w;
        if (hasUpper) w = std::max(w, up.w);
        if (hasLower) w = std::max(w, lo.w);

        /* Appendix G rule 13, in Knuth's own shape:
         *
         *     shift_up   = big_op_spacing3 - depth(upper)   floored at ...1
         *     shift_down = big_op_spacing4 - height(lower)  floored at ...2
         *
         * with a band of big_op_spacing5 outside each limit.  The floor is on
         * the SHIFT, not on the finished position, which stops being the same
         * thing as soon as a limit has any depth of its own -- the form this
         * replaces took the maximum after adding the depth in.
         *
         * The five spacings are MATH constants under names that do not match:
         * big_op_spacing1 is upperLimitBaselineRiseMin, ...2 is
         * lowerLimitGapMin, ...3 is upperLimitGapMin, ...4 is
         * lowerLimitBaselineDropMin -- crosswise, checked against
         * \the\fontdimen9..12\textfont3.
         *
         * big_op_spacing5 is 0.12 em.  \the\fontdimen13\textfont3 answers 0.1,
         * and that does not reproduce what TeX actually sets, so the number
         * here is measured from the boxes instead: asked for \sum^{N},
         * \sum_{n=1} and \lim_{x \to 0} separately, all three want 1.4397 pt
         * against a 12 point body, and with it the two shifts come out at
         * exactly big_op_spacing3 and big_op_spacing2.  Something in the
         * OpenType path is not reading fontdimen 13; three independent boxes
         * are the better witness than one parameter that disagrees with them. */
        const double up1 = mc.upperLimitBaselineRiseMin * sizePt;
        const double lo2 = mc.lowerLimitGapMin * sizePt;
        const double up3 = mc.upperLimitGapMin * sizePt;
        const double lo4 = mc.lowerLimitBaselineDropMin * sizePt;
        const double extra = 0.12 * sizePt;

        Layout out;
        out.absorb(op, (w - op.w) / 2.0, 0);
        out.w = w;
        out.asc = op.asc;
        out.desc = op.desc;
        if (hasUpper) {
            const double shift = std::max(up3 - up.desc, up1);
            const double lift = op.asc + shift + up.desc;
            out.absorb(up, (w - up.w) / 2.0, -lift);
            out.asc = std::max(out.asc, lift + up.asc + extra);
        }
        if (hasLower) {
            const double shift = std::max(lo4 - lo.asc, lo2);
            const double drop = op.desc + shift + lo.asc;
            out.absorb(lo, (w - lo.w) / 2.0, drop);
            out.desc = std::max(out.desc, drop + lo.desc + extra);
        }
        return out;
    }

    Layout layout_script(const ScriptNode& s, double sizePt,
                         const std::string& lp, int c) {
        const int subSlot = 1;
        const int supSlot = s.hasSub ? 2 : 1;

        /* \lim and its family put their scripts UNDER and OVER in display
         * style, not beside.  The name is already in the tree as a run of
         * TF_FUNCTION characters, so both this and the writers read it back
         * through the same two functions -- there is no node for it, and no
         * way for the picture and the paste to disagree about where the limit
         * went. */
        if (fracDepth_ == 0 &&
            mtef_name_takes_limits(function_name_of(s.base))) {
            Layout base = layout_list(s.base, sizePt, slot_path(lp, c, 0));
            return stack_limits(base, s.sub, s.sup, s.hasSub, s.hasSup,
                                sizePt, lp, c, subSlot, supSlot);
        }

        Layout base = layout_list(s.base, sizePt, slot_path(lp, c, 0));
        double ss = script_size(sizePt);

        /* The base's italic correction belongs to the SUPERSCRIPT when there
         * is a subscript, and to the width when there is not (TeX 755):
         *
         *     if subscript is empty and delta != 0 then append kern(delta)
         *     ... make_scripts(q, delta)   { delta shifts the superscript }
         *
         * A leaning letter needs the room above its lean, not below it.
         * Adding it to the width regardless left c_{i} 0.3 pt wide -- exactly
         * the italic correction of a c -- and that turned up again at every
         * level of anything built on it. */
        const mtef::MathFont& mfont = mtef::MathFont::math();
        double delta = 0;
        if (mfont.ok() && !base.glyphs.empty() && base.glyphs.back().symbol) {
            /* Only a glyph actually taken from the MATHS font has an italic
             * correction there.  Looking one up for a text-face letter finds
             * the maths font's own drawing of the same character and takes
             * ITS correction, which is not the same letter: that quietly
             * shortened \log_{2} by the 0.156 pt correction of a maths g. */
            const Glyph& last = base.glyphs.back();
            const uint16_t gid = last.glyph_id
                               ? last.glyph_id
                               : mfont.glyph_for(first_code_of(last.text));
            if (gid) delta = mfont.italics_correction(gid) * last.size;
        }
        if (s.hasSub && delta > 0) base.w -= delta;

        Layout out = base;
        double x = base.w;
        /* Where a script sits is stated by the font, and TeX reads the same
         * two numbers: 0.363 em up and 0.247 em down here.  Measured against
         * TeX, its superscript sits at exactly superscriptShiftUp and its
         * subscript at exactly subscriptShiftDown -- so the guesses this
         * replaces (0.45 of the type size, and the base's own extents scaled
         * by 0.35) were near enough to look right and wrong by a tenth.
         *
         * The rest is clearance, in the order MathML Core applies it: a
         * script starts at the font's shift, is pushed further by a tall or
         * deep base, and is pushed further again if it would otherwise reach
         * past the limits the font sets. */
        const mtef::MathFont& mf = mtef::MathFont::math();
        const mtef::MathConstants& mc = mf.constants();

        Layout sup, sub;
        if (s.hasSup) sup = layout_list(s.sup, ss, slot_path(lp, c, supSlot));
        if (s.hasSub) sub = layout_list(s.sub, ss, slot_path(lp, c, subSlot));

        double supShift = 0, subShift = 0;
        /* The two DROPS are measured at the SCRIPT size, not the type size:
         * Knuth writes sup_drop(t) and sub_drop(t) with t already stepped down
         * to script_size (Appendix G rule 18).  Scaling them by the full size
         * instead put a limit beside an integral nine tenths of a point too
         * high and seven tenths too low -- invisible on a letter, because for
         * a single-character nucleus the drop never binds, and plain on an
         * integral, whose box is tall. */
        if (s.hasSup) {
            supShift = std::max(mc.superscriptShiftUp * sizePt,
                                mc.superscriptBottomMin * sizePt + sup.desc);
            supShift = std::max(supShift,
                                base.asc - mc.superscriptBaselineDropMax * ss);
        }
        if (s.hasSub) {
            subShift = std::max(mc.subscriptShiftDown * sizePt,
                                sub.asc - mc.subscriptTopMax * sizePt);
            subShift = std::max(subShift,
                                base.desc + mc.subscriptBaselineDropMin * ss);
        }
        if (s.hasSup && s.hasSub) {
            /* The two must not close up on each other, and the superscript
             * must not ride so high that the pair looks unattached. */
            const double gap = (supShift - sup.desc) - (sub.asc - subShift);
            const double want = mc.subSuperscriptGapMin * sizePt;
            if (gap < want) {
                subShift += want - gap;
                /* Opening that gap by lowering the subscript alone drops the
                 * pair away from the base.  The font caps how low the
                 * superscript's foot may sit, so whatever room is left under
                 * that cap is taken by raising BOTH -- the gap is kept and the
                 * pair stays attached. */
                const double room = mc.superscriptBottomMaxWithSubscript * sizePt
                                  - (supShift - sup.desc);
                if (room > 0) {
                    supShift += room;
                    subShift -= room;
                }
            }
        }

        double wsub = 0, wsup = 0;
        if (s.hasSup) {
            const double dx = s.hasSub ? delta : 0.0;
            out.absorb(sup, x + dx, -supShift);
            out.asc = std::max(out.asc, supShift + sup.asc);
            wsup = dx + sup.w;
        }
        if (s.hasSub) {
            out.absorb(sub, x, subShift);
            out.desc = std::max(out.desc, subShift + sub.desc);
            wsub = sub.w;
        }
        /* TeX pads a script by \scriptspace so the next thing along does not
         * touch it -- 0.5 pt, added to the script box itself (Appendix G,
         * make_scripts).  Every x^2 in the comparison was exactly this much
         * narrower than TeX's. */
        out.w = x + std::max(wsub, wsup) + (s.hasSub || s.hasSup ? kScriptSpace : 0.0);
        return out;
    }

    /* node_slots(kFence) = { content } */
    /* A delimiter at the height the content needs.
     *
     * The font draws parentheses at a series of ready-made sizes and, past the
     * largest, ships pieces to assemble one of any height.  Scaling the small
     * one instead -- which is what this did -- thickens the stem in proportion
     * and, worse, keeps the SMALL glyph's advance: a parenthesis round a
     * fraction was drawn 25 pt tall and 4.5 pt wide, where TeX gives the same
     * delimiter 7.8 pt of width because the tall drawing IS wider.  That is
     * the same mistake the summation had, and the lesson the radical had
     * already learned. */
    Layout stretched_glyph(uint32_t cp, double needPt, double sizePt) {
        const mtef::MathFont& mf = mtef::MathFont::math();
        uint16_t chosen = 0;
        double gotEm = 0;
        if (mf.ok()) {
            const uint16_t base = mf.glyph_for(cp);
            if (const mtef::Stretch* st = mf.vertical(base)) {
                for (const auto& v : st->variants) {
                    const MetricCache::Box b = metrics().glyph_box_index(v.first);
                    chosen = v.first;
                    gotEm = b.asc + b.desc;      /* measured, not the record */
                    if (gotEm * sizePt >= needPt) break;
                }
            } else {
                chosen = base;
            }
        }
        if (!chosen)
            return glyph_layout(cp, sizePt, false, needs_math_face(cp));

        const MetricCache::Box b = metrics().glyph_box_index(chosen);
        Layout g;
        Glyph gg;
        gg.size = sizePt;
        gg.symbol = true;
        gg.glyph_id = chosen;
        gg.text = utf8_of(cp);
        g.glyphs.push_back(gg);
        g.w = b.ink_w * sizePt;
        g.asc = b.asc * sizePt;
        g.desc = b.desc * sizePt;
        /* Past the tallest drawing the font offers, the font ships PARTS to
         * build one of any height, and that is what is used.  Scaling the
         * largest instead -- which is what happened here -- stretches the
         * whole drawing including the hook and the stem, and it showed: a
         * radical over a fraction with a summation in it came out twice as
         * deep as TeX sets it. */
        if (gotEm > 0 && gotEm * sizePt < needPt) {
            if (Layout built = assembled_glyph(cp, needPt, sizePt); !built.glyphs.empty())
                return built;
            /* No parts either: scaling is all that is left. */
            const double k = needPt / std::max(gotEm * sizePt, 1e-6);
            for (auto& x : g.glyphs) x.stretchY = k;
            g.asc *= k;
            g.desc *= k;
        }
        return g;
    }

    /* A delimiter built from the font's own pieces.
     *
     * The font names a bottom, a top, sometimes a middle, and one or more
     * EXTENDERS that repeat as many times as the height needs.  Consecutive
     * pieces overlap, by at least minConnectorOverlap and at most what the
     * two connectors allow, so the joins do not show.
     *
     * Sizing follows MathML Core: take the fewest repeats that reach the
     * target with minimum overlap, then open the overlaps evenly to land ON
     * the target rather than past it. */
    Layout assembled_glyph(uint32_t cp, double needPt, double sizePt) {
        Layout out;
        const mtef::MathFont& mf = mtef::MathFont::math();
        if (!mf.ok()) return out;
        const mtef::Stretch* st = mf.vertical(mf.glyph_for(cp));
        if (!st || st->assembly.empty()) return out;

        const double minOv = mf.min_connector_overlap() * sizePt;
        std::vector<mtef::GlyphPart> parts;
        int repeats = 0;
        for (; repeats < 64; ++repeats) {
            parts.clear();
            for (const mtef::GlyphPart& p : st->assembly) {
                const int n = p.extender ? repeats : 1;
                for (int i = 0; i < n; ++i) parts.push_back(p);
            }
            if (parts.size() < 2) continue;
            double total = 0;
            for (const mtef::GlyphPart& p : parts) total += p.fullAdvance * sizePt;
            total -= minOv * double(parts.size() - 1);
            if (total >= needPt) break;
        }
        if (parts.size() < 2) return out;

        /* Open every join by the same amount, within what the connectors
         * allow, so the assembly is exactly as tall as it was asked for. */
        double total = 0;
        for (const mtef::GlyphPart& p : parts) total += p.fullAdvance * sizePt;
        const size_t joins = parts.size() - 1;
        double overlap = (total - needPt) / double(joins);
        double maxOv = 1e9;
        for (size_t i = 0; i + 1 < parts.size(); ++i)
            maxOv = std::min(maxOv, std::min(parts[i].endConnector,
                                             parts[i + 1].startConnector) * sizePt);
        overlap = std::max(minOv, std::min(overlap, maxOv));

        /* The assembly list runs bottom to top, so it is laid out upwards
         * from the baseline and reported as all ascent. */
        double y = 0;
        for (const mtef::GlyphPart& p : parts) {
            const double adv = p.fullAdvance * sizePt;
            Glyph g;
            g.size = sizePt;
            g.symbol = true;
            g.glyph_id = p.glyph;
            g.text = utf8_of(cp);
            g.x = 0;
            /* Each piece is drawn on its own baseline; the metric cache knows
             * where its ink sits relative to that. */
            const MetricCache::Box b = metrics().glyph_box_index(p.glyph);
            g.y = -(y + b.desc * sizePt);
            out.glyphs.push_back(g);
            y += adv - overlap;
        }
        out.asc = y + overlap;          /* the last piece added no join */
        out.desc = 0;
        out.w = 0;
        for (const mtef::GlyphPart& p : parts)
            out.w = std::max(out.w, metrics().glyph_box_index(p.glyph).ink_w * sizePt);
        return out;
    }

    Layout layout_fence(const FenceNode& f, double sizePt,
                        const std::string& lp, int c) {
        Layout inner = layout_list(f.content, sizePt, slot_path(lp, c, 0));
        auto gl = fence_glyphs(f.selector);

        const mtef::MathFont& mf = mtef::MathFont::math();
        const mtef::MathConstants& mc = mf.constants();
        const double axis = mc.axisHeight * sizePt;

        /* TeX's make_left_right decides the size from how far the content
         * reaches past the AXIS, not from its total height, because a
         * delimiter is set symmetrically about the axis and has to cover the
         * worse of the two sides on both:
         *
         *     delta2 = max(height - axis, depth + axis)
         *     delta  = max(delta2 * delimiterfactor / 500,
         *                  delta2 * 2 - delimitershortfall)
         *
         * with \delimiterfactor 901 and \delimitershortfall 5 pt.  The second
         * term is what stops a very tall content asking for a delimiter far
         * bigger than the font has: it may fall short, by up to 5 pt. */
        const double delta2 = std::max(inner.asc - axis, inner.desc + axis);
        const double need = std::max(delta2 * (901.0 / 500.0), delta2 * 2.0 - 5.0);

        Layout out;
        double x = 0;
        const bool left = (f.variation == 0 || f.variation == 1);
        const bool right = (f.variation == 0 || f.variation == 2);

        /* Centred on the axis, which is what makes a bracket look upright
         * beside a fraction instead of sitting a little low. */
        auto place = [&](uint32_t cp) {
            Layout g = stretched_glyph(cp, need, sizePt);
            const double half = (g.asc + g.desc) / 2.0;
            const double shift = g.asc - (half + axis);
            out.absorb(g, x, shift);
            out.asc = std::max(out.asc, g.asc - shift);
            out.desc = std::max(out.desc, g.desc + shift);
            x += g.w;
        };

        if (left) place(gl.first);
        out.absorb(inner, x, 0);
        out.asc = std::max(out.asc, inner.asc);
        out.desc = std::max(out.desc, inner.desc);
        x += inner.w;
        if (right) place(gl.second);

        out.w = x;
        return out;
    }

    /* node_slots(kFrac) = { numer, denom } */
    Layout layout_frac(const FracNode& f, double sizePt,
                       const std::string& lp, int c) {
        const double kidPt = frac_child_size(sizePt);
        ++fracDepth_;
        Layout num = layout_list(f.numer, kidPt, slot_path(lp, c, 0));
        Layout den = layout_list(f.denom, kidPt, slot_path(lp, c, 1));
        --fracDepth_;

        const mtef::MathFont& mf = mtef::MathFont::math();
        const mtef::MathConstants& mc = mf.constants();

        /* The axis is where a fraction bar sits and where a minus sign is
         * centred, so everything else hangs off it.  The guess this replaces
         * was 0.28 of the type size against the font's 0.2856 -- close enough
         * that nothing moves, which is the point: it is now a fact. */
        const double axis  = mc.axisHeight * sizePt;
        const double thick = mc.fractionRuleThickness * sizePt;

        /* TeX puts a null delimiter on each side of a fraction --
         * 
ulldelimiterspace, 1.2 pt -- so the BOX is wider than the parts,
         * but make_fraction sets the rule itself to width(x), exactly as wide
         * as the wider part.  The two are separate numbers and were conflated
         * here: widening the rule as well drew a bar that stuck out a point
         * past its own numerator, which is a fifth more black than TeX puts on
         * the page and the single largest disagreement the ink comparison
         * found on a plain fraction.
         *
         * (Equation Editor 3.1 does overhang, by one point on each side, and
         * that is what this used to imitate.  Appearance follows TeX.) */
        const double sidebearing = 1.2;

        /* Then the two parts are pushed apart until they clear the bar by at
         * least the font's minimum gap.  Equation Editor states this as
         * "Numerator height 35%" and "Denominator depth 100%" against its own
         * reference; the font states it as a target shift plus a floor, which
         * is the same idea with the floor made explicit. */
        /* Display style for the fraction a reader meets first, text style for
         * any inside it.  TeX makes the same distinction, and the font carries
         * two sets of numbers for it: the display gap is nearly twice the text
         * one.  Using the display set at every level was what still left a
         * fraction inside a fraction a fifth taller than TeX sets it, even
         * after its contents were stepped down in size. */
        const bool display = (fracDepth_ == 0);
        double shiftUp = (display ? mc.fractionNumeratorDisplayStyleShiftUp
                                  : mc.fractionNumeratorShiftUp) * sizePt;
        double shiftDown = (display ? mc.fractionDenominatorDisplayStyleShiftDown
                                    : mc.fractionDenominatorShiftDown) * sizePt;
        const double gapNum = (display ? mc.fractionNumDisplayStyleGapMin
                                       : mc.fractionNumeratorGapMin) * sizePt;
        const double gapDen = (display ? mc.fractionDenomDisplayStyleGapMin
                                       : mc.fractionDenominatorGapMin) * sizePt;

        shiftUp   = std::max(shiftUp,   axis + thick / 2.0 + gapNum + num.desc);
        shiftDown = std::max(shiftDown, -axis + thick / 2.0 + gapDen + den.asc);

        const double inner = std::max(num.w, den.w);
        const double w = inner + 2.0 * sidebearing;

        Layout out;
        out.w = w;
        out.absorb(num, (w - num.w) / 2.0, -shiftUp);
        out.absorb(den, (w - den.w) / 2.0, shiftDown);

        Rule bar;
        bar.x = sidebearing;
        bar.y = -axis - thick / 2.0;
        bar.w = inner;
        bar.h = thick;
        out.rules.push_back(bar);

        out.asc  = shiftUp + num.asc;
        out.desc = shiftDown + den.desc;
        return out;
    }

    /* Rows on a fixed pitch, columns as wide as their widest cell, the whole
     * grid centred on the axis.
     *
     * Measured against TeX, which reports the same three numbers for a matrix
     * of any width: one row 14.5 pt tall, two 29.0, three 43.5 -- a pitch of
     * 14.5 pt against a 12 point body, which is the array strut, and the box
     * always centred so that height minus depth is twice the axis.  The
     * columns are separated by 10 pt, which is 2\arraycolsep, checked on
     * "a & b & c": 17.16 pt of letters plus two gaps makes TeX's 37.16.
     *
     * Both are written here as fractions of the type size rather than as the
     * fixed dimensions TeX uses, so a matrix inside a script does not carry
     * full-size gaps.  At 12 point they are TeX's numbers exactly.
     *
     * The pitch is a floor, not a fixed step: TeX lets a tall cell overflow
     * its row and collide with the next, which is a well-known wart of array
     * and not worth reproducing. */
    Layout layout_grid(const NodeList& cells, int rows, int cols,
                       double sizePt, const std::string& lp, int c,
                       int halign = 0) {
        Layout out;
        if (rows <= 0 || cols <= 0) return out;

        const double pitch  = (14.5 / 12.0) * sizePt;
        /* A matrix separates its columns; an alignment does not -- the & is a
         * seam, and "a" and "= b" must meet there or the equals sign floats. */
        const double colGap = (halign == 1) ? 0.0 : (10.0 / 12.0) * sizePt;
        const double axis =
            mtef::MathFont::math().constants().axisHeight * sizePt;

        std::vector<Layout> cell(size_t(rows) * cols);
        std::vector<double> colW(cols, 0.0);
        std::vector<double> rowAsc(rows, 0.0), rowDesc(rows, 0.0);
        for (int r = 0; r < rows; ++r) {
            for (int k = 0; k < cols; ++k) {
                const size_t i = size_t(r) * cols + k;
                if (i < cells.size() && cells[i]) {
                    NodeList one;                    /* layout_list wants a list */
                    const Node& n = *cells[i];
                    if (n.tag() == Node::kLine) {
                        /* A cell is set in TEXT style, not display: TeX's
                         * array does that, and it is why a fraction in a
                         * matrix is a size smaller than the same fraction
                         * standing on its own. */
                        ++fracDepth_;
                        cell[i] = layout_list(
                            static_cast<const LineNode&>(n).children, sizePt,
                            slot_path(lp, c, int(i)));
                        --fracDepth_;
                    }
                }
                colW[k] = std::max(colW[k], cell[i].w);
                rowAsc[r] = std::max(rowAsc[r], cell[i].asc);
                rowDesc[r] = std::max(rowDesc[r], cell[i].desc);
            }
        }

        double total = 0;
        std::vector<double> rowY(rows, 0.0);
        for (int r = 0; r < rows; ++r) {
            const double step = std::max(pitch, rowAsc[r] + rowDesc[r]);
            rowY[r] = total + std::max(rowAsc[r], step - rowDesc[r] -
                                                  (step - rowAsc[r] - rowDesc[r]) / 2.0);
            total += step;
        }
        /* Put the middle of the stack on the axis. */
        const double top = -(total / 2.0 + axis);

        double w = 0;
        for (int k = 0; k < cols; ++k) w += colW[k] + (k ? colGap : 0.0);

        /* Where a cell sits in its column.  A matrix centres; a stack of
         * lines can be asked for flush left or flush right, which is the
         * Format menu's Align Left / Center / Right. */
        auto offset = [&](int col, double colWidth, double cellWidth) {
            switch (halign) {
                case 1:
                    /* An alignment: odd columns flush right, even flush left,
                     * so the & falls on one vertical line and the = beside it
                     * lines up down the page.  That is what align does. */
                    return (col % 2 == 0) ? colWidth - cellWidth : 0.0;
                case 2:  return 0.0;                            /* left  */
                case 3:  return colWidth - cellWidth;           /* right */
                default: return (colWidth - cellWidth) / 2.0;   /* centre */
            }
        };

        for (int r = 0; r < rows; ++r) {
            double x = 0;
            for (int k = 0; k < cols; ++k) {
                const size_t i = size_t(r) * cols + k;
                out.absorb(cell[i], x + offset(k, colW[k], cell[i].w),
                           top + rowY[r]);
                x += colW[k] + colGap;
            }
        }
        out.w = w;
        out.asc = total / 2.0 + axis;
        out.desc = total / 2.0 - axis;
        return out;
    }

    /* node_slots(kSqrt) = { content } + (hasIndex ? index) */
    Layout layout_sqrt(const SqrtNode& s, double sizePt,
                       const std::string& lp, int c) {
        Layout inner = layout_list(s.content, sizePt, slot_path(lp, c, 0));

        /* The font states all four of these; they used to be 0.04 and 0.18 of
         * the type size, picked by eye. */
        const mtef::MathFont& mf = mtef::MathFont::math();
        const mtef::MathConstants& mc = mf.constants();
        /* Display style takes the wider clearance.  TeX: clr is
         * default_rule_thickness + x_height/4 in display and 1.25 *
         * default_rule_thickness otherwise, and the font's two radical gaps
         * are those same two quantities -- 0.148 em against 0.050.  Only the
         * narrow one was ever used, so a displayed root sat a tenth of an em
         * too close to what it covers. */
        const double gap   = (fracDepth_ == 0 ? mc.radicalDisplayStyleVerticalGap
                                              : mc.radicalVerticalGap) * sizePt;
        const double thick = mc.radicalRuleThickness * sizePt;
        const double extra = mc.radicalExtraAscender * sizePt;

        /* A radical is not a character scaled up.  The font draws one at six
         * heights, so ask for the one that reaches -- the stroke stays the
         * weight the designer drew and the hook still meets the bar. */
        const double need = inner.asc + inner.desc + gap + thick;

        /* Choose the variant by what it MEASURES, not by what it advertises.
         *
         * A variant record states an advance along the stretch axis, which is
         * not the same as the glyph's ink: measured, this font's radicals are
         * shorter than their records claim.  Selecting on the record picked a
         * size larger than needed every time -- a root over a fraction came
         * out a sixth taller than TeX's, and TeX picks from the same list. */
        /* One routine picks the drawing, for the radical and for every fence:
         * the smallest ready-made size whose measured ink reaches, and past
         * the largest, the font's own pieces assembled.  This used to have its
         * own copy that scaled instead, which is why a radical over a
         * fraction with a summation in it came out twice as deep as TeX's. */
        Layout sign = stretched_glyph(0x221A, need, sizePt);

        /* The index sits ON the radical's left arm, not beside it.
         *
         * plain.tex's \root builds its box at \scriptscriptstyle -- two steps
         * down, not one -- and then pulls it back over the sign with two
         * kerns:  \mkern-5mu <index> \mkern-4mu <radical>.  A font states the
         * same two kerns in its MATH table, and this one does
         * (radicalKernBeforeDegree, radicalKernAfterDegree); they were read
         * and then never used.
         *
         * Giving the index its own full width at script size instead, which is
         * what this did, put a cube root a quarter wider than TeX sets one and
         * left the 3 stranded in clear space to the left of the arm it belongs
         * on -- by far the worst equation in the ink comparison. */
        Layout idx;
        double idxX = 0, signX = 0;
        if (s.hasIndex) {
            idx = layout_list(s.index, script_size(script_size(sizePt)),
                              slot_path(lp, c, 1));
            const double before = mc.radicalKernBeforeDegree * sizePt;
            const double after  = mc.radicalKernAfterDegree * sizePt;
            idxX  = std::max(0.0, before);
            signX = std::max(0.0, idxX + idx.w + after);
        }

        Layout out;

        /* The sign's FOOT sits with the radicand's, and the bar follows its
         * top.
         *
         * A font offers a radical at a few ready-made heights -- 1.0, 1.2,
         * 1.8, 2.4 em here -- so the one that reaches is usually taller than
         * needed.  Placing its top on the bar sent all of that excess below
         * the line: over a fraction the box came out two thirds deeper than
         * TeX's, and the sign dangled past what it was covering.
         *
         * TeX puts the excess above instead (Appendix G rule 11 widens the
         * clearance by it), which is why its root over a fraction is only a
         * fifteenth of a point deeper than the fraction alone.  Measured on
         * the same equation: TeX 8.43 against the fraction's own 8.36. */
        /* TeX's make_radical, which settles what happens to the slack.
         *
         *     delta = depth(y) - (height(x) + depth(x) + clr)
         *     if delta > 0 then clr = clr + half(delta)
         *
         * The font offers a radical at a few fixed heights, so the one that
         * reaches is nearly always taller than needed.  HALF that excess goes
         * into the clearance, widening the gap above the radicand; the other
         * half is left to hang below the line.  Sending all of it one way --
         * either way -- was a guess, and both guesses were wrong. */
        const double slack = (sign.asc + sign.desc)
                           - (inner.asc + inner.desc + gap + thick);
        const double clr = gap + (slack > 0 ? slack / 2.0 : 0.0);

        const double barTop = -(inner.asc + clr + thick);
        const double signBaseline = barTop + sign.asc;

        out.absorb(sign, signX, signBaseline);
        out.absorb(inner, signX + sign.w, 0);
        out.asc = std::max(-barTop + extra, inner.asc);
        /* The sign is set with its baseline BELOW the equation's, so that its
         * top reaches the bar; how far below is exactly how much of it hangs
         * under the line, and that has to be added to its own depth rather
         * than subtracted from it.  Getting the sign wrong made the box a
         * sixth shallower than the radical drawn in it -- measured against
         * TeX, which puts 0.86 pt of a 12 pt root below the baseline. */
        out.desc = std::max(signBaseline + sign.desc, inner.desc);

        Rule bar;
        bar.x = signX + sign.w;
        bar.y = barTop;
        bar.w = inner.w;
        bar.h = thick;
        out.rules.push_back(bar);

        out.w = signX + sign.w + inner.w;
        if (s.hasIndex) {
            /* The font says how far up the degree sits, as a fraction of the
             * radical's own height. */
            const double lift = mc.radicalDegreeBottomRaisePercent * (-barTop);
            out.absorb(idx, idxX, -lift);
            out.asc = std::max(out.asc, lift + idx.asc);
            /* A long index can reach past the sign it is kerned onto. */
            out.w = std::max(out.w, idxX + idx.w);
        }
        return out;
    }

    /* node_slots(kIntegral / kBigOp) =
     *     (hasLower ? lower) + (hasUpper ? upper) + { body } */
    Layout layout_bigop(uint32_t glyph, const NodeList& body,
                        const NodeList& lower, const NodeList& upper,
                        bool hasLower, bool hasUpper, bool stacked, double sizePt,
                        const std::string& lp, int c) {
        const int lowerSlot = 0;
        const int upperSlot = hasLower ? 1 : 0;
        const int bodySlot = (hasLower ? 1 : 0) + (hasUpper ? 1 : 0);
        const double opSize = st_.sym * (sizePt / std::max(st_.full, 1e-6));
        double ss = script_size(sizePt);
        Layout out;
        double x = 0;

        const mtef::MathFont& mf = mtef::MathFont::math();
        const mtef::MathConstants& mc = mf.constants();

        /* A summation or an integral set for display is not the text glyph at
         * a larger size: the font draws a taller one, and states in
         * displayOperatorMinHeight how tall it must be.  Asking for it is the
         * same move the radical makes, and it is why an integral came out a
         * seventh shorter than TeX sets one. */
        Layout op;
        const uint16_t baseGid = mf.ok() ? mf.glyph_for(glyph) : 0;
        double gotEm = 0;
        /* Which drawing was used, and the size it was drawn at.  The italic
         * correction below belongs to THAT glyph: an integral's display
         * variant leans further than the text one, and states 7.09 pt of
         * correction against the text glyph's 3.98. */
        uint16_t drawnGid = 0;
        double drawnPt = sizePt;
        /* How tall the font says a display operator must be, and nothing
         * else.  This was floored at Equation Editor's own ratio -- its Size
         * dialog sets symbols at 18 point against a 12 point body -- which
         * overrode the font with 1.5 where it asks for 1.3, and picked a
         * variant two sizes too large.  An integral came out a seventh taller
         * than TeX sets one.  Appearance follows TeX; the font is TeX's. */
        /* displayOperatorMinHeight is what its name says: the size a DISPLAY
         * operator must reach.  In text style the operator stays the size the
         * font draws it at, which is why a summation in a denominator is small
         * -- asking for the display size everywhere made one four and a half
         * points too wide inside a radical. */
        const double opTargetEm = (fracDepth_ == 0) ? mc.displayOperatorMinHeight : 0.0;
        const uint16_t bigGid =
            baseGid ? mf.vertical_variant(baseGid, opTargetEm, &gotEm) : 0;
        if (bigGid && bigGid != baseGid) {
            const MetricCache::Box b = metrics().glyph_box_index(bigGid);
            Glyph g;
            g.size = sizePt;
            g.symbol = true;
            g.glyph_id = bigGid;
            g.text = utf8_of(glyph);
            op.glyphs.push_back(g);
            op.w = b.ink_w * sizePt;
            op.asc = b.asc * sizePt;
            op.desc = b.desc * sizePt;
            drawnGid = bigGid;
            drawnPt = sizePt;
        } else {
            /* At the TYPE size, not Equation Editor's symbol size.
             *
             * st_.sym is 18 point against a 12 point body -- its Sizes dialog
             * -- and it only ever showed in text style, because in display the
             * font's own variant is picked instead and this branch is not
             * taken.  In a denominator it was, and a summation there came out
             * 6.3 pt wide of TeX: 1.056 em at 18 point rather than at 12.
             * Appearance follows TeX, and TeX draws the text-style operator at
             * the type size. */
            op = glyph_layout(glyph, sizePt, false, true);
            op.w = std::max(op.w, glyph_ink_width(glyph, sizePt, false, true));
            drawnGid = baseGid;
            drawnPt = sizePt;
        }

        /* Centre it on the maths axis, which is TeX's make_op:
         *
         *     shift_amount(x) = half(height(x) - depth(x)) - axis_height
         *
         * A summation or an integral is drawn about its own middle, not about
         * the baseline, and the axis is the line a fraction bar and a minus
         * sign already sit on -- so centring there is what puts an operator in
         * line with everything beside it.  Left on the baseline it rides high,
         * and its limits inherit the error.
         *
         * It applies to the SIDE-limit case only.  TeX sets shift_amount on
         * the operator's box, and a shift_amount is vertical in an hlist and
         * HORIZONTAL in a vlist -- so when the limits go above and below and
         * the operator is packed into a vlist instead, the centring is simply
         * discarded and the vlist takes the glyph's own height and depth.
         * Applying it there too made a summation with limits 0.24 pt short at
         * the top, which is exactly this shift. */
        if (!stacked) {
            const double lift = (op.asc - op.desc) / 2.0 - mc.axisHeight * sizePt;
            if (std::fabs(lift) > 1e-9) {
                Layout centred;
                centred.absorb(op, 0, lift);
                centred.w = op.w;
                centred.asc = op.asc - lift;
                centred.desc = op.desc + lift;
                op = centred;
            }
        }

        /* Limits go above and below in DISPLAY style only.  In text style --
         * inside a fraction, say -- TeX sets even a summation's limits beside
         * it, because a stacked one would blow the line apart.  Taking the
         * template's own flag as the whole answer made a summation in a
         * denominator twice as deep as TeX sets it. */
        stacked = stacked && (fracDepth_ == 0);

        if (stacked) {
            /* Limits above and below.  One routine, shared with \lim and its
             * family, so a limit under a summation and a limit under \lim are
             * placed by the same arithmetic and cannot drift apart. */
            Layout st = stack_limits(op, lower, upper, hasLower, hasUpper,
                                     sizePt, lp, c, lowerSlot, upperSlot);
            out.absorb(st, 0, 0);
            out.asc = st.asc;
            out.desc = st.desc;
            x = st.w;
        } else {
            /* Inline: limits sit beside the operator like ordinary scripts. */
            out.absorb(op, 0, 0);
            out.asc = op.asc;
            out.desc = op.desc;
            x = op.w;
            /* An integral is tall and slanted, and both of those matter here.
             *
             * Tall: shifting the limits by a fixed fraction of the TYPE size
             * put them inside an operator set at the symbol size, so they came
             * out lying across the integral's own stroke.  TeX hangs a script
             * off the box it belongs to -- the operator's height less a drop
             * the font states -- and takes whichever is further.
             *
             * Slanted: the top of the integral leans right of its foot, so the
             * upper limit must follow it.  That is what a font's italic
             * correction is for, and the font gives a large one here. */
            /* The plain script shifts, NOT hung off the operator's height.
             *
             * TeX's rule 18a: when the nucleus is a single character -- which
             * a large operator glyph is -- the scripts are not pushed out to
             * clear it.  They sit at the ordinary superscript and subscript
             * heights and overlap its vertical extent, which is exactly how
             * an integral with limits looks when TeX sets one.
             *
             * Hanging them off op.asc and op.desc instead made the box a
             * sixth taller than TeX's.  The thing that had looked wrong
             * before was never the height: it was that both limits sat in a
             * column against a sign that leans, for want of the italic
             * correction below. */
            /* The drops are at the SCRIPT size -- Knuth's sup_drop(t) with t
             * already stepped down.  On an operator this is the term that
             * binds, so having it at the wrong size is worth 0.9 pt. */
            const double supShift = std::max(
                mc.superscriptShiftUp * sizePt,
                op.asc - mc.superscriptBaselineDropMax * ss);
            const double subShift = std::max(
                mc.subscriptShiftDown * sizePt,
                op.desc + mc.subscriptBaselineDropMin * ss);
            const double ic =
                drawnGid ? mf.italics_correction(drawnGid) * drawnPt : 0.0;

            /* TeX's make_op, for an operator whose limits go at its side:
             *
             *     if (subscript is not empty) and (not \limits) then
             *         width(x) <- width(x) - delta   { remove ic }
             *     ...
             *     shift_amount(superscript) <- delta
             *
             * The correction is taken OUT of the operator's own width and
             * spent shifting the superscript right instead -- it is the lean
             * of the integral, not extra advance.  Keeping it in as well put
             * the limits a whole correction too far out: an integral came out
             * half again as wide as TeX's, a contour integral nearly twice.
             *
             * The check is exact, in points:  \oint_C is
             * (11.988 - 7.092) + 7.3164 + 0.5 = 12.7124, which is what TeX
             * reports to the last digit. */
            const bool sideScripts = hasLower || hasUpper;
            if (sideScripts) x -= ic;

            double wsub = 0, wsup = 0;
            if (hasUpper) {
                Layout up = layout_list(upper, ss, slot_path(lp, c, upperSlot));
                out.absorb(up, x + ic, -supShift);
                out.asc = std::max(out.asc, supShift + up.asc);
                wsup = ic + up.w;
            }
            if (hasLower) {
                Layout lo = layout_list(lower, ss, slot_path(lp, c, lowerSlot));
                out.absorb(lo, x, subShift);
                out.desc = std::max(out.desc, subShift + lo.desc);
                wsub = lo.w;
            }
            x += std::max(wsub, wsup);
            if (sideScripts) x += kScriptSpace;
        }

        /* The operand joins the operator here rather than through the atom
         * loop, so the Op-to-whatever space has to be applied by hand -- a
         * sigma otherwise touches the symbol that follows it. */
        for (const auto& n : body) {
            if (!n || n->tag() == Node::kSize) continue;
            x += space_mu(kOp, class_of(*n)) * sizePt / 18.0;
            break;
        }

        Layout bodyL = layout_list(body, sizePt, slot_path(lp, c, bodySlot));
        out.absorb(bodyL, x, 0);
        out.asc = std::max(out.asc, bodyL.asc);
        out.desc = std::max(out.desc, bodyL.desc);
        out.w = x + bodyL.w;
        return out;
    }

    /* Unhandled templates still show their content rather than vanishing. */
    Layout layout_fallback(const Node& n, double sizePt,
                           const std::string& lp, int c) {
        /* Each of these has { content } as its first slot. */
        const std::string p = slot_path(lp, c, 0);
        switch (n.tag()) {
            case Node::kDecoration:
                return layout_list(static_cast<const DecorationNode&>(n).content,
                                   sizePt, p);
            case Node::kBraceDeco:
                return layout_list(static_cast<const BraceDecoNode&>(n).content,
                                   sizePt, p);
            case Node::kEmbell:
                return layout_list(static_cast<const EmbellNode&>(n).content,
                                   sizePt, p);
            default:
                return Layout();
        }
    }
};

}  // namespace

/* ------------------------------------------------------------------ */
/* Public API                                                          */
/* ------------------------------------------------------------------ */
Layout layout_math(const LineNode& root, const SvgStyle& style) {
    Renderer r(style);
    return r.run(root);
}

const CaretStop* find_stop(const Layout& layout, const std::string& path,
                           int index) {
    for (const CaretStop& s : layout.stops)
        if (s.index == index && s.path == path) return &s;
    return nullptr;
}

const CaretStop* nearest_stop(const Layout& layout, double x, double y) {
    const CaretStop* best = nullptr;
    double best_d = 0;
    int best_depth = -1;
    for (const CaretStop& s : layout.stops) {
        /* Vertical distance dominates: a click in a fraction's denominator
         * means the denominator, however close the numerator is horizontally. */
        double dy = 0;
        if (y < s.top) dy = s.top - y;
        else if (y > s.bottom) dy = y - s.bottom;
        double d = dy * 4.0 + std::fabs(x - s.x);

        /* Positions genuinely coincide -- the end of an integrand and the end
         * of the integral are the same point -- so the tie is broken towards
         * the innermost, where typing continues what is already there rather
         * than starting something after it. */
        int depth = 0;
        for (char ch : s.path) if (ch == '/') ++depth;
        if (!s.path.empty()) ++depth;

        const bool better = !best || d < best_d - 1e-9 ||
                            (std::fabs(d - best_d) <= 1e-9 && depth > best_depth);
        if (better) { best = &s; best_d = d; best_depth = depth; }
    }
    return best;
}

std::string render_svg(const LineNode& root, const SvgStyle& style) {
    Layout L = layout_math(root, style);

    const double pad = style.padding;
    const double width = L.w + 2 * pad;
    const double height = L.asc + L.desc + 2 * pad;
    const double baseline = pad + L.asc;

    std::ostringstream o;
    o.setf(std::ios::fixed);
    o.precision(3);
    o << "<svg xmlns=\"http://www.w3.org/2000/svg\" "
      << "width=\"" << width << "pt\" height=\"" << height << "pt\" "
      << "viewBox=\"0 0 " << width << ' ' << height << "\">\n";
    for (const auto& r2 : L.rules) {
        o << "  <rect x=\"" << (r2.x + pad) << "\" y=\"" << (r2.y + baseline)
          << "\" width=\"" << r2.w << "\" height=\"" << r2.h
          << "\" fill=\"currentColor\"/>\n";
    }
    for (const auto& g : L.glyphs) {
        const std::string& family = g.cjk    ? style.cjk
                                  : g.symbol ? style.symbol
                                             : style.serif;
        o << "  <text x=\"" << (g.x + pad) << "\" y=\"" << (g.y + baseline)
          << "\" font-family=\"" << xml_escape(family) << "\""
          << " font-size=\"" << g.size << "\"";
        if (g.italic) o << " font-style=\"italic\"";
        if (std::fabs(g.stretchY - 1.0) > 1e-6) {
            o << " transform=\"translate(" << (g.x + pad) << ',' << (g.y + baseline)
              << ") scale(1," << g.stretchY << ") translate("
              << -(g.x + pad) << ',' << -(g.y + baseline) << ")\"";
        }
        o << ">" << xml_escape(g.text) << "</text>\n";
    }
    o << "</svg>\n";
    return o.str();
}

std::string mtef_to_svg(const uint8_t* data, size_t len, const SvgStyle& style) {
    MtefParser::Result res = MtefParser::parse(data, len);
    if (!res.root) return std::string();
    return render_svg(*res.root, style);
}

std::string tex_to_svg(const std::string& latex, const SvgStyle& style) {
    std::unique_ptr<LineNode> root = parse_latex(latex);
    if (!root) return std::string();
    return render_svg(*root, style);
}


bool tex_box(const std::string& latex, const SvgStyle& style,
             double& w, double& asc, double& desc) {
    std::unique_ptr<LineNode> root = parse_latex(latex);
    if (!root) return false;
    Layout L = layout_math(*root, style);
    w = L.w; asc = L.asc; desc = L.desc;
    return true;
}

AtomKind atom_kind(uint32_t cp) {
    return AtomKind(int(class_of_char(cp)));
}

int atom_space_mu(AtomKind l, AtomKind r) {
    return space_mu(AtomClass(int(l)), AtomClass(int(r)));
}

}  // namespace mtef
