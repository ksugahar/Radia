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
std::string utf8_of(uint32_t cp) {
    std::string s;
    if (cp < 0x80) {
        s += char(cp);
    } else if (cp < 0x800) {
        s += char(0xC0 | (cp >> 6));
        s += char(0x80 | (cp & 0x3F));
    } else {
        s += char(0xE0 | (cp >> 12));
        s += char(0x80 | ((cp >> 6) & 0x3F));
        s += char(0x80 | (cp & 0x3F));
    }
    return s;
}

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

HFONT make_font(bool italic, bool symbol, bool cjk) {
    LOGFONTW lf = {};
    lf.lfHeight = -kEm;
    lf.lfItalic = italic ? TRUE : FALSE;
    /* DEFAULT_CHARSET, never SYMBOL_CHARSET.  Cambria Math is a Unicode font;
     * asking for the symbol charset makes GDI apply the legacy Symbol code
     * page, so U+0028 is measured as whatever sits at 0x28 in that page and
     * U+2264 is not found at all.  That single flag was behind the spurious
     * gap after "(" and the relation overlapping the fraction bar. */
    lf.lfCharSet = DEFAULT_CHARSET;
    /* Yu Mincho is the serif that sits beside Times without looking borrowed;
     * neither Times nor Cambria Math has a single kana. */
    const wchar_t* face = cjk      ? L"Yu Mincho"
                        : symbol   ? L"Cambria Math"
                                   : L"Times New Roman";
    wcscpy_s(lf.lfFaceName, face);
    return CreateFontIndirectW(&lf);
}

struct MetricCache {
    HDC hdc = nullptr;
    std::map<int, HFONT> fonts;              /* key: italic*2 + symbol */
    std::map<std::pair<int, uint32_t>, double> widths;
    std::map<int, std::pair<double, double>> vmetrics;  /* asc, desc per em */

    MetricCache() { hdc = CreateCompatibleDC(nullptr); }
    ~MetricCache() {
        for (auto& kv : fonts) DeleteObject(kv.second);
        if (hdc) DeleteDC(hdc);
    }
    HFONT font(int key) {
        auto it = fonts.find(key);
        if (it != fonts.end()) return it->second;
        HFONT f = make_font((key & 1) != 0, (key & 2) != 0, (key & 4) != 0);
        fonts[key] = f;
        return f;
    }
    /* The face follows the character, so a caller cannot measure with one font
     * and have the renderer draw with another. */
    static int face_key(uint32_t cp, bool italic, bool symbol) {
        return (italic ? 1 : 0) | (symbol ? 2 : 0) | (is_cjk(cp) ? 4 : 0);
    }
    double width_em(uint32_t cp, bool italic, bool symbol) {
        int key = face_key(cp, italic, symbol);
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
    std::pair<double, double> vmetric(bool italic, bool symbol) {
        int key = (italic ? 1 : 0) | (symbol ? 2 : 0);
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
     * ascent and descent are the wrong measure: Cambria Math reserves room for
     * extensible brackets and integral signs, so using its global descent puts
     * the subscript of a sigma a third of a line too low, while the subscript
     * of a Times "B" sits correctly.  TeX has always used per-glyph height and
     * depth for exactly this reason. */
    struct Box { double asc = 0, desc = 0, ink_w = 0; };

    Box glyph_box(uint32_t cp, bool italic, bool symbol) {
        int key = face_key(cp, italic, symbol);
        auto k = std::make_pair(key, cp);
        auto it = boxes.find(k);
        if (it != boxes.end()) return it->second;

        auto v = vmetric(italic, symbol);
        Box b;
        b.asc = v.first;
        b.desc = v.second;
        b.ink_w = width_em(cp, italic, symbol);
        if (cp < 0x10000) {
            HGDIOBJ old = SelectObject(hdc, font(key));
            GLYPHMETRICS gm = {};
            MAT2 id = {{0, 1}, {0, 0}, {0, 0}, {0, 1}};
            DWORD r = GetGlyphOutlineW(hdc, cp, GGO_METRICS, &gm, 0, nullptr, &id);
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
        if (GetCharWidthI(hdc, gid, 1, nullptr, &adv)) b.ink_w = adv / kEm;
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
    return metrics().width_em(cp, italic, symbol) * sizePt;
}
/* The font's own extent, for things sized against the face rather than
 * against one glyph (fence stretching, the fallback line height). */
void char_vmetrics(double sizePt, bool italic, bool symbol,
                   double& asc, double& desc) {
    auto v = metrics().vmetric(italic, symbol);
    asc = v.first * sizePt;
    desc = v.second * sizePt;
}
/* The ink box of one glyph, which is what neighbours and scripts must clear. */
void glyph_vmetrics(uint32_t cp, double sizePt, bool italic, bool symbol,
                    double& asc, double& desc) {
    auto b = metrics().glyph_box(cp, italic, symbol);
    asc = b.asc * sizePt;
    desc = b.desc * sizePt;
}
/* How far the drawing actually reaches.  A large operator is drawn wider than
 * it advances, so laying the next atom out at the advance alone lets a sigma
 * touch the symbol after it. */
double glyph_ink_width(uint32_t cp, double sizePt, bool italic, bool symbol) {
    return metrics().glyph_box(cp, italic, symbol).ink_w * sizePt;
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

/* ------------------------------------------------------------------ */
/* TeX atom classes and the spacing between them                       */
/* ------------------------------------------------------------------ */
enum AtomClass { kOrd, kOp, kBin, kRel, kOpen, kClose, kPunct, kInner };

AtomClass class_of_char(uint32_t cp) {
    switch (cp) {
        case '+': case '-': case 0x2212: case 0x00B1: case 0x2213:
        case 0x00D7: case 0x00F7: case 0x22C5: case 0x2217: case 0x2218:
        case 0x2229: case 0x222A: case 0x2227: case 0x2228: case 0x2295:
        case 0x2297: case 0x2299: case 0x2296: case 0x228E: case 0x2216:
            return kBin;
        case '=': case '<': case '>': case 0x2260: case 0x2264: case 0x2265:
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

    double size_of(int sizeType) const {
        switch (sizeType) {
            case SIZETYPE_SUB:    return st_.sub;
            case SIZETYPE_SUB2:   return st_.sub2;
            case SIZETYPE_SYM:    return st_.sym;
            case SIZETYPE_SUBSYM: return st_.subsym;
            default:              return st_.full;
        }
    }
    /* One step down for scripts, floored at the sub-subscript size. */
    double script_size(double cur) const {
        if (cur > st_.sub + 1e-9) return st_.sub;
        return st_.sub2;
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
        Layout L;
        L.w = char_width(cp, sizePt, italic, symbol);
        glyph_vmetrics(cp, sizePt, italic, symbol, L.asc, L.desc);
        Glyph g;
        g.x = 0; g.y = 0; g.size = sizePt;
        g.italic = italic; g.symbol = symbol; g.cjk = is_cjk(cp);
        g.text = utf8_of(cp);
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

        /* A bar is a rule the width of what it covers; a fixed-width macron
         * over a wide expression would look wrong. */
        const double gap = sizePt * 0.08;
        if (!acc) {
            Rule r;
            r.x = 0;
            r.y = -(body.asc + gap + sizePt * 0.05);
            r.w = body.w;
            r.h = std::max(0.5, sizePt * 0.05);
            out.rules.push_back(r);
            out.asc = body.asc + gap + r.h + sizePt * 0.05;
            return out;
        }

        /* An accent is set smaller than what it sits over.  At full size an
         * arrow is as heavy as the letter and reads as a second symbol rather
         * than as a mark on the first. */
        const double markPt = sizePt * 0.62;
        Layout mark = glyph_layout(acc, markPt, false, needs_math_face(acc));
        const double dx = (body.w - mark.w) / 2;
        const double dy = -(body.asc + gap + mark.desc);
        out.absorb(mark, dx, dy);
        out.asc = std::max(out.asc, body.asc + gap + mark.asc + mark.desc);
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

        const double gap = sizePt * 0.10;
        Rule r;
        r.x = 0;
        r.w = body.w;
        r.h = std::max(0.5, sizePt * 0.05);
        if (above) {
            r.y = -(body.asc + gap + r.h);
            out.asc = body.asc + gap + r.h;
        } else {
            r.y = body.desc + gap;
            out.desc = body.desc + gap + r.h;
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
            case Node::kFence:
                return kInner;
            case Node::kFrac:
                return kInner;
            case Node::kLine: {
                /* A group takes the class of its first atom, so \sin(x)
                 * spaces like a function and not like a bare group. */
                const auto& l = static_cast<const LineNode&>(n);
                for (const auto& c : l.children)
                    if (c && c->tag() != Node::kSize) return class_of(*c);
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

            if (have_prev)
                x += space_mu(prev, cls) * cur / 18.0;

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
            out.w = 0.55 * sizePt;
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
                                    needs_math_face(cp));
            }
            case Node::kEmbell:
                return layout_embell(static_cast<const EmbellNode&>(n), sizePt,
                                     listPath, child);
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
    Layout layout_script(const ScriptNode& s, double sizePt,
                         const std::string& lp, int c) {
        const int subSlot = 1;
        const int supSlot = s.hasSub ? 2 : 1;
        Layout base = layout_list(s.base, sizePt, slot_path(lp, c, 0));
        double ss = script_size(sizePt);
        Layout out = base;
        double x = base.w;
        /* Scripts hang off the base's own extents, not off a fixed offset:
         * the exponent of a tall base has to clear that base, and the
         * subscript of a deep one has to sit below it. */
        const double supShift = std::max(0.45 * sizePt, base.asc - 0.35 * ss);
        const double subShift = std::max(0.22 * sizePt, base.desc + 0.12 * ss);
        double wsub = 0, wsup = 0;
        if (s.hasSup) {
            Layout sup = layout_list(s.sup, ss, slot_path(lp, c, supSlot));
            out.absorb(sup, x, -supShift);
            out.asc = std::max(out.asc, supShift + sup.asc);
            wsup = sup.w;
        }
        if (s.hasSub) {
            Layout sub = layout_list(s.sub, ss, slot_path(lp, c, subSlot));
            out.absorb(sub, x, subShift);
            out.desc = std::max(out.desc, subShift + sub.desc);
            wsub = sub.w;
        }
        out.w = x + std::max(wsub, wsup);
        return out;
    }

    /* node_slots(kFence) = { content } */
    Layout layout_fence(const FenceNode& f, double sizePt,
                        const std::string& lp, int c) {
        Layout inner = layout_list(f.content, sizePt, slot_path(lp, c, 0));
        auto gl = fence_glyphs(f.selector);
        /* Grow the fence to the content, keeping a plain glyph when it fits. */
        double need = inner.asc + inner.desc;
        double plainAsc, plainDesc;
        char_vmetrics(sizePt, false, false, plainAsc, plainDesc);
        double stretch = std::max(1.0, need / std::max(plainAsc + plainDesc, 1e-6));

        Layout out;
        double x = 0;
        bool left = (f.variation == 0 || f.variation == 1);
        bool right = (f.variation == 0 || f.variation == 2);
        if (left) {
            Layout g = glyph_layout(gl.first, sizePt, false, needs_math_face(gl.first));
            for (auto& gg : g.glyphs) gg.stretchY = stretch;
            g.asc *= stretch; g.desc *= stretch;
            out.absorb(g, x, 0);
            out.asc = std::max(out.asc, g.asc);
            out.desc = std::max(out.desc, g.desc);
            x += g.w;
        }
        out.absorb(inner, x, 0);
        out.asc = std::max(out.asc, inner.asc);
        out.desc = std::max(out.desc, inner.desc);
        x += inner.w;
        if (right) {
            Layout g = glyph_layout(gl.second, sizePt, false, needs_math_face(gl.second));
            for (auto& gg : g.glyphs) gg.stretchY = stretch;
            g.asc *= stretch; g.desc *= stretch;
            out.absorb(g, x, 0);
            out.asc = std::max(out.asc, g.asc);
            out.desc = std::max(out.desc, g.desc);
            x += g.w;
        }
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

        const mtef::MathFont& mf = mtef::MathFont::cambria();
        const mtef::MathConstants& mc = mf.constants();

        /* The axis is where a fraction bar sits and where a minus sign is
         * centred, so everything else hangs off it.  The guess this replaces
         * was 0.28 of the type size against the font's 0.2856 -- close enough
         * that nothing moves, which is the point: it is now a fact. */
        const double axis  = mc.axisHeight * sizePt;
        const double thick = mc.fractionRuleThickness * sizePt;

        /* Equation Editor 3.1 extends the bar past the wider of the two parts
         * by ONE POINT on each side (Format > Spacing, "Fraction bar
         * overhang"); the OpenType MATH table has no such notion and MathML
         * Core draws the bar exactly as wide as the parts.  The overhang is
         * part of how a fraction reads, so keep Equation Editor's -- but note
         * that this replaces 0.4 of the type size, which at 12 pt was two and
         * a half times wider on each side than the editor being imitated. */
        const double overhang = 1.0;

        /* Then the two parts are pushed apart until they clear the bar by at
         * least the font's minimum gap.  Equation Editor states this as
         * "Numerator height 35%" and "Denominator depth 100%" against its own
         * reference; the font states it as a target shift plus a floor, which
         * is the same idea with the floor made explicit. */
        double shiftUp   = mc.fractionNumeratorDisplayStyleShiftUp * sizePt;
        double shiftDown = mc.fractionDenominatorDisplayStyleShiftDown * sizePt;
        const double gapNum = mc.fractionNumDisplayStyleGapMin * sizePt;
        const double gapDen = mc.fractionDenomDisplayStyleGapMin * sizePt;

        shiftUp   = std::max(shiftUp,   axis + thick / 2.0 + gapNum + num.desc);
        shiftDown = std::max(shiftDown, -axis + thick / 2.0 + gapDen + den.asc);

        const double w = std::max(num.w, den.w) + 2.0 * overhang;

        Layout out;
        out.w = w;
        out.absorb(num, (w - num.w) / 2.0, -shiftUp);
        out.absorb(den, (w - den.w) / 2.0, shiftDown);

        Rule bar;
        bar.x = 0;
        bar.y = -axis - thick / 2.0;
        bar.w = w;
        bar.h = thick;
        out.rules.push_back(bar);

        out.asc  = shiftUp + num.asc;
        out.desc = shiftDown + den.desc;
        return out;
    }

    /* node_slots(kSqrt) = { content } + (hasIndex ? index) */
    Layout layout_sqrt(const SqrtNode& s, double sizePt,
                       const std::string& lp, int c) {
        Layout inner = layout_list(s.content, sizePt, slot_path(lp, c, 0));

        /* The font states all four of these; they used to be 0.04 and 0.18 of
         * the type size, picked by eye. */
        const mtef::MathFont& mf = mtef::MathFont::cambria();
        const mtef::MathConstants& mc = mf.constants();
        const double gap   = mc.radicalVerticalGap  * sizePt;
        const double thick = mc.radicalRuleThickness * sizePt;
        const double extra = mc.radicalExtraAscender * sizePt;

        /* A radical is not a character scaled up.  The font draws one at six
         * heights, so ask for the one that reaches -- the stroke stays the
         * weight the designer drew and the hook still meets the bar. */
        const double need = inner.asc + inner.desc + gap + thick;
        Layout sign;
        double got_em = 0;
        const uint16_t variant =
            mf.ok() ? mf.vertical_variant(mf.glyph_for(0x221A), need / sizePt,
                                          &got_em)
                    : 0;
        if (variant) {
            const MetricCache::Box b = metrics().glyph_box_index(variant);
            Glyph g;
            g.size = sizePt;
            g.symbol = true;
            g.glyph_id = variant;
            g.text = utf8_of(0x221A);
            sign.glyphs.push_back(g);
            sign.w = b.ink_w * sizePt;
            sign.asc = b.asc * sizePt;
            sign.desc = b.desc * sizePt;
        } else {
            sign = glyph_layout(0x221A, sizePt, false, true);
        }

        /* Past the tallest drawing the font offers, the remainder is taken by
         * scaling.  The font also ships parts to assemble one of any height;
         * until that is wired, this at least covers the radicand rather than
         * stopping short of it. */
        if (variant && got_em * sizePt < need) {
            const double k = need / std::max(got_em * sizePt, 1e-6);
            for (auto& g : sign.glyphs) g.stretchY = k;
            sign.asc *= k;
            sign.desc *= k;
        }

        /* The index sits above the radical's left arm, at script size.  It is
         * given its own width rather than being tucked under the sign, which
         * keeps a two-digit index from colliding with the radicand. */
        Layout idx;
        double idxW = 0;
        if (s.hasIndex) {
            idx = layout_list(s.index, script_size(sizePt), slot_path(lp, c, 1));
            idxW = idx.w;
        }

        Layout out;
        const double signX = idxW;

        /* Everything hangs off the bar: it sits one gap above the radicand,
         * and the sign is placed so its own top lands on it. */
        const double barTop = -(inner.asc + gap + thick);
        const double signBaseline = barTop + sign.asc;

        out.absorb(sign, signX, signBaseline);
        out.absorb(inner, signX + sign.w, 0);
        out.asc = std::max(-barTop + extra, inner.asc);
        out.desc = std::max(sign.desc - signBaseline, inner.desc);

        Rule bar;
        bar.x = signX + sign.w;
        bar.y = barTop;
        bar.w = inner.w;
        bar.h = thick;
        out.rules.push_back(bar);

        if (s.hasIndex) {
            /* The font says how far up the degree sits, as a fraction of the
             * radical's own height. */
            const double lift = mc.radicalDegreeBottomRaisePercent * (-barTop);
            out.absorb(idx, 0, -lift);
            out.asc = std::max(out.asc, lift + idx.asc);
        }
        out.w = signX + sign.w + inner.w;
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
        Layout op = glyph_layout(glyph, opSize, false, true);
        op.w = std::max(op.w, glyph_ink_width(glyph, opSize, false, true));
        double ss = script_size(sizePt);
        Layout out;
        double x = 0;

        if (stacked) {
            /* Limits above and below the operator, everything centred on the
             * widest of the three.  Centring the limits on the operator alone
             * lets a wide limit hang left of the origin and overlap whatever
             * precedes the operator, because the reported width never sees it. */
            Layout up, lo;
            if (hasUpper) up = layout_list(upper, ss, slot_path(lp, c, upperSlot));
            if (hasLower) lo = layout_list(lower, ss, slot_path(lp, c, lowerSlot));
            double w = op.w;
            if (hasUpper) w = std::max(w, up.w);
            if (hasLower) w = std::max(w, lo.w);

            out.absorb(op, (w - op.w) / 2.0, 0);
            out.asc = op.asc;
            out.desc = op.desc;
            if (hasUpper) {
                out.absorb(up, (w - up.w) / 2.0, -(op.asc + 0.15 * sizePt + up.desc));
                out.asc = std::max(out.asc, op.asc + 0.15 * sizePt + up.desc + up.asc);
            }
            if (hasLower) {
                out.absorb(lo, (w - lo.w) / 2.0, op.desc + 0.15 * sizePt + lo.asc);
                out.desc = std::max(out.desc, op.desc + 0.15 * sizePt + lo.asc + lo.desc);
            }
            x = w;
        } else {
            /* Inline: limits sit beside the operator like ordinary scripts. */
            out.absorb(op, 0, 0);
            out.asc = op.asc;
            out.desc = op.desc;
            x = op.w;
            double wsub = 0, wsup = 0;
            if (hasUpper) {
                Layout up = layout_list(upper, ss, slot_path(lp, c, upperSlot));
                out.absorb(up, x, -0.45 * sizePt);
                out.asc = std::max(out.asc, 0.45 * sizePt + up.asc);
                wsup = up.w;
            }
            if (hasLower) {
                Layout lo = layout_list(lower, ss, slot_path(lp, c, lowerSlot));
                out.absorb(lo, x, 0.22 * sizePt);
                out.desc = std::max(out.desc, 0.22 * sizePt + lo.desc);
                wsub = lo.w;
            }
            x += std::max(wsub, wsup);
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

}  // namespace mtef
