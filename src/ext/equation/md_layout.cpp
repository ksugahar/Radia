/*
 * md_layout.cpp -- laying a Markdown document out for viewing
 *
 * Text metrics come from GDI, measured once per font at a large reference em
 * and scaled.  That keeps the display list in points rather than pixels, so the
 * same layout serves the screen, an EMF, and a test that never opens a window.
 *
 * Line breaking handles Japanese, which has no spaces: a break is allowed
 * between two characters when either is CJK, minus the usual kinsoku rules
 * (nothing starts a line with a closing bracket or a full stop, nothing ends
 * one with an opening bracket).  Breaking only at spaces would put a whole
 * Japanese paragraph on one line, which is the failure this exists to avoid.
 */
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>

#include "md_layout.h"
#include "md_doc.h"
#include "tex_parser.h"

#include <algorithm>
#include <map>
#include <memory>

namespace mtef {
namespace {

/* Fonts are measured at this em and the widths scaled, so a caller may lay out
 * at any point size without a GDI object per size. */
const int kRefEm = 1000;

std::wstring to_utf16(const std::string& s) {
    if (s.empty()) return std::wstring();
    int n = MultiByteToWideChar(CP_UTF8, 0, s.c_str(), (int)s.size(), nullptr, 0);
    std::wstring w(n, L'\0');
    MultiByteToWideChar(CP_UTF8, 0, s.c_str(), (int)s.size(), &w[0], n);
    return w;
}

/* ---- code points -------------------------------------------------------- */

/* Decode one UTF-8 code point; `i` advances past it. */
uint32_t next_cp(const std::string& s, size_t& i) {
    unsigned char c = (unsigned char)s[i];
    if (c < 0x80) { ++i; return c; }
    int extra = (c >= 0xF0) ? 3 : (c >= 0xE0) ? 2 : 1;
    uint32_t cp = c & (0x3F >> extra);
    ++i;
    for (int k = 0; k < extra && i < s.size(); ++k, ++i)
        cp = (cp << 6) | ((unsigned char)s[i] & 0x3F);
    return cp;
}

bool is_cjk(uint32_t c) {
    return (c >= 0x3000 && c <= 0x30FF) ||    /* punctuation, kana        */
           (c >= 0x3400 && c <= 0x4DBF) ||    /* ideographs ext A         */
           (c >= 0x4E00 && c <= 0x9FFF) ||    /* ideographs               */
           (c >= 0xF900 && c <= 0xFAFF) ||    /* compatibility            */
           (c >= 0xFF00 && c <= 0xFF60) ||    /* fullwidth forms          */
           (c >= 0xFFE0 && c <= 0xFFE6);
}

bool in_set(uint32_t c, const uint32_t* set, size_t n) {
    for (size_t i = 0; i < n; ++i) if (set[i] == c) return true;
    return false;
}

/* Kinsoku: characters that may not begin a line.
 *
 * Written as code points rather than literals so the build does not depend on
 * the compiler guessing this file's encoding -- MSVC reads an unmarked source
 * in the ANSI code page, which on a Japanese machine mangles them. */
bool no_break_before(uint32_t c) {
    static const uint32_t set[] = {
        0x3001, 0x3002,                          /* ideographic comma, full stop */
        0xFF0C, 0xFF0E, 0xFF01, 0xFF1F,          /* fullwidth , . ! ?            */
        0xFF1A, 0xFF1B,                          /* fullwidth : ;                */
        0x300D, 0x300F, 0xFF09, 0xFF3D, 0xFF5D,  /* closing brackets             */
        0x3015, 0x3009, 0x300B,
        0x201D, 0x2019,                          /* closing quotes               */
        0x30FB, 0x3005, 0x309D, 0x309E,          /* middle dot, iteration marks  */
        0x30FC, 0x30FD, 0x30FE,                  /* prolonged sound, katakana it.*/
        0x3041, 0x3043, 0x3045, 0x3047, 0x3049,  /* small hiragana vowels        */
        0x3063, 0x3083, 0x3085, 0x3087, 0x308E,  /* small tsu, ya, yu, yo, wa    */
        0x30A1, 0x30A3, 0x30A5, 0x30A7, 0x30A9,  /* small katakana vowels        */
        0x30C3, 0x30E3, 0x30E5, 0x30E7, 0x30EE,
        ',', '.', '!', '?', ':', ';', ')', ']', '}',
    };
    return in_set(c, set, sizeof(set) / sizeof(set[0]));
}

/* Kinsoku: characters that may not end a line. */
bool no_break_after(uint32_t c) {
    static const uint32_t set[] = {
        0x300C, 0x300E, 0xFF08, 0xFF3B, 0xFF5B,  /* opening brackets */
        0x3014, 0x3008, 0x300A,
        0x201C, 0x2018,                          /* opening quotes   */
        '(', '[', '{',
    };
    return in_set(c, set, sizeof(set) / sizeof(set[0]));
}

/* ---- GDI metrics -------------------------------------------------------- */

struct FontFace {
    HFONT font = nullptr;
    double ascent = 0, descent = 0;   /* in reference-em units */
};

class Metrics {
public:
    Metrics() : dc_(CreateCompatibleDC(nullptr)) {}
    ~Metrics() {
        for (auto& kv : faces_) if (kv.second.font) DeleteObject(kv.second.font);
        if (dc_) DeleteDC(dc_);
    }
    Metrics(const Metrics&) = delete;
    Metrics& operator=(const Metrics&) = delete;

    /* Width of a string, in points, set at `size` points. */
    double width(const std::string& utf8, const std::string& face,
                 double size, bool bold, bool italic) {
        if (utf8.empty()) return 0.0;
        const FontFace& f = get(face, bold, italic);
        std::wstring w = to_utf16(utf8);
        HGDIOBJ old = SelectObject(dc_, f.font);
        SIZE sz{};
        GetTextExtentPoint32W(dc_, w.c_str(), (int)w.size(), &sz);
        SelectObject(dc_, old);
        return sz.cx * size / kRefEm;
    }

    void vmetrics(const std::string& face, double size, bool bold, bool italic,
                  double& asc, double& desc) {
        const FontFace& f = get(face, bold, italic);
        asc = f.ascent * size / kRefEm;
        desc = f.descent * size / kRefEm;
    }

private:
    const FontFace& get(const std::string& face, bool bold, bool italic) {
        std::string key = face + (bold ? "|b" : "|") + (italic ? "i" : "");
        auto it = faces_.find(key);
        if (it != faces_.end()) return it->second;

        LOGFONTW lf{};
        lf.lfHeight = -kRefEm;
        lf.lfWeight = bold ? FW_BOLD : FW_NORMAL;
        lf.lfItalic = italic ? TRUE : FALSE;
        /* DEFAULT_CHARSET, never SYMBOL_CHARSET: a symbol charset measures the
         * font through the legacy Symbol code page and loses most glyphs. */
        lf.lfCharSet = DEFAULT_CHARSET;
        lf.lfQuality = ANTIALIASED_QUALITY;
        std::wstring wf = to_utf16(face);
        wcsncpy_s(lf.lfFaceName, LF_FACESIZE, wf.c_str(), _TRUNCATE);

        FontFace f;
        f.font = CreateFontIndirectW(&lf);
        HGDIOBJ old = SelectObject(dc_, f.font);
        TEXTMETRICW tm{};
        GetTextMetricsW(dc_, &tm);
        f.ascent = tm.tmAscent;
        f.descent = tm.tmDescent;
        SelectObject(dc_, old);
        return faces_.emplace(key, f).first->second;
    }

    HDC dc_;
    std::map<std::string, FontFace> faces_;
};

/* ---- items on a line ---------------------------------------------------- */

struct Item {
    std::string text;         /* for a text run                              */
    bool math = false;
    Layout layout;            /* for an equation                             */
    std::string latex;
    bool display = false;
    double width = 0;
    double asc = 0, desc = 0;
    bool breakable = false;   /* a line may start here                       */
    bool space = false;       /* a space that vanishes at a line end         */

    /* Carried per item, not per block: an inline code span is set in the
     * monospace font inside an otherwise ordinary paragraph, and it must be
     * DRAWN in the font it was MEASURED in or the two disagree. */
    double size = 0;
    bool bold = false, italic = false, mono = false;
};

/* EQNEDT32's size ratios, scaled so the maths matches its surrounding text. */
SvgStyle math_style(double size, const DocStyle& st) {
    SvgStyle s;
    const double k = size * st.math_scale / 12.0;
    s.full = 12.0 * k;
    s.sub = 7.0 * k;
    s.sub2 = 5.0 * k;
    s.sym = 18.0 * k;
    s.subsym = 12.0 * k;
    s.padding = 0.0;          /* the document owns the spacing around it */
    return s;
}

/* Split a text run into pieces that must stay together, marking where a line
 * may begin. */
void push_text(std::vector<Item>& out, const std::string& text,
               Metrics& m, const std::string& face, double size,
               bool bold, bool italic, bool mono) {
    double asc = 0, desc = 0;
    m.vmetrics(face, size, bold, italic, asc, desc);

    size_t i = 0;
    std::string token;
    uint32_t prev = 0;

    auto make = [&](const std::string& s) {
        Item it;
        it.text = s;
        it.width = m.width(s, face, size, bold, italic);
        it.asc = asc;
        it.desc = desc;
        it.size = size;
        it.bold = bold;
        it.italic = italic;
        it.mono = mono;
        return it;
    };

    auto flush = [&](bool breakable) {
        if (token.empty()) return;
        Item it = make(token);
        it.breakable = breakable;
        out.push_back(it);
        token.clear();
    };

    bool next_breakable = out.empty();
    while (i < text.size()) {
        size_t start = i;
        uint32_t cp = next_cp(text, i);
        std::string ch = text.substr(start, i - start);

        if (cp == '\n' || cp == '\r') cp = ' ';   /* a wrapped paragraph */

        if (cp == ' ' || cp == '\t') {
            flush(next_breakable);
            next_breakable = true;
            /* A trailing space is kept on the line it ends but does not make
             * the line too long, so it gets its own zero-cost item. */
            Item sp = make(" ");
            sp.space = true;
            out.push_back(sp);
            prev = cp;
            continue;
        }

        bool can_break = prev != 0 && (is_cjk(prev) || is_cjk(cp)) &&
                         !no_break_before(cp) && !no_break_after(prev);
        if (can_break) {
            flush(next_breakable);
            next_breakable = true;
        }
        token += ch;
        prev = cp;
    }
    flush(next_breakable);
}

void push_math(std::vector<Item>& out, const std::string& latex, bool display,
               double size, const DocStyle& st) {
    std::unique_ptr<LineNode> root = parse_latex(latex);
    Item it;
    it.math = true;
    it.latex = latex;
    it.display = display;
    it.layout = layout_math(*root, math_style(size, st));
    it.width = it.layout.w;
    it.asc = it.layout.asc;
    it.desc = it.layout.desc;
    it.breakable = !out.empty();
    out.push_back(it);
}

}  // namespace

/* ------------------------------------------------------------------------- */

DocLayout layout_markdown(const std::string& markdown, double width,
                          const DocStyle& st) {
    DocLayout out;
    out.width = width;

    Metrics m;
    const double left = st.margin;
    const double right = std::max(left + 1.0, width - st.margin);
    const double gap = st.body * st.para_gap;

    double y = st.margin;
    int math_index = 0;

    std::vector<MdBlock> blocks = md_blocks(markdown);
    for (size_t bi = 0; bi < blocks.size(); ++bi) {
        const MdBlock& b = blocks[bi];
        /* Blank runs exist so the round trip is exact; the gap between blocks
         * is uniform, so they carry no layout of their own. */
        if (b.kind == MdBlock::kBlank) continue;

        const double top = y;
        double size = st.body;
        bool bold = false;
        std::string face = st.text_font;
        double indent = 0;

        switch (b.kind) {
            case MdBlock::kHeading: {
                int lv = std::min(std::max(b.level, 1), 6);
                size = st.heading[lv - 1];
                bold = true;
                break;
            }
            case MdBlock::kCode:
                face = st.mono_font;
                size = st.mono;
                break;
            case MdBlock::kBullet:
            case MdBlock::kNumbered:
                indent = st.body * st.list_indent * std::max(1, b.level);
                break;
            default:
                break;
        }

        if (b.kind == MdBlock::kCode) {
            /* Code is shown as written: one source line per line, no wrapping,
             * because a wrapped line of code says something the file does not. */
            double asc = 0, desc = 0;
            m.vmetrics(face, size, false, false, asc, desc);
            const double lh = (asc + desc) * st.line_spacing;
            double widest = 0;

            size_t pos = 0;
            while (pos <= b.text.size()) {
                size_t nl = b.text.find('\n', pos);
                std::string line = b.text.substr(
                    pos, nl == std::string::npos ? std::string::npos : nl - pos);
                DocRun r;
                r.text = line;
                r.x = left + st.body;
                r.baseline = y + asc;
                r.size = size;
                r.mono = true;
                widest = std::max(widest, m.width(line, face, size, false, false));
                out.runs.push_back(r);
                y += lh;
                if (nl == std::string::npos) break;
                pos = nl + 1;
            }
            Rule bg;
            bg.x = left;
            bg.y = top;
            bg.w = std::min(right - left, widest + 2 * st.body);
            bg.h = y - top;
            out.rules.push_back(bg);

            out.blocks.push_back({(int)bi, b.kind, top, y});
            y += gap;
            continue;
        }

        /* ---- a paragraph, heading, or list item -------------------------- */
        std::vector<Item> items;

        if (b.kind == MdBlock::kBullet || b.kind == MdBlock::kNumbered) {
            /* Take the marker from the source, so a numbered list shows the
             * number the file actually has rather than one we invented. */
            std::string marker;
            if (b.kind == MdBlock::kBullet) {
                marker = "\xE2\x80\xA2";           /* bullet */
            } else {
                size_t i = 0;
                while (i < b.source.size() &&
                       (b.source[i] == ' ' || b.source[i] == '\t')) ++i;
                size_t s = i;
                while (i < b.source.size() && b.source[i] >= '0' &&
                       b.source[i] <= '9') ++i;
                if (i < b.source.size() &&
                    (b.source[i] == '.' || b.source[i] == ')')) ++i;
                marker = b.source.substr(s, i - s);
            }
            DocRun r;
            r.text = marker;
            r.size = size;
            double asc = 0, desc = 0;
            m.vmetrics(face, size, bold, false, asc, desc);
            r.x = left + indent - st.body * 0.9;
            r.baseline = y + asc;   /* fixed up below once the first line lands */
            out.runs.push_back(r);
        }
        const size_t marker_run = out.runs.empty() ? 0 : out.runs.size() - 1;
        const bool has_marker =
            (b.kind == MdBlock::kBullet || b.kind == MdBlock::kNumbered);

        MarkdownDoc doc;
        doc.load(b.text);
        for (const MdSegment& seg : doc.segments()) {
            if (seg.is_math()) {
                push_math(items, seg.body, seg.kind == MdSegment::kDisplayMath,
                          size, st);
            } else if (seg.kind == MdSegment::kCodeSpan) {
                push_text(items, seg.body, m, st.mono_font, st.mono,
                          false, false, true);
            } else {
                push_text(items, seg.source(), m, face, size, bold, false, false);
            }
        }

        /* ---- greedy line breaking ---------------------------------------- */
        const double avail = right - (left + indent);
        size_t i = 0;
        bool first_line = true;
        while (i < items.size()) {
            size_t end = i;
            double w = 0, trailing = 0;
            while (end < items.size()) {
                const Item& it = items[end];
                if (end > i && it.breakable && w + it.width > avail && w > 0)
                    break;
                if (it.space) {
                    trailing += it.width;
                } else {
                    w += trailing + it.width;
                    trailing = 0;
                }
                ++end;
            }
            if (end == i) ++end;             /* one item is wider than the line */

            double asc = 0, desc = 0;
            for (size_t k = i; k < end; ++k) {
                asc = std::max(asc, items[k].asc);
                desc = std::max(desc, items[k].desc);
            }
            const double baseline = y + asc;

            /* A lone display equation is centred on its own line. */
            double x = left + indent;
            if (end == i + 1 && items[i].math && items[i].display)
                x = left + (right - left - items[i].width) / 2;

            for (size_t k = i; k < end; ++k) {
                Item& it = items[k];
                if (it.math) {
                    DocMath dm;
                    dm.layout = it.layout;
                    dm.x = x;
                    dm.baseline = baseline;
                    dm.latex = it.latex;
                    dm.display = it.display;
                    dm.block = (int)bi;
                    dm.index = math_index++;
                    out.maths.push_back(dm);
                } else {
                    DocRun r;
                    r.text = it.text;
                    r.x = x;
                    r.baseline = baseline;
                    r.size = it.size;
                    r.bold = it.bold;
                    r.italic = it.italic;
                    r.mono = it.mono;
                    out.runs.push_back(r);
                }
                x += it.width;
            }

            if (first_line && has_marker) {
                out.runs[marker_run].baseline = baseline;
                first_line = false;
            }
            y = baseline + desc * st.line_spacing;
            i = end;
        }
        if (items.empty()) y += (st.body * st.line_spacing);

        out.blocks.push_back({(int)bi, b.kind, top, y});
        y += gap;
    }

    out.height = y - gap + st.margin;
    if (out.height < st.margin * 2) out.height = st.margin * 2;
    return out;
}

int block_at(const DocLayout& doc, double x, double y) {
    (void)x;
    for (const DocBlockBox& b : doc.blocks)
        if (y >= b.top && y < b.bottom) return b.block;
    return -1;
}

int math_at(const DocLayout& doc, double x, double y) {
    for (const DocMath& m : doc.maths) {
        if (x >= m.x && x <= m.x + m.layout.w &&
            y >= m.baseline - m.layout.asc && y <= m.baseline + m.layout.desc)
            return m.index;
    }
    return -1;
}

}  // namespace mtef
