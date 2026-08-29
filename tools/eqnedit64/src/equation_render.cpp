/*
 * equation_render.cpp -- TeX equation tree -> native canvas / SVG
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
#include "equation_render.h"
#include "tex_parser.h"

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <iomanip>
#include <sstream>
#include <vector>

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX          /* else windows.h's max/min macros eat std::max */
#include <windows.h>
#endif

namespace eqnedit {
namespace {

/* ------------------------------------------------------------------ */
/* Display list                                                        */
/* ------------------------------------------------------------------ */
struct Glyph {
    double x = 0, y = 0;        /* y is the baseline */
    double size = 0;            /* pt */
    bool italic = false;
    bool bold = false;
    bool symbol = false;        /* pick the symbol font family */
    bool cjk = false;           /* use a real CJK face, never math-font linking */
    double stretchY = 1.0;      /* vertical scale for grown fences */
    /* Draw only the part of the glyph left of this x, in layout units.
     * Used for the radical: Cambria's surd ends in a flat flag thicker than
     * TeX's rule, so drawing it whole left a stub at the left of the bar.
     * Zero means draw the whole glyph. */
    double clipRight = 0.0;
    /* A designed size variant is reachable only by glyph index -- it has no
     * character.  When this is set, `text` is only a label for tests and the
     * drawing goes through the index; `outline` carries the same shape as an
     * SVG path, in em units, because an SVG cannot name the glyph either. */
    unsigned short glyphIndex = 0;
    std::string outline;
    std::string text;           /* UTF-8 */
};

struct Rule { double x = 0, y = 0, w = 0, h = 0; };   /* y is the TOP edge */

struct Placeholder { double x = 0, y = 0, w = 0, h = 0; };

enum class EditMarkKind { Space, Alignment };
struct EditMark {
    EditMarkKind kind = EditMarkKind::Space;
    double x = 0;
    double top = 0;
    double bottom = 0;
};

struct CaretSite {
    const NodeList* slot = nullptr;
    int index = 0;
    double x = 0;
    double top = 0;
    double bottom = 0;
};

struct Layout {
    double w = 0, asc = 0, desc = 0;
    std::vector<Glyph> glyphs;
    std::vector<Rule> rules;
    std::vector<Placeholder> placeholders;
    std::vector<EditMark> editMarks;
    std::vector<CaretSite> carets;

    void translate(double dx, double dy) {
        /* clipRight is an x in the same space as glyph.x, so it has to move
         * with it.  Leaving it behind clipped the radical sign away entirely
         * as soon as it sat anywhere but the origin -- inside `aligned`, or
         * after anything else on the line. */
        for (auto& g : glyphs) {
            g.x += dx; g.y += dy;
            if (g.clipRight != 0.0) g.clipRight += dx;
        }
        for (auto& r : rules)  { r.x += dx; r.y += dy; }
        for (auto& p : placeholders) { p.x += dx; p.y += dy; }
        for (auto& m : editMarks) {
            m.x += dx; m.top += dy; m.bottom += dy;
        }
        for (auto& c : carets) { c.x += dx; c.top += dy; c.bottom += dy; }
    }
    void absorb(const Layout& other, double dx, double dy) {
        Layout t = other;
        t.translate(dx, dy);
        glyphs.insert(glyphs.end(), t.glyphs.begin(), t.glyphs.end());
        rules.insert(rules.end(), t.rules.begin(), t.rules.end());
        placeholders.insert(placeholders.end(), t.placeholders.begin(),
                            t.placeholders.end());
        editMarks.insert(editMarks.end(), t.editMarks.begin(),
                         t.editMarks.end());
        carets.insert(carets.end(), t.carets.begin(), t.carets.end());
    }
};

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
    } else if (cp < 0x10000) {
        s += char(0xE0 | (cp >> 12));
        s += char(0x80 | ((cp >> 6) & 0x3F));
        s += char(0x80 | (cp & 0x3F));
    } else {
        s += char(0xF0 | (cp >> 18));
        s += char(0x80 | ((cp >> 12) & 0x3F));
        s += char(0x80 | ((cp >> 6) & 0x3F));
        s += char(0x80 | (cp & 0x3F));
    }
    return s;
}

std::vector<uint32_t> utf8_codes(const std::string& s) {
    std::vector<uint32_t> out;
    for (size_t i = 0; i < s.size();) {
        unsigned char c = static_cast<unsigned char>(s[i]);
        uint32_t cp = 0;
        size_t n = 1;
        if (c < 0x80) cp = c;
        else if ((c & 0xE0) == 0xC0 && i + 1 < s.size()) {
            cp = c & 0x1F; n = 2;
        } else if ((c & 0xF0) == 0xE0 && i + 2 < s.size()) {
            cp = c & 0x0F; n = 3;
        } else if ((c & 0xF8) == 0xF0 && i + 3 < s.size()) {
            cp = c & 0x07; n = 4;
        } else {
            cp = 0xFFFD;
        }
        if (n > 1) {
            bool valid = true;
            for (size_t j = 1; j < n; ++j) {
                unsigned char d = static_cast<unsigned char>(s[i + j]);
                if ((d & 0xC0) != 0x80) { valid = false; break; }
                cp = (cp << 6) | (d & 0x3F);
            }
            if (!valid) { cp = 0xFFFD; n = 1; }
        }
        out.push_back(cp);
        i += n;
    }
    return out;
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

/* An address inside this module, so the loader below can find which module
 * it is running in without hard-coding a name. */
void make_font_marker() {}

/* Make Latin Modern Math available to this process, once.
 *
 * The executable carries the .otf as RCDATA and writes only a verified,
 * content-addressed per-user cache: no font install or registry entry is
 * created, and the portable single-file input rule still holds.  The
 * canonical Python test module embeds the same RCDATA resource.  A target
 * that omits it fails the explicit math_font_loaded() health check. */
bool g_mathFontLoaded = false;
bool g_mathFontRegistered = false;


/* Is the face usable, or did it merely register?  Measured, not asked:
 * AddFontMemResourceEx can report success while GDI still has nothing to draw
 * with, and every glyph then measures zero wide.  That is what produced a
 * parenthesis 0.00 pt across and a radical with no sign -- three unrelated
 * checks failing at once, in a fresh process, only under process churn. */
bool math_face_measures() {
    LOGFONTW lf = {};
    lf.lfHeight = -kEm;
    lf.lfCharSet = DEFAULT_CHARSET;
    wcscpy_s(lf.lfFaceName, L"Latin Modern Math");
    HFONT probe = CreateFontIndirectW(&lf);
    if (!probe) return false;
    HDC dc = CreateCompatibleDC(nullptr);
    bool ok = false;
    if (dc) {
        HGDIOBJ old = SelectObject(dc, probe);
        const wchar_t sample = L'(';
        SIZE size = {};
        ok = GetTextExtentPoint32W(dc, &sample, 1, &size) && size.cx > 0;
        SelectObject(dc, old);
        DeleteDC(dc);
    }
    DeleteObject(probe);
    return ok;
}

std::filesystem::path font_cache_root() {
    const DWORD needed = GetEnvironmentVariableW(L"LOCALAPPDATA", nullptr, 0);
    if (needed > 1) {
        std::vector<wchar_t> value(needed);
        if (GetEnvironmentVariableW(L"LOCALAPPDATA", value.data(), needed))
            return std::filesystem::path(value.data()) / L"Eqnedit64" /
                   L"fonts";
    }
    wchar_t temporary[MAX_PATH] = {};
    if (GetTempPathW(_countof(temporary), temporary))
        return std::filesystem::path(temporary) / L"Eqnedit64" / L"fonts";
    return {};
}

uint64_t font_bytes_hash(const unsigned char* bytes, size_t size) {
    uint64_t hash = UINT64_C(14695981039346656037);
    for (size_t i = 0; i < size; ++i) {
        hash ^= bytes[i];
        hash *= UINT64_C(1099511628211);
    }
    return hash;
}

bool file_matches_bytes(const std::filesystem::path& path,
                        const unsigned char* bytes, size_t size) {
    std::error_code error;
    if (!std::filesystem::is_regular_file(path, error) || error ||
        std::filesystem::file_size(path, error) != size || error)
        return false;
    std::ifstream file(path, std::ios::binary);
    if (!file) return false;
    std::vector<unsigned char> buffer(64 * 1024);
    size_t offset = 0;
    while (offset < size) {
        const size_t count = std::min(buffer.size(), size - offset);
        file.read(reinterpret_cast<char*>(buffer.data()),
                  std::streamsize(count));
        if (size_t(file.gcount()) != count ||
            memcmp(buffer.data(), bytes + offset, count) != 0)
            return false;
        offset += count;
    }
    return true;
}

std::filesystem::path cache_embedded_math_font(const unsigned char* bytes,
                                               size_t size) {
    const std::filesystem::path directory = font_cache_root();
    if (directory.empty()) return {};
    std::error_code error;
    std::filesystem::create_directories(directory, error);
    if (error) return {};

    std::wostringstream filename;
    filename << L"latinmodern-math-" << std::hex << std::setw(16)
             << std::setfill(L'0') << font_bytes_hash(bytes, size) << L".otf";
    const std::filesystem::path target = directory / filename.str();
    if (file_matches_bytes(target, bytes, size)) return target;

    std::wostringstream temporaryName;
    temporaryName << filename.str() << L"." << GetCurrentProcessId() << L"."
                  << GetTickCount64() << L".new";
    const std::filesystem::path temporary = directory / temporaryName.str();
    {
        std::ofstream file(temporary, std::ios::binary | std::ios::trunc);
        if (!file) return {};
        file.write(reinterpret_cast<const char*>(bytes),
                   std::streamsize(size));
        file.flush();
        if (!file) {
            file.close();
            std::filesystem::remove(temporary, error);
            return {};
        }
    }
    if (!file_matches_bytes(temporary, bytes, size)) {
        std::filesystem::remove(temporary, error);
        return {};
    }
    if (!MoveFileExW(temporary.c_str(), target.c_str(),
                     MOVEFILE_REPLACE_EXISTING | MOVEFILE_WRITE_THROUGH)) {
        std::filesystem::remove(temporary, error);
        if (!file_matches_bytes(target, bytes, size)) return {};
    }
    return file_matches_bytes(target, bytes, size) ? target
                                                    : std::filesystem::path{};
}

bool load_math_font() {
    HMODULE self = nullptr;
    GetModuleHandleExW(GET_MODULE_HANDLE_EX_FLAG_FROM_ADDRESS |
                           GET_MODULE_HANDLE_EX_FLAG_UNCHANGED_REFCOUNT,
                       reinterpret_cast<LPCWSTR>(&make_font_marker), &self);
    if (!self) return false;
    /* RT_RCDATA is MAKEINTRESOURCE(10), which is the narrow macro unless
     * UNICODE is defined -- and the test module is not built with it. */
    HRSRC found = FindResourceW(self, MAKEINTRESOURCEW(200),
                               MAKEINTRESOURCEW(10));
    if (!found) return false;
    const DWORD size = SizeofResource(self, found);
    HGLOBAL block = LoadResource(self, found);
    const auto* bytes = block
        ? static_cast<const unsigned char*>(LockResource(block)) : nullptr;
    if (!bytes || !size) return false;

    /* AddFontMemResourceEx repeatedly crashed Server 2022's per-session
     * fontdrvhost.exe (0xc0000005) while every old hidden test still passed.
     * Extract the verified embedded bytes once and use the file-backed,
     * process-private API instead.  This is a cache, not an installation:
     * there is no registry entry and the EXE remains the only input. */
    const std::filesystem::path path = cache_embedded_math_font(bytes, size);
    if (path.empty()) return false;
    return AddFontResourceExW(path.c_str(), FR_PRIVATE | FR_NOT_ENUM,
                              nullptr) > 0;
}

/* Keep trying until the face actually measures.  Registering it is not the
 * same as being able to draw with it, and the gap between the two is real on
 * this machine when processes start while others are exiting. */
void ensure_math_font() {
    if (g_mathFontLoaded) return;
    /* Experiment kill switch (2026-08-24).  With EQNEDIT64_NO_FONT_REG set,
     * never register anything: the canvas then draws with whatever GDI
     * substitutes, which is unusable for real work but exactly what the
     * session-font-poisoning experiment needs -- hundreds of process starts
     * that are identical except for the registrations.  If the fonts die
     * under the normal build and survive under this switch, the
     * registrations are the vector; if they die under both, the churn is. */
    if (GetEnvironmentVariableW(L"EQNEDIT64_NO_FONT_REG", nullptr, 0)) return;
    /* Register ONCE, then retry only the measurement.  The earlier loop
     * re-registered the font on every attempt -- up to twenty
     * AddFontMemResourceEx calls per process, times the hundreds of test
     * processes a working day starts.  Session font tables do not enjoy
     * that: this machine''s GDI reached a state where healthy faces drew
     * blanks or the wrong glyphs, and our own churn is the leading suspect.
     * Registration is per-process and freed at exit; there is never a
     * reason to do it more than once. */
    /* A font host finishing another short-lived test process can reject
     * registration transiently. Retry only failed API calls. Once a handle
     * succeeds, remember it permanently for this process so rendering can
     * never register the same font twice. If every attempt fails, a later
     * health check may try again because no font was added. */
    for (int attempt = 0; attempt < 40 && !g_mathFontRegistered; ++attempt) {
        g_mathFontRegistered = load_math_font();
        if (!g_mathFontRegistered) Sleep(50);
    }
    if (!g_mathFontRegistered) return;
    for (int attempt = 0; attempt < 30; ++attempt) {
        g_mathFontLoaded = math_face_measures();
        if (g_mathFontLoaded) return;
        Sleep(25);
    }
}

const std::wstring& cjk_face_name() {
    static const std::wstring face = []() {
        HDC dc = CreateCompatibleDC(nullptr);
        if (!dc) return std::wstring(L"Yu Mincho");
        const wchar_t* candidates[] = {
            L"Yu Mincho", L"Yu Gothic UI", L"Meiryo", L"MS Mincho"
        };
        const wchar_t sample[] = L"を代入する";
        for (const wchar_t* candidate : candidates) {
            HFONT font = CreateFontW(-64, 0, 0, 0, FW_NORMAL, FALSE, FALSE,
                                     FALSE, DEFAULT_CHARSET, OUT_TT_PRECIS,
                                     CLIP_DEFAULT_PRECIS, CLEARTYPE_QUALITY,
                                     DEFAULT_PITCH | FF_DONTCARE, candidate);
            if (!font) continue;
            HGDIOBJ old = SelectObject(dc, font);
            WORD indices[_countof(sample) - 1] = {};
            const DWORD got = GetGlyphIndicesW(
                dc, sample, int(_countof(sample) - 1), indices,
                GGI_MARK_NONEXISTING_GLYPHS);
            bool owns = got != GDI_ERROR;
            for (WORD index : indices) owns = owns && index != 0xFFFF;
            wchar_t actual[LF_FACESIZE] = {};
            if (owns) GetTextFaceW(dc, _countof(actual), actual);
            SelectObject(dc, old);
            DeleteObject(font);
            if (owns && actual[0]) {
                DeleteDC(dc);
                return std::wstring(actual);
            }
        }
        DeleteDC(dc);
        /* A system text face may still reach CJK through normal font linking.
         * Unlike the process-private math font, its linked glyphs retain the
         * requested em size. */
        return std::wstring(L"Segoe UI");
    }();
    return face;
}

HFONT make_font(bool italic, bool symbol, bool cjk) {
    ensure_math_font();
    LOGFONTW lf = {};
    lf.lfHeight = -kEm;
    lf.lfItalic = italic ? TRUE : FALSE;
    /* DEFAULT_CHARSET, never SYMBOL_CHARSET.  Cambria Math is a Unicode font;
     * asking for the symbol charset makes GDI apply the legacy Symbol code
     * page, so U+0028 is measured as whatever sits at 0x28 in that page and
     * U+2264 is not found at all.  That single flag was behind the spurious
     * gap after "(" and the relation overlapping the fraction bar. */
    lf.lfCharSet = DEFAULT_CHARSET;
    /* One face for everything.  Latin Modern Math is the OpenType Computer
     * Modern -- the typeface TeX itself sets with -- and it carries the text
     * letters, the math alphabets, and every symbol, so the canvas and a
     * pdflatex run finally draw the same shapes and not merely the same
     * boxes.  Its MATH table agrees with the constants taken from tex.web:
     * RadicalRuleThickness and FractionRuleThickness are both 0.040 em and
     * AxisHeight is 0.250 em. */
    const wchar_t* face = cjk ? cjk_face_name().c_str()
                              : L"Latin Modern Math";
    (void)symbol;
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
    /* ink_bottom is where the glyph's ink STOPS above the baseline, and
     * it is deliberately not clamped: an accent such as U+02DC is drawn
     * entirely above the baseline, so its ink bottom is a large positive
     * number and `desc` is zero.  Placing such a glyph by its baseline puts
     * it that whole distance too high. */
    struct Box { double asc = 0, desc = 0, ink_w = 0, ink_left = 0,
                 ink_bottom = 0; };

    HFONT font(int key) {
        auto it = fonts.find(key);
        if (it != fonts.end()) return it->second;
        HFONT f = make_font((key & 1) != 0, (key & 2) != 0,
                            (key & 4) != 0);
        fonts[key] = f;
        return f;
    }

    /* The font's own design units, which is what the MATH table measures in.
     * GDI reports metrics in a 1000-unit em because the face is created at
     * that size, but the table is in the font's units and the two only
     * happen to agree for Latin Modern. */
    double upem(int key) {
        auto it = upems.find(key);
        if (it != upems.end()) return it->second;
        double value = kEm;
        HGDIOBJ old = SelectObject(hdc, font(key));
        const UINT bytes = GetOutlineTextMetricsW(hdc, 0, nullptr);
        if (bytes) {
            std::vector<uint8_t> buffer(bytes);
            auto* otm = reinterpret_cast<OUTLINETEXTMETRICW*>(buffer.data());
            otm->otmSize = bytes;
            if (GetOutlineTextMetricsW(hdc, bytes, otm) && otm->otmEMSquare)
                value = double(otm->otmEMSquare);
        }
        SelectObject(hdc, old);
        upems[key] = value;
        return value;
    }

    /* ---- OpenType MATH: the font's own larger sizes -------------------
     *
     * A display integral is not the text one made bigger.  Scaling the base
     * glyph widens it exactly as much as it heightens it, and Latin Modern's
     * base integral is 0.498 wide per unit of height where TeX's display one
     * is 0.400 -- so the sign grew fat and its curls reached into the limits.
     * The font carries the designed glyph (integral.v1, 0.399) in the MATH
     * table's vertical variants; this reads it out and uses it.
     *
     * The table comes from GDI rather than from the file, so it is whatever
     * the selected face actually provides, and a font without a MATH table
     * simply yields nothing and the caller falls back to scaling. */
    struct Variant { uint16_t glyph = 0; double heightEm = 0; bool valid = false; };

    const std::vector<uint8_t>& math_table(int key) {
        auto it = mathTables.find(key);
        if (it != mathTables.end()) return it->second;
        std::vector<uint8_t> data;
        HGDIOBJ old = SelectObject(hdc, font(key));
        const DWORD tag = 0x4854414D;              /* 'HTAM' == "MATH" LE */
        const DWORD size = GetFontData(hdc, tag, 0, nullptr, 0);
        if (size != GDI_ERROR && size > 0) {
            data.resize(size);
            if (GetFontData(hdc, tag, 0, data.data(), size) == GDI_ERROR)
                data.clear();
        }
        SelectObject(hdc, old);
        return mathTables.emplace(key, std::move(data)).first->second;
    }

    static uint16_t be16(const std::vector<uint8_t>& d, size_t at) {
        return size_t(at + 1) < d.size()
             ? uint16_t((d[at] << 8) | d[at + 1]) : uint16_t(0);
    }

    /* Index of `glyph` in a coverage table, or -1. */
    static int coverage_index(const std::vector<uint8_t>& d, size_t at,
                              uint16_t glyph) {
        const uint16_t format = be16(d, at);
        if (format == 1) {
            const uint16_t count = be16(d, at + 2);
            for (uint16_t i = 0; i < count; ++i)
                if (be16(d, at + 4 + size_t(i) * 2) == glyph) return int(i);
        } else if (format == 2) {
            const uint16_t count = be16(d, at + 2);
            for (uint16_t i = 0; i < count; ++i) {
                const size_t rec = at + 4 + size_t(i) * 6;
                const uint16_t first = be16(d, rec), last = be16(d, rec + 2);
                if (glyph >= first && glyph <= last)
                    return int(be16(d, rec + 4) + (glyph - first));
            }
        }
        return -1;
    }

    /* The smallest designed size at least `wantEm` tall, else the largest. */
    Variant vertical_variant(uint32_t cp, int key, double wantEm) {
        auto k = std::make_pair(std::make_pair(key, cp),
                                int(std::lround(wantEm * 1000.0)));
        auto cached = variants.find(k);
        if (cached != variants.end()) return cached->second;

        Variant best;
        const std::vector<uint8_t>& d = math_table(key);
        uint16_t base = 0;
        {
            HGDIOBJ old = SelectObject(hdc, font(key));
            wchar_t ch = wchar_t(cp);
            if (cp <= 0xFFFF)
                GetGlyphIndicesW(hdc, &ch, 1, &base, GGI_MARK_NONEXISTING_GLYPHS);
            SelectObject(hdc, old);
        }
        if (d.size() > 10 && base && base != 0xFFFF) {
            const size_t variantsAt = be16(d, 8);          /* MathVariants */
            if (variantsAt && variantsAt + 10 <= d.size()) {
                const size_t covAt = variantsAt + be16(d, variantsAt + 2);
                const uint16_t count = be16(d, variantsAt + 6);
                const int at = coverage_index(d, covAt, base);
                if (at >= 0 && at < int(count)) {
                    const size_t conAt = variantsAt +
                        be16(d, variantsAt + 10 + size_t(at) * 2);
                    const uint16_t records = be16(d, conAt + 2);
                    for (uint16_t i = 0; i < records; ++i) {
                        const size_t rec = conAt + 4 + size_t(i) * 4;
                        const uint16_t glyph = be16(d, rec);
                        const double height = double(be16(d, rec + 2)) / upem(key);
                        if (!glyph) continue;
                        const bool better = !best.valid ||
                            (best.heightEm < wantEm && height > best.heightEm) ||
                            (height >= wantEm && height < best.heightEm);
                        if (better) { best.glyph = glyph; best.heightEm = height;
                                      best.valid = true; }
                    }
                }
            }
        }
        variants[k] = best;
        return best;
    }

    /* The font's own italic correction for a character, in em.
     *
     * TeX puts a superscript at the base's width PLUS this, which is why f^2
     * and V^2 need noticeably more room than x^2: the letter leans out past
     * the box it advances by.  Without it the 2 sat on the f's hook -- 0.097
     * em short for f, 0.202 em for V, while x and a looked fine, so it read
     * as a problem with one letter rather than a missing rule. */
    double italic_correction(uint32_t cp, int key) {
        auto k = std::make_pair(key, cp);
        auto it = italics.find(k);
        if (it != italics.end()) return it->second;
        double value = 0.0;
        const std::vector<uint8_t>& d = math_table(key);
        const uint16_t glyph = glyph_id(cp, key);
        if (d.size() > 10 && glyph) {
            const size_t infoAt = be16(d, 6);            /* MathGlyphInfo */
            if (infoAt && infoAt + 2 <= d.size()) {
                const size_t italAt = infoAt + be16(d, infoAt);
                if (italAt > infoAt && italAt + 4 <= d.size()) {
                    const size_t covAt = italAt + be16(d, italAt);
                    const uint16_t count = be16(d, italAt + 2);
                    const int at = coverage_index(d, covAt, glyph);
                    if (at >= 0 && at < int(count)) {
                        const int16_t raw =
                            int16_t(be16(d, italAt + 4 + size_t(at) * 4));
                        value = double(raw) / upem(key);
                    }
                }
            }
        }
        italics[k] = value;
        return value;
    }

    std::map<std::pair<int, uint32_t>, double> italics;

    /* Ink box of a glyph named by index rather than by character. */
    Box glyph_index_box(uint16_t glyph, int key) {
        auto k = std::make_pair(key, uint32_t(glyph) | 0x80000000u);
        auto it = boxes.find(k);
        if (it != boxes.end()) return it->second;
        Box b;
        HGDIOBJ old = SelectObject(hdc, font(key));
        GLYPHMETRICS gm = {};
        MAT2 id = {{0, 1}, {0, 0}, {0, 0}, {0, 1}};
        const DWORD r = GetGlyphOutlineW(hdc, glyph,
            GGO_METRICS | GGO_GLYPH_INDEX, &gm, 0, nullptr, &id);
        SelectObject(hdc, old);
        if (r != GDI_ERROR && gm.gmBlackBoxY > 0) {
            b.asc = std::max(double(gm.gmptGlyphOrigin.y) / kEm, 0.0);
            b.desc = std::max(
                double(int(gm.gmBlackBoxY) - gm.gmptGlyphOrigin.y) / kEm, 0.0);
            b.ink_w = double(gm.gmptGlyphOrigin.x + int(gm.gmBlackBoxX)) / kEm;
            b.ink_left = std::min(0.0, double(gm.gmptGlyphOrigin.x) / kEm);
        }
        boxes[k] = b;
        return b;
    }

    static uint32_t be32(const std::vector<uint8_t>& d, size_t at) {
        return size_t(at + 3) < d.size()
             ? (uint32_t(d[at]) << 24) | (uint32_t(d[at + 1]) << 16) |
               (uint32_t(d[at + 2]) << 8) | uint32_t(d[at + 3])
             : 0u;
    }

    /* Glyph index for any code point, including the math alphabets.
     *
     * GetGlyphIndicesW takes a UTF-16 code unit, so it cannot name U+1D453 --
     * and every italic letter lives up there.  Reading the font's own cmap is
     * the way to reach them, and an index is what GetGlyphOutlineW needs to
     * hand back an outline. */
    uint16_t glyph_id(uint32_t cp, int key) {
        auto k = std::make_pair(key, cp);
        auto it = ids.find(k);
        if (it != ids.end()) return it->second;
        uint16_t id = 0;
        if (cp <= 0xFFFF) {
            HGDIOBJ old = SelectObject(hdc, font(key));
            wchar_t ch = wchar_t(cp);
            GetGlyphIndicesW(hdc, &ch, 1, &id, GGI_MARK_NONEXISTING_GLYPHS);
            SelectObject(hdc, old);
            if (id == 0xFFFF) id = 0;
        } else {
            const std::vector<uint8_t>& d = cmap_table(key);
            const uint16_t tables = be16(d, 2);
            size_t best = 0;
            for (uint16_t i = 0; i < tables; ++i) {
                const size_t rec = 4 + size_t(i) * 8;
                const uint16_t platform = be16(d, rec);
                const uint16_t encoding = be16(d, rec + 2);
                const size_t sub = be32(d, rec + 4);
                /* Only a format-12 subtable covers beyond the BMP. */
                if (sub && be16(d, sub) == 12 &&
                    ((platform == 3 && encoding == 10) || platform == 0))
                    best = sub;
            }
            if (best) {
                const uint32_t groups = be32(d, best + 12);
                for (uint32_t gi = 0; gi < groups; ++gi) {
                    const size_t rec = best + 16 + size_t(gi) * 12;
                    const uint32_t first = be32(d, rec);
                    const uint32_t last = be32(d, rec + 4);
                    if (cp >= first && cp <= last) {
                        id = uint16_t(be32(d, rec + 8) + (cp - first));
                        break;
                    }
                }
            }
        }
        ids[k] = id;
        return id;
    }

    const std::vector<uint8_t>& cmap_table(int key) {
        auto it = cmaps.find(key);
        if (it != cmaps.end()) return it->second;
        std::vector<uint8_t> data;
        HGDIOBJ old = SelectObject(hdc, font(key));
        const DWORD tag = 0x70616D63;             /* 'pamc' == "cmap" LE */
        const DWORD size = GetFontData(hdc, tag, 0, nullptr, 0);
        if (size != GDI_ERROR && size > 0) {
            data.resize(size);
            if (GetFontData(hdc, tag, 0, data.data(), size) == GDI_ERROR)
                data.clear();
        }
        SelectObject(hdc, old);
        return cmaps.emplace(key, std::move(data)).first->second;
    }

    std::map<std::pair<int, uint32_t>, uint16_t> ids;
    std::map<int, std::vector<uint8_t>> cmaps;

    /* The glyph's outline as an SVG path, in a 1-em coordinate space with y
     * downward, so an SVG can draw a glyph the reader's viewer has no way to
     * name.  A size variant has no cmap entry -- it is reachable only by
     * index -- so this is the only way to put integral.v1 into an SVG. */
    const std::string& outline_path(uint16_t glyph, int key) {
        auto k = std::make_pair(key, uint32_t(glyph));
        auto it = outlines.find(k);
        if (it != outlines.end()) return it->second;
        std::string path;
        HGDIOBJ old = SelectObject(hdc, font(key));
        GLYPHMETRICS gm = {};
        MAT2 id = {{0, 1}, {0, 0}, {0, 0}, {0, 1}};
        const DWORD need = GetGlyphOutlineW(hdc, glyph,
            GGO_NATIVE | GGO_GLYPH_INDEX, &gm, 0, nullptr, &id);
        if (need != GDI_ERROR && need > 0) {
            std::vector<uint8_t> buffer(need);
            if (GetGlyphOutlineW(hdc, glyph, GGO_NATIVE | GGO_GLYPH_INDEX,
                                 &gm, need, buffer.data(), &id) != GDI_ERROR) {
                std::ostringstream out;
                out << std::fixed << std::setprecision(4);
                const double s = 1.0 / kEm;
                auto fx = [&](const FIXED& f) {
                    return double(f.value) + double(f.fract) / 65536.0;
                };
                auto pt = [&](const POINTFX& p, char sep) {
                    out << sep << (fx(p.x) * s) << ' ' << (-fx(p.y) * s);
                };
                size_t at = 0;
                while (at + sizeof(TTPOLYGONHEADER) <= buffer.size()) {
                    const auto* head =
                        reinterpret_cast<const TTPOLYGONHEADER*>(&buffer[at]);
                    const size_t end = at + head->cb;
                    if (head->cb == 0 || end > buffer.size()) break;
                    pt(head->pfxStart, 'M');
                    size_t p = at + sizeof(TTPOLYGONHEADER);
                    while (p + sizeof(TTPOLYCURVE) - sizeof(POINTFX) <= end) {
                        const auto* curve =
                            reinterpret_cast<const TTPOLYCURVE*>(&buffer[p]);
                        const WORD n = curve->cpfx;
                        if (curve->wType == TT_PRIM_LINE) {
                            for (WORD i = 0; i < n; ++i) pt(curve->apfx[i], 'L');
                        } else {
                            /* Quadratic B-spline: every point but the last is
                             * a control point, and consecutive controls imply
                             * an on-curve midpoint between them. */
                            for (WORD i = 0; i + 1 < n; ++i) {
                                POINTFX endPoint = curve->apfx[i + 1];
                                if (i + 2 < n) {
                                    endPoint.x.value = SHORT((fx(curve->apfx[i].x) +
                                        fx(curve->apfx[i + 1].x)) / 2.0);
                                    endPoint.x.fract = 0;
                                    endPoint.y.value = SHORT((fx(curve->apfx[i].y) +
                                        fx(curve->apfx[i + 1].y)) / 2.0);
                                    endPoint.y.fract = 0;
                                }
                                pt(curve->apfx[i], 'Q');
                                pt(endPoint, ' ');
                            }
                        }
                        p += sizeof(TTPOLYCURVE) - sizeof(POINTFX) +
                             sizeof(POINTFX) * n;
                    }
                    out << 'Z';
                    at = end;
                }
                path = out.str();
            }
        }
        SelectObject(hdc, old);
        return outlines.emplace(k, std::move(path)).first->second;
    }

    std::map<std::pair<int, uint32_t>, std::string> outlines;
    std::map<int, double> upems;
    std::map<int, std::vector<uint8_t>> mathTables;
    std::map<std::pair<std::pair<int, uint32_t>, int>, Variant> variants;
    double width_em(uint32_t cp, bool italic, bool symbol, bool cjk) {
        int key = (italic ? 1 : 0) | (symbol ? 2 : 0) | (cjk ? 4 : 0);
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
    std::pair<double, double> vmetric(bool italic, bool symbol, bool cjk) {
        int key = (italic ? 1 : 0) | (symbol ? 2 : 0) | (cjk ? 4 : 0);
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
    /* ink_left is the glyph's left side bearing: where its ink begins
     * relative to the pen.  It is negative for a glyph that leans left of its
     * origin -- an italic x or f -- which is why a fraction's denominator
     * could poke a few pixels past the left end of its rule. */

    Box glyph_box(uint32_t cp, bool italic, bool symbol, bool cjk) {
        int key = (italic ? 1 : 0) | (symbol ? 2 : 0) | (cjk ? 4 : 0);
        auto k = std::make_pair(key, cp);
        auto it = boxes.find(k);
        if (it != boxes.end()) return it->second;

        auto v = vmetric(italic, symbol, cjk);
        Box b;
        b.asc = v.first;
        b.desc = v.second;
        b.ink_w = width_em(cp, italic, symbol, cjk);
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
                b.ink_left = std::min(0.0, double(gm.gmptGlyphOrigin.x) / kEm);
                b.ink_bottom =
                    double(gm.gmptGlyphOrigin.y - int(gm.gmBlackBoxY)) / kEm;
            }
        } else if (ink_box_by_drawing(cp, key, b)) {
            /* handled */
        }
        boxes[k] = b;
        return b;
    }

    /* Ink box for a character above the BMP, by drawing it and looking.
     *
     * GetGlyphOutlineW takes a UTF-16 code unit, so it cannot see U+1D44E and
     * friends -- and the math alphabets, which is where TeX's italic letters
     * actually live, are all up there.  Falling back to the font's own ascent
     * and descent is badly wrong for a math font: Latin Modern reserves room
     * for extensible integrals, so a one-letter fraction came out 164 pt tall.
     * TextOutW handles the surrogate pair fine, so draw it once and measure
     * the ink.  Cached per glyph like every other metric here. */
    bool ink_box_by_drawing(uint32_t cp, int key, Box& b) {
        const int box = int(kEm) * 3;
        const int originX = int(kEm), originY = int(kEm) * 2;
        BITMAPINFO bi = {};
        bi.bmiHeader.biSize = sizeof(bi.bmiHeader);
        bi.bmiHeader.biWidth = box;
        bi.bmiHeader.biHeight = -box;          /* top-down */
        bi.bmiHeader.biPlanes = 1;
        bi.bmiHeader.biBitCount = 32;
        bi.bmiHeader.biCompression = BI_RGB;
        void* bits = nullptr;
        HBITMAP bmp = CreateDIBSection(hdc, &bi, DIB_RGB_COLORS, &bits, nullptr, 0);
        if (!bmp || !bits) { if (bmp) DeleteObject(bmp); return false; }
        HDC mem = CreateCompatibleDC(hdc);
        HGDIOBJ oldBmp = SelectObject(mem, bmp);
        RECT all{0, 0, box, box};
        FillRect(mem, &all, HBRUSH(GetStockObject(WHITE_BRUSH)));
        HGDIOBJ oldFont = SelectObject(mem, font(key));
        SetBkMode(mem, TRANSPARENT);
        SetTextColor(mem, RGB(0, 0, 0));
        SetTextAlign(mem, TA_LEFT | TA_BASELINE | TA_NOUPDATECP);
        wchar_t pair[2];
        const uint32_t v = cp - 0x10000;
        pair[0] = wchar_t(0xD800 + (v >> 10));
        pair[1] = wchar_t(0xDC00 + (v & 0x3FF));
        TextOutW(mem, originX, originY, pair, 2);
        GdiFlush();
        const auto* px = static_cast<const unsigned char*>(bits);
        int top = box, bottom = -1, leftMost = box, rightMost = -1;
        for (int y = 0; y < box; ++y) {
            for (int x = 0; x < box; ++x) {
                if (px[(size_t(y) * box + x) * 4] < 200) {   /* blue channel */
                    if (y < top) top = y;
                    if (y > bottom) bottom = y;
                    if (x < leftMost) leftMost = x;
                    if (x > rightMost) rightMost = x;
                }
            }
        }
        SelectObject(mem, oldFont);
        SelectObject(mem, oldBmp);
        DeleteDC(mem);
        DeleteObject(bmp);
        if (bottom < 0) return false;
        b.asc = std::max(double(originY - top) / kEm, 0.0);
        b.desc = std::max(double(bottom - originY + 1) / kEm, 0.0);
        b.ink_w = std::max(b.ink_w, double(rightMost - originX + 1) / kEm);
        b.ink_left = std::min(0.0, double(leftMost - originX) / kEm);
        b.ink_bottom = double(originY - bottom) / kEm;
        return true;
    }

    std::map<std::pair<int, uint32_t>, Box> boxes;

    /* The flat top of a glyph that a rule is meant to continue -- the flag of
     * a radical.  A vinculum drawn at the font's nominal rule thickness and
     * starting at the advance width does not meet it: in Cambria Math the
     * flag is 0.065 em thick and ends at 0.710 em while the advance is
     * 0.655 em, so the bar came out thinner than the flag and started past
     * its end.  Measuring the glyph is the only way to join them exactly. */
    struct Plateau { double thickness = 0, left = 0, right = 0; bool valid = false; };

    /* `byIndex` selects a glyph the font has no character for -- a designed
     * size variant -- so a grown radical can be measured the same way as the
     * base one. */
    Plateau plateau(uint32_t cp, bool italic, bool symbol,
                    bool byIndex = false) {
        int key = (italic ? 1 : 0) | (symbol ? 2 : 0);
        auto k = std::make_pair(byIndex ? key | 4 : key, cp);
        auto it = plateaus.find(k);
        if (it != plateaus.end()) return it->second;

        const UINT format = GGO_GRAY8_BITMAP | (byIndex ? GGO_GLYPH_INDEX : 0);
        Plateau p;
        HGDIOBJ old = SelectObject(hdc, font(key));
        GLYPHMETRICS gm = {};
        MAT2 id = {{0, 1}, {0, 0}, {0, 0}, {0, 1}};
        const DWORD need =
            GetGlyphOutlineW(hdc, cp, format, &gm, 0, nullptr, &id);
        if (need != GDI_ERROR && need > 0 && gm.gmBlackBoxX > 0 &&
            gm.gmBlackBoxY > 0) {
            std::vector<unsigned char> bits(need, 0);
            if (GetGlyphOutlineW(hdc, cp, format, &gm, need,
                                 bits.data(), &id) != GDI_ERROR) {
                const int w = int(gm.gmBlackBoxX), h = int(gm.gmBlackBoxY);
                const int stride = (w + 3) & ~3;
                auto row_extent = [&](int row, int& first, int& last) {
                    first = -1; last = -1;
                    for (int c = 0; c < w; ++c) {
                        const size_t at = size_t(row) * size_t(stride) + size_t(c);
                        if (at < bits.size() && bits[at] > 8) {
                            if (first < 0) first = c;
                            last = c;
                        }
                    }
                };
                int f0 = -1, l0 = -1;
                row_extent(0, f0, l0);
                if (f0 >= 0) {
                    const double reference = double(l0 - f0 + 1);
                    /* The flat top starts at the TOP row's left edge.  Taking
                     * the minimum over every plateau row walked down the
                     * diagonal instead, which clipped off the part of the
                     * stroke that rises to meet the bar. */
                    int rows = 1, right = l0;
                    const int left = f0;
                    for (int r = 1; r < h; ++r) {
                        int f = -1, l = -1;
                        row_extent(r, f, l);
                        if (f < 0 || double(l - f + 1) < reference * 0.6) break;
                        right = std::max(right, l);
                        ++rows;
                    }
                    p.thickness = double(rows) / kEm;
                    p.left = double(gm.gmptGlyphOrigin.x + left) / kEm;
                    p.right = double(gm.gmptGlyphOrigin.x + right + 1) / kEm;
                    p.valid = true;
                }
            }
        }
        SelectObject(hdc, old);
        plateaus[k] = p;
        return p;
    }

    std::map<std::pair<int, uint32_t>, Plateau> plateaus;
};

MetricCache& metrics() {
    static MetricCache c;
    return c;
}

double char_width(uint32_t cp, double sizePt, bool italic, bool symbol,
                  bool cjk = false) {
    return metrics().width_em(cp, italic, symbol, cjk) * sizePt;
}
/* The font's own extent, for things sized against the face rather than
 * against one glyph (fence stretching, the fallback line height). */
void char_vmetrics(double sizePt, bool italic, bool symbol, bool cjk,
                   double& asc, double& desc) {
    auto v = metrics().vmetric(italic, symbol, cjk);
    asc = v.first * sizePt;
    desc = v.second * sizePt;
}
/* The ink box of one glyph, which is what neighbours and scripts must clear. */
void glyph_vmetrics(uint32_t cp, double sizePt, bool italic, bool symbol,
                    bool cjk,
                    double& asc, double& desc) {
    auto b = metrics().glyph_box(cp, italic, symbol, cjk);
    asc = b.asc * sizePt;
    desc = b.desc * sizePt;
}
/* How far the drawing actually reaches.  A large operator is drawn wider than
 * it advances, so laying the next atom out at the advance alone lets a sigma
 * touch the symbol after it. */
double glyph_ink_width(uint32_t cp, double sizePt, bool italic, bool symbol,
                       bool cjk = false) {
    return metrics().glyph_box(cp, italic, symbol, cjk).ink_w * sizePt;
}

/* Where a glyph's ink stops, measured up from its baseline.  Positive for
 * an accent, which floats entirely above it. */
double glyph_ink_bottom(uint32_t cp, double sizePt, bool italic, bool symbol,
                        bool cjk = false) {
    return metrics().glyph_box(cp, italic, symbol, cjk).ink_bottom * sizePt;
}

/* Left side bearing in points, <= 0.  Zero for a glyph whose ink starts at or
 * after its origin. */
double glyph_ink_left(uint32_t cp, double sizePt, bool italic, bool symbol,
                      bool cjk = false) {
    return metrics().glyph_box(cp, italic, symbol, cjk).ink_left * sizePt;
}

/* The font's designed size at least `wantEm` tall for this character, if it
 * has one.  Returns false when the font carries no MATH variants. */
bool glyph_size_variant(uint32_t cp, double wantEm, unsigned short& glyph,
                        double& heightEm) {
    auto v = metrics().vertical_variant(cp, 2 /* symbol face */, wantEm);
    if (!v.valid) return false;
    glyph = v.glyph;
    heightEm = v.heightEm;
    return true;
}

void glyph_index_metrics(unsigned short glyph, double sizePt,
                         double& asc, double& desc, double& inkWidth) {
    auto b = metrics().glyph_index_box(glyph, 2);
    asc = b.asc * sizePt;
    desc = b.desc * sizePt;
    inkWidth = b.ink_w * sizePt;
}

std::string glyph_index_outline(unsigned short glyph) {
    return metrics().outline_path(glyph, 2);
}


void ensure_math_font_public() { ensure_math_font(); }
bool math_font_is_loaded() { return g_mathFontLoaded; }

/* The font's italic correction for a character, in points at this size.  Zero
 * for an upright glyph, and for a font with no MATH table. */
double glyph_italic_correction(uint32_t cp, double sizePt) {
    return metrics().italic_correction(cp, 2) * sizePt;
}

/* Glyph index for a character in the face that will actually draw it. */
unsigned short glyph_id_of(uint32_t cp, bool cjk = false) {
    return metrics().glyph_id(cp, cjk ? 4 : 2);
}

/* Thickness and right edge of a glyph's flat top, in points, so a rule can
 * continue it seamlessly.  False when the glyph has no such top. */
bool glyph_top_plateau(uint32_t cp, double sizePt, bool italic, bool symbol,
                       double& thickness, double& left, double& right,
                       bool byIndex = false) {
    auto p = metrics().plateau(cp, italic, symbol, byIndex);
    if (!p.valid) return false;
    thickness = p.thickness * sizePt;
    left = p.left * sizePt;
    right = p.right * sizePt;
    return true;
}
#else
/* Portable fallback: enough to keep the tree walking and the tests honest
 * about being unmeasured, not enough for production output. */
double char_width(uint32_t, double sizePt, bool, bool, bool = false) {
    return 0.5 * sizePt;
}
void char_vmetrics(double sizePt, bool, bool, bool,
                   double& asc, double& desc) {
    asc = 0.75 * sizePt;
    desc = 0.25 * sizePt;
}
void glyph_vmetrics(uint32_t, double sizePt, bool, bool, bool,
                    double& asc, double& desc) {
    asc = 0.70 * sizePt;
    desc = 0.05 * sizePt;
}
double glyph_ink_width(uint32_t, double sizePt, bool, bool, bool = false) {
    return 0.5 * sizePt;
}
double glyph_ink_left(uint32_t, double, bool, bool, bool = false) { return 0.0; }
#endif

/* ------------------------------------------------------------------ */
/* Typeface mapping                                                    */
/* ------------------------------------------------------------------ */
bool typeface_is_italic(int tf) {
    return tf == TF_VARIABLE || tf == TF_LCGREEK ||
           tf == TF_MATH_ITALIC;
}

/* Japanese must never reach the process-private math font through GDI font
 * linking.  The linked CJK glyph can inherit the math face's 1000-unit metric
 * scale twice: a nominal 12 pt character then advances more than 50 pt and is
 * painted correspondingly huge.  Give CJK its own measured system face. */
bool needs_cjk_face(uint32_t cp) {
    return (cp >= 0x2E80 && cp <= 0x2FFF) ||
           (cp >= 0x3000 && cp <= 0x30FF) ||
           (cp >= 0x31F0 && cp <= 0x31FF) ||
           (cp >= 0x3400 && cp <= 0x4DBF) ||
           (cp >= 0x4E00 && cp <= 0x9FFF) ||
           (cp >= 0xF900 && cp <= 0xFAFF) ||
           (cp >= 0xFF00 && cp <= 0xFFEF) ||
           (cp >= 0x20000 && cp <= 0x323AF);
}

/* Which family draws a code point.  This is a property of the character, not
 * of the stored typeface: a TF_SYMBOL "(" is still an ordinary parenthesis and
 * belongs in the text face. */
bool needs_math_face(uint32_t cp) {
    return cp > 0xFF && !needs_cjk_face(cp);
}

/* The math-italic code point for a letter TeX would set in math italic, or 0
 * when the character has none and should be drawn as it is.  Uppercase Greek
 * is deliberately absent: LaTeX sets it upright. */
uint32_t math_italic_of(uint32_t cp) {
    if (cp >= 'a' && cp <= 'z') {
        /* Unicode has no MATHEMATICAL ITALIC SMALL H: the slot is reserved
         * because U+210E PLANCK CONSTANT already is that glyph. */
        if (cp == 'h') return 0x210E;
        return 0x1D44E + (cp - 'a');
    }
    if (cp >= 'A' && cp <= 'Z') return 0x1D434 + (cp - 'A');
    if (cp >= 0x03B1 && cp <= 0x03C9) return 0x1D6FC + (cp - 0x03B1);
    switch (cp) {                       /* the variant Greek letters */
        case 0x03D1: return 0x1D717;    /* vartheta  */
        case 0x03D5: return 0x1D719;    /* phi       */
        case 0x03D6: return 0x1D71B;    /* varpi     */
        case 0x03F0: return 0x1D718;    /* varkappa  */
        case 0x03F1: return 0x1D71A;    /* varrho    */
        case 0x03F5: return 0x1D716;    /* epsilon   */
        default: return 0;
    }
}

/* The designed Unicode mathematical-bold glyph for a character, or 0 when
 * the font has no dedicated bold alphabet entry.  `CreateFontW(FW_BOLD)` on
 * the process-private regular math face is only a synthetic request and was
 * visibly ignored on the native canvas: \mathbf{E} looked identical to E.
 * Latin Modern Math ships these real glyphs, so use the same kind of explicit
 * math-alphabet mapping already used for italic variables above. */
uint32_t math_bold_of(uint32_t cp) {
    if (cp >= 'A' && cp <= 'Z') return 0x1D400 + (cp - 'A');
    if (cp >= 'a' && cp <= 'z') return 0x1D41A + (cp - 'a');
    if (cp >= '0' && cp <= '9') return 0x1D7CE + (cp - '0');
    if (cp >= 0x0391 && cp <= 0x03A9 && cp != 0x03A2)
        return 0x1D6A8 + (cp - 0x0391);
    if (cp >= 0x03B1 && cp <= 0x03C9)
        return 0x1D6C2 + (cp - 0x03B1);
    switch (cp) {
        case 0x2202: return 0x1D6DB;    /* partial differential */
        case 0x03F5: return 0x1D6DC;    /* epsilon symbol */
        case 0x03D1: return 0x1D6DD;    /* theta symbol */
        case 0x03F0: return 0x1D6DE;    /* kappa symbol */
        case 0x03D5: return 0x1D6DF;    /* phi symbol */
        case 0x03F1: return 0x1D6E0;    /* rho symbol */
        case 0x03D6: return 0x1D6E1;    /* pi symbol */
        default: return 0;
    }
}

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

/* Computer Modern font dimensions, in em, as TeX uses them for display
 * style.  These are not taste: each was checked against pdfLaTeX geometry
 * by tools/tex_geometry.py, and tests/test_tex_geometry.py holds the
 * comparison so they cannot drift back to being guesses. */
constexpr double kAxisHeight = 0.250;      /* axis_height   */
constexpr double kRuleThickness = 0.040;   /* default_rule_thickness */

/* Every horizontal rule -- fraction bar, vinculum, overline, brace -- is
 * default_rule_thickness in TeX, and none of them varies with what it sits
 * over.  These were all floored at 0.6 pt, which at 12 pt overrode TeX's
 * 0.48 pt and made every bar 25 % too thick; the radical also scaled its bar
 * by the glyph stretch on top of that.  No floor is needed for visibility:
 * the painter floors the top and ceils the bottom, so a rule always covers
 * at least one pixel.
 */
inline double rule_thickness(double sizePt) {
    return std::max(0.05, kRuleThickness * sizePt);
}
constexpr double kSupShiftUp = 0.413;      /* sup1  */
constexpr double kSubShiftDown = 0.150;    /* sub1  */

uint32_t first_codepoint(const std::string& utf8) {
    if (utf8.empty()) return 0;
    const unsigned char lead = static_cast<unsigned char>(utf8[0]);
    size_t count = lead < 0x80 ? 1
                 : (lead & 0xE0) == 0xC0 ? 2
                 : (lead & 0xF0) == 0xE0 ? 3
                 : (lead & 0xF8) == 0xF0 ? 4 : 1;
    if (count == 1 || utf8.size() < count) return lead;
    uint32_t cp = lead & (0xFF >> (count + 1));
    for (size_t i = 1; i < count; ++i)
        cp = (cp << 6) | (static_cast<unsigned char>(utf8[i]) & 0x3F);
    return cp;
}

/* The rightmost ink, which is not the same as the advance width: an italic
 * `f` leans past its own advance, so a box measured from advances alone
 * clipped it off the right edge of every rendered image. */
double layout_ink_right(const Layout& L) {
    double right = L.w;
    for (const auto& g : L.glyphs) {
        if (g.text.empty()) continue;
        const uint32_t cp = first_codepoint(g.text);
        right = std::max(right,
                         g.x + glyph_ink_width(cp, g.size, g.italic, g.symbol,
                                               g.cjk));
    }
    for (const auto& rule : L.rules) right = std::max(right, rule.x + rule.w);
    return right;
}

/* The leftmost ink, which can be left of 0: an italic x or f leans past its
 * origin, so content centred on advances alone pokes past a fraction rule
 * that starts at 0.  Returns <= 0. */
double layout_ink_left(const Layout& L) {
    double left = 0.0;
    for (const auto& g : L.glyphs) {
        if (g.text.empty()) continue;
        const uint32_t cp = first_codepoint(g.text);
        left = std::min(left,
                        g.x + glyph_ink_left(cp, g.size, g.italic, g.symbol,
                                             g.cjk));
    }
    for (const auto& rule : L.rules) left = std::min(left, rule.x);
    return left;
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

    Layout run(const LineNode& root) { return layout_list(root.children, st_.full); }

private:
    const SvgStyle& st_;

    /* TeX's style chain: display, text, script, scriptscript.  Display and
     * text share a size, so the size alone cannot say which one a construct
     * is in -- and a fraction's parts drop a style whether or not that
     * changes the size.  Without tracking it, nested fractions stayed at full
     * size where pdfLaTeX sets 12, 8 and 6 pt, so they came out a third too
     * tall and dragged the fences around them with them. */
    enum Style { kDisplay = 0, kText = 1, kScript = 2, kScriptScript = 3 };
    int style_ = kDisplay;

    /* The style a fraction's numerator and denominator are set in, and the
     * size that goes with it.  tex.web 702: D and T give T, S and SS give SS. */
    int fraction_part_style() const {
        return style_ == kDisplay ? kText
             : style_ == kText    ? kScript : kScriptScript;
    }
    double size_for_child_style(double sizePt, int child) const {
        if (child <= kText) return sizePt;              /* same size as D/T */
        if (style_ <= kText) return script_size(sizePt);
        return script_size(sizePt);
    }

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

    Layout glyph_layout(uint32_t cp, double sizePt, bool italic, bool symbol) {
        /* Latin Modern Math has no italic face.  Asking GDI for one makes it
         * slant the upright glyphs, which is not the same shape as TeX's math
         * italic at all.  The real italic letters live in the Unicode math
         * alphanumeric block, so a variable is drawn as its math-italic code
         * point, upright, out of the math font. */
        if (italic) {
            const uint32_t mathItalic = math_italic_of(cp);
            if (mathItalic) { cp = mathItalic; italic = false; symbol = true; }
        }
        const bool cjk = needs_cjk_face(cp);
        if (cjk) {
            italic = false;
            symbol = false;
        }
        Layout L;
        L.w = char_width(cp, sizePt, italic, symbol, cjk);
        glyph_vmetrics(cp, sizePt, italic, symbol, cjk, L.asc, L.desc);
        Glyph g;
        g.x = 0; g.y = 0; g.size = sizePt;
        g.italic = italic; g.symbol = symbol; g.cjk = cjk;
        g.text = utf8_of(cp);
        L.glyphs.push_back(g);
        return L;
    }

    Layout text_layout(const std::string& text, double sizePt,
                       bool italic = false, bool bold = false) {
        Layout out;
        double x = 0;
        for (uint32_t cp : utf8_codes(text)) {
            Layout g = glyph_layout(cp, sizePt, italic, needs_math_face(cp));
            for (auto& gg : g.glyphs) gg.bold = bold;
            out.absorb(g, x, 0);
            x += g.w;
            out.asc = std::max(out.asc, g.asc);
            out.desc = std::max(out.desc, g.desc);
        }
        out.w = x;
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
            case Node::kFunction:
            case Node::kLim:
                return kOp;
            case Node::kFence:
                return kInner;
            case Node::kFrac:
            case Node::kMatrix:
            case Node::kPile:
                return kInner;
            case Node::kLine: {
                /* A group takes the class of its first atom, so \sin(x)
                 * spaces like a function and not like a bare group.
                 *
                 * A function name is the exception: \sin arrives as a run of
                 * Function-styled letters, whose first atom is an ordinary
                 * "s", so the whole word was spaced as ordinary and TeX's
                 * thin space after an operator went missing -- "sin\omega t"
                 * set as sinωt where TeX sets sin ωt.  The letters inside
                 * keep their own ordinary class, so the word does not fall
                 * apart into s i n. */
                const auto& l = static_cast<const LineNode&>(n);
                for (const auto& c : l.children) {
                    if (!c || c->tag() == Node::kSize) continue;
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

    Layout layout_list(const NodeList& list, double sizePt) {
        Layout out;
        double x = 0, cur = sizePt;
        bool have_prev = false;
        AtomClass prev = kOrd;

        if (list.empty()) {
            /* Empty template holes are real editing positions.  Giving them
             * geometry keeps a blank numerator or matrix cell clickable and
             * lets the native editor draw the familiar dotted placeholder. */
            out.w = 0.55 * sizePt;
            out.asc = 0.68 * sizePt;
            out.desc = 0.18 * sizePt;
            out.placeholders.push_back(
                {0, -out.asc, out.w, out.asc + out.desc});
            out.carets.push_back({&list, 0, 0.5 * out.w,
                                  -out.asc, out.desc});
            return out;
        }

        std::vector<double> boundary(list.size() + 1, 0.0);

        for (size_t ni = 0; ni < list.size(); ++ni) {
            const auto& n = list[ni];
            if (!n) continue;
            if (n->tag() == Node::kSize) {
                cur = size_of(static_cast<const SizeNode*>(n.get())->sizeType);
                boundary[ni + 1] = x;
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
                double gap = space_mu(prev, cls) * cur / 18.0;
                boundary[ni] = x + gap * 0.5;
                x += gap;
            } else {
                boundary[ni] = x;
            }

            Layout piece = layout_node(*n, cur);
            out.absorb(piece, x, 0);
            x += piece.w;
            boundary[ni + 1] = x;
            out.asc = std::max(out.asc, piece.asc);
            out.desc = std::max(out.desc, piece.desc);
            prev = cls;
            have_prev = true;
        }
        out.w = x;
        double caretAsc = out.asc, caretDesc = out.desc;
        if (caretAsc + caretDesc < 0.25 * sizePt) {
            caretAsc = 0.68 * sizePt;
            caretDesc = 0.18 * sizePt;
            out.asc = std::max(out.asc, caretAsc);
            out.desc = std::max(out.desc, caretDesc);
        }
        for (size_t i = 0; i < boundary.size(); ++i)
            out.carets.push_back({&list, int(i), boundary[i],
                                  -caretAsc, caretDesc});
        return out;
    }

    Layout layout_node(const Node& n, double sizePt) {
        switch (n.tag()) {
            case Node::kLine: {
                const auto& ln = static_cast<const LineNode&>(n);
                if (ln.isNull) return Layout();
                return layout_list(ln.children, sizePt);
            }
            case Node::kChar: {
                const auto& c = static_cast<const CharNode&>(n);
                /* Explicit spacing is an advance with no ink.  Its height
                 * stays zero so it never stretches a fence or a radical. */
                if (c.typeface == TF_SPACE) {
                    Layout space;
                    space.w = space_width_em(c.latex.c_str()) * sizePt;
                    space.editMarks.push_back({EditMarkKind::Space,
                        space.w * 0.5, -0.08 * sizePt, 0.08 * sizePt});
                    return space;
                }
                uint32_t cp = c.charCode ? c.charCode : uint32_t(uint8_t(c.ch));
                if (!cp) return Layout();
                if (c.typeface == TF_VECTOR) {
                    if (const uint32_t mathBold = math_bold_of(cp))
                        return glyph_layout(mathBold, sizePt, false, true);
                    Layout fallback = glyph_layout(cp, sizePt, false,
                                                   needs_math_face(cp));
                    for (auto& glyph : fallback.glyphs) glyph.bold = true;
                    return fallback;
                }
                return glyph_layout(cp, sizePt, typeface_is_italic(c.typeface),
                                    needs_math_face(cp));
            }
            case Node::kFunction: {
                const auto& f = static_cast<const FunctionNode&>(n);
                return text_layout(f.name, sizePt, false);
            }
            case Node::kText: {
                const auto& t = static_cast<const TextNode&>(n);
                return text_layout(t.text, sizePt, false);
            }
            case Node::kScript:   return layout_script(static_cast<const ScriptNode&>(n), sizePt);
            case Node::kFence:    return layout_fence(static_cast<const FenceNode&>(n), sizePt);
            case Node::kFrac:     return layout_frac(static_cast<const FracNode&>(n), sizePt);
            case Node::kSqrt:     return layout_sqrt(static_cast<const SqrtNode&>(n), sizePt);
            case Node::kIntegral: {
                const auto& i = static_cast<const IntegralNode&>(n);
                return layout_bigop(bigop_glyph(i.selector), i.body, i.lower, i.upper,
                                    i.hasLower, i.hasUpper, i.hasLimits, sizePt);
            }
            case Node::kBigOp: {
                const auto& b = static_cast<const BigOpNode&>(n);
                return layout_bigop(bigop_glyph(b.selector), b.body, b.lower, b.upper,
                                    b.hasLower, b.hasUpper, b.hasLimits, sizePt);
            }
            case Node::kMatrix:
                return layout_matrix(static_cast<const MatrixNode&>(n), sizePt);
            case Node::kPile:
                return layout_pile(static_cast<const PileNode&>(n), sizePt);
            case Node::kEmbell:
                return layout_embell(static_cast<const EmbellNode&>(n), sizePt);
            case Node::kDecoration:
                return layout_decoration(static_cast<const DecorationNode&>(n), sizePt);
            case Node::kBraceDeco:
                return layout_brace(static_cast<const BraceDecoNode&>(n), sizePt);
            case Node::kDirac:
                return layout_dirac(static_cast<const DiracNode&>(n), sizePt);
            case Node::kLim:
                return layout_lim(static_cast<const LimNode&>(n), sizePt);
            case Node::kMathbf: {
                Layout b = layout_list(static_cast<const MathbfNode&>(n).content, sizePt);
                for (auto& g : b.glyphs) g.bold = true;
                return b;
            }
            case Node::kGroup:
                return layout_list(static_cast<const GroupNode&>(n).children, sizePt);
            case Node::kPrime: {
                Layout out;
                int count = std::max(1, static_cast<const PrimeNode&>(n).count);
                for (int i = 0; i < count; ++i) {
                    Layout p = glyph_layout(0x2032, script_size(sizePt), false, true);
                    out.absorb(p, out.w, -0.35 * sizePt);
                    out.w += p.w;
                    out.asc = std::max(out.asc, 0.35 * sizePt + p.asc);
                    out.desc = std::max(out.desc, p.desc);
                }
                return out;
            }
            case Node::kDegree:
                return glyph_layout(0x00B0, script_size(sizePt), false, false);
            case Node::kOverset:
                return layout_overset(static_cast<const OversetNode&>(n), sizePt);
            default:
                return layout_fallback(n, sizePt);
        }
    }

    /* The italic correction of the last character in a base, which is the one
     * a script has to clear.  Zero unless the base ends in a plain character:
     * a fraction or a fence has no lean to correct for. */
    double last_italic_correction(const NodeList& base, double sizePt) const {
        for (auto it = base.rbegin(); it != base.rend(); ++it) {
            if (!*it) continue;
            if ((*it)->tag() != Node::kChar) return 0.0;
            const auto& ch = static_cast<const CharNode&>(**it);
            if (!typeface_is_italic(ch.typeface)) return 0.0;
            uint32_t cp = ch.charCode ? ch.charCode : uint32_t(uint8_t(ch.ch));
            const uint32_t mathItalic = math_italic_of(cp);
            if (mathItalic) cp = mathItalic;
            return glyph_italic_correction(cp, sizePt);
        }
        return 0.0;
    }

    Layout layout_script(const ScriptNode& s, double sizePt) {
        Layout base = layout_list(s.base, sizePt);
        double ss = script_size(sizePt);
        /* A script is in script style, and one inside a script is in
         * scriptscript, so anything nested in an exponent -- a fraction, say
         * -- reduces from there rather than starting over at display. */
        const int outerStyle = style_;
        const int scriptStyle = style_ <= kText ? kScript : kScriptScript;
        Layout out = base;
        double x = base.w;
        /* Scripts hang off the base's own extents, not off a fixed offset:
         * the exponent of a tall base has to clear that base, and the
         * subscript of a deep one has to sit below it. */
        /* sup1 and sub1 are the Computer Modern shifts TeX uses in display
         * style; measured against pdfLaTeX they are 0.411 and 0.149 em, and
         * the old 0.45/0.22 sat a superscript 0.04 em high and a subscript
         * 0.07 em low.  The max() keeps a tall or deep base clear. */
        const double supShift = std::max(kSupShiftUp * sizePt,
                                         base.asc - 0.35 * ss);
        const double subShift = std::max(kSubShiftDown * sizePt,
                                         base.desc + 0.12 * ss);
        /* A superscript starts at the base's width PLUS the base's italic
         * correction; a subscript starts at the width alone.  That asymmetry
         * is TeX's (tex.web 756) and it is what keeps the 2 of f^2 clear of
         * the f's hook -- 0.097 em of clearance for f, 0.202 for V, and next
         * to nothing for x, which is why leaving it out looked like a problem
         * with one letter rather than a missing rule. */
        const double italic = last_italic_correction(s.base, sizePt);
        double wsub = 0, wsup = 0;
        style_ = scriptStyle;
        if (s.hasSup) {
            Layout sup = layout_list(s.sup, ss);
            out.absorb(sup, x + italic, -supShift);
            out.asc = std::max(out.asc, supShift + sup.asc);
            wsup = italic + sup.w;
        }
        if (s.hasSub) {
            Layout sub = layout_list(s.sub, ss);
            out.absorb(sub, x, subShift);
            out.desc = std::max(out.desc, subShift + sub.desc);
            wsub = sub.w;
        }
        style_ = outerStyle;
        out.w = x + std::max(wsub, wsup);
        return out;
    }

    /* One delimiter, grown to `needHeight` and centred on the math axis.
     *
     * Prefers the font's designed larger glyph, exactly as a big operator
     * does.  Stretching the base glyph widens nothing while the height grows,
     * so a tall parenthesis came out 0.378 wide per unit of height where TeX
     * draws 0.762 -- a thin, rubbery bracket instead of a designed one.
     * Falls back to stretching when the font offers no variant tall enough,
     * which is what has to happen for a delimiter around a very deep box. */
    Layout sized_delimiter(uint32_t cp, double sizePt, double needHeight) {
        const bool symbol = needs_math_face(cp);
        const double axis = kAxisHeight * sizePt;
        Layout g;
        unsigned short variant = 0;
        double variantHeight = 0;
        if (glyph_size_variant(cp, needHeight / sizePt, variant,
                               variantHeight) && variantHeight > 1e-6) {
            double vAsc = 0, vDesc = 0, vInk = 0;
            glyph_index_metrics(variant, sizePt, vAsc, vDesc, vInk);
            Glyph glyph;
            glyph.x = 0; glyph.y = 0; glyph.size = sizePt; glyph.symbol = true;
            glyph.glyphIndex = variant;
            glyph.outline = glyph_index_outline(variant);
            glyph.text = utf8_of(cp);
            g.glyphs.push_back(glyph);
            g.asc = vAsc; g.desc = vDesc;
            g.w = std::max(vInk, char_width(cp, sizePt, false, symbol));
            /* The largest designed size can still fall short of a very deep
             * box; stretch that one the rest of the way. */
            const double have = g.asc + g.desc;
            if (have > 1e-6 && needHeight > have) {
                const double extra = needHeight / have;
                for (auto& gg : g.glyphs) gg.stretchY = extra;
                g.asc *= extra; g.desc *= extra;
            }
        } else {
            g = glyph_layout(cp, sizePt, false, symbol);
            const double natural = std::max(g.asc + g.desc, 1e-6);
            const double stretch = std::max(1.0, needHeight / natural);
            if (stretch > 1.0) {
                for (auto& gg : g.glyphs) gg.stretchY = stretch;
                g.asc *= stretch; g.desc *= stretch;
            }
        }
        /* TeX centres a delimiter on the axis, like a display operator. */
        const double lift = axis - (g.asc - g.desc) / 2.0;
        if (std::fabs(lift) > 1e-9) {
            Layout centred;
            centred.absorb(g, 0, -lift);
            centred.w = g.w;
            centred.asc = g.asc + lift;
            centred.desc = g.desc - lift;
            return centred;
        }
        return g;
    }

    Layout layout_fence(const FenceNode& f, double sizePt) {
        Layout inner = layout_list(f.content, sizePt);
        auto gl = fence_glyphs(f.selector);
        /* \left( x \right] keeps its own closing delimiter. */
        gl.second = fence_glyphs(f.right_selector()).second;
        /* How tall the delimiter has to be, measured the way TeX measures it:
         * symmetrically about the math axis, since that is where a delimiter
         * is centred.  Content that reaches further above the axis than below
         * it still needs a delimiter tall enough for both halves. */
        const double axis = kAxisHeight * sizePt;
        const double half = std::max(inner.asc - axis, inner.desc + axis);
        const double need = 2.0 * half;

        Layout out;
        double x = 0;
        bool left = (f.variation == 0 || f.variation == 1);
        bool right = (f.variation == 0 || f.variation == 2);
        if (left) {
            Layout g = sized_delimiter(gl.first, sizePt, need);
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
            Layout g = sized_delimiter(gl.second, sizePt, need);
            out.absorb(g, x, 0);
            out.asc = std::max(out.asc, g.asc);
            out.desc = std::max(out.desc, g.desc);
            x += g.w;
        }
        out.w = x;
        return out;
    }

    Layout layout_frac(const FracNode& f, double sizePt) {
        /* Set the parts one style down, as TeX does.  The size only changes
         * once the chain reaches script, which is why the first level of
         * nesting stays at full size and the second drops to 8 pt. */
        const int outer = style_;
        const int part = fraction_part_style();
        const double partSize = size_for_child_style(sizePt, part);
        style_ = part;
        Layout num = layout_list(f.numer, partSize);
        Layout den = layout_list(f.denom, partSize);
        style_ = outer;
        /* sizePt stays the fraction's own: TeX takes the bar thickness and
         * the num/denom shifts from the style the fraction is in, not from
         * the smaller one its parts are set in. */
        /* TeX's rule (tex.web 704): the parts sit at fixed shifts from the
         * baseline, not at a fixed gap from the bar, and the rule spans the
         * fraction exactly.  Measured against pdfLaTeX the old constants put
         * the numerator 0.24 em too close to the bar and made the bar 0.39 em
         * wider than the fraction.  num1/denom1/axis are the Computer Modern
         * font dimensions for display style. */
        const double axis = kAxisHeight * sizePt;
        const double thick = rule_thickness(sizePt);
        /* A fraction outside display style is set tighter: TeX swaps num1 and
         * denom1 for num2 and denom2 and drops the clearance from three rule
         * thicknesses to one (tex.web 704).  Using the display numbers
         * everywhere left a nested fraction far taller than its width -- and
         * matching the font sizes alone made that worse, because the parts
         * shrank while the gaps around them did not. */
        const bool display = (style_ == kDisplay);
        const double clr = (display ? 3.0 : 1.0) * thick;
        double up = (display ? 0.677 : 0.394) * sizePt;    /* num1  / num2  */
        double down = (display ? 0.686 : 0.345) * sizePt;  /* denom1/denom2 */

        const double gapUp = (up - num.desc) - (axis + thick / 2.0);
        if (gapUp < clr) up += clr - gapUp;
        const double gapDown = (axis - thick / 2.0) - (den.asc - down);
        if (gapDown < clr) down += clr - gapDown;

        /* TeX's rule is as wide as the wider part's box (its advance width),
         * which is what the 0.512-em reference for \frac{a}{b} is.  The parts
         * are centred by that box width.  The only correction is for ink that
         * leans outside the box -- an italic x's ink starts left of its
         * origin -- so the rule is then extended just far enough to cover it,
         * rather than padded symmetrically (which widened \frac{a}{b} past
         * TeX). */
        const double boxW = std::max(num.w, den.w);
        const double numX = (boxW - num.w) / 2.0;
        const double denX = (boxW - den.w) / 2.0;
        const double inkLeft = std::min({0.0,
            numX + layout_ink_left(num), denX + layout_ink_left(den)});
        const double inkRight = std::max({boxW,
            numX + layout_ink_right(num), denX + layout_ink_right(den)});

        /* A layout must start at x >= 0, or it misaligns when something wraps
         * it -- a radical drew its vinculum in two pieces over a fraction
         * whose italic ink leaned to negative x.  Shift everything right so
         * the leftmost ink sits at 0; for upright content the shift is 0. */
        const double shift = -inkLeft;
        Layout out;
        out.w = inkRight + shift;
        out.absorb(num, numX + shift, -up);
        out.absorb(den, denX + shift, down);
        Rule bar;
        bar.x = inkLeft + shift; bar.y = -axis - thick / 2.0;
        bar.w = inkRight - inkLeft; bar.h = thick;
        out.rules.push_back(bar);
        out.asc = up + num.asc;
        out.desc = down + den.desc;
        return out;
    }

    Layout layout_sqrt(const SqrtNode& s, double sizePt) {
        Layout inner = layout_list(s.content, sizePt);
        const double thick = rule_thickness(sizePt);
        /* The vinculum floats a fixed clearance above the content, not at the
         * glyph's own top -- placing it at the glyph top stood \sqrt{x}'s bar
         * 0.135 em too high, because Cambria's radical is 0.96 em tall and
         * small content does not need all of it.  The clearance is measured
         * from a display \sqrt in pdfLaTeX: the bar sits 0.344 em above the
         * content box.  That is larger than the textbook rule+x_height/4
         * because it folds in the radical's own characteristic space, and it
         * is a constant gap, which is how TeX floats the rule above the
         * tallest content regardless of its height. */
        /* tex.web 737.  The surd is asked for content + clearance + rule, the
         * clearance being one rule thickness plus a quarter of the maths
         * x-height in display style; whatever the chosen size overshoots by
         * is then split above and below.  A flat 0.344 em was used for both
         * the clearance and the sizing, which is nearly twice TeX's clearance
         * -- harmless while the glyph was simply stretched to fit, but once
         * designed sizes are chosen it asked for one size too large and the
         * sign swallowed the equation. */
        const double theta = rule_thickness(sizePt);
        const double xHeight = 0.430 * sizePt;
        double clearance = theta + xHeight / 4.0;
        const double content = inner.asc + inner.desc;
        const double wantHeight = content + clearance + theta;

        /* Take the font's designed larger radical rather than stretching the
         * base one, exactly as fences and big operators do.  Stretching left
         * the sign thin and rubbery over tall content -- a radical over a
         * fraction came out 30% off TeX's proportions -- and thickened its
         * strokes with the height.  Latin Modern carries five designed sizes.
         * The largest can still fall short of a very deep box, and that one
         * is stretched the rest of the way. */
        Layout sign;
        unsigned short signVariant = 0;
        {
            const double want = wantHeight;
            double variantHeight = 0;
            unsigned short variant = 0;
            if (glyph_size_variant(0x221A, want / sizePt, variant,
                                   variantHeight) && variantHeight > 1e-6) {
                double vAsc = 0, vDesc = 0, vInk = 0;
                glyph_index_metrics(variant, sizePt, vAsc, vDesc, vInk);
                Glyph glyph;
                glyph.x = 0; glyph.y = 0; glyph.size = sizePt;
                glyph.symbol = true;
                glyph.glyphIndex = variant;
                glyph.outline = glyph_index_outline(variant);
                glyph.text = utf8_of(0x221A);
                sign.glyphs.push_back(glyph);
                sign.asc = vAsc; sign.desc = vDesc;
                sign.w = std::max(vInk, char_width(0x221A, sizePt, false, true));
                signVariant = variant;
            } else {
                sign = glyph_layout(0x221A, sizePt, false, true);
            }
        }
        /* The sign hangs from the bar: its flag is at barTop and it reaches
         * down toward the content's foot.  For content that fits inside the
         * natural glyph the sign is only shifted down, exactly as TeX uses
         * the base radical and lets it dip below the baseline; only taller
         * content stretches it (a size-variant glyph would be sharper here,
         * noted in tests/test_tex_geometry.py). */
        const double ascNat = sign.asc;
        const double descNat = sign.desc;
        const double naturalTotal = std::max(ascNat + descNat, 1e-6);
        /* Keep the chosen size and split its surplus, as TeX does.  Easing
         * the glyph down to exactly `wantHeight` was tried and is wrong: the
         * surplus is what lifts the bar clear of the content, and removing it
         * squashed the sign to 0.609 em where TeX has 0.806.
         *
         * For heights between two designed sizes TeX assembles the radical
         * from extensible pieces and we take the next size up, so a wide but
         * shallow radical still comes out one size tall.  That is the
         * remaining difference tools/tex_sweep.py reports for it. */
        double stretch = 1.0;
        if (wantHeight > naturalTotal) stretch = wantHeight / naturalTotal;
        if (stretch > 1.0) {
            for (auto& glyph : sign.glyphs) glyph.stretchY = stretch;
        }
        /* Split whatever the chosen size overshoots by, so the bar rises with
         * the sign instead of the content hanging off its bottom. */
        const double signHeight = naturalTotal * stretch;
        const double surplus = signHeight - (content + clearance);
        if (surplus > 0) clearance += surplus / 2.0;
        const double barTop = inner.asc + clearance;
        /* Place the glyph baseline so its stretched flag lands exactly on the
         * bar.  The renderer scales the glyph about its baseline, so the flag
         * (ascNat above the baseline) draws ascNat*stretch above it; for the
         * flag to sit at -barTop the baseline goes to (ascNat*stretch - barTop)
         * below the equation baseline.  The previous form negated this, which
         * for a large stretch -- \sqrt over a fraction -- drove the flag above
         * the box and clipped it into two pieces. */
        const double signBaseline = ascNat * stretch - barTop;
        const double signBottom = signBaseline + descNat * stretch;

        /* The bar starts where the radical's flag ends, but its THICKNESS is
         * TeX's default_rule_thickness -- a constant.  It used to be the
         * flag's own thickness times the vertical stretch, which made it
         * 0.065 em for \sqrt{x} against TeX's 0.040, and 0.085 em over taller
         * content: a bar that grew as the content grew, which TeX never does. */
        const double barThick = thick;
        double signW = sign.w;
        double flagThick = 0, flagLeft = 0, flagRight = 0;
        /* Measure the glyph actually being drawn.  A designed variant has a
         * different flag from the base radical, so measuring the character
         * would join the bar to the wrong place. */
        const bool haveFlag = signVariant
            ? glyph_top_plateau(signVariant, sizePt, false, true, flagThick,
                                flagLeft, flagRight, /*byIndex=*/true)
            : glyph_top_plateau(0x221A, sizePt, false, true, flagThick,
                                flagLeft, flagRight);
        if (haveFlag) {
            /* Only clip when the font's own flag is thicker than the rule.
             *
             * Latin Modern's RadicalRuleThickness is 0.040 em -- the same
             * default_rule_thickness the bar uses -- so its flag continues
             * the bar exactly and the whole glyph should be drawn, shoulder
             * and all, which is what Computer Modern looks like.  Cambria's
             * flag is 0.065 em, and drawing that whole left a short deeper
             * stub under the left end of the bar that reads as a chipped
             * bar; there the glyph is clipped where its shoulder stops
             * tapering and the bar starts from that point instead. */
            if (flagThick > barThick * 1.15) {
                signW = flagLeft;
                for (auto& glyph : sign.glyphs) glyph.clipRight = flagLeft;
            } else {
                signW = flagRight;
            }
        } else {
            signW = std::max(sign.w,
                             glyph_ink_width(0x221A, sizePt, false, true));
        }

        Layout out;
        out.absorb(sign, 0, signBaseline);  /* flag now rests on the bar */
        out.absorb(inner, signW, 0);        /* content keeps its own baseline */
        out.asc = barTop;
        out.desc = std::max(signBottom, inner.desc);
        if (s.hasIndex) {
            Layout idx = layout_list(s.index, script_size(sizePt));
            double ix = std::max(0.0, signW * 0.55 - idx.w);
            double iy = -(barTop * 0.55 + idx.desc);
            out.absorb(idx, ix, iy);
            out.asc = std::max(out.asc, -iy + idx.asc);
        }
        Rule bar;
        bar.x = signW;
        bar.y = -barTop;
        /* A hair of overhang past the content keeps the right end from
         * looking clipped, as TeX's \sqrt does. */
        bar.w = inner.w + 0.08 * sizePt;
        bar.h = barThick;
        out.rules.push_back(bar);
        out.w = signW + bar.w;
        return out;
    }

    Layout layout_matrix(const MatrixNode& m, double sizePt) {
        const int rows = std::max(1, m.rows);
        const int cols = std::max(1, m.cols);
        std::vector<Layout> cells(size_t(rows * cols));
        std::vector<double> cw(size_t(cols), 0.0);
        std::vector<double> ra(size_t(rows), 0.0), rd(size_t(rows), 0.0);
        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                size_t k = size_t(r * cols + c);
                if (k < m.elements.size() && m.elements[k] &&
                    m.elements[k]->tag() == Node::kLine) {
                    const auto& line = static_cast<const LineNode&>(*m.elements[k]);
                    cells[k] = layout_list(line.children, sizePt);
                } else {
                    static const NodeList empty;
                    cells[k] = layout_list(empty, sizePt);
                }
                cw[size_t(c)] = std::max(cw[size_t(c)], cells[k].w);
                ra[size_t(r)] = std::max(ra[size_t(r)], cells[k].asc);
                rd[size_t(r)] = std::max(rd[size_t(r)], cells[k].desc);
            }
        }
        const double colGap = 0.75 * sizePt;
        const double rowGap = 0.35 * sizePt;
        double totalH = 0;
        for (int r = 0; r < rows; ++r) totalH += ra[size_t(r)] + rd[size_t(r)];
        totalH += rowGap * (rows - 1);
        double top = -0.5 * totalH + 0.12 * sizePt;
        double y = top;

        Layout out;
        for (int r = 0; r < rows; ++r) {
            double baseline = y + ra[size_t(r)];
            double x = 0;
            for (int c = 0; c < cols; ++c) {
                size_t k = size_t(r * cols + c);
                double inset = (cw[size_t(c)] - cells[k].w) * 0.5;
                if (m.layoutKind == MatrixNode::kAlignedLayout)
                    inset = (c % 2 == 0) ? (cw[size_t(c)] - cells[k].w) : 0.0;
                else if (m.layoutKind == MatrixNode::kCasesLayout)
                    inset = 0.0;        /* cases left-aligns every column */
                out.absorb(cells[k], x + inset, baseline);
                x += cw[size_t(c)] + (c + 1 < cols ? colGap : 0);
            }
            out.w = std::max(out.w, x);
            y += ra[size_t(r)] + rd[size_t(r)] +
                 (r + 1 < rows ? rowGap : 0);
        }
        out.asc = std::max(0.0, -top);
        out.desc = std::max(0.0, top + totalH);
        if (m.layoutKind == MatrixNode::kAlignedLayout && cols > 1) {
            double boundary = 0;
            for (int c = 0; c + 1 < cols; ++c) {
                boundary += cw[size_t(c)];
                out.editMarks.push_back({EditMarkKind::Alignment,
                    boundary + 0.5 * colGap, top, top + totalH});
                boundary += colGap;
            }
        }
        return out;
    }

    Layout layout_pile(const PileNode& p, double sizePt) {
        std::vector<Layout> lines;
        for (const auto& n : p.lines) {
            if (n && n->tag() == Node::kLine)
                lines.push_back(layout_list(static_cast<const LineNode&>(*n).children,
                                            sizePt));
        }
        if (lines.empty()) {
            static const NodeList empty;
            return layout_list(empty, sizePt);
        }
        const double gap = 0.35 * sizePt;
        double width = 0, total = gap * (lines.size() - 1);
        for (const auto& l : lines) {
            width = std::max(width, l.w);
            total += l.asc + l.desc;
        }
        double top = -0.5 * total + 0.12 * sizePt;
        double y = top;
        Layout out;
        for (const auto& l : lines) {
            double baseline = y + l.asc;
            double x = (p.halign == 1) ? 0 : (width - l.w) * 0.5;
            out.absorb(l, x, baseline);
            y += l.asc + l.desc + gap;
        }
        out.w = width;
        out.asc = std::max(0.0, -top);
        out.desc = std::max(0.0, top + total);
        return out;
    }

    Layout layout_embell(const EmbellNode& e, double sizePt) {
        Layout inner = layout_list(e.content, sizePt);
        /* Accent characters, not their ASCII lookalikes.  U+005E is a full
         * circumflex punctuation mark: at the base's size it drew 1.24 times
         * the letter's width where TeX's hat is 0.58. */
        uint32_t cp = 0x02C6;                     /* hat */
        if (e.embellType == 2) cp = 0x02D9;       /* dot */
        else if (e.embellType == 8) cp = 0x02DC;  /* tilde */
        else if (e.embellType == 11) cp = 0x2192; /* vector */
        else if (e.embellType == 17) cp = 0x00AF; /* bar */
        /* Accents are drawn at the base's own size, not at script size.  TeX
         * takes them from the math font at the current size; shrinking them
         * made every accent about two thirds of TeX's width -- a tilde 0.49
         * of the letter's width where TeX draws 0.73.  The arrow of \vec is
         * a full-width glyph rather than an accent, so it keeps the smaller
         * size or it would overhang the letter badly. */
        const double markSize =
            (e.embellType == EM_RARROW || e.embellType == EM_LARROW ||
             e.embellType == EM_BARROW) ? script_size(sizePt) : sizePt;
        Layout mark = glyph_layout(cp, markSize, false, needs_math_face(cp));
        Layout out = inner;
        double x = (inner.w - mark.w) * 0.5;
        /* Sit the accent's INK on the base, not its baseline.  U+02DC and the
         * other combining-accent characters are drawn entirely above their
         * baseline, so `desc` is zero and placing them by it left the whole
         * of that rise as empty space: the tilde floated 1.1 base-heights
         * above the letter where TeX puts it at 0.30. */
        const double inkBottom =
            glyph_ink_bottom(cp, markSize, false, needs_math_face(cp));
        double y = -(inner.asc + 0.10 * sizePt) + inkBottom;
        out.absorb(mark, x, y);
        out.asc = std::max(out.asc, -y + mark.asc);
        out.w = std::max(out.w, mark.w);
        return out;
    }

    Layout layout_decoration(const DecorationNode& d, double sizePt) {
        Layout out = layout_list(d.content, sizePt);
        const double thick = rule_thickness(sizePt);
        Rule bar;
        bar.x = 0;
        bar.w = out.w;
        bar.h = thick;
        if (d.selector == tmUBAR) {
            bar.y = out.desc + 0.10 * sizePt;
            out.desc = bar.y + thick;
        } else {
            bar.y = -(out.asc + 0.10 * sizePt + thick);
            out.asc = -bar.y;
        }
        out.rules.push_back(bar);
        return out;
    }

    Layout layout_brace(const BraceDecoNode& b, double sizePt) {
        Layout body = layout_list(b.content, sizePt);
        Layout label = layout_list(b.label, script_size(sizePt));
        Layout out = body;
        const bool over = b.selector == tmUHBRACE;
        uint32_t cp = over ? 0x23DE : 0x23DF;
        Layout brace = glyph_layout(cp, sizePt, false, true);
        double stretch = std::max(1.0, body.w / std::max(brace.w, 1e-6));
        /* Horizontal growth is represented by a rule when a single brace
         * glyph would become unreadably distorted. */
        double by = over ? -(body.asc + 0.18 * sizePt)
                         :  (body.desc + 0.60 * sizePt);
        if (stretch < 1.8) {
            out.absorb(brace, (body.w - brace.w) * 0.5, by);
        } else {
            Rule r{0, over ? by - 0.15 * sizePt : by,
                   body.w, rule_thickness(sizePt)};
            out.rules.push_back(r);
        }
        if (!label.glyphs.empty() || !label.placeholders.empty()) {
            double ly = over ? by - 0.35 * sizePt - label.desc
                             : by + 0.35 * sizePt + label.asc;
            out.absorb(label, (body.w - label.w) * 0.5, ly);
            if (over) out.asc = std::max(out.asc, -ly + label.asc);
            else out.desc = std::max(out.desc, ly + label.desc);
        }
        out.w = std::max(body.w, label.w);
        return out;
    }

    Layout layout_dirac(const DiracNode& d, double sizePt) {
        Layout bra = layout_list(d.bra, sizePt);
        Layout ket = layout_list(d.ket, sizePt);
        /* Grow the brackets and the separator to the taller side, the way
         * \left ... \middle ... \right does in TeX.  A bra-ket whose bar is
         * one line tall around a fraction is the tell that it was assembled
         * from plain characters. */
        const double need = std::max(bra.asc + bra.desc, ket.asc + ket.desc);
        /* The angle bracket's own ink, for the same reason as in
         * layout_fence: a math font's global ascent and descent reserve room
         * for extensible signs, so measuring against them made every bracket
         * believe it was already tall enough. */
        double plainAsc = 0, plainDesc = 0;
        glyph_vmetrics(0x27E8, sizePt, false, true, false,
                       plainAsc, plainDesc);
        const double stretch =
            std::max(1.0, need / std::max(plainAsc + plainDesc, 1e-6));

        Layout l = glyph_layout(0x27E8, sizePt, false, true);
        Layout mid = glyph_layout('|', sizePt, false, false);
        Layout r = glyph_layout(0x27E9, sizePt, false, true);
        for (Layout* part : {&l, &mid, &r}) {
            for (auto& gg : part->glyphs) gg.stretchY = stretch;
            part->asc *= stretch;
            part->desc *= stretch;
        }
        Layout out;
        double x = 0;
        for (Layout* part : std::vector<Layout*>{&l, &bra, &mid, &ket, &r}) {
            out.absorb(*part, x, 0);
            x += part->w;
            out.asc = std::max(out.asc, part->asc);
            out.desc = std::max(out.desc, part->desc);
        }
        out.w = x;
        return out;
    }

    Layout layout_lim(const LimNode& l, double sizePt) {
        Layout word = text_layout("lim", sizePt, false);
        Layout sub = layout_list(l.content, script_size(sizePt));
        Layout out = word;
        double x = (word.w - sub.w) * 0.5;
        double y = word.desc + 0.15 * sizePt + sub.asc;
        out.absorb(sub, x, y);
        out.desc = std::max(out.desc, y + sub.desc);
        out.w = std::max(word.w, sub.w);
        return out;
    }

    Layout layout_overset(const OversetNode& o, double sizePt) {
        Layout base = layout_list(o.base, sizePt);
        Layout over = layout_list(o.over, script_size(sizePt));
        Layout out = base;
        double x = (base.w - over.w) * 0.5;
        if (o.under) {
            double y = base.desc + 0.12 * sizePt + over.asc;
            out.absorb(over, x, y);
            out.desc = std::max(out.desc, y + over.desc);
        } else {
            double y = -(base.asc + 0.12 * sizePt + over.desc);
            out.absorb(over, x, y);
            out.asc = std::max(out.asc, -y + over.asc);
        }
        out.w = std::max(base.w, over.w);
        return out;
    }

    Layout layout_bigop(uint32_t glyph, const NodeList& body,
                        const NodeList& lower, const NodeList& upper,
                        bool hasLower, bool hasUpper, bool stacked, double sizePt) {
        /* Display operators are sized to the ink height TeX gives them, not
         * to one nominal point size.  Measured from pdfLaTeX at 12 pt, the
         * integral family stands 2.22 em tall and the sum-like ones 1.40 em.
         * Cambria Math's integral is proportionally shorter than cmex's, so a
         * single nominal size served neither: the integral came out 27 % too
         * short, which is why its correctly-spaced limits sat above and below
         * a sign too small to reach them. */
        const bool integralFamily = glyph >= 0x222B && glyph <= 0x2230;
        /* Rasterised from a 12 pt pdfLaTeX \displaystyle operator and
         * measured as ink on both sides -- comparing a PDF character box
         * against rendered ink is what once made the integral look too big
         * when it was in fact 27% too small. */
        const double targetInk = (integralFamily ? 2.215 : 1.395) * sizePt;
        double opSize = st_.sym * (sizePt / std::max(st_.full, 1e-6));
        {
            double signAsc = 0, signDesc = 0;
            glyph_vmetrics(glyph, opSize, false, true, false,
                           signAsc, signDesc);
            const double ink = signAsc + signDesc;
            if (ink > 1e-6) opSize *= targetInk / ink;
        }
        /* Prefer the font's designed larger glyph over the base one scaled up.
         * Scaling widens a sign exactly as much as it heightens it: Latin
         * Modern's base integral is 0.498 wide per unit of height where its
         * designed display integral (integral.v1) is 0.399 -- the same 0.400
         * pdfLaTeX produces.  With the scaled one the curls reached into the
         * limits, which is what a reader sees as the limits overlapping. */
        Layout op;
        unsigned short variant = 0;
        double variantHeight = 0;
        if (glyph_size_variant(glyph, targetInk / sizePt, variant,
                               variantHeight) && variantHeight > 1e-6) {
            double vAsc = 0, vDesc = 0, vInk = 0;
            glyph_index_metrics(variant, sizePt, vAsc, vDesc, vInk);
            const double natural = vAsc + vDesc;
            const double size = natural > 1e-6 ? sizePt * targetInk / natural
                                               : sizePt;
            glyph_index_metrics(variant, size, vAsc, vDesc, vInk);
            Glyph g;
            g.x = 0; g.y = 0; g.size = size; g.symbol = true;
            g.glyphIndex = variant;
            g.outline = glyph_index_outline(variant);
            g.text = utf8_of(glyph);          /* a label, for tests */
            op.glyphs.push_back(g);
            op.asc = vAsc; op.desc = vDesc;
            op.w = vInk;
        } else {
            op = glyph_layout(glyph, opSize, false, true);
            op.w = std::max(op.w, glyph_ink_width(glyph, opSize, false, true));
        }
        /* TeX centres a display operator on the math axis, so the sign
         * straddles the baseline.  Setting it on the baseline with its own
         * ink metrics -- which is what happened -- stood the integral about
         * 0.24 em too high: measured as rendered ink, 72.6% of the sign sat
         * above the baseline where TeX puts 61.5%.  That is why its top curl
         * ran into the upper limit even though the limit was in the right
         * place.  Centring on the axis predicts 61.3%. */
        const double lift = kAxisHeight * sizePt - (op.asc - op.desc) / 2.0;
        const double opAsc = op.asc + lift;
        const double opDesc = op.desc - lift;
        double ss = script_size(sizePt);
        Layout out;
        double x = 0;

        if (stacked) {
            /* Limits above and below the operator, everything centred on the
             * widest of the three.  Centring the limits on the operator alone
             * lets a wide limit hang left of the origin and overlap whatever
             * precedes the operator, because the reported width never sees it. */
            Layout up, lo;
            if (hasUpper) up = layout_list(upper, ss);
            if (hasLower) lo = layout_list(lower, ss);
            double w = op.w;
            if (hasUpper) w = std::max(w, up.w);
            if (hasLower) w = std::max(w, lo.w);

            /* TeX's big_op_spacing leaves more room below the operator than
             * above it, and the flat 0.15 em used here put the limits 0.39 em
             * closer together than a display sum in LaTeX.  Measured from
             * pdfLaTeX at 12 pt: 0.30 em from the operator's top to the upper
             * limit, 0.42 em from its bottom to the lower limit. */
            /* Re-derived against a pdfLaTeX display \sum once the operator
             * ink height was measured rather than guessed: the pair has to
             * put the two limit baselines 2.391 em apart at 12 pt. */
            constexpr double kOpSpacingAbove = 0.289;
            constexpr double kOpSpacingBelow = 0.406;
            out.absorb(op, (w - op.w) / 2.0, -lift);
            out.asc = opAsc;
            out.desc = opDesc;
            if (hasUpper) {
                const double gap = kOpSpacingAbove * sizePt;
                out.absorb(up, (w - up.w) / 2.0, -(opAsc + gap + up.desc));
                out.asc = std::max(out.asc, opAsc + gap + up.desc + up.asc);
            }
            if (hasLower) {
                const double gap = kOpSpacingBelow * sizePt;
                out.absorb(lo, (w - lo.w) / 2.0, opDesc + gap + lo.asc);
                out.desc = std::max(out.desc, opDesc + gap + lo.asc + lo.desc);
            }
            x = w;
        } else {
            /* TeX places integral limits much farther from the baseline than
             * ordinary letter scripts, and applies the integral's italic
             * correction to the upper limit.  The ratios below come from a
             * 12 pt TeX reference: upper/main baselines -1.097 em,
             * lower/main +0.903 em, upper/lower x difference +0.444 em. */
            /* Both limits are placed from the operator's ORIGIN, using the
             * offsets a 12 pt pdfLaTeX \int_a^b actually uses: lower +0.554
             * em, upper +0.996 em, the 0.442 em between them being the
             * italic correction.
             *
             * They used to be placed at op.w, and op.w had been forced up to
             * the glyph's full ink width just above.  For a slanted sign the
             * ink is widest near the top, so measuring from there threw both
             * limits 0.65 em right of where TeX puts them -- the visible gap
             * between the integral and its limits.  TeX tucks the lower limit
             * under the bottom curl, where the sign's ink has already swung
             * left, which is why it can sit so close. */
            constexpr double kIntLowerX = 0.554;
            constexpr double kIntUpperX = 0.996;
            out.absorb(op, 0, -lift);
            out.asc = opAsc;
            out.desc = opDesc;
            /* Only the integral family slants; a \sum\nolimits has an upright
             * sign whose advance is the right place for its scripts. */
            const double lowerX = integralFamily ? kIntLowerX * sizePt : op.w;
            const double upperX = integralFamily ? kIntUpperX * sizePt
                                                 : op.w + 0.44 * sizePt;
            x = op.w;
            if (hasUpper) {
                Layout up = layout_list(upper, ss);
                const double upperShift = 1.10 * sizePt;
                out.absorb(up, upperX, -upperShift);
                out.asc = std::max(out.asc, upperShift + up.asc);
                x = std::max(x, upperX + up.w);
            }
            if (hasLower) {
                Layout lo = layout_list(lower, ss);
                const double lowerShift = 0.90 * sizePt;
                out.absorb(lo, lowerX, lowerShift);
                out.desc = std::max(out.desc, lowerShift + lo.desc);
                x = std::max(x, lowerX + lo.w);
            }
        }

        /* The operand joins the operator here rather than through the atom
         * loop, so the Op-to-whatever space has to be applied by hand -- a
         * sigma otherwise touches the symbol that follows it. */
        for (const auto& n : body) {
            if (!n || n->tag() == Node::kSize) continue;
            x += space_mu(kOp, class_of(*n)) * sizePt / 18.0;
            break;
        }

        Layout bodyL = layout_list(body, sizePt);
        out.absorb(bodyL, x, 0);
        out.asc = std::max(out.asc, bodyL.asc);
        out.desc = std::max(out.desc, bodyL.desc);
        out.w = x + bodyL.w;
        return out;
    }

    /* Unhandled templates still show their content rather than vanishing. */
    Layout layout_fallback(const Node& n, double sizePt) {
        switch (n.tag()) {
            case Node::kDecoration:
                return layout_list(static_cast<const DecorationNode&>(n).content, sizePt);
            case Node::kBraceDeco:
                return layout_list(static_cast<const BraceDecoNode&>(n).content, sizePt);
            case Node::kEmbell:
                return layout_list(static_cast<const EmbellNode&>(n).content, sizePt);
            default:
                return Layout();
        }
    }
};

}  // namespace

/* ------------------------------------------------------------------ */
/* Public API                                                          */
/* ------------------------------------------------------------------ */
std::string render_svg(const LineNode& root, const SvgStyle& style) {
    Renderer r(style);
    Layout L = r.run(root);

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
        if (g.glyphIndex && !g.outline.empty()) {
            /* No viewer can name a size variant, so its shape travels as a
             * path.  data-char names the character it stands for, which is
             * what the geometry tool matches on. */
            o << "  <path data-char=\"" << xml_escape(g.text)
              << "\" data-size=\"" << g.size
              << "\" transform=\"translate(" << (g.x + pad) << ','
              << (g.y + baseline) << ") scale(" << g.size << ")\" d=\""
              << g.outline << "\" fill=\"currentColor\"/>\n";
            continue;
        }
        const std::string& family = g.cjk ? style.cjk
                                  : g.symbol ? style.symbol : style.serif;
        o << "  <text x=\"" << (g.x + pad) << "\" y=\"" << (g.y + baseline)
          << "\" font-family=\"" << xml_escape(family) << "\""
          << " font-size=\"" << g.size << "\"";
        if (g.italic) o << " font-style=\"italic\"";
        if (g.bold) o << " font-weight=\"bold\"";
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

RenderMetrics measure_equation(const LineNode& root, const SvgStyle& style) {
    Renderer r(style);
    Layout L = r.run(root);
    RenderMetrics m;
    m.width = layout_ink_right(L) + 2 * style.padding;
    m.height = L.asc + L.desc + 2 * style.padding;
    m.baseline = style.padding + L.asc;
    return m;
}

bool hit_test_equation(const LineNode& root, double x, double y,
                       const SvgStyle& style,
                       const NodeList** slot, int* index) {
    if (!slot || !index) return false;
    Renderer r(style);
    Layout L = r.run(root);
    const double lx = x - style.padding;
    const double ly = y - (style.padding + L.asc);
    const CaretSite* best = nullptr;
    double bestScore = std::numeric_limits<double>::infinity();
    for (const auto& c : L.carets) {
        double dx = std::fabs(lx - c.x);
        double dy = 0;
        if (ly < c.top) dy = c.top - ly;
        else if (ly > c.bottom) dy = ly - c.bottom;
        double score = dx + 1.8 * dy;
        if (score < bestScore) { bestScore = score; best = &c; }
    }
    if (!best) return false;
    *slot = best->slot;
    *index = best->index;
    return true;
}

#ifdef _WIN32
namespace {

std::wstring wide_utf8(const std::string& s) {
    if (s.empty()) return std::wstring();
    int n = MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS,
                                s.data(), int(s.size()), nullptr, 0);
    if (n <= 0) return L"\xFFFD";
    std::wstring w(size_t(n), L'\0');
    MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS,
                        s.data(), int(s.size()), w.data(), n);
    return w;
}

std::wstring first_family(const std::string& families) {
    size_t comma = families.find(',');
    std::string one = families.substr(0, comma);
    while (!one.empty() && one.front() == ' ') one.erase(one.begin());
    while (!one.empty() && one.back() == ' ') one.pop_back();
    return wide_utf8(one);
}

/* Drawing fonts, kept between repaints.
 *
 * The draw loop used to create and destroy an HFONT for every glyph, so a
 * thirty-glyph equation cost thirty CreateFontW/DeleteObject pairs on every
 * keystroke -- and CreateFontW is among the most expensive GDI calls there
 * is.  Layout, by contrast, measures at 0.017 ms for a typical equation, so
 * this loop was the whole of the difference in responsiveness.  Eqnedit32
 * creates fonts at twelve places in the entire program; it does not make one
 * per character either. */
struct DrawFontKey {
    std::wstring face;
    int height = 0;
    bool bold = false;
    bool italic = false;
    bool operator<(const DrawFontKey& o) const {
        if (height != o.height) return height < o.height;
        if (bold != o.bold) return bold < o.bold;
        if (italic != o.italic) return italic < o.italic;
        return face < o.face;
    }
};

struct DrawFontCache {
    std::map<DrawFontKey, HFONT> fonts;
    ~DrawFontCache() {
        for (auto& kv : fonts) DeleteObject(kv.second);
    }
    bool enabled = true;
    HFONT get(const DrawFontKey& key) {
        if (!enabled)                       /* benchmark comparison path */
            return CreateFontW(key.height, 0, 0, 0,
                key.bold ? FW_BOLD : FW_NORMAL,
                key.italic ? TRUE : FALSE, FALSE, FALSE, DEFAULT_CHARSET,
                OUT_TT_PRECIS, CLIP_DEFAULT_PRECIS, CLEARTYPE_QUALITY,
                DEFAULT_PITCH | FF_DONTCARE, key.face.c_str());
        auto it = fonts.find(key);
        if (it != fonts.end()) return it->second;
        /* Zoom is continuous, so the key space is unbounded in principle.
         * Drop everything rather than grow without limit; in practice a
         * session settles on a handful of sizes. */
        if (fonts.size() > 96) {
            for (auto& kv : fonts) DeleteObject(kv.second);
            fonts.clear();
        }
        HFONT f = CreateFontW(key.height, 0, 0, 0,
            key.bold ? FW_BOLD : FW_NORMAL,
            key.italic ? TRUE : FALSE, FALSE, FALSE, DEFAULT_CHARSET,
            OUT_TT_PRECIS, CLIP_DEFAULT_PRECIS, CLEARTYPE_QUALITY,
            DEFAULT_PITCH | FF_DONTCARE, key.face.c_str());
        fonts[key] = f;
        return f;
    }
};

DrawFontCache& draw_fonts() {
    static DrawFontCache cache;
    return cache;
}

}  // namespace

void set_draw_font_cache_enabled(bool enabled) { draw_fonts().enabled = enabled; }

namespace {

const CaretSite* find_site(const Layout& l, const NodeList* slot, int index) {
    if (!slot || index < 0) return nullptr;
    for (const auto& c : l.carets)
        if (c.slot == slot && c.index == index) return &c;
    return nullptr;
}

}  // namespace

bool caret_geometry_equation(const LineNode& root, const NodeList* slot,
                             int index, const SvgStyle& style,
                             CaretGeometry* geometry) {
    if (!geometry) return false;
    Renderer renderer(style);
    const Layout layout = renderer.run(root);
    const CaretSite* caret = find_site(layout, slot, index);
    if (!caret) return false;
    const double baseline = style.padding + layout.asc;
    geometry->x = style.padding + caret->x;
    geometry->top = baseline + caret->top;
    geometry->bottom = baseline + caret->bottom;
    return true;
}

/* Draw one glyph as its outline, recording geometry rather than text.
 *
 * A metafile stores the calls it was given, so a TextOut in it is still a
 * TextOut when PowerPoint replays it -- and the font it names is loaded into
 * this process alone, so PowerPoint substitutes and the equation stops
 * looking like TeX.  Bracketing the TextOut in a path does not help for the
 * same reason.  Fetching the outline here and recording MoveTo/PolylineTo
 * inside the path puts pure geometry in the metafile: no font is consulted at
 * playback, and it stays vector.
 *
 * Returns false when the glyph has no outline (a space, or a face that will
 * not give one), so the caller can fall back to drawing it as text.
 */
bool emit_glyph_outline(HDC hdc, HFONT face, const Glyph& g,
                        int gx, int gy, double scale) {
    /* Always by glyph index.  A character code only reaches the BMP, and the
     * italic letters are all above it, so naming glyphs by character left
     * every variable out of the metafile entirely. */
    const UINT format = GGO_NATIVE | GGO_GLYPH_INDEX;
    const UINT which = g.glyphIndex ? UINT(g.glyphIndex)
                                    : UINT(glyph_id_of(first_codepoint(g.text),
                                                       g.cjk));
    if (!which || !face) return false;
    /* Ask a real DC for the outline.  The target here is a metafile DC, which
     * records rather than renders and does not answer glyph queries. */
    HDC probe = CreateCompatibleDC(nullptr);
    if (!probe) return false;
    HGDIOBJ oldFace = SelectObject(probe, face);
    GLYPHMETRICS gm = {};
    MAT2 identity = {{0, 1}, {0, 0}, {0, 0}, {0, 1}};
    const DWORD need = GetGlyphOutlineW(probe, which, format, &gm, 0, nullptr,
                                        &identity);
    std::vector<uint8_t> buffer;
    bool got = (need != GDI_ERROR);
    if (got && need > 0) {
        buffer.resize(need);
        got = GetGlyphOutlineW(probe, which, format, &gm, need, buffer.data(),
                               &identity) != GDI_ERROR;
    }
    SelectObject(probe, oldFace);
    DeleteDC(probe);
    if (!got) return false;
    if (buffer.empty()) return true;            /* blank, nothing to draw */

    auto fx = [](const FIXED& f) {
        return double(f.value) + double(f.fract) / 65536.0;
    };
    const double sy = g.stretchY;
    auto place = [&](double px, double py) {
        return POINT{LONG(std::lround(gx + px)),
                     LONG(std::lround(gy - py * sy))};
    };

    /* One PolyPolygon record rather than a path.  BeginPath/EndPath/FillPath
     * round-trips through GDI but GDI+ -- which is what Office and most
     * viewers replay a metafile with -- does not honour it reliably, and the
     * equation came back empty.  A filled poly-polygon is a single primitive
     * every reader understands. */
    std::vector<POINT> points;
    std::vector<INT> counts;
    size_t at = 0;
    while (at + sizeof(TTPOLYGONHEADER) <= buffer.size()) {
        const auto* head =
            reinterpret_cast<const TTPOLYGONHEADER*>(&buffer[at]);
        const size_t end = at + head->cb;
        if (head->cb == 0 || end > buffer.size()) break;
        double cx = fx(head->pfxStart.x), cy = fx(head->pfxStart.y);
        const size_t contourStart = points.size();
        points.push_back(place(cx, cy));
        std::vector<POINT>& run = points;
        size_t p = at + sizeof(TTPOLYGONHEADER);
        while (p + sizeof(TTPOLYCURVE) - sizeof(POINTFX) <= end) {
            const auto* curve = reinterpret_cast<const TTPOLYCURVE*>(&buffer[p]);
            const WORD n = curve->cpfx;
            if (curve->wType == TT_PRIM_LINE) {
                for (WORD i = 0; i < n; ++i) {
                    cx = fx(curve->apfx[i].x); cy = fx(curve->apfx[i].y);
                    run.push_back(place(cx, cy));
                }
            } else {
                /* Quadratic B-spline.  Two consecutive control points imply an
                 * on-curve midpoint between them.  Flattened rather than sent
                 * as beziers: the outline is already in device units here, and
                 * sixteen segments is finer than any printer resolution. */
                for (WORD i = 0; i + 1 < n; ++i) {
                    const double c1x = fx(curve->apfx[i].x);
                    const double c1y = fx(curve->apfx[i].y);
                    double ex = fx(curve->apfx[i + 1].x);
                    double ey = fx(curve->apfx[i + 1].y);
                    if (i + 2 < n) {
                        ex = (c1x + ex) / 2.0;
                        ey = (c1y + ey) / 2.0;
                    }
                    constexpr int kSteps = 16;
                    for (int s = 1; s <= kSteps; ++s) {
                        const double t = double(s) / kSteps, u = 1.0 - t;
                        run.push_back(place(
                            u * u * cx + 2 * u * t * c1x + t * t * ex,
                            u * u * cy + 2 * u * t * c1y + t * t * ey));
                    }
                    cx = ex; cy = ey;
                }
            }
            p += sizeof(TTPOLYCURVE) - sizeof(POINTFX) + sizeof(POINTFX) * n;
        }
        const size_t n = points.size() - contourStart;
        if (n >= 3) counts.push_back(INT(n));
        else points.resize(contourStart);
        at = end;
    }
    if (counts.empty()) return true;            /* nothing to draw */
    /* The brush matters: the default one is white, and filling white on white
     * reads exactly like the metafile being empty. */
    const int mode = SetPolyFillMode(hdc, WINDING);
    HBRUSH ink = CreateSolidBrush(RGB(20, 20, 20));
    HGDIOBJ oldBrush = SelectObject(hdc, ink);
    HGDIOBJ oldPen = SelectObject(hdc, GetStockObject(NULL_PEN));
    PolyPolygon(hdc, points.data(), counts.data(), INT(counts.size()));
    SelectObject(hdc, oldPen);
    SelectObject(hdc, oldBrush);
    DeleteObject(ink);
    SetPolyFillMode(hdc, mode);
    (void)scale;
    return true;
}

void draw_equation_gdi(const LineNode& root, HDC hdc,
                       double left, double top, double scale,
                       const SvgStyle& style,
                       const NodeList* caret_slot, int caret_index,
                       const NodeList* selection_slot,
                       int selection_first, int selection_last,
                       bool show_placeholders, bool show_caret,
                       bool text_as_outlines) {
    if (!hdc || scale <= 0) return;
    Renderer renderer(style);
    Layout L = renderer.run(root);
    const double baseline = style.padding + L.asc;
    const int saved = SaveDC(hdc);
    SetBkMode(hdc, TRANSPARENT);

    if (selection_slot && selection_first >= 0 && selection_last >= 0 &&
        selection_first != selection_last) {
        if (selection_last < selection_first) std::swap(selection_first, selection_last);
        const CaretSite* a = find_site(L, selection_slot, selection_first);
        const CaretSite* b = find_site(L, selection_slot, selection_last);
        if (a && b) {
            RECT rr;
            rr.left = LONG(std::floor(left + (std::min(a->x, b->x) + style.padding) * scale));
            rr.right = LONG(std::ceil(left + (std::max(a->x, b->x) + style.padding) * scale));
            rr.top = LONG(std::floor(top + (baseline + std::min(a->top, b->top)) * scale));
            rr.bottom = LONG(std::ceil(top + (baseline + std::max(a->bottom, b->bottom)) * scale));
            if (rr.right <= rr.left) rr.right = rr.left + 2;
            HBRUSH sel = CreateSolidBrush(RGB(204, 229, 255));
            FillRect(hdc, &rr, sel);
            DeleteObject(sel);
        }
    }

    if (show_placeholders) {
        HPEN pen = CreatePen(PS_DOT, 1, RGB(150, 150, 150));
        HGDIOBJ oldPen = SelectObject(hdc, pen);
        HGDIOBJ oldBrush = SelectObject(hdc, GetStockObject(NULL_BRUSH));
        for (const auto& p : L.placeholders) {
            Rectangle(hdc,
                int(std::lround(left + (p.x + style.padding) * scale)),
                int(std::lround(top + (p.y + baseline) * scale)),
                int(std::lround(left + (p.x + p.w + style.padding) * scale)),
                int(std::lround(top + (p.y + p.h + baseline) * scale)));
        }
        SelectObject(hdc, oldBrush);
        SelectObject(hdc, oldPen);
        DeleteObject(pen);

        HPEN guide = CreatePen(PS_DOT, 1, RGB(75, 135, 155));
        oldPen = SelectObject(hdc, guide);
        for (const auto& mark : L.editMarks) {
            if (mark.kind != EditMarkKind::Alignment) continue;
            const int x = int(std::lround(
                left + (mark.x + style.padding) * scale));
            MoveToEx(hdc, x, int(std::lround(
                top + (mark.top + baseline) * scale)), nullptr);
            LineTo(hdc, x, int(std::lround(
                top + (mark.bottom + baseline) * scale)));
        }
        SelectObject(hdc, oldPen);
        DeleteObject(guide);

        HBRUSH dot = CreateSolidBrush(RGB(75, 135, 155));
        HGDIOBJ oldDot = SelectObject(hdc, dot);
        oldPen = SelectObject(hdc, GetStockObject(NULL_PEN));
        for (const auto& mark : L.editMarks) {
            if (mark.kind != EditMarkKind::Space) continue;
            const int x = int(std::lround(
                left + (mark.x + style.padding) * scale));
            const int y = int(std::lround(
                top + (0.5 * (mark.top + mark.bottom) + baseline) * scale));
            const int radius = std::max(1, int(std::lround(scale)));
            Ellipse(hdc, x - radius, y - radius,
                    x + radius + 1, y + radius + 1);
        }
        SelectObject(hdc, oldPen);
        SelectObject(hdc, oldDot);
        DeleteObject(dot);
    }

    HBRUSH ink = CreateSolidBrush(RGB(20, 20, 20));
    for (const auto& rule : L.rules) {
        RECT rr{
            LONG(std::floor(left + (rule.x + style.padding) * scale)),
            LONG(std::floor(top + (rule.y + baseline) * scale)),
            LONG(std::ceil(left + (rule.x + rule.w + style.padding) * scale)),
            LONG(std::ceil(top + (rule.y + rule.h + baseline) * scale))};
        FillRect(hdc, &rr, ink);
    }
    DeleteObject(ink);

    SetTextColor(hdc, RGB(20, 20, 20));
    SetTextAlign(hdc, TA_LEFT | TA_BASELINE | TA_NOUPDATECP);
    /* The two faces are fixed for the whole run, so resolve them once instead
     * of parsing the family list per glyph. */
    const std::wstring symbolFace = first_family(style.symbol);
    const std::wstring serifFace = first_family(style.serif);
    const std::wstring cjkFace = cjk_face_name();
    HGDIOBJ previousFont = nullptr;
    HFONT selected = nullptr;
    std::wstring text;
    for (const auto& g : L.glyphs) {
        DrawFontKey key;
        key.face = g.cjk ? cjkFace : g.symbol ? symbolFace : serifFace;
        key.height = -std::max(1, int(std::lround(g.size * scale)));
        key.bold = g.bold;
        key.italic = g.italic;
        HFONT font = draw_fonts().get(key);
        const bool ownFont = !draw_fonts().enabled;
        if (ownFont || font != selected) {  /* consecutive glyphs share it */
            HGDIOBJ old = SelectObject(hdc, font);
            if (!previousFont) previousFont = old;
            selected = font;
        }
        const int gx = int(std::lround(left + (g.x + style.padding) * scale));
        const int gy = int(std::lround(top + (g.y + baseline) * scale));
        int glyphSaved = 0;
        const bool clipped = g.clipRight > 0.0;
        if (clipped || std::fabs(g.stretchY - 1.0) > 1e-6) {
            glyphSaved = SaveDC(hdc);
        }
        /* Clip before the world transform, so the boundary is in device
         * space and does not move when the glyph is stretched. */
        if (clipped) {
            const double boxH = (L.asc + L.desc + 2 * style.padding) * scale;
            IntersectClipRect(hdc,
                LONG(std::floor(left + (g.x + style.padding) * scale)) - 1,
                LONG(std::floor(top)) - 1,
                LONG(std::ceil(left + (g.clipRight + style.padding) * scale)),
                LONG(std::ceil(top + boxH)) + 1);
        }
        if (std::fabs(g.stretchY - 1.0) > 1e-6) {
            SetGraphicsMode(hdc, GM_ADVANCED);
            XFORM xf{1.0f, 0.0f, 0.0f, FLOAT(g.stretchY),
                     0.0f, FLOAT((1.0 - g.stretchY) * gy)};
            SetWorldTransform(hdc, &xf);
        }
        if (text_as_outlines &&
            emit_glyph_outline(hdc, font, g, gx, gy, scale)) {
            /* drawn as geometry */
        } else if (g.glyphIndex) {
            /* A designed size variant has no character; it is drawn by index. */
            const WORD index = g.glyphIndex;
            ExtTextOutW(hdc, gx, gy, ETO_GLYPH_INDEX, nullptr,
                        reinterpret_cast<LPCWSTR>(&index), 1, nullptr);
        } else {
            text = wide_utf8(g.text);
            TextOutW(hdc, gx, gy, text.data(), int(text.size()));
        }
        if (glyphSaved) {
            RestoreDC(hdc, glyphSaved);
            selected = nullptr;          /* RestoreDC put the old font back */
        }
        if (ownFont) {
            SelectObject(hdc, previousFont);
            selected = nullptr;
            DeleteObject(font);
        }
    }
    if (previousFont) SelectObject(hdc, previousFont);

    if (show_caret) {
        const CaretSite* c = find_site(L, caret_slot, caret_index);
        if (c) {
            int x = int(std::lround(left + (c->x + style.padding) * scale));
            int y1 = int(std::lround(top + (baseline + c->top) * scale));
            int y2 = int(std::lround(top + (baseline + c->bottom) * scale));
            HPEN pen = CreatePen(PS_SOLID, std::max(1, int(std::lround(scale))),
                                 RGB(0, 85, 170));
            HGDIOBJ old = SelectObject(hdc, pen);
            MoveToEx(hdc, x, y1, nullptr);
            LineTo(hdc, x, y2);
            SelectObject(hdc, old);
            DeleteObject(pen);
        }
    }
    RestoreDC(hdc, saved);
}
#endif

/* Did Latin Modern Math actually load?  When it does not, GDI substitutes
 * silently and every measurement shifts by a little -- which reads as a
 * flaky test rather than as a missing font.  Checks assert this so the
 * failure names itself instead of looking like noise. */
bool math_font_loaded() {
    ensure_math_font_public();
    return math_font_is_loaded();
}

void ensure_math_font_ready() { ensure_math_font_public(); }
std::string tex_to_svg(const std::string& latex, const SvgStyle& style) {
    std::unique_ptr<LineNode> root = parse_latex(latex);
    if (!root) return std::string();
    return render_svg(*root, style);
}

}  // namespace eqnedit
