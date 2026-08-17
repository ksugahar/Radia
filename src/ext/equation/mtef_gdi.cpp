/*
 * mtef_gdi.cpp -- draw a laid-out equation with GDI
 *
 * The drawing routine takes a device context and knows nothing about what it
 * is: a metafile DC records it, a memory DC rasterises it, and a window DC
 * would put it on screen.  That is the whole reason to draw with GDI rather
 * than emit a picture format directly.
 */
#include "mtef_gdi.h"
#include "tex_parser.h"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>
#include <objidl.h>
#include <gdiplus.h>
#endif

namespace mtef {
namespace {

#ifdef _WIN32

/* GDI works in device units; the layout is in points.  One unit per 0.01 mm is
 * what an enhanced metafile's frame wants, and drawing at a fixed multiple of a
 * point keeps the arithmetic exact enough that a rule lands on a pixel. */
constexpr double kUnitsPerPt = 20.0;      /* 20 device units to the point */

std::wstring widen(const std::string& utf8) {
    if (utf8.empty()) return std::wstring();
    int n = MultiByteToWideChar(CP_UTF8, 0, utf8.data(), int(utf8.size()),
                                nullptr, 0);
    std::wstring w(size_t(n), L'\0');
    MultiByteToWideChar(CP_UTF8, 0, utf8.data(), int(utf8.size()), &w[0], n);
    return w;
}

HFONT make_font(double sizePt, bool italic, bool symbol, bool cjk,
                double units_per_pt) {
    LOGFONTW lf = {};
    lf.lfHeight = -LONG(std::lround(sizePt * units_per_pt));
    lf.lfItalic = italic ? TRUE : FALSE;
    /* Never SYMBOL_CHARSET: Cambria Math is a Unicode font, and the legacy
     * Symbol code page mangles ASCII and loses anything past it. */
    lf.lfCharSet = DEFAULT_CHARSET;
    lf.lfQuality = ANTIALIASED_QUALITY;
    /* The same three faces the metric layer measures with, or the drawing
     * would not match the widths the layout was built from. */
    wcscpy_s(lf.lfFaceName, cjk    ? L"Yu Mincho"
                          : symbol ? L"Cambria Math"
                                   : L"Times New Roman");
    return CreateFontIndirectW(&lf);
}

/* The metafile and bitmap writers both put the equation one padding in from
 * the top left; a window places it wherever it likes. */
void draw(HDC hdc, const Layout& L, const SvgStyle& style, double scale) {
    const double upp = kUnitsPerPt * scale;
    const double pad = style.padding * upp;
    draw_layout(hdc, L, style, upp,
                int(std::lround(pad)),
                int(std::lround(pad + L.asc * upp)), RGB(0, 0, 0));
}

int gdiplus_encoder(const wchar_t* mime, CLSID* out) {
    UINT num = 0, size = 0;
    Gdiplus::GetImageEncodersSize(&num, &size);
    if (!size) return 0;
    std::vector<char> buf(size);
    auto* info = reinterpret_cast<Gdiplus::ImageCodecInfo*>(buf.data());
    Gdiplus::GetImageEncoders(num, size, info);
    for (UINT i = 0; i < num; ++i)
        if (wcscmp(info[i].MimeType, mime) == 0) { *out = info[i].Clsid; return 1; }
    return 0;
}

/* GDI+ has to be started once per process; the equation writers are called
 * from a Python process that may or may not already have done so. */
struct GdiPlusOnce {
    ULONG_PTR token = 0;
    GdiPlusOnce() {
        Gdiplus::GdiplusStartupInput in;
        Gdiplus::GdiplusStartup(&token, &in, nullptr);
    }
    ~GdiPlusOnce() { if (token) Gdiplus::GdiplusShutdown(token); }
};

#endif  /* _WIN32 */

}  // namespace

#ifdef _WIN32

void draw_layout(HDC hdc, const Layout& L, const SvgStyle& style,
                 double units_per_pt, int originX, int originY,
                 COLORREF colour) {
    (void)style;
    SetBkMode(hdc, TRANSPARENT);
    SetTextColor(hdc, colour);
    SetTextAlign(hdc, TA_LEFT | TA_BASELINE);

    HBRUSH brush = CreateSolidBrush(colour);
    for (const Rule& r : L.rules) {
        RECT rc;
        rc.left   = originX + LONG(std::lround(r.x * units_per_pt));
        rc.top    = originY + LONG(std::lround(r.y * units_per_pt));
        rc.right  = originX + LONG(std::lround((r.x + r.w) * units_per_pt));
        rc.bottom = originY + LONG(std::lround((r.y + r.h) * units_per_pt));
        if (rc.bottom <= rc.top) rc.bottom = rc.top + 1;   /* never vanish */
        FillRect(hdc, &rc, brush);
    }
    DeleteObject(brush);

    for (const Glyph& g : L.glyphs) {
        /* A stretched fence is the one thing GDI cannot do with a font size
         * alone, so it gets a taller font rather than a transform -- keeping
         * the metafile free of world transforms that some readers ignore. */
        HFONT f = make_font(g.size * (g.stretchY > 1.0 ? g.stretchY : 1.0),
                            g.italic, g.symbol, g.cjk, units_per_pt);
        HGDIOBJ old = SelectObject(hdc, f);
        std::wstring w = widen(g.text);
        TextOutW(hdc,
                 originX + int(std::lround(g.x * units_per_pt)),
                 originY + int(std::lround(g.y * units_per_pt)),
                 w.c_str(), int(w.size()));
        SelectObject(hdc, old);
        DeleteObject(f);
    }
}

std::string render_emf(const Layout& layout, const SvgStyle& style) {
    const double w_pt = layout.w + 2 * style.padding;
    const double h_pt = layout.asc + layout.desc + 2 * style.padding;
    if (w_pt <= 0 || h_pt <= 0) return std::string();

    /* The frame is in 0.01 mm; a point is 25.4/72 mm. */
    const double mm100_per_pt = 2540.0 / 72.0;
    RECT frame = {0, 0,
                  LONG(std::lround(w_pt * mm100_per_pt)),
                  LONG(std::lround(h_pt * mm100_per_pt))};

    HDC ref = GetDC(nullptr);
    HDC meta = CreateEnhMetaFileW(ref, nullptr, &frame,
                                  L"radia.equation\0equation\0\0");
    ReleaseDC(nullptr, ref);
    if (!meta) return std::string();

    /* Map the device units the drawing uses onto the frame. */
    SetMapMode(meta, MM_ANISOTROPIC);
    SetWindowExtEx(meta, int(std::lround(w_pt * kUnitsPerPt)),
                         int(std::lround(h_pt * kUnitsPerPt)), nullptr);
    SetViewportExtEx(meta, frame.right, frame.bottom, nullptr);

    draw(meta, layout, style, 1.0);

    HENHMETAFILE emf = CloseEnhMetaFile(meta);
    if (!emf) return std::string();
    UINT n = GetEnhMetaFileBits(emf, 0, nullptr);
    std::string out(n, '\0');
    if (n) GetEnhMetaFileBits(emf, n, reinterpret_cast<LPBYTE>(&out[0]));
    DeleteEnhMetaFile(emf);
    return out;
}

/* One rasterisation, two payloads: the PNG an application asks for by name and
 * the packed DIB it reads when it pastes a picture. */
static std::string rasterize(const Layout& layout, const SvgStyle& style,
                             double scale, bool as_dib) {
    const double w_pt = layout.w + 2 * style.padding;
    const double h_pt = layout.asc + layout.desc + 2 * style.padding;
    if (w_pt <= 0 || h_pt <= 0) return std::string();

    const int W = std::max(1, int(std::lround(w_pt * kUnitsPerPt * scale / kUnitsPerPt)));
    const int H = std::max(1, int(std::lround(h_pt * kUnitsPerPt * scale / kUnitsPerPt)));

    GdiPlusOnce gdip;

    HDC screen = GetDC(nullptr);
    HDC mem = CreateCompatibleDC(screen);
    BITMAPINFO bi = {};
    bi.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
    bi.bmiHeader.biWidth = W;
    bi.bmiHeader.biHeight = -H;              /* top-down */
    bi.bmiHeader.biPlanes = 1;
    bi.bmiHeader.biBitCount = 32;
    bi.bmiHeader.biCompression = BI_RGB;
    void* bits = nullptr;
    HBITMAP bmp = CreateDIBSection(mem, &bi, DIB_RGB_COLORS, &bits, nullptr, 0);
    ReleaseDC(nullptr, screen);
    if (!bmp) { DeleteDC(mem); return std::string(); }
    HGDIOBJ oldBmp = SelectObject(mem, bmp);

    RECT all = {0, 0, W, H};
    HBRUSH white = CreateSolidBrush(RGB(255, 255, 255));
    FillRect(mem, &all, white);
    DeleteObject(white);

    /* The drawing is in device units at kUnitsPerPt per point; the bitmap is
     * one pixel per point times the scale. */
    SetMapMode(mem, MM_ANISOTROPIC);
    SetWindowExtEx(mem, int(std::lround(w_pt * kUnitsPerPt * scale)),
                        int(std::lround(h_pt * kUnitsPerPt * scale)), nullptr);
    SetViewportExtEx(mem, W, H, nullptr);
    draw(mem, layout, style, scale);

    std::string out;
    if (as_dib) {
        /* A packed DIB: header, then the pixels bottom-up.  CF_DIB is what an
         * application reads when it pastes an image, and Windows does NOT
         * synthesise it from a metafile -- without it a paste into a browser
         * lands as text, which is what Google Slides was doing. */
        BITMAPINFOHEADER h = bi.bmiHeader;
        h.biHeight = H;                        /* positive: bottom-up */
        h.biSizeImage = DWORD(W) * DWORD(H) * 4;
        out.resize(sizeof(h) + h.biSizeImage);
        memcpy(&out[0], &h, sizeof(h));
        const uint8_t* src = static_cast<const uint8_t*>(bits);
        uint8_t* dst = reinterpret_cast<uint8_t*>(&out[sizeof(h)]);
        const size_t stride = size_t(W) * 4;
        for (int y = 0; y < H; ++y)            /* flip: our DIB is top-down */
            memcpy(dst + size_t(y) * stride,
                   src + size_t(H - 1 - y) * stride, stride);
    } else {
        Gdiplus::Bitmap image(bmp, nullptr);
        CLSID png;
        if (gdiplus_encoder(L"image/png", &png)) {
            IStream* stream = nullptr;
            if (CreateStreamOnHGlobal(nullptr, TRUE, &stream) == S_OK) {
                if (image.Save(stream, &png, nullptr) == Gdiplus::Ok) {
                    HGLOBAL h = nullptr;
                    GetHGlobalFromStream(stream, &h);
                    SIZE_T n = GlobalSize(h);
                    void* p = GlobalLock(h);
                    out.assign(static_cast<char*>(p), n);
                    GlobalUnlock(h);
                }
                stream->Release();
            }
        }
    }

    SelectObject(mem, oldBmp);
    DeleteObject(bmp);
    DeleteDC(mem);
    return out;
}

std::string render_png(const Layout& layout, const SvgStyle& style,
                       double scale) {
    return rasterize(layout, style, scale, false);
}

std::string render_dib(const Layout& layout, const SvgStyle& style,
                       double scale) {
    return rasterize(layout, style, scale, true);
}

#else   /* not Windows */

std::string render_emf(const Layout&, const SvgStyle&) { return std::string(); }
std::string render_png(const Layout&, const SvgStyle&, double) { return std::string(); }
std::string render_dib(const Layout&, const SvgStyle&, double) { return std::string(); }

#endif

std::string tex_to_emf(const std::string& latex, const SvgStyle& style) {
    std::unique_ptr<LineNode> root = parse_latex(latex);
    if (!root) return std::string();
    return render_emf(layout_math(*root, style), style);
}

std::string tex_to_png(const std::string& latex, const SvgStyle& style,
                       double scale) {
    std::unique_ptr<LineNode> root = parse_latex(latex);
    if (!root) return std::string();
    return render_png(layout_math(*root, style), style, scale);
}

std::string tex_to_dib(const std::string& latex, const SvgStyle& style,
                       double scale) {
    std::unique_ptr<LineNode> root = parse_latex(latex);
    if (!root) return std::string();
    return render_dib(layout_math(*root, style), style, scale);
}

}  // namespace mtef
