/*
 * eq_window.cpp -- the editor window
 *
 * Three things here are deliberate.
 *
 * The key chords are BUILT FROM Equation::shortcuts(), not written out again.
 * That table is what says the editor keeps Equation Editor 3.0's feel; a second
 * copy in the window would drift from it and the feel would quietly rot.
 *
 * Every paint recomputes the layout from the tree and draws the whole thing
 * into a back buffer.  Equation Editor needed a Redraw command because its
 * incremental display could go stale; nothing here can, because nothing is
 * kept between paints.
 *
 * And the window is per-monitor DPI aware, because the single most dated thing
 * about the old editor on a modern screen is that it is blurry.
 */
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>
#include <windowsx.h>          /* GET_X_LPARAM */
#include <commdlg.h>           /* the open / save dialogs */

#include "eq_window.h"

#include "eq_chords.h"
#include "eq_edit.h"
#include "gvml_clip.h"
#include "mtef_gdi.h"
#include "mtef_mathml.h"
#include "mtef_omml.h"
#include "mtef_rtf.h"
#include "mtef_svg.h"

#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

namespace mtef {
namespace {

const wchar_t* kClassName = L"Eqnedt64Window";

/* The raster the clipboard carries.  Equation Editor's era could not have
 * afforded this; an equation at 600 dpi is a few tens of kilobytes and takes
 * milliseconds, so the picture is print quality rather than screen quality. */
const double kPasteDpi = 600.0;
/* The equation is shown larger than life, as the old editor did at 200%.
 * Equation Editor offers 100/200/400 and a custom value; the wheel covers the
 * same ground continuously and the two ends just keep it usable. */
const double kZoomDefault = 2.0;
const double kZoomMin = 0.5;
const double kZoomMax = 8.0;
const int kMargin = 24;            /* device-independent pixels */

std::wstring widen(const std::string& s) {
    if (s.empty()) return std::wstring();
    int n = MultiByteToWideChar(CP_UTF8, 0, s.data(), int(s.size()), nullptr, 0);
    std::wstring w(size_t(n), L'\0');
    MultiByteToWideChar(CP_UTF8, 0, s.data(), int(s.size()), &w[0], n);
    return w;
}

std::string narrow(const std::wstring& w) {
    if (w.empty()) return std::string();
    int n = WideCharToMultiByte(CP_UTF8, 0, w.data(), int(w.size()),
                                nullptr, 0, nullptr, nullptr);
    std::string s(size_t(n), '\0');
    WideCharToMultiByte(CP_UTF8, 0, w.data(), int(w.size()), &s[0], n,
                        nullptr, nullptr);
    return s;
}

/* ---- the palette bar ---------------------------------------------------- */

const wchar_t* kPopupClass = L"Eqnedt64Palette";

/* Half again the size Equation Editor drew its bar at.
 *
 * Its palette was designed for a 96 dpi CRT and a 14-inch screen, and the
 * samples in it are EQUATIONS -- a fraction with two slots, a matrix with
 * four -- drawn at 9 and 12 point in a 34-pixel cell.  On a modern display
 * that is a smudge you have to lean towards.  Everything about the bar scales
 * together, the type inside the buttons with the buttons, so the proportions
 * are the ones that were there and only the reading distance changes.
 *
 * This multiplies the DPI scaling rather than replacing it: a 150 % display
 * still gets 150 % of this. */
const double kPaletteScale = 1.5;

const int kBtnW = int(46 * kPaletteScale), kBtnH = int(24 * kPaletteScale);
const int kCellW = int(34 * kPaletteScale), kCellH = int(30 * kPaletteScale);
const int kBarPad = int(3 * kPaletteScale);
const double kBarPt = 9.0 * kPaletteScale;    /* type size inside a button */
const double kCellPt = 12.0 * kPaletteScale;  /* type size inside a cell   */

/* How much room an empty slot takes on screen.  A template nobody has typed
 * into has no extent of its own, so without this a fresh fraction is a bar
 * floating in space and Tab appears to go nowhere. */
const double kEmptySlotEm = 0.55;

SvgStyle editing_style() {
    SvgStyle st;
    st.empty_slot_em = kEmptySlotEm;
    return st;
}

/* What a cell shows is what inserting it produces: the sample is rendered by
 * actually performing the insertion into a scratch equation.  A hand-written
 * table of sample LaTeX would drift from what the templates really are. */
Layout sample_layout(const Equation::PaletteItem& item, double sizePt) {
    SvgStyle st;
    const double k = sizePt / 12.0;
    st.full = 12.0 * k; st.sub = 7.0 * k; st.sub2 = 5.0 * k;
    st.sym = 18.0 * k;  st.subsym = 12.0 * k;
    st.padding = 0.0;
    st.empty_slot_em = kEmptySlotEm;   /* on screen: show there is a slot */

    Equation e;
    if (item.is_template) e.insert_template(item.command);
    else                  e.insert_symbol(item.command);
    return e.layout(st);
}

struct Cell {
    Equation::PaletteItem item;
    Layout layout;
    RECT rc{};
};

struct Button {
    const Equation::PaletteGroup* group = nullptr;
    Layout sample;                      /* a few of its members, drawn small */
    RECT rc{};
};

/* ---- the window --------------------------------------------------------- */

struct Editor {
    Equation eq;
    /* The editing view, which is the one that shows empty slots.  What gets
     * pasted out is laid out with a plain SvgStyle, so a picture on a slide is
     * the equation and nothing else. */
    SvgStyle style = editing_style();
    bool copied = false;
    bool caret_on = true;
    int dpi = 96;
    KeyState keys;                  /* a chord waiting for its second key */
    bool swallow_char = false;      /* this press was a chord, not typing */

    std::wstring path;              /* where Ctrl+S writes; empty until asked */
    std::string saved;              /* the equation as last written or read */
    std::wstring title_shown;       /* so the bar is not rewritten per key */

    bool dragging = false;          /* the mouse is selecting a range */
    std::vector<Button> bar;
    HWND popup = nullptr;
    std::vector<Cell> cells;        /* what the open popup is showing */
    int hot_cell = -1;

    double zoom = kZoomDefault;
    double units_per_pt() const { return dpi / 72.0 * zoom; }
    int scaled(int dip) const { return MulDiv(dip, dpi, 96); }
    int bar_height() const { return scaled(kBtnH) * 2 + scaled(kBarPad) * 3; }
    int status_height() const { return scaled(22); }
};

/* Two rows, symbols above templates -- Equation Editor's arrangement. */
void build_bar(Editor& ed) {
    ed.bar.clear();
    const int w = ed.scaled(kBtnW), h = ed.scaled(kBtnH), pad = ed.scaled(kBarPad);
    const std::vector<Equation::PaletteGroup>* rows[2] = {
        &Equation::symbol_palettes(), &Equation::template_palettes()};

    for (int r = 0; r < 2; ++r) {
        int x = pad;
        const int y = pad + r * (h + pad);
        for (const Equation::PaletteGroup& g : *rows[r]) {
            Button b;
            b.group = &g;
            b.rc = {x, y, x + w, y + h};
            /* The button wears its own contents: the first few members, drawn
             * by the same routine that will draw them once inserted. */
            SvgStyle st;
            const double k = kBarPt / 12.0;
            st.full = 12.0 * k; st.sub = 7.0 * k; st.sub2 = 5.0 * k;
            st.sym = 18.0 * k;  st.subsym = 12.0 * k;
            st.padding = 0.0;
            st.empty_slot_em = kEmptySlotEm;
            Equation e;
            for (size_t i = 0; i < g.items.size() && i < 3; ++i) {
                if (g.items[i].is_template) e.insert_template(g.items[i].command);
                else                        e.insert_symbol(g.items[i].command);
            }
            b.sample = e.layout(st);
            ed.bar.push_back(b);
            x += w + pad;
        }
    }
}

Editor* editor_of(HWND hwnd) {
    return reinterpret_cast<Editor*>(GetWindowLongPtrW(hwnd, GWLP_USERDATA));
}

/* Centre a laid-out equation in a rectangle and draw it. */
void draw_centred(HDC dc, const Layout& L, const RECT& rc, double upp,
                  COLORREF colour) {
    if (L.glyphs.empty() && L.rules.empty() && L.empty_slots.empty()) return;
    SvgStyle st;

    /* Shrink to fit before centring.  A palette cell is a fixed rectangle and
     * what goes in it is a real equation laid out by the real layout, so its
     * size is not ours to choose -- a radical that now asks the font for a
     * taller drawing, or a matrix with more rows, will simply be bigger than
     * the button.  Centring alone let it spill over the neighbours. */
    const double availW = std::max(1.0, double(rc.right - rc.left) - 2.0);
    const double availH = std::max(1.0, double(rc.bottom - rc.top) - 2.0);
    const double needW = L.w * upp;
    const double needH = (L.asc + L.desc) * upp;
    double fit = 1.0;
    if (needW > availW) fit = std::min(fit, availW / needW);
    if (needH > availH) fit = std::min(fit, availH / needH);
    const double u = upp * fit;

    const int x = rc.left + int(((rc.right - rc.left) - L.w * u) / 2);
    const int y = rc.top + int(((rc.bottom - rc.top) + (L.asc - L.desc) * u) / 2);
    /* Slot boxes on: a palette cell for a template is ALL empty slots, so
     * without them the button for scripts, matrices and accents is blank --
     * which is exactly how those three buttons came to look broken. */
    draw_layout(dc, L, st, u, x, y, colour, true);
}

/* The strip along the bottom, the way Equation Editor has one.
 *
 * A mode you cannot see is a trap: after Ctrl+Shift+G everything typed comes
 * out Greek, and without somewhere saying so the only way to find out is to
 * type something and be surprised.  Equation Editor solved this by reading
 * "Style: Math" along the bottom at all times, so this does too. */
std::wstring pretty_style(const std::string& s) {
    if (s == "math")     return L"Math";
    if (s == "text")     return L"Text";
    if (s == "function") return L"Function";
    if (s == "variable") return L"Variable";
    if (s == "vector")   return L"Matrix-Vector";
    if (s == "greek")    return L"Greek";
    return widen(s);
}

void paint_status(HDC dc, Editor& ed, const RECT& rc) {
    FillRect(dc, &rc, GetSysColorBrush(COLOR_BTNFACE));
    HPEN pen = CreatePen(PS_SOLID, 1, GetSysColor(COLOR_BTNSHADOW));
    HGDIOBJ oldPen = SelectObject(dc, pen);
    MoveToEx(dc, rc.left, rc.top, nullptr);
    LineTo(dc, rc.right, rc.top);
    SelectObject(dc, oldPen);
    DeleteObject(pen);

    LOGFONTW lf = {};
    lf.lfHeight = -ed.scaled(11);
    lf.lfCharSet = DEFAULT_CHARSET;
    lf.lfQuality = ANTIALIASED_QUALITY;
    wcscpy_s(lf.lfFaceName, L"Segoe UI");
    HFONT f = CreateFontIndirectW(&lf);
    HGDIOBJ oldFont = SelectObject(dc, f);
    SetBkMode(dc, TRANSPARENT);
    SetTextColor(dc, GetSysColor(COLOR_BTNTEXT));

    const std::wstring cells[3] = {
        L"Style: " + pretty_style(ed.eq.style()),
        L"Size: Full",
        L"Zoom: " + std::to_wstring(int(std::lround(ed.zoom * 100))) + L"%",
    };
    const int w = (rc.right - rc.left) / 4;
    for (int i = 0; i < 3; ++i) {
        RECT cell = {rc.left + ed.scaled(6) + i * w, rc.top,
                     rc.left + (i + 1) * w, rc.bottom};
        DrawTextW(dc, cells[i].c_str(), -1, &cell,
                  DT_LEFT | DT_SINGLELINE | DT_VCENTER);
    }
    SelectObject(dc, oldFont);
    DeleteObject(f);
}

void close_popup(Editor& ed) {
    if (ed.popup) { DestroyWindow(ed.popup); ed.popup = nullptr; }
    ed.cells.clear();
    ed.hot_cell = -1;
}

LRESULT CALLBACK popup_proc(HWND hwnd, UINT msg, WPARAM wp, LPARAM lp);

/* Open the grid for one button, laid out under it. */
void open_popup(HWND owner, Editor& ed, const Button& b) {
    close_popup(ed);
    const int n = int(b.group->items.size());
    if (n == 0) return;

    int cols = 1;
    while (cols * cols < n) ++cols;
    if (cols > 10) cols = 10;
    const int rows = (n + cols - 1) / cols;

    const int cw = ed.scaled(kCellW), ch = ed.scaled(kCellH);
    const double upp = ed.dpi / 72.0;

    ed.cells.clear();
    for (int i = 0; i < n; ++i) {
        Cell c;
        c.item = b.group->items[size_t(i)];
        c.layout = sample_layout(c.item, kCellPt);
        const int cx = (i % cols) * cw, cy = (i / cols) * ch;
        c.rc = {cx, cy, cx + cw, cy + ch};
        ed.cells.push_back(c);
    }
    (void)upp;

    POINT at = {b.rc.left, b.rc.bottom};
    ClientToScreen(owner, &at);
    ed.popup = CreateWindowExW(
        WS_EX_TOOLWINDOW | WS_EX_TOPMOST, kPopupClass, L"",
        WS_POPUP | WS_BORDER, at.x, at.y, cols * cw + 2, rows * ch + 2,
        owner, nullptr, GetModuleHandleW(nullptr), &ed);
    if (ed.popup) {
        ShowWindow(ed.popup, SW_SHOWNOACTIVATE);
        SetCapture(ed.popup);
    }
}

int cell_at(const Editor& ed, POINT p) {
    for (size_t i = 0; i < ed.cells.size(); ++i)
        if (PtInRect(&ed.cells[i].rc, p)) return int(i);
    return -1;
}

LRESULT CALLBACK popup_proc(HWND hwnd, UINT msg, WPARAM wp, LPARAM lp) {
    Editor* ed = editor_of(hwnd);
    switch (msg) {
        case WM_CREATE: {
            auto* cs = reinterpret_cast<CREATESTRUCTW*>(lp);
            SetWindowLongPtrW(hwnd, GWLP_USERDATA,
                              reinterpret_cast<LONG_PTR>(cs->lpCreateParams));
            return 0;
        }
        case WM_PAINT: {
            if (!ed) return 0;
            PAINTSTRUCT ps;
            HDC hdc = BeginPaint(hwnd, &ps);
            RECT all;
            GetClientRect(hwnd, &all);
            HDC mem = CreateCompatibleDC(hdc);
            HBITMAP bmp = CreateCompatibleBitmap(hdc, all.right, all.bottom);
            HGDIOBJ old = SelectObject(mem, bmp);

            HBRUSH bg = CreateSolidBrush(GetSysColor(COLOR_WINDOW));
            FillRect(mem, &all, bg);
            DeleteObject(bg);

            for (size_t i = 0; i < ed->cells.size(); ++i) {
                if (int(i) == ed->hot_cell) {
                    HBRUSH hot = CreateSolidBrush(GetSysColor(COLOR_HIGHLIGHT));
                    FillRect(mem, &ed->cells[i].rc, hot);
                    DeleteObject(hot);
                }
                draw_centred(mem, ed->cells[i].layout, ed->cells[i].rc,
                             ed->dpi / 72.0,
                             int(i) == ed->hot_cell
                                 ? GetSysColor(COLOR_HIGHLIGHTTEXT)
                                 : GetSysColor(COLOR_WINDOWTEXT));
            }
            BitBlt(hdc, 0, 0, all.right, all.bottom, mem, 0, 0, SRCCOPY);
            SelectObject(mem, old);
            DeleteObject(bmp);
            DeleteDC(mem);
            EndPaint(hwnd, &ps);
            return 0;
        }
        case WM_ERASEBKGND:
            return 1;
        case WM_MOUSEMOVE: {
            if (!ed) return 0;
            POINT p = {GET_X_LPARAM(lp), GET_Y_LPARAM(lp)};
            const int hit = cell_at(*ed, p);
            if (hit != ed->hot_cell) {
                ed->hot_cell = hit;
                InvalidateRect(hwnd, nullptr, FALSE);
            }
            return 0;
        }
        case WM_LBUTTONDOWN: {
            if (!ed) return 0;
            POINT p = {GET_X_LPARAM(lp), GET_Y_LPARAM(lp)};
            const int hit = cell_at(*ed, p);
            HWND owner = GetWindow(hwnd, GW_OWNER);
            if (hit >= 0) {
                const Equation::PaletteItem item = ed->cells[size_t(hit)].item;
                if (item.is_template) ed->eq.insert_template(item.command);
                else                  ed->eq.insert_symbol(item.command);
            }
            ReleaseCapture();
            close_popup(*ed);
            if (owner) {
                InvalidateRect(owner, nullptr, FALSE);
                SetFocus(owner);      /* typing continues where it left off */
            }
            return 0;
        }
        case WM_CAPTURECHANGED:
            return 0;
        case WM_DESTROY:
            if (GetCapture() == hwnd) ReleaseCapture();
            return 0;
    }
    return DefWindowProcW(hwnd, msg, wp, lp);
}

void paint(HWND hwnd, Editor& ed) {
    PAINTSTRUCT ps;
    HDC hdc = BeginPaint(hwnd, &ps);

    RECT rc;
    GetClientRect(hwnd, &rc);
    const int W = rc.right - rc.left, H = rc.bottom - rc.top;

    /* Draw into a back buffer: a half-drawn equation on screen is exactly the
     * flicker the old editor's Redraw button existed to paper over. */
    HDC mem = CreateCompatibleDC(hdc);
    HBITMAP bmp = CreateCompatibleBitmap(hdc, W, H);
    HGDIOBJ oldBmp = SelectObject(mem, bmp);

    HBRUSH bg = CreateSolidBrush(GetSysColor(COLOR_WINDOW));
    FillRect(mem, &rc, bg);
    DeleteObject(bg);

    /* ---- the palette bar ---- */
    const int bar_h = ed.bar_height();
    RECT bar_rc = {0, 0, W, bar_h};
    HBRUSH bar_bg = CreateSolidBrush(GetSysColor(COLOR_BTNFACE));
    FillRect(mem, &bar_rc, bar_bg);
    DeleteObject(bar_bg);
    for (const Button& b : ed.bar) {
        HBRUSH face = CreateSolidBrush(GetSysColor(COLOR_WINDOW));
        FillRect(mem, &b.rc, face);
        DeleteObject(face);
        FrameRect(mem, &b.rc, HBRUSH(GetStockObject(GRAY_BRUSH)));
        draw_centred(mem, b.sample, b.rc, ed.dpi / 72.0,
                     GetSysColor(COLOR_WINDOWTEXT));
    }

    const Layout L = ed.eq.layout(ed.style);
    const double upp = ed.units_per_pt();
    const int originX = ed.scaled(kMargin);
    const int originY = bar_h +
        int(((H - bar_h - ed.status_height()) + L.asc * upp - L.desc * upp) / 2);

    /* The highlight goes down first so the equation draws over it, which keeps
     * the glyphs their own colour instead of inverting them. */
    Equation::SelectionBox sel = ed.eq.selection_geometry(ed.style);
    if (sel.found) {
        RECT r;
        r.left   = originX + int(std::min(sel.x0, sel.x1) * upp);
        r.right  = originX + int(std::max(sel.x0, sel.x1) * upp);
        r.top    = originY + int(sel.top * upp);
        r.bottom = originY + int(sel.bottom * upp);
        if (r.right <= r.left) r.right = r.left + 1;
        HBRUSH hl = CreateSolidBrush(GetSysColor(COLOR_HIGHLIGHT));
        FillRect(mem, &r, hl);
        DeleteObject(hl);
    }

    draw_layout(mem, L, ed.style, upp, originX, originY,
                GetSysColor(COLOR_WINDOWTEXT), true);

    if (ed.caret_on && !sel.found) {
        Equation::CaretGeometry g = ed.eq.caret_geometry(ed.style);
        if (g.found) {
            RECT c;
            c.left   = originX + int(g.x * upp);
            c.right  = c.left + std::max(1, ed.scaled(1));
            c.top    = originY + int(g.top * upp);
            c.bottom = originY + int(g.bottom * upp);
            if (c.bottom <= c.top) c.bottom = c.top + ed.scaled(8);
            HBRUSH caret = CreateSolidBrush(GetSysColor(COLOR_WINDOWTEXT));
            FillRect(mem, &c, caret);
            DeleteObject(caret);
        }
    }

    /* ---- the status strip ---- */
    RECT st_rc = {0, H - ed.status_height(), W, H};
    paint_status(mem, ed, st_rc);

    BitBlt(hdc, 0, 0, W, H, mem, 0, 0, SRCCOPY);
    SelectObject(mem, oldBmp);
    DeleteObject(bmp);
    DeleteDC(mem);
    EndPaint(hwnd, &ps);
}

void redraw(HWND hwnd) { InvalidateRect(hwnd, nullptr, FALSE); }

/* ---- the file ----------------------------------------------------------- */

/* UTF-8 with a BOM.  The equation may hold Japanese, and without the mark
 * every editor on this machine reads the file as cp932 and shows nothing but
 * mojibake.  It is stripped again on load. */
const char kBom[] = "\xEF\xBB\xBF";

bool modified(const Editor& ed) { return ed.eq.latex() != ed.saved; }

std::wstring file_name(const std::wstring& path) {
    const size_t cut = path.find_last_of(L"\\/");
    return cut == std::wstring::npos ? path : path.substr(cut + 1);
}

std::wstring ask_for_path(HWND hwnd, bool saving) {
    wchar_t buf[MAX_PATH] = L"";
    OPENFILENAMEW ofn{};
    ofn.lStructSize = sizeof(ofn);
    ofn.hwndOwner   = hwnd;
    ofn.lpstrFilter = L"LaTeX equation\0*.tex\0All files\0*.*\0";
    ofn.lpstrFile   = buf;
    ofn.nMaxFile    = MAX_PATH;
    ofn.lpstrDefExt = L"tex";
    ofn.Flags = OFN_PATHMUSTEXIST | OFN_NOCHANGEDIR
              | (saving ? OFN_OVERWRITEPROMPT : OFN_FILEMUSTEXIST);
    const BOOL ok = saving ? GetSaveFileNameW(&ofn) : GetOpenFileNameW(&ofn);
    return ok ? std::wstring(buf) : std::wstring();
}

bool write_file(const std::wstring& path, const std::string& text) {
    HANDLE h = CreateFileW(path.c_str(), GENERIC_WRITE, 0, nullptr,
                           CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, nullptr);
    if (h == INVALID_HANDLE_VALUE) return false;
    DWORD wrote = 0;
    const BOOL ok = WriteFile(h, text.data(), DWORD(text.size()), &wrote,
                              nullptr);
    CloseHandle(h);
    return ok && wrote == text.size();
}

bool read_file(const std::wstring& path, std::string& out) {
    HANDLE h = CreateFileW(path.c_str(), GENERIC_READ, FILE_SHARE_READ, nullptr,
                           OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, nullptr);
    if (h == INVALID_HANDLE_VALUE) return false;
    LARGE_INTEGER size{};
    if (!GetFileSizeEx(h, &size) || size.QuadPart > (1 << 20)) {
        CloseHandle(h);
        return false;
    }
    out.resize(size_t(size.QuadPart));
    DWORD got = 0;
    const BOOL ok = out.empty() ? TRUE
                  : ReadFile(h, &out[0], DWORD(out.size()), &got, nullptr);
    CloseHandle(h);
    if (!ok) return false;
    out.resize(got);
    if (out.compare(0, 3, kBom) == 0) out.erase(0, 3);
    while (!out.empty() && (out.back() == '\n' || out.back() == '\r'))
        out.pop_back();
    return true;
}

/* The title carries the build, the file and the zoom.  The build because an
 * install pointing at a stale tree has bitten this repository repeatedly; the
 * file and its unsaved mark because an editor that cannot tell you whether
 * your work is on disk is one you cannot trust it to; the zoom because the
 * wheel is continuous, and a magnification you cannot read off is one you
 * cannot get back to.
 *
 * Called from WM_PAINT, so the mark tracks every edit without a flag threaded
 * through each of them.  It writes only on a change: the title bar flickers
 * otherwise, once per keystroke. */
void update_title(HWND hwnd, Editor& ed) {
#ifdef EQNEDT64_VERSION
    std::wstring t = L"EQNEDT64 " + widen(EQNEDT64_VERSION);
#else
    std::wstring t = L"EQNEDT64 (development build)";
#endif
    t += L"  -  " + (ed.path.empty() ? std::wstring(L"untitled")
                                     : file_name(ed.path));
    if (modified(ed)) t += L" *";
    t += L"  -  " + std::to_wstring(int(std::lround(ed.zoom * 100))) + L"%";
    if (t == ed.title_shown) return;
    ed.title_shown = t;
    SetWindowTextW(hwnd, t.c_str());
}

/* Ctrl+S, and Ctrl+Shift+S to choose a new name.  False when the user backed
 * out of the dialog, so closing can be called off too. */
bool save_equation(HWND hwnd, Editor& ed, bool ask_where) {
    if (ask_where || ed.path.empty()) {
        std::wstring chosen = ask_for_path(hwnd, true);
        if (chosen.empty()) return false;
        ed.path = chosen;
    }
    if (!write_file(ed.path, kBom + ed.eq.latex() + "\n")) {
        MessageBoxW(hwnd, L"The file could not be written.", L"EQNEDT64",
                    MB_OK | MB_ICONERROR);
        return false;
    }
    ed.saved = ed.eq.latex();
    update_title(hwnd, ed);
    return true;
}

void open_equation(HWND hwnd, Editor& ed) {
    const std::wstring chosen = ask_for_path(hwnd, false);
    if (chosen.empty()) return;
    std::string text;
    if (!read_file(chosen, text)) {
        MessageBoxW(hwnd, L"The file could not be read.", L"EQNEDT64",
                    MB_OK | MB_ICONERROR);
        return;
    }
    ed.eq.load_latex(text);
    ed.eq.move_end();
    ed.path  = chosen;
    ed.saved = ed.eq.latex();
    update_title(hwnd, ed);
    InvalidateRect(hwnd, nullptr, TRUE);
}

/* Ctrl+V.  The editor was a one-way door without this: an equation could be
 * written and sent out, and nothing could be brought back in.
 *
 * The clipboard's plain text is taken as LaTeX, which is the format this
 * editor stores anyway -- so what it copied out is exactly what it reads back,
 * and so is anything from a paper, a note, or another tool. */
void paste_from_clipboard(HWND hwnd, Editor& ed) {
    if (!OpenClipboard(hwnd)) return;
    std::string latex;
    if (HANDLE h = GetClipboardData(CF_UNICODETEXT)) {
        if (const wchar_t* w = static_cast<const wchar_t*>(GlobalLock(h))) {
            latex = narrow(w);
            GlobalUnlock(h);
        }
    }
    CloseClipboard();
    if (latex.empty()) return;
    if (ed.eq.insert_latex(latex)) redraw(hwnd);
}

/* The window's whole part in this: read the modifier state, which only a
 * window can, and hand the press to the table.  Everything that DECIDES
 * anything lives in press_key, where it can be tested by calling it. */
bool handle_key(HWND hwnd, Editor& ed, UINT vk) {
    const bool ctrl  = (GetKeyState(VK_CONTROL) & 0x8000) != 0;
    const bool shift = (GetKeyState(VK_SHIFT)   & 0x8000) != 0;
    const bool alt   = (GetKeyState(VK_MENU)    & 0x8000) != 0;

    const KeyResult r = press_key(ed.eq, ed.keys, vk, ctrl, shift, alt);
    ed.swallow_char = ed.keys.swallow_char;
    if (r == KeyResult::Ignored) return false;
    if (r == KeyResult::Consumed) redraw(hwnd);
    return true;
}

LRESULT CALLBACK proc(HWND hwnd, UINT msg, WPARAM wp, LPARAM lp) {
    Editor* ed = editor_of(hwnd);

    switch (msg) {
        case WM_CREATE: {
            auto* cs = reinterpret_cast<CREATESTRUCTW*>(lp);
            SetWindowLongPtrW(hwnd, GWLP_USERDATA,
                              reinterpret_cast<LONG_PTR>(cs->lpCreateParams));
            SetTimer(hwnd, 1, GetCaretBlinkTime(), nullptr);
            return 0;
        }
        case WM_TIMER:
            if (ed) { ed->caret_on = !ed->caret_on; redraw(hwnd); }
            return 0;

        case WM_PAINT:
            if (ed) { update_title(hwnd, *ed); paint(hwnd, *ed); }
            return 0;

        case WM_ERASEBKGND:
            return 1;              /* the back buffer already covers it */

        case WM_DPICHANGED: {
            if (ed) { ed->dpi = HIWORD(wp); build_bar(*ed); }
            const RECT* want = reinterpret_cast<const RECT*>(lp);
            SetWindowPos(hwnd, nullptr, want->left, want->top,
                         want->right - want->left, want->bottom - want->top,
                         SWP_NOZORDER | SWP_NOACTIVATE);
            redraw(hwnd);
            return 0;
        }

        case WM_LBUTTONDOWN: {
            if (!ed) return 0;
            const POINT p = {GET_X_LPARAM(lp), GET_Y_LPARAM(lp)};

            for (const Button& b : ed->bar) {
                if (PtInRect(&b.rc, p)) { open_popup(hwnd, *ed, b); return 0; }
            }
            if (p.y < ed->bar_height()) return 0;   /* bar background */

            RECT rc;
            GetClientRect(hwnd, &rc);
            const int bar_h = ed->bar_height();
            const Layout L = ed->eq.layout(ed->style);
            const double upp = ed->units_per_pt();
            const int originX = ed->scaled(kMargin);
            const int originY = bar_h +
                int(((rc.bottom - bar_h - ed->status_height())
                     + L.asc * upp - L.desc * upp) / 2);
            const double ex = (p.x - originX) / upp;
            const double ey = (p.y - originY) / upp;

            /* Shift+click extends from where the caret already is, the way it
             * does in every text editor. */
            if (GetKeyState(VK_SHIFT) & 0x8000) ed->eq.extend_to_point(ex, ey, ed->style);
            else                                ed->eq.move_to_point(ex, ey, ed->style);

            ed->dragging = true;
            SetCapture(hwnd);
            ed->caret_on = true;
            redraw(hwnd);
            return 0;
        }

        case WM_MOUSEMOVE: {
            if (!ed || !ed->dragging) return 0;
            RECT rc;
            GetClientRect(hwnd, &rc);
            const int bar_h = ed->bar_height();
            const Layout L = ed->eq.layout(ed->style);
            const double upp = ed->units_per_pt();
            const int originX = ed->scaled(kMargin);
            const int originY = bar_h +
                int(((rc.bottom - bar_h - ed->status_height())
                     + L.asc * upp - L.desc * upp) / 2);
            ed->eq.extend_to_point((GET_X_LPARAM(lp) - originX) / upp,
                                   (GET_Y_LPARAM(lp) - originY) / upp, ed->style);
            redraw(hwnd);
            return 0;
        }

        case WM_LBUTTONUP:
            if (ed && ed->dragging) { ed->dragging = false; ReleaseCapture(); }
            return 0;

        case WM_SIZE:
            if (ed) { build_bar(*ed); redraw(hwnd); }
            return 0;

        case WM_MOUSEWHEEL: {
            if (!ed) return 0;
            /* There is nothing to scroll -- one equation, always fully
             * visible -- so the wheel is free to do the useful thing. */
            const int notches = GET_WHEEL_DELTA_WPARAM(wp) / WHEEL_DELTA;
            ed->zoom *= std::pow(1.15, notches);
            ed->zoom = std::max(kZoomMin, std::min(kZoomMax, ed->zoom));
            update_title(hwnd, *ed);
            redraw(hwnd);
            return 0;
        }

        case WM_KEYDOWN: {
            if (!ed) return 0;
            /* These are PLAIN ctrl bindings, so they must not answer to a
             * chord that adds shift: Ctrl+Shift+X is the matrix-vector style,
             * and cutting the selection instead would delete the very thing
             * the user asked to make bold. */
            const bool ctrl = (GetKeyState(VK_CONTROL) & 0x8000) != 0
                           && (GetKeyState(VK_SHIFT)   & 0x8000) == 0;
            const bool ctrl_shift = (GetKeyState(VK_CONTROL) & 0x8000) != 0
                                 && (GetKeyState(VK_SHIFT)   & 0x8000) != 0;
            if (wp == VK_ESCAPE) { PostMessageW(hwnd, WM_CLOSE, 0, 0); return 0; }

            if (ctrl       && wp == 'S') { save_equation(hwnd, *ed, false); return 0; }
            if (ctrl_shift && wp == 'S') { save_equation(hwnd, *ed, true);  return 0; }
            if (ctrl       && wp == 'O') { open_equation(hwnd, *ed);        return 0; }

            /* Copy takes the selection when there is one and the whole
             * equation otherwise, which is what every editor does and what
             * makes "copy this bit into a slide" possible at all. */
            if (ctrl && (wp == 'C' || wp == 'X' || wp == VK_RETURN)) {
                const bool cut = (wp == 'X');
                const std::string what = ed->eq.has_selection()
                                       ? ed->eq.selected_latex()
                                       : ed->eq.latex();
                ed->copied = copy_equation_to_clipboard(what);
                if (cut && ed->eq.has_selection()) {
                    ed->eq.delete_selection();
                    redraw(hwnd);
                }
                return 0;
            }
            if (ctrl && wp == 'V') { paste_from_clipboard(hwnd, *ed); return 0; }

            /* Equation Editor's View menu offers 100 / 200 / 400; the digits
             * are the same, without a menu to put them in. */
            if (ctrl && (wp == '1' || wp == '2' || wp == '4')) {
                ed->zoom = (wp == '1') ? 1.0 : (wp == '2') ? 2.0 : 4.0;
                update_title(hwnd, *ed);
                redraw(hwnd);
                return 0;
            }
            /* WM_CHAR follows its own WM_KEYDOWN, so clearing here and setting
             * it below cannot leak into the next key. */
            ed->swallow_char = false;
            handle_key(hwnd, *ed, UINT(wp));
            return 0;
        }

        /* Alt-modified keys arrive HERE, not at WM_KEYDOWN -- Windows routes
         * anything pressed with Alt held as a system key.  The published
         * table has such a chord (Ctrl+Alt+Space, the quad), and without this
         * case it was dead: the press went to DefWindowProc and became menu
         * activation instead of ever reaching the table. */
        case WM_SYSKEYDOWN:
            if (ed) {
                ed->swallow_char = false;
                if (handle_key(hwnd, *ed, UINT(wp))) return 0;
            }
            break;   /* DefWindowProc keeps Alt+F4 and Alt+Space working */

        /* The char half of a consumed system key.  Left to DefWindowProc it
         * would open the system menu -- Alt+Space's other meaning. */
        case WM_SYSCHAR:
            if (ed && ed->swallow_char) { ed->swallow_char = false; return 0; }
            break;

        case WM_CHAR: {
            if (!ed) return 0;
            const wchar_t ch = wchar_t(wp);
            if (ch < 0x20) return 0;                       /* control keys */
            if (GetKeyState(VK_CONTROL) & 0x8000) return 0; /* a chord, not text */
            /* The second key of "Ctrl+T, S" is pressed WITHOUT ctrl, so the
             * test above lets it through and the S would be typed on top of
             * the summation the chord just built. */
            if (ed->swallow_char) { ed->swallow_char = false; return 0; }
            ed->eq.insert_text(narrow(std::wstring(1, ch)));
            ed->caret_on = true;
            redraw(hwnd);
            return 0;
        }

        /* Escape used to destroy the window outright, so an afternoon's
         * equation left with it.  Ask first. */
        case WM_CLOSE:
            if (ed && modified(*ed)) {
                const int answer = MessageBoxW(
                    hwnd, L"Save the changes to this equation?", L"EQNEDT64",
                    MB_YESNOCANCEL | MB_ICONWARNING);
                if (answer == IDCANCEL) return 0;
                if (answer == IDYES && !save_equation(hwnd, *ed, false))
                    return 0;          /* the save dialog was cancelled */
            }
            DestroyWindow(hwnd);
            return 0;

        case WM_DESTROY:
            KillTimer(hwnd, 1);
            PostQuitMessage(0);
            return 0;
    }
    return DefWindowProcW(hwnd, msg, wp, lp);
}

/* ---- the clipboard ------------------------------------------------------ */

bool put(UINT format, const std::string& bytes) {
    if (bytes.empty()) return false;
    HGLOBAL h = GlobalAlloc(GMEM_MOVEABLE, bytes.size() + 1);
    if (!h) return false;
    void* p = GlobalLock(h);
    memcpy(p, bytes.data(), bytes.size());
    static_cast<char*>(p)[bytes.size()] = '\0';
    GlobalUnlock(h);
    if (!SetClipboardData(format, h)) { GlobalFree(h); return false; }
    return true;
}

/* Register both window classes once.  The editor and the self-test create the
 * same windows, so they must register the same classes -- a second copy of
 * these two blocks would drift.  Registering twice in one process is normal
 * (Python can open the editor repeatedly), so "already exists" is success. */
bool register_window_classes() {
    WNDCLASSEXW wc = {};
    wc.cbSize = sizeof(wc);
    wc.lpfnWndProc = proc;
    wc.hInstance = GetModuleHandleW(nullptr);
    wc.hCursor = LoadCursorW(nullptr, IDC_IBEAM);
    wc.lpszClassName = kClassName;
    if (!RegisterClassExW(&wc) &&
        GetLastError() != ERROR_CLASS_ALREADY_EXISTS) return false;

    WNDCLASSEXW pc = {};
    pc.cbSize = sizeof(pc);
    pc.lpfnWndProc = popup_proc;
    pc.hInstance = wc.hInstance;
    pc.hCursor = LoadCursorW(nullptr, IDC_ARROW);
    pc.lpszClassName = kPopupClass;
    if (!RegisterClassExW(&pc) &&
        GetLastError() != ERROR_CLASS_ALREADY_EXISTS) return false;
    return true;
}

}  // namespace

bool copy_equation_to_clipboard(const std::string& latex, bool display,
                                bool pictures) {
    RtfOptions rtf;
    rtf.display = display;
    MathMLOptions mml;
    mml.display = display;

    const std::string rtf_bytes = tex_to_rtf(latex, rtf);
    const std::string mml_bytes = tex_to_mathml(latex, mml);
    /* PowerPoint's own shape format, which is the only one that can state a
     * SIZE -- see gvml_clip.h.  Offered before the MathML, which PowerPoint
     * would otherwise take and paste at whatever the destination box is. */
    const std::string gvml_bytes = tex_to_gvml(latex, kPasteSizePt, display);

    SvgStyle style;
    std::string emf_bytes, png_bytes, dib_bytes;
    if (pictures) {
        emf_bytes = tex_to_emf(latex, style);
        png_bytes = tex_to_png(latex, style, kPasteDpi / 72.0);
        dib_bytes = tex_to_dib(latex, style, kPasteDpi / 72.0);
    }

    if (!OpenClipboard(nullptr)) return false;
    EmptyClipboard();

    bool ok = put(RegisterClipboardFormatW(L"Rich Text Format"), rtf_bytes);
    /* GVML first: PowerPoint takes the richest format it recognises, and this
     * is the one that arrives at the size we asked for rather than the size
     * of whatever box it landed in. */
    put(RegisterClipboardFormatW(L"Art::GVML ClipFormat"), gvml_bytes);
    ok = put(RegisterClipboardFormatW(L"MathML"), mml_bytes) && ok;

    if (!emf_bytes.empty()) {
        /* CF_ENHMETAFILE wants a metafile handle, and the system owns it once
         * it is on the clipboard -- it must not be deleted here. */
        HENHMETAFILE emf = SetEnhMetaFileBits(UINT(emf_bytes.size()),
            reinterpret_cast<const BYTE*>(emf_bytes.data()));
        if (emf && !SetClipboardData(CF_ENHMETAFILE, emf)) DeleteEnhMetaFile(emf);
    }
    if (!png_bytes.empty())
        put(RegisterClipboardFormatW(L"PNG"), png_bytes);
    /* CF_DIB is what an application takes when it pastes a picture, and it is
     * the only image format a browser will find here -- Windows synthesises it
     * from a bitmap but not from a metafile. */
    if (!dib_bytes.empty())
        put(CF_DIB, dib_bytes);

    /* The LaTeX itself, for Markdown, Jupyter, and any plain editor. */
    std::wstring text = widen(latex);
    HGLOBAL h = GlobalAlloc(GMEM_MOVEABLE, (text.size() + 1) * sizeof(wchar_t));
    if (h) {
        void* p = GlobalLock(h);
        memcpy(p, text.c_str(), (text.size() + 1) * sizeof(wchar_t));
        GlobalUnlock(h);
        if (!SetClipboardData(CF_UNICODETEXT, h)) GlobalFree(h);
    }

    CloseClipboard();
    return ok;
}

EditorResult run_equation_window(const std::string& latex,
                                 const std::wstring& path) {
    /* Per-monitor aware, so the equation is sharp on a 4K screen.  It fails
     * harmlessly when the host process already chose its awareness -- which is
     * the normal case when this is called from Python. */
    SetProcessDpiAwarenessContext(DPI_AWARENESS_CONTEXT_PER_MONITOR_AWARE_V2);

    Editor ed;
    std::string start = latex;
    if (!path.empty() && read_file(path, start)) ed.path = path;
    ed.eq.load_latex(start);
    ed.eq.move_end();
    ed.saved = ed.eq.latex();   /* a file just opened is not modified */

    if (!register_window_classes()) return EditorResult{};

    /* The title says which build this is.  An install pointing at a stale
     * tree has bitten this repository repeatedly, and "is the fix actually in
     * the thing on screen" has to be answerable by looking at it. */
#ifdef EQNEDT64_VERSION
    const std::wstring title = L"EQNEDT64 " + widen(EQNEDT64_VERSION);
#else
    const std::wstring title = L"EQNEDT64 (development build)";
#endif

    HWND hwnd = CreateWindowExW(
        0, kClassName, title.c_str(),
        WS_OVERLAPPEDWINDOW, CW_USEDEFAULT, CW_USEDEFAULT, 760, 260,
        nullptr, nullptr, GetModuleHandleW(nullptr), &ed);
    if (!hwnd) return EditorResult{};

    ed.dpi = int(GetDpiForWindow(hwnd));
    build_bar(ed);
    update_title(hwnd, ed);
    ShowWindow(hwnd, SW_SHOW);
    SetFocus(hwnd);

    MSG msg;
    while (GetMessageW(&msg, nullptr, 0, 0) > 0) {
        TranslateMessage(&msg);
        DispatchMessageW(&msg);
    }

    EditorResult r;
    r.copied = ed.copied;
    r.latex = ed.eq.latex();
    return r;
}

/* ==== the interaction-layer self-test ==================================== */
/*
 * Drives the real window through the real WndProc by injecting messages.
 * eq_window.h says what this is for; the short of it is that this file is the
 * one layer the Python tests cannot reach, and the editor's first crash in
 * real use happened here with nobody able to say which operation did it.
 *
 * Nothing below touches the real keyboard, mouse or foreground.  Modifier
 * keys are set on THIS THREAD's keyboard state (SetKeyboardState), which is
 * what the WndProc's GetKeyState reads for messages sent on the same thread,
 * and every message is aimed at the window by handle.  A person can keep
 * typing in another window while this runs, and it runs the same on a
 * headless CI desktop.
 */
namespace {

/* xorshift64* -- deterministic and seedable, so a crashing walk replays. */
struct Fuzz {
    unsigned long long s;
    explicit Fuzz(unsigned long long seed)
        : s(seed ? seed : 0x9E3779B97F4A7C15ull) {}
    unsigned long long next() {
        s ^= s >> 12;
        s ^= s << 25;
        s ^= s >> 27;
        return s * 0x2545F4914F6CDD1Dull;
    }
    unsigned pick(unsigned n) { return n ? unsigned(next() % n) : 0; }
};

/* The journal.  Each step is written and FLUSHED before it runs, so the last
 * line of a crashed run names the step that killed it -- which is the whole
 * difference between this harness and "it crashed after a minute". */
struct Journal {
    FILE* f = nullptr;
    ~Journal() {
        if (f) fclose(f);
    }
    void line(const char* fmt, ...) {
        if (!f) return;
        va_list ap;
        va_start(ap, fmt);
        vfprintf(f, fmt, ap);
        va_end(ap);
        fputc('\n', f);
        fflush(f);
    }
};

struct Driver {
    HWND hwnd = nullptr;
    Editor* ed = nullptr;
    Journal* log = nullptr;
    const SelftestOptions* opt = nullptr;
    int failures = 0;

    void fail(const char* what) {
        ++failures;
        log->line("FAIL %s", what);
    }

    /* ---- message plumbing ------------------------------------------- */

    void pump() {
        MSG m;
        while (PeekMessageW(&m, nullptr, 0, 0, PM_REMOVE)) {
            if (m.message == WM_QUIT) return;
            TranslateMessage(&m);
            DispatchMessageW(&m);
        }
    }

    /* A synchronous repaint through the real WM_PAINT, popup included.  The
     * invalidation gives BeginPaint a region even on a desktop nobody sees. */
    void paint() {
        InvalidateRect(hwnd, nullptr, FALSE);
        SendMessageW(hwnd, WM_PAINT, 0, 0);
        if (ed->popup) {
            InvalidateRect(ed->popup, nullptr, FALSE);
            SendMessageW(ed->popup, WM_PAINT, 0, 0);
        }
        pump();
    }

    /* Modifiers for the WndProc's GetKeyState, set on this thread only. */
    static void mods(bool ctrl, bool shift, bool alt) {
        BYTE ks[256] = {};
        GetKeyboardState(ks);
        ks[VK_CONTROL] = BYTE(ctrl ? 0x80 : 0);
        ks[VK_LCONTROL] = ks[VK_CONTROL];
        ks[VK_SHIFT] = BYTE(shift ? 0x80 : 0);
        ks[VK_LSHIFT] = ks[VK_SHIFT];
        ks[VK_MENU] = BYTE(alt ? 0x80 : 0);
        ks[VK_LMENU] = ks[VK_MENU];
        SetKeyboardState(ks);
    }

    /* What TranslateMessage would deliver after this key, for the alphabet
     * this test uses.  Letters, digits and space are enough: they are what
     * the two-step chords' second keys are, which is the path with the
     * history of typing the chord's own key. */
    static wchar_t char_of(unsigned vk, bool shift) {
        if (vk >= 'A' && vk <= 'Z') return wchar_t(shift ? vk : vk + 32);
        if (vk >= '0' && vk <= '9' && !shift) return wchar_t(vk);
        if (vk == VK_SPACE) return L' ';
        return 0;
    }

    /* One key press the way Windows delivers it: the down, the WM_CHAR that
     * TranslateMessage would synthesize, the up -- and WM_SYSKEYDOWN when Alt
     * is held, which is the routing that once made the quad-space chord
     * unreachable. */
    void key(unsigned vk, bool ctrl, bool shift, bool alt) {
        mods(ctrl, shift, alt);
        SendMessageW(hwnd, alt ? WM_SYSKEYDOWN : WM_KEYDOWN, vk, 0);
        const wchar_t c = char_of(vk, shift);
        if (c && !ctrl && !alt) SendMessageW(hwnd, WM_CHAR, WPARAM(c), 0);
        SendMessageW(hwnd, alt ? WM_SYSKEYUP : WM_KEYUP, vk, 0);
        mods(false, false, false);
    }

    void chr(unsigned cp) { SendMessageW(hwnd, WM_CHAR, WPARAM(cp), 0); }

    void click(int x, int y, bool shift = false) {
        mods(false, shift, false);
        SendMessageW(hwnd, WM_LBUTTONDOWN, MK_LBUTTON, MAKELPARAM(x, y));
        SendMessageW(hwnd, WM_LBUTTONUP, 0, MAKELPARAM(x, y));
        mods(false, false, false);
    }

    void drag(int x0, int y0, int x1, int y1) {
        SendMessageW(hwnd, WM_LBUTTONDOWN, MK_LBUTTON, MAKELPARAM(x0, y0));
        for (int i = 1; i <= 4; ++i) {
            const int x = x0 + (x1 - x0) * i / 4;
            const int y = y0 + (y1 - y0) * i / 4;
            SendMessageW(hwnd, WM_MOUSEMOVE, MK_LBUTTON, MAKELPARAM(x, y));
        }
        SendMessageW(hwnd, WM_LBUTTONUP, 0, MAKELPARAM(x1, y1));
    }

    void wheel(int notches) {
        const WORD delta = WORD(short(notches * WHEEL_DELTA));
        SendMessageW(hwnd, WM_MOUSEWHEEL, MAKEWPARAM(0, delta),
                     MAKELPARAM(0, 0));
    }

    void resize(int w, int h) {
        SetWindowPos(hwnd, nullptr, 0, 0, w, h,
                     SWP_NOMOVE | SWP_NOZORDER | SWP_NOACTIVATE);
        pump();
    }

    void dpi(unsigned d) {
        RECT want;
        GetWindowRect(hwnd, &want);
        SendMessageW(hwnd, WM_DPICHANGED, MAKEWPARAM(d, d), LPARAM(&want));
    }

    void blink() { SendMessageW(hwnd, WM_TIMER, 1, 0); }

    /* ---- fixture states ---------------------------------------------- */

    static const char* state_name(int s) {
        switch (s) {
            case 0: return "empty";
            case 1: return "text-end";
            case 2: return "in-slot";
            default: return "selection";
        }
    }

    /* The handover's four caret situations: an empty document, the end of a
     * line of text, inside a template's slot, and with a selection. */
    void state(int s) {
        close_popup(*ed);
        ed->dragging = false;
        ed->keys = KeyState{};
        ed->swallow_char = false;
        ed->zoom = kZoomDefault;
        switch (s) {
            case 0:
                ed->eq.load_latex("");
                break;
            case 1:
                ed->eq.load_latex("ab+c");
                ed->eq.move_end();
                break;
            case 2:
                ed->eq.load_latex("");
                ed->eq.insert_template("frac");
                ed->eq.insert_text("x");
                break;
            default:
                ed->eq.load_latex("ab");
                ed->eq.move_end();
                ed->eq.extend_left();
                break;
        }
    }

    /* Keys a headless run may not send: they open modal dialogs.  Escape is
     * here because its WM_CLOSE asks about unsaved changes. */
    bool modal(unsigned vk, bool ctrl, bool shift) const {
        if (vk == VK_ESCAPE) return true;
        if (ctrl && !shift && (vk == 'S' || vk == 'O')) return true;
        if (ctrl && shift && vk == 'S') return true;
        return false;
    }

    /* The clipboard is the user's; only a run that asked may write it. */
    bool clipboardy(unsigned vk, bool ctrl, bool shift) const {
        return ctrl && !shift &&
               (vk == 'C' || vk == 'X' || vk == 'V' || vk == VK_RETURN);
    }

    bool allowed(const Chord& c) const {
        for (const Step& st : c.steps) {
            if (modal(st.vk, st.ctrl, st.shift)) return false;
            if (!opt->clipboard && clipboardy(st.vk, st.ctrl, st.shift))
                return false;
        }
        return true;
    }

    /* ---- the wiring, checked end to end ------------------------------ */
    /* Four presses whose EFFECT is asserted, not just survived.  Each one is
     * a wiring bug this editor has actually had (or has just been given a
     * fix for): typing, a Ctrl chord, the two-step chord whose second key
     * must not also be typed, and the Alt chord that WM_SYSKEYDOWN routing
     * silently killed. */
    void sweep_wiring() {
        log->line("== wiring ==");

        state(0);
        log->line("[wiring] type 'a'");
        chr('a');
        paint();
        if (ed->eq.latex().find('a') == std::string::npos)
            fail("typing 'a' did not reach the equation");

        state(0);
        log->line("[wiring] Ctrl+F");
        key('F', true, false, false);
        paint();
        if (ed->eq.latex().find("\\dfrac") == std::string::npos)
            fail("Ctrl+F did not insert a fraction");

        state(0);
        log->line("[wiring] Ctrl+T, S with its WM_CHAR");
        key('T', true, false, false);
        key('S', false, false, false);
        paint();
        {
            std::string t = ed->eq.latex();
            if (t.find("\\sum") == std::string::npos)
                fail("Ctrl+T, S did not insert a summation");
            const size_t at = t.find("\\sum");
            if (at != std::string::npos) t.erase(at, 4);
            if (t.find('s') != std::string::npos)
                fail("the chord's second key was also typed");
        }

        state(0);
        log->line("[wiring] Ctrl+Alt+Space through WM_SYSKEYDOWN");
        key(VK_SPACE, true, false, true);
        paint();
        if (ed->eq.latex().empty())
            fail("Ctrl+Alt+Space never reached the chord table");
    }

    /* ---- every published chord from every caret state ----------------- */
    void sweep_chords() {
        const std::vector<Chord>& table = chords();
        log->line("== chords: %u x 4 states ==", unsigned(table.size()));
        for (int s = 0; s < 4; ++s) {
            for (size_t i = 0; i < table.size(); ++i) {
                const Chord& c = table[i];
                if (!allowed(c)) {
                    log->line("[chord %u/%u] state=%s %s (skipped: modal/clipboard)",
                              unsigned(i + 1), unsigned(table.size()),
                              state_name(s), c.command.c_str());
                    continue;
                }
                state(s);
                log->line("[chord %u/%u] state=%s %s",
                          unsigned(i + 1), unsigned(table.size()),
                          state_name(s), c.command.c_str());
                for (const Step& st : c.steps)
                    key(st.vk, st.ctrl, st.shift, st.alt);
                paint();
                (void)ed->eq.latex();   /* the serializer survives too */
            }
        }
    }

    /* ---- the window's own keys: zoom, and (opted in) the clipboard ---- */
    void sweep_window_keys() {
        log->line("== window keys ==");
        for (int s = 0; s < 4; ++s) {
            state(s);
            log->line("[winkey] state=%s zoom 1/2/4 and the wheel", state_name(s));
            key('1', true, false, false);
            paint();
            key('4', true, false, false);
            paint();
            key('2', true, false, false);
            paint();
            wheel(+3);
            paint();
            wheel(-30);   /* clamps at the floor */
            paint();
            wheel(+40);   /* clamps at the ceiling */
            paint();
            wheel(-5);
            paint();
            if (opt->clipboard) {
                log->line("[winkey] state=%s copy / paste / cut", state_name(s));
                key('C', true, false, false);
                paint();
                key('V', true, false, false);
                paint();
                key('X', true, false, false);
                paint();
                key(VK_RETURN, true, false, false);
                paint();
            }
        }
    }

    /* ---- every palette cell, through the real popup -------------------- */
    void sweep_palette() {
        /* Two states, not four: a cell inserts the same thing from any caret,
         * and every cell is already every group's worth of popup paints. */
        for (int s = 0; s < 4; s += 2) {
            state(s);
            const size_t buttons = ed->bar.size();
            log->line("== palette: %u buttons, state=%s ==",
                      unsigned(buttons), state_name(s));
            for (size_t b = 0; b < buttons; ++b) {
                const RECT rc = ed->bar[b].rc;
                const int bx = (rc.left + rc.right) / 2;
                const int by = (rc.top + rc.bottom) / 2;

                /* Open once to hover every cell and count them. */
                click(bx, by);
                if (!ed->popup) {
                    fail("a palette button opened no popup");
                    continue;
                }
                const size_t cells = ed->cells.size();
                for (size_t j = 0; j < cells; ++j) {
                    const RECT cc = ed->cells[j].rc;
                    SendMessageW(ed->popup, WM_MOUSEMOVE, 0,
                                 MAKELPARAM((cc.left + cc.right) / 2,
                                            (cc.top + cc.bottom) / 2));
                }
                paint();
                /* A click outside the grid closes it -- the escape route a
                 * person uses when the popup was the wrong one. */
                SendMessageW(ed->popup, WM_LBUTTONDOWN, MK_LBUTTON,
                             MAKELPARAM(-100, -100));

                for (size_t j = 0; j < cells; ++j) {
                    state(s);
                    log->line("[palette] state=%s button=%u cell=%u",
                              state_name(s), unsigned(b), unsigned(j));
                    click(bx, by);
                    if (!ed->popup || j >= ed->cells.size()) {
                        fail("the popup changed shape between opens");
                        break;
                    }
                    const RECT cc = ed->cells[j].rc;
                    SendMessageW(ed->popup, WM_LBUTTONDOWN, MK_LBUTTON,
                                 MAKELPARAM((cc.left + cc.right) / 2,
                                            (cc.top + cc.bottom) / 2));
                    paint();
                }
            }
        }
    }

    /* ---- the mouse on the canvas -------------------------------------- */
    void sweep_mouse() {
        log->line("== mouse ==");
        RECT rc;
        GetClientRect(hwnd, &rc);
        const int W = rc.right, H = rc.bottom;
        const int bar = ed->bar_height();
        state(1);
        const int xs[] = {-20, 0, 24, W / 2, W - 1, W + 40};
        const int ys[] = {-10, 0, bar - 1, bar + 1, (bar + H) / 2, H - 1, H + 30};
        for (int y : ys) {
            for (int x : xs) {
                log->line("[mouse] click (%d,%d)", x, y);
                click(x, y);
                paint();
                log->line("[mouse] shift-click (%d,%d)", x, y);
                click(x, y, true);
                paint();
            }
        }
        log->line("[mouse] drags");
        state(1);
        drag(24, bar + 10, W - 30, H - 30);
        paint();
        drag(W - 30, H - 30, 24, bar + 10);
        paint();
        drag(W / 2, (bar + H) / 2, W / 2, (bar + H) / 2);
        paint();
        state(2);
        drag(-30, -30, W + 30, H + 30);
        paint();
    }

    /* ---- blink, resize, minimize, DPI ---------------------------------- */
    void sweep_environment() {
        log->line("== environment ==");
        state(1);
        log->line("[env] caret blink x3");
        blink();
        paint();
        blink();
        paint();
        blink();
        paint();
        static const int sizes[][2] = {
            {200, 120}, {90, 40}, {1400, 700}, {760, 260}};
        for (const auto& wh : sizes) {
            log->line("[env] resize %dx%d", wh[0], wh[1]);
            resize(wh[0], wh[1]);
            paint();
        }
        log->line("[env] minimize, paint while iconic, restore");
        ShowWindow(hwnd, SW_MINIMIZE);
        pump();
        paint();
        ShowWindow(hwnd, SW_SHOWNOACTIVATE);
        pump();
        paint();
        static const unsigned dpis[] = {144, 192, 96};
        for (unsigned d : dpis) {
            log->line("[env] dpi %u", d);
            dpi(d);
            paint();
        }
        resize(760, 260);
    }

    /* ---- a seeded editing session -------------------------------------- */
    void walk(unsigned seed, unsigned steps) {
        Fuzz rng(seed);
        state(int(seed % 4));
        const unsigned long gdi0 =
            GetGuiResources(GetCurrentProcess(), GR_GDIOBJECTS);
        log->line("== walk seed=%u steps=%u gdi=%lu ==", seed, steps, gdi0);

        static const unsigned navs[] = {
            VK_LEFT, VK_RIGHT, VK_UP,     VK_DOWN,   VK_HOME,  VK_END,
            VK_TAB,  VK_BACK,  VK_DELETE, VK_INSERT, VK_RETURN, VK_PRIOR,
            VK_NEXT};
        static const unsigned cps[] = {
            'a', 'b', 'c', 'x', 'y', 'z', '0', '1', '2', '+', '-', '=',
            '(', ')', '!', '?', ',', '.',
            0x3042,   /* HIRAGANA A: the Japanese path */
            0x03B1,   /* alpha */
            0x2202};  /* partial */

        const std::vector<Chord>& table = chords();
        std::vector<const Chord*> ok;
        ok.reserve(table.size());
        for (const Chord& c : table)
            if (allowed(c)) ok.push_back(&c);

        RECT rc;
        for (unsigned i = 0; i < steps; ++i) {
            GetClientRect(hwnd, &rc);
            const unsigned r = rng.pick(100);
            if (r < 30) {
                const unsigned cp =
                    cps[rng.pick(unsigned(sizeof(cps) / sizeof(cps[0])))];
                log->line("[walk %u:%u] char U+%04X", seed, i, cp);
                chr(cp);
            } else if (r < 50) {
                const unsigned vk =
                    navs[rng.pick(unsigned(sizeof(navs) / sizeof(navs[0])))];
                const bool shift = rng.pick(3) == 0;
                const bool ctrl = rng.pick(4) == 0;
                if (modal(vk, ctrl, shift) ||
                    (!opt->clipboard && clipboardy(vk, ctrl, shift))) {
                    log->line("[walk %u:%u] key skipped", seed, i);
                    continue;
                }
                log->line("[walk %u:%u] key vk=0x%02X ctrl=%d shift=%d",
                          seed, i, vk, int(ctrl), int(shift));
                key(vk, ctrl, shift, false);
            } else if (r < 70) {
                const Chord& c = *ok[rng.pick(unsigned(ok.size()))];
                log->line("[walk %u:%u] chord %s", seed, i, c.command.c_str());
                for (const Step& st : c.steps)
                    key(st.vk, st.ctrl, st.shift, st.alt);
            } else if (r < 80) {
                if (ed->bar.empty()) continue;
                const size_t b = rng.pick(unsigned(ed->bar.size()));
                const RECT brc = ed->bar[b].rc;
                log->line("[walk %u:%u] palette button %u", seed, i,
                          unsigned(b));
                click((brc.left + brc.right) / 2, (brc.top + brc.bottom) / 2);
                if (ed->popup && !ed->cells.empty()) {
                    const size_t j = rng.pick(unsigned(ed->cells.size()));
                    const RECT cc = ed->cells[j].rc;
                    log->line("[walk %u:%u] palette cell %u", seed, i,
                              unsigned(j));
                    SendMessageW(ed->popup, WM_LBUTTONDOWN, MK_LBUTTON,
                                 MAKELPARAM((cc.left + cc.right) / 2,
                                            (cc.top + cc.bottom) / 2));
                }
            } else if (r < 90) {
                const int x = int(rng.pick(unsigned(rc.right + 100))) - 50;
                const int y = int(rng.pick(unsigned(rc.bottom + 100))) - 50;
                if (rng.pick(2)) {
                    log->line("[walk %u:%u] click (%d,%d)", seed, i, x, y);
                    click(x, y, rng.pick(3) == 0);
                } else {
                    const int x2 = int(rng.pick(unsigned(rc.right + 100))) - 50;
                    const int y2 = int(rng.pick(unsigned(rc.bottom + 100))) - 50;
                    log->line("[walk %u:%u] drag (%d,%d)->(%d,%d)", seed, i, x,
                              y, x2, y2);
                    drag(x, y, x2, y2);
                }
            } else if (r < 94) {
                const int n = int(rng.pick(7)) - 3;
                log->line("[walk %u:%u] wheel %d", seed, i, n);
                wheel(n);
            } else if (r < 97) {
                log->line("[walk %u:%u] blink", seed, i);
                blink();
            } else {
                static const int sz[][2] = {
                    {760, 260}, {420, 200}, {1100, 520}, {240, 140}};
                const auto& wh = sz[rng.pick(4)];
                log->line("[walk %u:%u] resize %dx%d", seed, i, wh[0], wh[1]);
                resize(wh[0], wh[1]);
            }
            paint();

            if (i % 100 == 99) {
                const unsigned long gdi =
                    GetGuiResources(GetCurrentProcess(), GR_GDIOBJECTS);
                const unsigned long usr =
                    GetGuiResources(GetCurrentProcess(), GR_USEROBJECTS);
                log->line("[walk %u:%u] gdi=%lu user=%lu latex=%u", seed, i,
                          gdi, usr, unsigned(ed->eq.latex().size()));
                /* The slow death: a paint path that leaks a handle per frame
                 * hits the 10,000-object GDI limit after a few minutes of
                 * blinking, and dies somewhere unrelated. */
                if (gdi > gdi0 + 512)
                    fail("GDI object count grew by over 512: a paint leak");
                if (ed->eq.latex().size() > 4000) {
                    log->line("[walk %u:%u] equation reset (too large)", seed,
                              i);
                    state(int(seed % 4));
                }
            }
        }
        SendMessageW(hwnd, WM_LBUTTONUP, 0, MAKELPARAM(0, 0));
        close_popup(*ed);
    }
};

}  // namespace

int run_window_selftest(const SelftestOptions& opt) {
    SetProcessDpiAwarenessContext(DPI_AWARENESS_CONTEXT_PER_MONITOR_AWARE_V2);

    std::wstring log_path = opt.log_path;
    if (log_path.empty()) {
        wchar_t tmp[MAX_PATH] = L"";
        GetTempPathW(MAX_PATH, tmp);
        log_path = std::wstring(tmp) + L"eqnedt64_selftest.log";
    }
    Journal log;
    /* No journal, no run: a crash with no last line is the old situation. */
    if (_wfopen_s(&log.f, log_path.c_str(), L"w") != 0 || !log.f) return 2;

#ifdef EQNEDT64_VERSION
    log.line("eqnedt64 --selftest  build %s", EQNEDT64_VERSION);
#else
    log.line("eqnedt64 --selftest  (development build)");
#endif
    log.line("walks=%u walk_steps=%u clipboard=%d", opt.walks, opt.walk_steps,
             int(opt.clipboard));

    if (!register_window_classes()) {
        log.line("FAIL RegisterClassExW");
        return 2;
    }

    Editor ed;
    ed.eq.load_latex("");
    HWND hwnd = CreateWindowExW(0, kClassName, L"EQNEDT64 selftest",
                                WS_OVERLAPPEDWINDOW, CW_USEDEFAULT,
                                CW_USEDEFAULT, 760, 260, nullptr, nullptr,
                                GetModuleHandleW(nullptr), &ed);
    if (!hwnd) {
        log.line("FAIL CreateWindowExW");
        return 2;
    }
    ed.dpi = int(GetDpiForWindow(hwnd));
    build_bar(ed);
    /* Visible so painting is honest, NOACTIVATE so nobody's typing is
     * stolen: every message is injected by handle, so focus is never
     * needed -- which is what lets this run beside a working user. */
    ShowWindow(hwnd, SW_SHOWNOACTIVATE);

    Driver d;
    d.hwnd = hwnd;
    d.ed = &ed;
    d.log = &log;
    d.opt = &opt;

    const unsigned long gdi0 =
        GetGuiResources(GetCurrentProcess(), GR_GDIOBJECTS);

    d.pump();
    d.paint();
    d.sweep_wiring();
    d.sweep_chords();
    d.sweep_window_keys();
    d.sweep_palette();
    d.sweep_mouse();
    d.sweep_environment();
    for (unsigned s = 1; s <= opt.walks; ++s) d.walk(s, opt.walk_steps);

    /* DestroyWindow, not WM_CLOSE: the close path asks about unsaved
     * changes, and a message box has no place in a headless run. */
    log.line("== shutdown ==");
    DestroyWindow(hwnd);
    d.pump();

    const unsigned long gdi1 =
        GetGuiResources(GetCurrentProcess(), GR_GDIOBJECTS);
    log.line("%s failures=%d gdi_start=%lu gdi_end=%lu",
             d.failures == 0 ? "PASS" : "FAIL", d.failures, gdi0, gdi1);
    return d.failures;
}

}  // namespace mtef
