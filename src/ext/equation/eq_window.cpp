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
#include "mtef_gdi.h"
#include "mtef_mathml.h"
#include "mtef_omml.h"
#include "mtef_rtf.h"
#include "mtef_svg.h"

#include <algorithm>
#include <cmath>
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
const int kBtnW = 46, kBtnH = 24;      /* device-independent pixels */
const int kCellW = 34, kCellH = 30;
const int kBarPad = 3;
const double kBarPt = 9.0;             /* type size inside a button   */
const double kCellPt = 12.0;           /* type size inside a cell     */

/* What a cell shows is what inserting it produces: the sample is rendered by
 * actually performing the insertion into a scratch equation.  A hand-written
 * table of sample LaTeX would drift from what the templates really are. */
Layout sample_layout(const Equation::PaletteItem& item, double sizePt) {
    SvgStyle st;
    const double k = sizePt / 12.0;
    st.full = 12.0 * k; st.sub = 7.0 * k; st.sub2 = 5.0 * k;
    st.sym = 18.0 * k;  st.subsym = 12.0 * k;
    st.padding = 0.0;

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
    SvgStyle style;
    bool copied = false;
    bool caret_on = true;
    int dpi = 96;
    std::vector<Step> pending;      /* a chord waiting for its second key */
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
        int(((H - bar_h) + L.asc * upp - L.desc * upp) / 2);

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

bool handle_key(HWND hwnd, Editor& ed, UINT vk) {
    Step cur;
    cur.vk = vk;
    cur.ctrl  = (GetKeyState(VK_CONTROL) & 0x8000) != 0;
    cur.shift = (GetKeyState(VK_SHIFT)   & 0x8000) != 0;
    cur.alt   = (GetKeyState(VK_MENU)    & 0x8000) != 0;

    if (!ed.pending.empty()) {
        std::vector<Step> want = ed.pending;
        want.push_back(cur);
        ed.pending.clear();
        for (const Chord& c : chords()) {
            if (c.steps == want) { ed.eq.command(c.command); redraw(hwnd); return true; }
        }
        return true;         /* the prefix was consumed either way */
    }

    for (const Chord& c : chords()) {
        if (c.steps.empty() || !(c.steps[0] == cur)) continue;
        if (c.steps.size() == 1) {
            ed.eq.command(c.command);
            redraw(hwnd);
        } else {
            ed.pending.push_back(cur);
        }
        return true;
    }
    return false;
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
                int(((rc.bottom - bar_h) + L.asc * upp - L.desc * upp) / 2);
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
                int(((rc.bottom - bar_h) + L.asc * upp - L.desc * upp) / 2);
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
            if (handle_key(hwnd, *ed, UINT(wp))) ed->swallow_char = true;
            return 0;
        }

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

}  // namespace

bool copy_equation_to_clipboard(const std::string& latex, bool display,
                                bool pictures) {
    RtfOptions rtf;
    rtf.display = display;
    MathMLOptions mml;
    mml.display = display;

    const std::string rtf_bytes = tex_to_rtf(latex, rtf);
    const std::string mml_bytes = tex_to_mathml(latex, mml);

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

    WNDCLASSEXW wc = {};
    wc.cbSize = sizeof(wc);
    wc.lpfnWndProc = proc;
    wc.hInstance = GetModuleHandleW(nullptr);
    wc.hCursor = LoadCursorW(nullptr, IDC_IBEAM);
    wc.lpszClassName = kClassName;
    RegisterClassExW(&wc);

    WNDCLASSEXW pc = {};
    pc.cbSize = sizeof(pc);
    pc.lpfnWndProc = popup_proc;
    pc.hInstance = wc.hInstance;
    pc.hCursor = LoadCursorW(nullptr, IDC_ARROW);
    pc.lpszClassName = kPopupClass;
    RegisterClassExW(&pc);

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
        nullptr, nullptr, wc.hInstance, &ed);
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

}  // namespace mtef
