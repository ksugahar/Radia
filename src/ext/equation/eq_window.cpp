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

#include "eq_window.h"

#include "eq_edit.h"
#include "mtef_gdi.h"
#include "mtef_mathml.h"
#include "mtef_omml.h"
#include "mtef_rtf.h"
#include "mtef_svg.h"

#include <algorithm>
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
const double kZoom = 2.0;          /* the equation is shown larger than life,
                                    * as the old editor did at 200% */
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

/* ---- chords ------------------------------------------------------------- */

struct Step {
    UINT vk = 0;
    bool ctrl = false, shift = false, alt = false;
    bool operator==(const Step& o) const {
        return vk == o.vk && ctrl == o.ctrl && shift == o.shift && alt == o.alt;
    }
};

struct Chord {
    std::vector<Step> steps;
    std::string command;
};

UINT named_key(const std::string& name) {
    if (name == "Tab")       return VK_TAB;
    if (name == "Left")      return VK_LEFT;
    if (name == "Right")     return VK_RIGHT;
    if (name == "Up")        return VK_UP;
    if (name == "Down")      return VK_DOWN;
    if (name == "Home")      return VK_HOME;
    if (name == "End")       return VK_END;
    if (name == "Backspace") return VK_BACK;
    if (name == "Delete")    return VK_DELETE;
    if (name == "Escape")    return VK_ESCAPE;
    if (name == "Enter")     return VK_RETURN;
    return 0;
}

/* "Ctrl+T, S" -> two steps; "Shift+Tab" -> one.  Parsing the published table
 * rather than restating it keeps one source for the chords. */
bool parse_chord(const std::string& text, Chord& out) {
    size_t pos = 0;
    while (pos <= text.size()) {
        size_t comma = text.find(", ", pos);
        std::string part = text.substr(
            pos, comma == std::string::npos ? std::string::npos : comma - pos);

        Step st;
        size_t i = 0;
        for (;;) {
            size_t plus = part.find('+', i);
            /* A lone '+' at the end is the key itself, not a separator. */
            if (plus == std::string::npos || plus + 1 >= part.size()) break;
            std::string mod = part.substr(i, plus - i);
            if (mod == "Ctrl")       st.ctrl = true;
            else if (mod == "Shift") st.shift = true;
            else if (mod == "Alt")   st.alt = true;
            else return false;
            i = plus + 1;
        }
        std::string key = part.substr(i);
        if (key.empty()) return false;

        if (UINT vk = named_key(key)) {
            st.vk = vk;
        } else if (key.size() == 1) {
            /* Ask the layout, so "Ctrl+[" and "Ctrl+{" resolve to the same
             * physical key with and without shift. */
            SHORT r = VkKeyScanW(wchar_t((unsigned char)key[0]));
            if (r == -1) return false;
            st.vk = LOBYTE(r);
            if (HIBYTE(r) & 1) st.shift = true;
        } else {
            return false;
        }
        out.steps.push_back(st);

        if (comma == std::string::npos) break;
        pos = comma + 2;
    }
    return !out.steps.empty();
}

const std::vector<Chord>& chords() {
    static const std::vector<Chord> table = [] {
        std::vector<Chord> v;
        for (const Equation::Binding& b : Equation::shortcuts()) {
            Chord c;
            c.command = b.command;
            if (parse_chord(b.chord, c)) v.push_back(c);
        }
        return v;
    }();
    return table;
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

    std::vector<Button> bar;
    HWND popup = nullptr;
    std::vector<Cell> cells;        /* what the open popup is showing */
    int hot_cell = -1;

    double units_per_pt() const { return dpi / 72.0 * kZoom; }
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
    if (L.glyphs.empty() && L.rules.empty()) return;
    SvgStyle st;
    const int x = rc.left + int(((rc.right - rc.left) - L.w * upp) / 2);
    const int y = rc.top + int(((rc.bottom - rc.top) +
                                (L.asc - L.desc) * upp) / 2);
    draw_layout(dc, L, st, upp, x, y, colour);
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

    draw_layout(mem, L, ed.style, upp, originX, originY,
                GetSysColor(COLOR_WINDOWTEXT));

    if (ed.caret_on) {
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
            if (ed) paint(hwnd, *ed);
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
            ed->eq.move_to_point((p.x - originX) / upp,
                                 (p.y - originY) / upp, ed->style);
            ed->caret_on = true;
            redraw(hwnd);
            return 0;
        }

        case WM_SIZE:
            if (ed) { build_bar(*ed); redraw(hwnd); }
            return 0;

        case WM_KEYDOWN: {
            if (!ed) return 0;
            const bool ctrl = (GetKeyState(VK_CONTROL) & 0x8000) != 0;
            if (wp == VK_ESCAPE) { DestroyWindow(hwnd); return 0; }
            if (ctrl && (wp == 'C' || wp == VK_RETURN)) {
                ed->copied = copy_equation_to_clipboard(ed->eq.latex());
                return 0;
            }
            if (handle_key(hwnd, *ed, UINT(wp))) return 0;
            return 0;
        }

        case WM_CHAR: {
            if (!ed) return 0;
            const wchar_t ch = wchar_t(wp);
            if (ch < 0x20) return 0;                       /* control keys */
            if (GetKeyState(VK_CONTROL) & 0x8000) return 0; /* a chord, not text */
            ed->eq.insert_text(narrow(std::wstring(1, ch)));
            ed->caret_on = true;
            redraw(hwnd);
            return 0;
        }

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

EditorResult run_equation_window(const std::string& latex) {
    /* Per-monitor aware, so the equation is sharp on a 4K screen.  It fails
     * harmlessly when the host process already chose its awareness -- which is
     * the normal case when this is called from Python. */
    SetProcessDpiAwarenessContext(DPI_AWARENESS_CONTEXT_PER_MONITOR_AWARE_V2);

    Editor ed;
    ed.eq.load_latex(latex);
    ed.eq.move_end();

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

    HWND hwnd = CreateWindowExW(
        0, kClassName, L"EQNEDT64",
        WS_OVERLAPPEDWINDOW, CW_USEDEFAULT, CW_USEDEFAULT, 760, 260,
        nullptr, nullptr, wc.hInstance, &ed);
    if (!hwnd) return EditorResult{};

    ed.dpi = int(GetDpiForWindow(hwnd));
    build_bar(ed);
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
