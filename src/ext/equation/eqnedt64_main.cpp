/*
 * eqnedt64_main.cpp -- the application
 *
 * Nothing but an entry point.  The window is a library function so the same
 * editor opens from a Python call, from a Jupyter button, or from here, and
 * there is one editor rather than one per host.
 *
 * Usage:  eqnedt64.exe [file.tex | initial LaTeX]
 * An argument naming a file that exists is opened; anything else is taken as
 * the equation itself.  So a .tex opens by double-click and a one-liner still
 * works from a shell.
 * On exit the equation is printed, so a shell or a script can take it back.
 */
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>

#include "eq_window.h"

#include <cstdio>
#include <string>

static std::wstring widen_utf8(const std::string& s) {
    int n = MultiByteToWideChar(CP_UTF8, 0, s.data(), int(s.size()), nullptr, 0);
    std::wstring w;
    w.resize(size_t(n));
    MultiByteToWideChar(CP_UTF8, 0, s.data(), int(s.size()), &w[0], n);
    return w;
}

int WINAPI wWinMain(HINSTANCE, HINSTANCE, PWSTR cmdline, int) {
    std::string latex;
    if (cmdline && *cmdline) {
        int n = WideCharToMultiByte(CP_UTF8, 0, cmdline, -1, nullptr, 0,
                                    nullptr, nullptr);
        if (n > 1) {
            latex.resize(size_t(n) - 1);
            WideCharToMultiByte(CP_UTF8, 0, cmdline, -1, &latex[0], n,
                                nullptr, nullptr);
        }
    }

    /* Strip the quotes Explorer puts around a dropped path. */
    if (latex.size() > 1 && latex.front() == '"' && latex.back() == '"')
        latex = latex.substr(1, latex.size() - 2);

    std::wstring open_path;
    if (!latex.empty()) {
        const std::wstring w = widen_utf8(latex);
        const DWORD attr = GetFileAttributesW(w.c_str());
        if (attr != INVALID_FILE_ATTRIBUTES &&
            !(attr & FILE_ATTRIBUTE_DIRECTORY)) {
            open_path = w;
            latex.clear();
        }
    }

    mtef::EditorResult r = mtef::run_equation_window(latex, open_path);

    /* Attach to the launching console when there is one, so running this from
     * a shell shows the result and running it from Explorer stays silent. */
    if (AttachConsole(ATTACH_PARENT_PROCESS)) {
        FILE* out = nullptr;
        if (freopen_s(&out, "CONOUT$", "w", stdout) == 0 && out) {
            std::fputs(r.latex.c_str(), out);
            std::fputc('\n', out);
        }
        FreeConsole();
    }
    return r.copied ? 0 : 1;
}
