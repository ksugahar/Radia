/*
 * eqnedt64_main.cpp -- the application
 *
 * Nothing but an entry point.  The window is a library function so the same
 * editor opens from a Python call, from a Jupyter button, or from here, and
 * there is one editor rather than one per host.
 *
 * Usage:  eqnedt64.exe [initial LaTeX]
 * On exit the equation is printed, so a shell or a script can take it back.
 */
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>

#include "eq_window.h"

#include <cstdio>
#include <string>

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

    mtef::EditorResult r = mtef::run_equation_window(latex);

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
