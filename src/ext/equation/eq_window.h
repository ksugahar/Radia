/*
 * eq_window.h -- the editor window: one equation, and a way out to Office
 *
 * Equation Editor 3.0 was never a document editor.  It was a small window that
 * edited one equation and handed it to whatever was waiting.  That is the shape
 * worth keeping, because everything else -- the document, the file, the list --
 * is already someone else's job and better done by them.
 *
 * What is deliberately absent is a Redraw command.  Equation Editor needed one
 * because its display could go stale; here every paint recomputes the layout
 * from the tree, so there is no stale state for a button to repair.
 *
 * Windows only.
 */
#ifndef EQ_WINDOW_H
#define EQ_WINDOW_H

#include <string>

namespace mtef {

struct EditorResult {
    bool copied = false;      /* the equation reached the clipboard */
    std::string latex;        /* what the equation was when the window closed */
};

/* Open the window on `latex` and run it until the user closes it.  Blocks. */
/* `path` opens that file instead, and becomes where Ctrl+S writes. */
EditorResult run_equation_window(const std::string& latex,
                                 const std::wstring& path = std::wstring());

/* Put one equation on the clipboard for every target at once.
 *
 * Which format each target needs was measured, not assumed:
 *
 *   Rich Text Format   Word reads it as maths
 *   MathML             PowerPoint reads it as maths
 *   CF_ENHMETAFILE     a vector picture, for anywhere with no equation object
 *   PNG                the picture that always works
 *   CF_UNICODETEXT     the LaTeX itself, for Markdown, Jupyter, any editor
 *
 * The pictures go on after the equation formats deliberately: Word and
 * PowerPoint were re-checked with them present and still produce native
 * equations, but an application that prefers a picture would silently
 * downgrade.  `pictures` is there to take them away.
 */
bool copy_equation_to_clipboard(const std::string& latex, bool display = false,
                                bool pictures = true);

/* The interaction-layer self-test: `eqnedt64.exe --selftest`.
 *
 * The window is the one layer the Python tests cannot reach -- the model
 * behind it is ordinary functions, but WM_KEYDOWN wiring, painting, the
 * palette popups and the mouse all live in a WndProc.  The first crash this
 * editor had in real use happened in exactly that gap: some fifty seconds of
 * ordinary editing, an access violation, and nobody able to say which
 * operation did it.
 *
 * So this drives the REAL window through the REAL WndProc by injecting window
 * messages -- no keyboard is touched, no foreground is stolen, so it can run
 * beside a working user and on a headless CI desktop.  Every published chord
 * and every palette cell is applied from several caret states, then seeded
 * random walks mix keys, mouse, resizes and repaints the way an editing
 * session does.  Each step is journalled and flushed BEFORE it runs, so when
 * a step crashes the process, the journal's last line names it, and the WER
 * LocalDumps entry for eqnedt64.exe (see the handover) holds the dump.
 *
 * Returns 0 when every step survived; the failure count otherwise.  A crash
 * of course returns nothing -- the exit code is the exception code, which is
 * what the pytest wrapper asserts on. */
struct SelftestOptions {
    std::wstring log_path;      /* the journal; empty puts it in %TEMP% */
    unsigned walks = 2;         /* random walks, seeded 1..walks */
    unsigned walk_steps = 1200; /* injected messages per walk */
    /* The clipboard is the user's; a test that clobbers it may not do so by
     * default.  Enabling this adds Ctrl+C / Ctrl+X / Ctrl+V round-trips. */
    bool clipboard = false;
};
int run_window_selftest(const SelftestOptions& opt);

}  // namespace mtef

#endif /* EQ_WINDOW_H */
