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
EditorResult run_equation_window(const std::string& latex);

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

}  // namespace mtef

#endif /* EQ_WINDOW_H */
