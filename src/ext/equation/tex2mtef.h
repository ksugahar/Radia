/*
 * tex2mtef.h -- LaTeX -> MTEF v3 binary
 *
 *   int outLen;
 *   uint8_t *mtef = tex_to_mtef(latex, &outLen);
 *   if (mtef) { ... free(mtef); }
 *
 * Only needed to hand an equation back to Equation Editor as a .eqn file;
 * everything else in this library works from LaTeX directly.
 */

#ifndef TEX2MTEF_H
#define TEX2MTEF_H

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Convert a LaTeX string to MTEF v3 binary.
 *
 * latex:  UTF-8 LaTeX, in the subset mtef2tex.cpp emits
 * outLen: receives the length of the result; may be NULL
 *
 * Returns a malloc'd MTEF buffer, or NULL on failure.  The caller frees it.
 */
uint8_t *tex_to_mtef(const char *latex, int *outLen);

/*
 * Guess whether a piece of text is LaTeX math, for deciding what to do with
 * something read off the clipboard.
 *
 * Returns 1 for probable LaTeX, 0 for ordinary text.
 */
int looks_like_latex(const char *text);

/*
 * Strip dollar delimiters in place: $$...$$ -> ..., $...$ -> ...
 * `text` must be a writable buffer.
 */
void strip_dollar_delimiters(char *text);

/*
 * Look up a LaTeX command (with leading backslash, e.g. "\nabla") in the
 * shared command->Unicode table.  Returns the code point, or -1 when the
 * command is not a symbol.  Exposed so the LaTeX parser and the MTEF writer
 * share one table instead of drifting apart.
 */
int tex_command_to_unicode(const char *cmd);

#ifdef __cplusplus
}
#endif

#endif /* TEX2MTEF_H */
