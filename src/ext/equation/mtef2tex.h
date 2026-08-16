/*
 * mtef2tex.h -- MTEF v3 binary -> LaTeX
 *
 *   char *latex = mtef_to_latex_c(data, len);
 *   if (latex) { ... free(latex); }
 *
 * The returned string is malloc'd UTF-8 and is the caller's to free.
 */

#ifndef MTEF2TEX_H
#define MTEF2TEX_H

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Convert MTEF binary data to LaTeX.
 *
 * data: MTEF bytes, with or without the Equation Native header
 * len:  length of that buffer
 *
 * Returns a malloc'd UTF-8 LaTeX string, or NULL on failure.  The caller
 * frees it.
 */
char *mtef_to_latex_c(const uint8_t *data, size_t len);

#ifdef __cplusplus
}
#endif

#endif /* MTEF2TEX_H */
