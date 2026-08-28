/*
 * tex_document.h -- TeX file envelope and paste normalization
 *
 * Eqnedit64 stores TeX.  This module deliberately knows nothing about the
 * editor window or the clipboard so opening, saving, and pasting can be tested
 * in a background process.
 */
#ifndef TEX_DOCUMENT_H
#define TEX_DOCUMENT_H

#include <string>

namespace eqnedit {

struct TexDocument {
    std::string body;
    std::string prefix;
    std::string suffix;
    bool numbered = true;
    bool hadEquationEnvironment = false;
};

/* Read the first equation/equation* environment in a TeX file.  When there is
 * no such environment, the complete input is treated as a math fragment. */
TexDocument parse_tex_document(const std::string& text);

/* Remove only clipboard-level display wrappers.  An aligned environment is
 * retained because it carries the equation's row/column structure. */
std::string normalize_tex_paste(const std::string& text);

/* Write one equation environment while retaining any surrounding TeX text
 * that was present when the file was opened. */
std::string compose_tex_document(const std::string& body, bool numbered,
                                 const std::string& prefix = {},
                                 const std::string& suffix = {});

}  // namespace eqnedit

#endif
