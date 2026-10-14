#ifndef MATH_SYMBOLS_H
#define MATH_SYMBOLS_H

#include <cstdint>
#include <string>
#include <vector>

namespace eqnedit {

/* Unicode code point for a TeX symbol command including its leading slash;
 * -1 when the command names a structure rather than a single glyph. */
int latex_symbol_codepoint(const std::string& command);

/* True when the Unicode scalar has a named TeX symbol in the shared table. */
bool is_latex_symbol_codepoint(uint32_t codepoint);

/* Every command in the shared table, so a test can sweep the whole symbol
 * set instead of the handful a hand-written corpus happens to mention. */
std::vector<std::string> latex_symbol_commands();

/* Logical typeface for a symbol code point: TF_UCGREEK, TF_LCGREEK, or
 * TF_SYMBOL.  The variant Greek letters live outside the contiguous Greek
 * block, so they have to be named individually or they lose the italic that
 * every other lowercase Greek letter gets. */
int typeface_for_code(uint32_t codepoint);

}  // namespace eqnedit

#endif
