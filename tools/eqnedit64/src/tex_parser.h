/* TeX math -> structural equation tree.  Every slot is filled directly. */
#ifndef TEX_PARSER_H
#define TEX_PARSER_H

#include "equation_node.h"

#include <memory>
#include <string>

namespace eqnedit {

/* Parse LaTeX math into a node tree.  Surrounding $...$ / $$...$$ / \[...\]
 * delimiters are accepted and stripped.  Never returns null: unknown commands
 * become literal text rather than aborting the parse, so a typo costs one
 * wrong glyph instead of the whole equation. */
std::unique_ptr<LineNode> parse_latex(const std::string& latex,
                                      bool* depthExceeded = nullptr);

}  // namespace eqnedit

#endif /* TEX_PARSER_H */
