/* TeX math -> structural equation tree.  Every slot is filled directly. */
#ifndef TEX_PARSER_H
#define TEX_PARSER_H

#include "equation_node.h"

#include <memory>
#include <string>

namespace eqnedit {

/* Parse LaTeX math into a node tree.  Surrounding $...$ / $$...$$ / \[...\]
 * delimiters are accepted and stripped. Unsupported colour names/models return
 * null and an optional diagnostic. Other unknown control words still follow
 * the historical permissive fallback; do not assume strict TeX validation. */
std::unique_ptr<LineNode> parse_latex(const std::string& latex,
                                      bool* depthExceeded = nullptr,
                                      std::string* error = nullptr);

}  // namespace eqnedit

#endif /* TEX_PARSER_H */
