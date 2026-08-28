/*
 * tex_parser.h -- LaTeX math -> shared node tree
 *
 * LaTeX is the working format: an equation is stored as LaTeX inside a
 * Markdown file, and every output (OMML for Office, SVG for display) is
 * rendered from the tree this parser builds.  Retired binary equation formats
 * are not accepted by this path.
 *
 * The tree is the one mtef_node.h was designed for -- its header already
 * names "LaTeXParser (text->tree)" as the second producer -- so the OMML
 * emitter and the SVG renderer serve both sources.
 *
 * Every slot is filled directly: a script carries its own base, and a big
 * operator carries its own limits and body.  Historical repair passes must not
 * be run on a tree from here.
 */
#ifndef TEX_PARSER_H
#define TEX_PARSER_H

#include "mtef_node.h"

#include <memory>
#include <string>

namespace mtef {

/* Parse LaTeX math into a node tree.  Surrounding $...$ / $$...$$ / \[...\]
 * delimiters are accepted and stripped.  Never returns null: unknown commands
 * become literal text rather than aborting the parse, so a typo costs one
 * wrong glyph instead of the whole equation. */
std::unique_ptr<LineNode> parse_latex(const std::string& latex);

}  // namespace mtef

#endif /* TEX_PARSER_H */
