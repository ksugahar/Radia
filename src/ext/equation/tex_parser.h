/*
 * tex_parser.h -- LaTeX math -> shared node tree
 *
 * LaTeX is the working format: an equation is stored as LaTeX inside a
 * Markdown file, and every output (OMML for Office, SVG for display) is
 * rendered from the tree this parser builds.  MTEF is not on this path at all;
 * it survives only to read legacy Equation Editor objects, and that direction
 * goes MTEF -> LaTeX through the corpus-validated mtef2tex converter.
 *
 * The tree is the one mtef_node.h was designed for -- its header already
 * names "LaTeXParser (text->tree)" as the second producer -- so the OMML
 * emitter and the SVG renderer serve both sources.
 *
 * Unlike the MTEF parser this fills every slot directly: a script carries its
 * own base, a big operator carries its own limits and body.  The LINE pass
 * pipeline reconstructs those from EQNEDT32's sibling layout and must NOT be
 * run on a tree from here.
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
