/*
 * mtef_svg.h -- MTEF node tree -> SVG renderer
 *
 * The third consumer of the unified node tree (mtef_node.h), next to
 * latex_emitter (tree -> LaTeX) and tex2mtef (tree -> MTEF).  Its purpose is
 * to make EQNEDT32 unnecessary for PRODUCING pictures: the editor stays as the
 * GUI a human uses to author an equation, while every programmatic path --
 * LaTeX -> MTEF -> picture -> slide -- runs in our own code.
 *
 * Sizes follow EQNEDT32's Sizes dialog (Full 12 / Subscript 7 /
 * Sub-subscript 5 / Symbol 18 / Sub-symbol 12 pt) and are parameters, not
 * constants, so a deck can be typeset at its own scale.
 */
#ifndef MTEF_SVG_H
#define MTEF_SVG_H

#include "mtef_node.h"
#include <string>

namespace mtef {

struct SvgStyle {
    /* EQNEDT32 "Sizes" dialog, in points. */
    double full = 12.0;
    double sub = 7.0;
    double sub2 = 5.0;
    double sym = 18.0;
    double subsym = 12.0;

    /* Font families emitted into the SVG.  A list is used so a viewer that
     * lacks the first still finds the glyph. */
    std::string serif = "Times New Roman, Times, serif";
    std::string symbol = "Cambria Math, Symbol, Segoe UI Symbol, serif";
    /* Neither of the two above contains a single kana, so Japanese needs a
     * face of its own rather than whatever a viewer happens to substitute. */
    std::string cjk = "Yu Mincho, MS Mincho, serif";

    double padding = 1.0;   /* pt of white space around the equation */
};

/* Render a parsed equation to a standalone SVG document. */
std::string render_svg(const LineNode& root, const SvgStyle& style = SvgStyle());

/* Convenience: MTEF binary -> SVG.  Returns an empty string when the MTEF
 * cannot be parsed. */
std::string mtef_to_svg(const uint8_t* data, size_t len,
                        const SvgStyle& style = SvgStyle());

/* LaTeX -> SVG.  The working path: it parses LaTeX straight into the tree
 * rather than going through MTEF, so the renderer sees filled base/body/limit
 * slots instead of EQNEDT32's sibling layout. */
std::string tex_to_svg(const std::string& latex,
                       const SvgStyle& style = SvgStyle());

}  // namespace mtef

#endif  /* MTEF_SVG_H */
