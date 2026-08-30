/* Structural equation tree -> Presentation MathML for Microsoft Office. */
#ifndef MATHML_EMITTER_H
#define MATHML_EMITTER_H

#include "equation_node.h"

#include <string>

namespace eqnedit {

/* The root mathsize is deliberately explicit. The Office clipboard caller
 * uses the accepted 18 pt left-aligned paste contract; render/export callers
 * may continue to request a 24 pt source size. */
std::string tree_to_mathml(const LineNode& root, double pointSize = 24.0);
std::string latex_to_mathml(const std::string& latex,
                            double pointSize = 24.0);

/* Office's HTML MathML importer renders MathML alignment elements as visible
 * ampersands and centres one-column mtables. For an outer, unanchored aligned
 * environment, recursively flatten nested unanchored aligned wrappers and
 * publish one editable inline MathML root per leaf row separated by HTML line
 * breaks. Other input remains one MathML root. */
std::string latex_to_office_mathml_fragment(
    const std::string& latex, double pointSize = 18.0);

}  // namespace eqnedit

#endif  // MATHML_EMITTER_H
