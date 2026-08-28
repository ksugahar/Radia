/* Structural equation tree -> Presentation MathML for Microsoft Office. */
#ifndef MATHML_EMITTER_H
#define MATHML_EMITTER_H

#include "equation_node.h"

#include <string>

namespace eqnedit {

/* The root mathsize is deliberately explicit: PowerPoint otherwise applies
 * its 18 pt insertion default even when Eqnedit64 is rendering at 24 pt. */
std::string tree_to_mathml(const LineNode& root, double pointSize = 24.0);
std::string latex_to_mathml(const std::string& latex,
                            double pointSize = 24.0);

}  // namespace eqnedit

#endif  // MATHML_EMITTER_H
