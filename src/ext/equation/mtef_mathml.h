/*
 * mtef_mathml.h -- MathML, the one clipboard format the whole of Office takes
 *
 * Measured: with only the `MathML` clipboard format present, PowerPoint turns
 * a paste into a native equation.  RTF does not -- it arrives as a text box --
 * and RTF is therefore Word's route only.  MathML is also a W3C standard, so
 * it costs 142 bytes where Word's HTML costs 42 KB and DrawingML would have to
 * be reverse-engineered.
 *
 * What Word itself writes, for reference:
 *
 *   \frac{a}{b}      <mfrac><mi>a</mi><mi>b</mi></mfrac>
 *   \sqrt[3]{x}      <mroot><mi>x</mi><mn>3</mn></mroot>     (body first)
 *   \left(x\right)   <mrow><mo>(</mo><mi>x</mi><mo>)</mo></mrow>
 *   \sum_{i}^{n} a   <mrow><munderover><mo stretchy="false">S</mo>
 *                    <mi>i</mi><mi>n</mi></munderover><mi>a</mi></mrow>
 */
#ifndef MTEF_MATHML_H
#define MTEF_MATHML_H

#include "mtef_node.h"

#include <cstddef>
#include <cstdint>
#include <string>

namespace mtef {

struct MathMLOptions {
    /* display="block" for a standalone equation, "inline" within text. */
    bool display = true;
    /* Emit the MathML namespace on <math>.  Needed for a standalone payload;
     * a host that already declares it can turn it off. */
    bool declare_namespace = true;
};

std::string render_mathml(const LineNode& root,
                          const MathMLOptions& opt = MathMLOptions(),
                          bool run_passes = true);

std::string tex_to_mathml(const std::string& latex,
                          const MathMLOptions& opt = MathMLOptions());

}  // namespace mtef

#endif /* MTEF_MATHML_H */
