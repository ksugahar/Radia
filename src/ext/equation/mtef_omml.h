/*
 * mtef_omml.h -- node tree -> OMML (Office Math Markup Language)
 *
 * The Office-native equation format (Word/PowerPoint 2007+).  Unlike the SVG
 * renderer this needs no layout at all: Office does the typesetting, so the
 * emitter is a structural mapping and inherits Office's spacing, fonts and
 * colour-follows-text behaviour for free.
 *
 * Equations pasted this way are NOT pictures and NOT OLE objects -- they are
 * editable with PowerPoint's own equation tools on any machine, with no
 * Equation Editor installed.
 *
 * The reverse direction already exists as python/eqnedt32/omml2tex.py, which
 * makes a round-trip the natural test.
 */
#ifndef MTEF_OMML_H
#define MTEF_OMML_H

#include "mtef_node.h"
#include <string>

namespace mtef {

struct OmmlOptions {
    /* Wrap the result in <m:oMathPara> (a display equation on its own line)
     * rather than a bare inline <m:oMath>. */
    bool display = false;
    /* Emit the m: namespace declaration on the root element.  Needed when the
     * fragment is inserted standalone; PowerPoint's a14:m wrapper already
     * declares it. */
    bool declare_namespace = true;
    /* Style of the run properties: OMML marks upright/italic per run. */
    bool italic_variables = true;
};

/* `run_passes` repairs EQNEDT32's sibling layout and belongs only to trees
 * that came from MTEF.  A tree from the LaTeX parser already has every slot
 * filled and must be rendered with the passes off. */
std::string render_omml(const LineNode& root, const OmmlOptions& opt = OmmlOptions(),
                        bool run_passes = true);

/* LaTeX -> OMML.  This is the working path: LaTeX is the stored form of an
 * equation and Office receives a native, editable equation.  No MTEF involved. */
std::string tex_to_omml(const std::string& latex,
                        const OmmlOptions& opt = OmmlOptions());

/* MTEF binary -> OMML.  Empty string when the MTEF cannot be parsed. */
std::string mtef_to_omml(const uint8_t* data, size_t len,
                         const OmmlOptions& opt = OmmlOptions());

}  // namespace mtef

#endif /* MTEF_OMML_H */
