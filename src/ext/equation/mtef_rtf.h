/*
 * mtef_rtf.h -- an equation as RTF, which is how it reaches Office by paste
 *
 * Word offers "Rich Text Format" on the clipboard when an equation is copied,
 * and the maths inside it is OMML spelled as control words:
 *
 *     {\mf{\mfPr{\mctrlPr}}{\mnum{\mr\mscr0\msty2 a}}{\mden{\mr\mscr0\msty2 b}}}
 *
 * That correspondence was measured by copying each construct out of Word
 * rather than inferred, so the writer transcribes rather than guesses.  The
 * structure walk is shared with the OMML output (see math_writer.h).
 *
 * `rtf_document` wraps the maths in the smallest complete RTF file, which is
 * what a clipboard payload has to be: Word will not accept a bare fragment.
 */
#ifndef MTEF_RTF_H
#define MTEF_RTF_H

#include "mtef_node.h"

#include <cstddef>
#include <cstdint>
#include <string>

namespace mtef {

struct RtfOptions {
    /* A display equation stands on its own line; an inline one sits in the
     * surrounding text. */
    bool display = false;
    /* Point size of the equation, doubled as RTF's \fsN wants half-points. */
    double font_size_pt = 11.0;
};

/* Just the maths destination, for embedding in a larger RTF document. */
std::string render_rtf_math(const LineNode& root,
                            const RtfOptions& opt = RtfOptions(),
                            bool run_passes = true);

/* A complete, minimal RTF document containing one equation -- the form the
 * clipboard needs. */
std::string tex_to_rtf(const std::string& latex,
                       const RtfOptions& opt = RtfOptions());

}  // namespace mtef

#endif /* MTEF_RTF_H */
