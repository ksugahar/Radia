/*
 * math_writer.h -- one tree walk, several serialisations
 *
 * Office writes an equation in two places with the same structure and
 * different spelling.  In a .docx it is OMML:
 *
 *     <m:f><m:fPr/><m:num>...</m:num><m:den>...</m:den></m:f>
 *
 * and on the clipboard it is RTF, which is the same element names as control
 * words -- measured, not guessed, by copying an equation out of Word:
 *
 *     {\mf{\mfPr{\mctrlPr}}{\mnum ...}{\mden ...}}
 *
 * Because the structures coincide, the interesting part -- deciding what the
 * tree means -- is written once here, and each output supplies only a
 * `MathSyntax` that knows how to spell an element, a property and a run.
 * Duplicating the walk instead would guarantee the two drift apart.
 *
 * The walk is not a plain recursion: MTEF stores a script's base as the
 * preceding sibling and a big operator's operand as the following run, so
 * both have to be absorbed or Office draws a placeholder box where the
 * content should be.
 */
#ifndef MATH_WRITER_H
#define MATH_WRITER_H

#include "mtef_node.h"

#include <string>

namespace mtef {

/* How one output spells the pieces.  Element names are passed without a
 * prefix ("f", "num", "sSubSup"); the syntax adds its own. */
class MathSyntax {
public:
    virtual ~MathSyntax() = default;

    /* An element with content: <m:num>inner</m:num> / {\mnum inner} */
    virtual std::string group(const char* name, const std::string& inner) const = 0;

    /* A property carrying a value: <m:chr m:val="v"/> / {\mchr v} */
    virtual std::string prop(const char* name, const std::string& value) const = 0;

    /* A property that is simply on: <m:degHide m:val="1"/> / {\mdegHide on} */
    virtual std::string flag(const char* name) const = 0;

    /* The formatting placeholder Word writes inside every property group. */
    virtual std::string ctrl() const = 0;

    /* One run of characters.  style: 0 upright, 1 italic, 2 bold italic. */
    virtual std::string run(const std::string& utf8, int style) const = 0;

    /* Wrap the finished body as a complete equation. */
    virtual std::string document(const std::string& inner, bool display) const = 0;
};

/* Walk a tree and serialise it through `syntax`.
 *
 * `run_passes` repairs EQNEDT32's sibling layout and belongs only to trees
 * that came from MTEF; a tree from the LaTeX parser already has every slot
 * filled and must be walked with the passes off. */
std::string write_math(const LineNode& root, const MathSyntax& syntax,
                       bool display, bool run_passes);

}  // namespace mtef

#endif /* MATH_WRITER_H */
