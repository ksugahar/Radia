/*
 * math_writer.h -- one tree walk, several serialisations
 *
 * Office takes an equation in three different spellings, and which one is
 * needed depends on where it is going.  All three were measured rather than
 * inferred, by copying an equation out of Word and reading the clipboard:
 *
 *   OMML   inside a .docx / .pptx      <m:f><m:num/><m:den/></m:f>
 *   RTF    clipboard, Word only        {\mf{\mfPr{\mctrlPr}}{\mnum ...}{\mden ...}}
 *   MathML clipboard, all of Office    <mfrac><mi>a</mi><mi>b</mi></mfrac>
 *
 * The structure is the same in every case; only the spelling differs.  So the
 * interesting part -- deciding what the tree means -- is written once here and
 * each output supplies a `MathSyntax`.  Duplicating the walk would guarantee
 * they drift apart, and the walk is where the judgement lives: MTEF stores a
 * script's base as the *preceding* sibling and a big operator's operand as the
 * *following* run, so both have to be absorbed or Office draws a placeholder
 * box where the content should be.
 *
 * The interface is deliberately semantic (fraction, radical, script) rather
 * than a set of element names.  OMML and RTF happen to share element names, but
 * MathML does not: its root takes its arguments the other way round, its
 * delimiters are ordinary operators inside a row, and its runs are split into
 * identifier / number / operator.  Naming the meaning instead of the tag keeps
 * that difference inside the one sink it belongs to.
 */
#ifndef MATH_WRITER_H
#define MATH_WRITER_H

#include "mtef_node.h"

#include <cstdint>
#include <string>
#include <vector>

namespace mtef {

class MathSyntax {
public:
    virtual ~MathSyntax() = default;

    /* One run of characters.  style: 0 upright, 1 italic, 2 bold italic. */
    virtual std::string run(const std::string& utf8, int style) const = 0;

    /* A slot's contents, grouped so the parent can treat it as one argument. */
    virtual std::string row(const std::string& inner) const = 0;

    virtual std::string fraction(const std::string& num, const std::string& den,
                                 bool slashed) const = 0;

    /* `index` is meaningful only when has_index. */
    virtual std::string radical(const std::string& body,
                                const std::string& index, bool has_index) const = 0;

    virtual std::string script(const std::string& base, const std::string& sub,
                               const std::string& sup,
                               bool has_sub, bool has_sup) const = 0;

    /* `beg` / `end` may be null for the default parentheses, or empty for a
     * one-sided fence. */
    virtual std::string fence(const std::string& body,
                              const char* beg, const char* end) const = 0;

    /* A big operator.  `stacked` puts the limits above and below rather than
     * beside. */
    virtual std::string nary(uint32_t chr, const std::string& lower,
                             const std::string& upper, const std::string& body,
                             bool stacked, bool has_lower, bool has_upper) const = 0;

    virtual std::string matrix(int rows, int cols,
                               const std::vector<std::string>& cells) const = 0;

    /* Stacked lines: cases, an aligned block. */
    virtual std::string stack(const std::vector<std::string>& lines) const = 0;

    /* `chr` is the combining accent, 0 for the default circumflex. */
    virtual std::string accent(uint32_t chr, const std::string& body) const = 0;

    virtual std::string bar(const std::string& body, bool over) const = 0;

    virtual std::string document(const std::string& inner, bool display) const = 0;
};

/* Walk a tree and serialise it through `syntax`.
 *
 * `run_passes` repairs EQNEDT32's sibling layout and belongs only to trees that
 * came from MTEF; a tree from the LaTeX parser already has every slot filled
 * and must be walked with the passes off. */
std::string write_math(const LineNode& root, const MathSyntax& syntax,
                       bool display, bool run_passes);

/* An embellishment's accent, in the two spellings it is needed in.
 *
 * Markup wants the COMBINING character, because Office and MathML position it
 * themselves.  Drawing wants the SPACING one, because a combining character
 * handed to TextOut has nowhere to combine with and lands as a stray mark or
 * as nothing.  Both live here so they cannot come to disagree about which
 * accent an embellishment is; 0 means "no glyph", which for a bar means draw a
 * rule instead. */
uint32_t accent_char(int embellType);          /* combining, for markup  */
uint32_t accent_drawing_char(int embellType);  /* spacing, for a picture */

}  // namespace mtef

#endif /* MATH_WRITER_H */
