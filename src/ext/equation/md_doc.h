/*
 * md_doc.h -- the Markdown document, seen as text with equations in it
 *
 * A Markdown editor that also edits equations needs exactly one thing from the
 * Markdown side: to know which spans of the file are math and which are not.
 * Everything else -- headings, emphasis, lists -- is the text editor's usual
 * business and is deliberately not modelled here.
 *
 * The scan is conservative rather than clever, because the cost of a false
 * positive is high: a `$` inside a fenced code block or an inline code span is
 * NOT math, and treating it as math would silently rewrite someone's code when
 * the file is saved.  Fences and code spans are therefore skipped wholesale,
 * and the inline `$...$` rule follows Pandoc's (no space just inside the
 * delimiters, no digit right after the closing one) so that "$5 and $6" stays
 * prose.
 *
 * Serialisation is exact.  Loading and saving an untouched file reproduces it
 * byte for byte -- the delimiters each equation actually used are kept, so a
 * file written with \( \) does not come back as $ $.
 */
#ifndef MD_DOC_H
#define MD_DOC_H

#include <string>
#include <vector>

namespace mtef {

struct MdSegment {
    enum Kind {
        kText,          /* ordinary Markdown, passed through untouched  */
        kInlineMath,    /* $...$  or \(...\)                            */
        kDisplayMath,   /* $$...$$ or \[...\]                           */
        kCodeSpan,      /* `...`                                        */
        kCodeBlock,     /* ```...``` or ~~~...~~~                       */
    };

    Kind kind = kText;
    std::string open;    /* the exact opening delimiter, or empty        */
    std::string body;    /* the content between the delimiters           */
    std::string close;   /* the exact closing delimiter, or empty        */

    bool is_math() const { return kind == kInlineMath || kind == kDisplayMath; }
    std::string source() const { return open + body + close; }
};

class MarkdownDoc {
public:
    void load(const std::string& markdown);
    std::string text() const;                 /* exact round trip */

    const std::vector<MdSegment>& segments() const { return segs_; }

    /* Equations, in document order. */
    int math_count() const;
    std::string math_latex(int i) const;
    bool set_math_latex(int i, const std::string& latex);
    bool math_is_display(int i) const;

    /* Segment index of the i-th equation, for an editor that needs to map a
     * caret position to a segment. */
    int math_segment_index(int i) const;

private:
    std::vector<MdSegment> segs_;
};

}  // namespace mtef

#endif /* MD_DOC_H */
