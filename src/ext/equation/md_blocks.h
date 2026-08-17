/*
 * md_blocks.h -- the block structure of a Markdown file
 *
 * `MarkdownDoc` answers which *spans* are maths.  This answers which *blocks*
 * the file is made of, which is what a viewer needs in order to lay it out and
 * what a writer needs in order to choose a paragraph style.
 *
 * Deliberately shallow.  A viewer for technical notes needs headings,
 * paragraphs, fenced code, and the two kinds of list; tables, block quotes,
 * footnotes and images are not modelled, and their source passes through as
 * ordinary text rather than being dropped.  Guessing at more structure than is
 * really understood produces a viewer that silently renders something other
 * than what the file says.
 *
 * This lives in C++ because the viewer cannot afford Python: importing the
 * package costs about 1.4 s, against 4 ms for this module itself, and an
 * equation editor that takes a second and a half to open is not the tool this
 * is trying to be.  Having it here also keeps one implementation -- the Office
 * writers call the same scanner rather than carrying a second copy that would
 * drift.
 */
#ifndef MD_BLOCKS_H
#define MD_BLOCKS_H

#include <string>
#include <vector>

namespace mtef {

struct MdBlock {
    enum Kind {
        kParagraph,
        kHeading,      /* `level` is 1..6                          */
        kBullet,       /* one item; `level` is the indent depth     */
        kNumbered,     /* one item; `level` is the indent depth     */
        kCode,         /* a fenced block, fences already removed    */
        kBlank,        /* a run of empty lines, kept so a round trip is exact */
    };

    Kind kind = kParagraph;
    int level = 0;
    std::string text;      /* the content, without the marker or the fences */
    std::string info;      /* a code fence's language, when it gave one     */
    std::string source;    /* the original lines, so the file can be rebuilt */
};

/* Split a Markdown file into blocks.  Concatenating every `source` reproduces
 * the input exactly. */
std::vector<MdBlock> md_blocks(const std::string& markdown);

}  // namespace mtef

#endif /* MD_BLOCKS_H */
