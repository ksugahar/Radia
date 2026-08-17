/*
 * md_layout.h -- a Markdown document, laid out for viewing
 *
 * The same separation as math_layout.h, one level up: laying the document out
 * and drawing it are different jobs, and only the first one has any judgement
 * in it.  The result is a display list of text runs and *equations*, where an
 * equation carries its own math layout rather than a picture of one.
 *
 * That is the point of the whole exercise: the equation on screen and the
 * equation pasted into Word come from the same layout, so they cannot disagree.
 * A viewer that renders maths through a browser engine has two renderers and
 * they always drift.
 *
 * It also gives the two things an editor needs beyond drawing: which block a
 * point is in (so a click can open that cell) and whether a point is inside an
 * equation (so a click can open the equation widget instead).
 */
#ifndef MD_LAYOUT_H
#define MD_LAYOUT_H

#include "math_layout.h"
#include "md_blocks.h"

#include <string>
#include <vector>

namespace mtef {

struct DocStyle {
    double body = 11.0;            /* pt */
    double heading[6] = {20.0, 17.0, 14.5, 12.5, 11.5, 11.0};
    double mono = 9.5;
    double line_spacing = 1.35;    /* multiple of the run's own height */
    double para_gap = 0.6;         /* blank space between blocks, in body ems */
    double list_indent = 1.6;      /* in body ems */
    double margin = 8.0;           /* pt around the whole document */

    std::string text_font = "Yu Gothic UI";
    std::string mono_font = "Consolas";

    /* The maths is set relative to the surrounding text, so a heading's
     * equation is as big as the heading. */
    double math_scale = 1.0;
};

struct DocRun {
    std::string text;              /* UTF-8 */
    double x = 0, baseline = 0;
    double size = 0;               /* pt */
    bool bold = false, italic = false, mono = false;
};

struct DocMath {
    Layout layout;                 /* the equation's own display list */
    double x = 0, baseline = 0;
    std::string latex;             /* so a click can open the widget on it */
    int block = -1;                /* which block it came from */
    int index = -1;                /* which equation in the document */
    bool display = false;
};

/* One block's extent, for turning a click into a cell. */
struct DocBlockBox {
    int block = -1;
    MdBlock::Kind kind = MdBlock::kParagraph;
    double top = 0, bottom = 0;
};

struct DocLayout {
    double width = 0, height = 0;
    std::vector<DocRun> runs;
    std::vector<DocMath> maths;
    std::vector<Rule> rules;       /* code backgrounds, heading rules */
    std::vector<DocBlockBox> blocks;
};

/* Lay a Markdown document out to `width` points.  Windows only, like the maths
 * layout it embeds: the metrics come from GDI. */
DocLayout layout_markdown(const std::string& markdown, double width,
                          const DocStyle& style = DocStyle());

/* Which block a point falls in, or -1. */
int block_at(const DocLayout& doc, double x, double y);

/* The equation a point falls in, or -1 -- the test for "open the equation
 * widget rather than the text editor". */
int math_at(const DocLayout& doc, double x, double y);

}  // namespace mtef

#endif /* MD_LAYOUT_H */
