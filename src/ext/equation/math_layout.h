/*
 * math_layout.h -- the laid-out equation, shared by every way of drawing it
 *
 * Laying an equation out and drawing it are separate jobs.  The layout is a
 * display list: glyphs and rules positioned relative to the equation's origin,
 * which sits ON THE BASELINE at its left edge, plus the box extents.  Nothing
 * about SVG, GDI or a window appears in it.
 *
 * That separation is what lets one layout serve four purposes: an SVG picture,
 * an enhanced metafile for the clipboard, a bitmap for the targets that take
 * nothing else, and -- when the editor exists -- the drawing on screen and the
 * caret position, which is the reason the layout has to be a first-class thing
 * rather than a step inside a renderer.
 */
#ifndef MATH_LAYOUT_H
#define MATH_LAYOUT_H

#include "mtef_node.h"
#include "mtef_svg.h"          /* SvgStyle: the type sizes, shared by all backends */

#include <cstdint>
#include <string>
#include <vector>

namespace mtef {

struct Glyph {
    double x = 0, y = 0;        /* y is the baseline */
    double size = 0;            /* pt */
    bool italic = false;
    bool symbol = false;        /* draw with the maths face rather than the text one */
    bool cjk = false;           /* draw with the Japanese face                       */
    double stretchY = 1.0;      /* vertical scale, for a fence grown to its content */
    std::string text;           /* UTF-8 */
    /* When set, the math font has a drawing of exactly this size and the
     * backend should ask for that glyph rather than scale the character: a
     * stretched radical thins its own stroke and bends its hook, which is how
     * a bar comes to float above a sign it no longer meets.  `text` still
     * holds the character, so a backend that cannot address glyphs directly
     * has something to draw and the geometry stays the same either way. */
    uint16_t glyph_id = 0;
};

/* True for the ranges a Latin font has no glyphs for.  `symbol` and `cjk` are
 * mutually exclusive and both are decided from the character's own code point,
 * so a glyph's face is a property of the character rather than something a
 * caller chooses -- and the width it is measured with is the width it is drawn
 * with.  Getting that wrong drew Japanese as correctly spaced blank paper. */
bool is_cjk(uint32_t cp);

struct Rule { double x = 0, y = 0, w = 0, h = 0; };   /* y is the TOP edge */

/* Where an insertion point sits.  The editing model addresses the caret as a
 * path of (child, slot) steps plus an index in the slot it lands in; the layout
 * knows where things are.  These join the two, which is what lets an editor
 * draw the caret and turn a click into a position.
 *
 * `path` uses the same spelling as Equation::caret() minus the ":index", and
 * `slot` numbering must agree with node_slots() -- a mismatch puts the cursor
 * somewhere other than where it was clicked, which is the kind of bug that is
 * very hard to see and very annoying to use. */
struct CaretStop {
    std::string path;      /* "" for the outermost line */
    int index = 0;         /* insertion index within that slot */
    double x = 0;          /* relative to the equation origin */
    double top = 0;        /* the caret's extent, baseline-relative */
    double bottom = 0;
};

/* An empty slot: the numerator of a fraction nobody has typed into yet.
 *
 * Such a slot has no extent of its own, so a renderer that draws only what is
 * there shows a fresh fraction as a bare bar floating in space, and gives no
 * clue where Tab is about to go.  Equation Editor drew these as dotted boxes,
 * and that is most of what made its structure legible.
 *
 * The layout REPORTS them and leaves the drawing to the caller, because the two
 * callers want opposite things: the editor must show them, and a picture on its
 * way to a slide must not. */
struct SlotBox { double x = 0, y = 0, w = 0, h = 0; };   /* y is the TOP edge */

struct Layout {
    double w = 0, asc = 0, desc = 0;
    std::vector<Glyph> glyphs;
    std::vector<Rule> rules;
    std::vector<SlotBox> empty_slots;

    std::vector<CaretStop> stops;

    void translate(double dx, double dy) {
        for (auto& g : glyphs) { g.x += dx; g.y += dy; }
        for (auto& r : rules)  { r.x += dx; r.y += dy; }
        for (auto& b : empty_slots) { b.x += dx; b.y += dy; }
        for (auto& s : stops)  { s.x += dx; s.top += dy; s.bottom += dy; }
    }
    void absorb(const Layout& other, double dx, double dy) {
        Layout t = other;
        t.translate(dx, dy);
        glyphs.insert(glyphs.end(), t.glyphs.begin(), t.glyphs.end());
        rules.insert(rules.end(), t.rules.begin(), t.rules.end());
        empty_slots.insert(empty_slots.end(),
                           t.empty_slots.begin(), t.empty_slots.end());
        stops.insert(stops.end(), t.stops.begin(), t.stops.end());
    }
};

/* Lay an equation out.  Font metrics come from GDI, so this is Windows-only;
 * the OMML, RTF and MathML paths are not. */
Layout layout_math(const LineNode& root, const SvgStyle& style);

/* The caret at `path`:`index`, or null when the layout has no such position --
 * which happens for a construct the layout does not draw yet. */
const CaretStop* find_stop(const Layout& layout, const std::string& path, int index);

/* The insertion point a click at (x, y) means.  Coordinates are relative to the
 * equation origin, so a caller subtracts the padding and the baseline first.
 * Null only for an empty layout. */
const CaretStop* nearest_stop(const Layout& layout, double x, double y);

}  // namespace mtef

#endif /* MATH_LAYOUT_H */
