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

#include <string>
#include <vector>

namespace mtef {

struct Glyph {
    double x = 0, y = 0;        /* y is the baseline */
    double size = 0;            /* pt */
    bool italic = false;
    bool symbol = false;        /* draw with the maths face rather than the text one */
    double stretchY = 1.0;      /* vertical scale, for a fence grown to its content */
    std::string text;           /* UTF-8 */
};

struct Rule { double x = 0, y = 0, w = 0, h = 0; };   /* y is the TOP edge */

struct Layout {
    double w = 0, asc = 0, desc = 0;
    std::vector<Glyph> glyphs;
    std::vector<Rule> rules;

    void translate(double dx, double dy) {
        for (auto& g : glyphs) { g.x += dx; g.y += dy; }
        for (auto& r : rules)  { r.x += dx; r.y += dy; }
    }
    void absorb(const Layout& other, double dx, double dy) {
        Layout t = other;
        t.translate(dx, dy);
        glyphs.insert(glyphs.end(), t.glyphs.begin(), t.glyphs.end());
        rules.insert(rules.end(), t.rules.begin(), t.rules.end());
    }
};

/* Lay an equation out.  Font metrics come from GDI, so this is Windows-only;
 * the OMML, RTF and MathML paths are not. */
Layout layout_math(const LineNode& root, const SvgStyle& style);

}  // namespace mtef

#endif /* MATH_LAYOUT_H */
