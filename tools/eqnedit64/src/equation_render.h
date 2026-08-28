/* Shared native-canvas and SVG renderer for the TeX equation tree. */
#ifndef EQUATION_RENDER_H
#define EQUATION_RENDER_H

#include "equation_node.h"
#include <string>

#ifdef _WIN32
#ifndef _WINDEF_
struct HDC__;
using HDC = HDC__*;
#endif
#endif

namespace eqnedit {

struct SvgStyle {
    /* Full, script, second-level script, operator and small-operator sizes.
     * The script sizes are TeX's: a 12 pt body sets scripts at 8 pt and
     * second-level scripts at 6 pt.  The 7/5 this used to ship made every
     * exponent noticeably smaller than the same equation typeset by LaTeX. */
    double full = 12.0;
    double sub = 8.0;
    double sub2 = 6.0;
    double sym = 18.0;
    double subsym = 12.0;

    /* Font families emitted into the SVG.  A list is used so a viewer that
     * lacks the first still finds the glyph. */
    std::string serif = "Latin Modern Math, Cambria Math, serif";
    std::string symbol = "Latin Modern Math, Cambria Math, Symbol, serif";
    std::string cjk = "Yu Mincho, Yu Gothic UI, Meiryo, MS Mincho, serif";

    double padding = 1.0;   /* pt of white space around the equation */
};

/* Geometry shared by SVG export and the native editor canvas.  Coordinates
 * are points relative to the top-left of the rendered equation. */
struct RenderMetrics {
    double width = 0;
    double height = 0;
    double baseline = 0;
};

struct CaretGeometry {
    double x = 0;
    double top = 0;
    double bottom = 0;
};

RenderMetrics measure_equation(const LineNode& root,
                               const SvgStyle& style = SvgStyle());

/* Find the structurally nearest insertion site.  `x` and `y` use the same
 * top-left point coordinate system as RenderMetrics. */
bool hit_test_equation(const LineNode& root, double x, double y,
                       const SvgStyle& style,
                       const NodeList** slot, int* index);

bool caret_geometry_equation(const LineNode& root, const NodeList* slot,
                             int index, const SvgStyle& style,
                             CaretGeometry* geometry);

#ifdef _WIN32
/* Paint the native editor surface with the exact display list used for SVG.
 * `scale` is device pixels per point.  Caret and selection pointers may be
 * null for an export/preview without editing chrome. */
/* Benchmark hook: turning the drawing-font cache off restores the old
 * behaviour of creating and destroying an HFONT per glyph, so the two paths
 * can be timed against each other in one binary instead of one being
 * asserted to be faster. */
void set_draw_font_cache_enabled(bool enabled);

void draw_equation_gdi(const LineNode& root, HDC hdc,
                       double left, double top, double scale,
                       const SvgStyle& style,
                       const NodeList* caret_slot = nullptr,
                       int caret_index = -1,
                       const NodeList* selection_slot = nullptr,
                       int selection_first = -1,
                       int selection_last = -1,
                       bool show_placeholders = true,
                       bool show_caret = true,
                       /* Record glyph outlines instead of text, so a metafile
                        * does not depend on the reader having the font. */
                       bool text_as_outlines = false);
#endif

/* False when the embedded math font could not be loaded and GDI is
 * substituting; every measurement is then subtly wrong. */
bool math_font_loaded();

/* Load the embedded math font now (idempotent); cheap when already loaded. */
void ensure_math_font_ready();

/* Render a parsed equation to a standalone SVG document. */
std::string render_svg(const LineNode& root, const SvgStyle& style = SvgStyle());

/* TeX fragment -> SVG. */
std::string tex_to_svg(const std::string& latex,
                       const SvgStyle& style = SvgStyle());

}  // namespace eqnedit

#endif  /* EQUATION_RENDER_H */
