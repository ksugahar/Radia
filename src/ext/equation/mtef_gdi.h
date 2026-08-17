/*
 * mtef_gdi.h -- draw a laid-out equation with GDI: metafile and bitmap
 *
 * One drawing routine, three uses.  The same calls that put an equation on a
 * window's device context also record an enhanced metafile and rasterise a
 * bitmap, because in GDI those are all just device contexts.
 *
 * Which one a target needs was measured:
 *
 *   Word            takes RTF as a native equation             (no picture)
 *   PowerPoint      takes MathML as a native equation          (no picture)
 *   Excel           takes a pasted Word equation as a picture, keeping the
 *                   metafile alongside the bitmap -- so EMF is not lossy there
 *   Google Slides   has no equation object at all, rejects SVG on upload, and
 *                   accepts raster; EMF gets in via Google Drawings if the
 *                   vector matters
 *
 * So the metafile is the better payload wherever a picture is unavoidable, and
 * the bitmap is the one that always works.  The bitmap is produced by playing
 * the metafile, so there is one drawing routine and no separate encoder.
 *
 * Windows only: this is GDI, as is the font metric layer the layout rests on.
 */
#ifndef MTEF_GDI_H
#define MTEF_GDI_H

#include "math_layout.h"

#include <cstddef>
#include <cstdint>
#include <string>

#ifdef _WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX      /* windows.h defines min/max as macros and breaks std:: */
#endif
#include <windows.h>
#endif

namespace mtef {

#ifdef _WIN32
/* Draw the display list onto a device context.
 *
 * `units_per_pt` says how many device units go to a point -- a metafile records
 * at its own fixed rate, a bitmap at that times its scale, and a window passes
 * its DPI over 72 to get pixels.  The equation's left edge, on its baseline,
 * lands at (originX, originY).
 *
 * This is the routine the metafile and bitmap writers use, so a window and a
 * paste are drawn by the same code and cannot disagree about the equation. */
void draw_layout(HDC hdc, const Layout& layout, const SvgStyle& style,
                 double units_per_pt, int originX, int originY,
                 COLORREF colour, bool show_empty_slots = false);
#endif

/* An enhanced metafile of the equation, as bytes.  Empty on failure. */
std::string render_emf(const Layout& layout, const SvgStyle& style);

/* A PNG of the equation, as bytes.  `scale` multiplies the point size, so 4.0
 * gives roughly 288 dpi -- enough that a slide does not show the pixels. */
std::string render_png(const Layout& layout, const SvgStyle& style,
                       double scale = 4.0);

/* The same raster as a packed device-independent bitmap: a BITMAPINFOHEADER
 * followed by bottom-up pixels, which is what CF_DIB wants.
 *
 * Needed because Windows synthesises CF_DIB from a bitmap but NOT from a
 * metafile, and an application pasting a picture reads CF_DIB.  Offering only
 * the named "PNG" format left a paste into a browser taking the plain text
 * instead -- so an equation dropped into Google Slides, which has no equation
 * object and takes raster only, arrived as its LaTeX. */
std::string render_dib(const Layout& layout, const SvgStyle& style,
                       double scale = 4.0);
std::string tex_to_dib(const std::string& latex,
                       const SvgStyle& style = SvgStyle(), double scale = 4.0);

std::string tex_to_emf(const std::string& latex,
                       const SvgStyle& style = SvgStyle());
std::string tex_to_png(const std::string& latex,
                       const SvgStyle& style = SvgStyle(), double scale = 4.0);

}  // namespace mtef

#endif /* MTEF_GDI_H */
