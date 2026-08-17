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

namespace mtef {

/* An enhanced metafile of the equation, as bytes.  Empty on failure. */
std::string render_emf(const Layout& layout, const SvgStyle& style);

/* A PNG of the equation, as bytes.  `scale` multiplies the point size, so 4.0
 * gives roughly 288 dpi -- enough that a slide does not show the pixels. */
std::string render_png(const Layout& layout, const SvgStyle& style,
                       double scale = 4.0);

std::string tex_to_emf(const std::string& latex,
                       const SvgStyle& style = SvgStyle());
std::string tex_to_png(const std::string& latex,
                       const SvgStyle& style = SvgStyle(), double scale = 4.0);

}  // namespace mtef

#endif /* MTEF_GDI_H */
