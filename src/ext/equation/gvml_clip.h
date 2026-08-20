/*
 * gvml_clip.h -- the clipboard format PowerPoint uses for its own shapes
 *
 * An equation pasted into PowerPoint used to come out at whatever size the
 * destination box happened to be -- 18 pt in a body placeholder, which is too
 * small to read on a slide.  MathML cannot say otherwise: `mathsize` on the
 * <math> element and on an <mstyle> were both measured and both ignored, and
 * the size that arrived was the placeholder's own.
 *
 * PowerPoint's own copy says it, though.  Copying a shape out of PowerPoint
 * puts "Art::GVML ClipFormat" on the clipboard: an OPC package -- a ZIP --
 * holding a DrawingML `lockedCanvas`, whose text carries an explicit
 * `sz="2400"`.  So this writes that, the same way the RTF here was written by
 * transcribing what Word puts on the clipboard.
 *
 * What arrives is a real, editable PowerPoint equation at the size asked for,
 * not a picture: the OMML goes inside the <a14:m> wrapper a slide uses, which
 * is the same wrapper markdown_to_pptx writes into a deck.
 *
 * The package is written with stored (uncompressed) entries -- valid ZIP,
 * accepted by PowerPoint (measured), and no deflate implementation to carry
 * for a payload of a few kilobytes.
 */
#ifndef MTEF_GVML_CLIP_H
#define MTEF_GVML_CLIP_H

#include <string>

namespace mtef {

/* The size a pasted equation lands at, in points.
 *
 * 24 pt because 18 is too small to read from the back of a room, which is
 * what a slide is for.  PowerPoint's own default body text is 18 pt at the
 * first outline level and shrinks from there, so an equation that merely
 * inherits it is always at the small end of the deck. */
extern const double kPasteSizePt;

/* One equation as an Art::GVML ClipFormat package.  Empty if the LaTeX does
 * not parse into anything. */
std::string tex_to_gvml(const std::string& latex,
                        double size_pt = kPasteSizePt,
                        bool display = false);

}  // namespace mtef

#endif /* MTEF_GVML_CLIP_H */
