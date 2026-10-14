# Third-party notices

Eqnedit64's application icon incorporates the **Superscript** glyph from
[Bootstrap Icons](https://icons.getbootstrap.com/icons/superscript/).

Bootstrap Icons is Copyright (c) 2019-2025 The Bootstrap Authors and is
licensed under the MIT License. The full license is available at
<https://github.com/twbs/icons/blob/main/LICENSE>.

The icon's background, colors, composition, and rasterization are specific to
Eqnedit64.

## Latin Modern Math

Eqnedit64 embeds **Eqnedit Math** (`eqnedit-math.ttf`), a modified TrueType-outline
derivative of Latin Modern Math, maintained upstream by GUST, the Polish TeX
Users Group. The original `latinmodern-math.otf` remains in the source tree
for reproducible conversion and diagnosis; neither native resource embeds it.

Latin Modern is licensed under the **GUST Font License** (a LaTeX Project
Public License variant); the full text is in `assets/GUST-FONT-LICENSE.txt`.
The derivative converts cubic CFF outlines to quadratic glyf outlines with
fontTools 4.60.1, corrects TrueType bearings, and changes the family name.
Glyph order, advances, cmap, MATH, GSUB and GPOS are preserved. Curve conversion
and integer rounding mean rasterization is not claimed to be pixel-identical.
See `MANIFEST-eqnedit-math.txt` for provenance and regeneration.
The executable verifies and extracts the derivative to a
content-addressed per-user cache, then loads it for the running process only
with `AddFontResourceExW(FR_PRIVATE | FR_NOT_ENUM)`.  It is never registered
as an installed Windows font and writes no font-registry entry.
