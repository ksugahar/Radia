# Third-party notices

Eqnedit64's application icon incorporates the **Superscript** glyph from
[Bootstrap Icons](https://icons.getbootstrap.com/icons/superscript/).

Bootstrap Icons is Copyright (c) 2019-2025 The Bootstrap Authors and is
licensed under the MIT License. The full license is available at
<https://github.com/twbs/icons/blob/main/LICENSE>.

The icon's background, colors, composition, and rasterization are specific to
Eqnedit64.

## Latin Modern Math

Eqnedit64 embeds **Latin Modern Math** (`latinmodern-math.otf`), the OpenType
Computer Modern maintained by GUST, the Polish TeX Users Group.  It is the
typeface TeX itself sets with, which is why the canvas draws the same shapes
as a pdflatex run rather than merely the same geometry.

Latin Modern is licensed under the **GUST Font License** (a LaTeX Project
Public License variant); the full text is in `assets/GUST-FONT-LICENSE.txt`.
The font is unmodified.  It is loaded for the running process only, with
`AddFontMemResourceEx`, and is never installed on the machine.
