"""Equations as LaTeX: edit them, render them, put them into Office.

An equation is stored as LaTeX, usually inside a Markdown file.  Everything
here works from that:

    >>> import radia.equation as eq
    >>> eq.tex_to_omml(r"\\frac{a+b}{c}")        # Office-native math
    >>> eq.tex_to_svg(r"\\frac{a+b}{c}")         # a picture
    >>> doc = eq.MarkdownDoc(); doc.load(text)   # which spans are math
    >>> e = eq.Equation(); e.insert_template("frac")

`Equation` is the editing model, not a widget: an insertion point that lives
inside the structure, templates whose empty slots you tab through, and a
backspace that unwraps a template rather than swallowing it.  A front end binds
keys to the command names in `Equation.shortcuts()`; the model itself knows
nothing about keyboards.

Office receives OMML -- a native, editable equation, not a picture and not an
OLE object -- so nothing needs to be installed on the reader's machine.

MTEF and `.eqn` are retired formats.  This package deliberately exposes no
reader, writer, or converter for them; preserved documents must be migrated to
TeX before entering the supported workflow.
"""

from radia._equation import (  # noqa: F401
    DocBlockBox,
    DocLayout,
    DocMath,
    DocRun,
    DocStyle,
    Equation,
    MarkdownDoc,
    MathMLOptions,
    MdBlock,
    MdSegment,
    OmmlOptions,
    PaletteGroup,
    PaletteItem,
    RtfOptions,
    SvgStyle,
    layout_markdown,
    md_blocks,
    symbol_palettes,
    template_palettes,
    tex_empty_slots,
    tex_metrics,
    AtomKind,
    atom_kind,
    atom_space_mu,
    math_constants,
    math_glyph,
    math_stretch,
    math_variant_for_height,
    tex_to_dib,
    tex_normalize,
    tex_to_emf,
    tex_to_mathml,
    tex_to_omml,
    tex_to_png,
    tex_to_rtf,
    tex_to_svg,
)

from radia.equation.office import (  # noqa: F401
    markdown_to_docx,
    markdown_to_pptx,
    omml_paragraph,
    split_math,
)

__all__ = [
    "DocBlockBox",
    "DocLayout",
    "DocMath",
    "DocRun",
    "DocStyle",
    "Equation",
    "MarkdownDoc",
    "MathMLOptions",
    "MdBlock",
    "MdSegment",
    "OmmlOptions",
    "PaletteGroup",
    "PaletteItem",
    "RtfOptions",
    "SvgStyle",
    "layout_markdown",
    "markdown_to_docx",
    "markdown_to_pptx",
    "md_blocks",
    "symbol_palettes",
    "template_palettes",
    "tex_empty_slots",
    "tex_metrics",
    "AtomKind",
    "atom_kind",
    "atom_space_mu",
    "math_constants",
    "math_glyph",
    "math_stretch",
    "math_variant_for_height",
    "tex_to_dib",
    "omml_paragraph",
    "split_math",
    "tex_normalize",
    "tex_to_emf",
    "tex_to_mathml",
    "tex_to_omml",
    "tex_to_png",
    "tex_to_rtf",
    "tex_to_svg",
]
