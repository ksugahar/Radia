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

MTEF, the Equation Editor 3.x / MathType binary format, is supported in one
direction for a reason: to read equations out of documents that already
contain them.  `mtef_to_tex` and `mtef_to_omml` are the import path;
`tex_to_mtef` exists so an equation can be handed back to Equation Editor as a
.eqn file.  Nothing else in this package goes through MTEF.
"""

from radia._equation import (  # noqa: F401
    Equation,
    MarkdownDoc,
    MdSegment,
    OmmlOptions,
    SvgStyle,
    dump_tree,
    mtef_to_omml,
    mtef_to_svg,
    mtef_to_tex,
    read_eqn,
    tex_dump_tree,
    tex_normalize,
    tex_to_mtef,
    tex_to_omml,
    tex_to_svg,
    write_eqn,
)

from radia.equation.office import (  # noqa: F401
    markdown_to_docx,
    omml_paragraph,
    split_math,
)

__all__ = [
    "Equation",
    "MarkdownDoc",
    "MdSegment",
    "OmmlOptions",
    "SvgStyle",
    "dump_tree",
    "markdown_to_docx",
    "mtef_to_omml",
    "mtef_to_svg",
    "mtef_to_tex",
    "omml_paragraph",
    "read_eqn",
    "split_math",
    "tex_dump_tree",
    "tex_normalize",
    "tex_to_mtef",
    "tex_to_omml",
    "tex_to_svg",
    "write_eqn",
]
